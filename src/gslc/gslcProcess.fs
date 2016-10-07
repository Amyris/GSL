/// Top-level compiler operations.
module gslcProcess
open System
open System.IO

// lexing, parsing, expanding
open GSLParser
open GSLLexer
open Microsoft.FSharp.Text.Lexing
open Expansion  // Serial GSL expansion logic
open docstrings // support for /// docstrings
open constants // global constants
open pragmaProcess
open parseTypes
open commonTypes

open sgdrefformat
open shared // common routines used in several modules

// output generation
open thumper // Thumper output formats, dumping etc
open ryse    // RYSE architecture
open cloneManager
open ape
open jsonAssembly
open DnaCreation
open PrimerCreation
open dumpFlat

// Helper libs for oligo design, sequence parsing all in Amyris.Bio.dll
open Amyris.Bio
open utils

/// Starting design parameters for construction
let initialDesignParams =
   {pp = primercore.defaultParams;
    targetTm = ryseLinkerTargetDefault;
    seamlessTm = seamlessTargetDefault;
    seamlessOverlapTm = seamlessTargetDefault;
    overlapParams = primercore.defaultParams;
    overlapMinLen = overlapMinLenDefault}

/// Helpful wrapper type for handing around GSLC's static assets and caches.
type GlobalAssets =
    {seqLibrary: Map<string, char []>;
     codonTableCache: plugins.CodonTableCache;
     rgs: Map<string, GenomeDef>;}

/// generate list of available reference genome folders
let enumerateLibs (opts:ParsedOptions) =
    Directory.EnumerateDirectories(opts.libDir)
        |> Seq.map (Amyris.Bio.utils.baseName) |> List.ofSeq

/// Load static assets and initialize global caches.
let loadGlobalAssets (opts:ParsedOptions) =
    let lib = opj opts.libDir "lib.fa"

    // Crude sequence library for misc pieces
    let library =
        if File.Exists lib then
            Amyris.Bio.biolib.readReference lib
            |> Seq.map (fun kv -> (kv.Key.ToUpper(),kv.Value |> basesUpper ))
            |> Map.ofSeq
        else Map.empty

    // Get ready for codon tables we will load on demand
    let codonTableCache = plugins.CodonTableCache(opts.libDir)

    if opts.verbose then printfn "opts.libDir=%s" opts.libDir

    let availRefs = enumerateLibs opts
        
    if opts.verbose then printfn "availrefs=%A" availRefs

    let rgs =
        seq {
            for s in availRefs do
                let p = opj opts.libDir s
                if not (Directory.Exists(p)) then
                    failwithf "ERROR: unable to find genome reference dir %s\n" p
                if File.Exists(opj p (sprintf "%s.fsa" s)) then
                    yield (s, new GenomeDef(p))
            } |> Map.ofSeq

    // Debugging - dump list of available genomes
    if opts.verbose then
        printf "loadedgenomes %A\n"
            (rgs |> Seq.map (fun kv -> kv.Key) |> List.ofSeq)

    {seqLibrary = library;
     codonTableCache = codonTableCache;
     rgs = rgs}


/// Validate the names of any provided genes against the reference genome and
/// provided gene library.  Raises an exception if any genes do not resolve.
let checkGeneNames
    (rgs:GenomeDefs)
    (library:Map<string,char array>)
    (aList : Assembly list) =
    let problems =
        seq {
            for a in aList do
                for ppp in a.parts do
                    yield
                        match ppp.part with
                        | GENEPART(gp) ->
                            let g = gp.part.gene.[1..].ToUpper()
                            let rg' = getRG a rgs ppp.pr
                            let ok = (rg'.IsValid(g)) || library.ContainsKey(g)
                            if not ok then Some(g) else None
                        | _ -> None
        }
        |> Seq.choose id
        |> Seq.groupBy id
        |> Seq.map (fun (g,s) ->
            let count = (Seq.length s)
            sprintf "ERROR: unknown gene '%s' found %d time%s"
                g count (if count>1 then "s" else ""))
        |> List.ofSeq

    match problems with
    | [] -> ()
    | x -> failwith (String.concat "\n" x)

let expandTree
    (rgs:GenomeDefs)
    (opts:ParsedOptions)
    (library:Map<string,char array>)
    (aList:seq<Assembly>)
    (proxyURL:string option) =

    if opts.verbose then printf "Processing %d assemblies\n" (Seq.length aList)

    let assemblies =
        aList
        |> Seq.mapi
            (fun i a ->
                {id = Some(i)
                 dnaParts = expandAssembly opts.verbose rgs library a proxyURL
                 name =
                    match a.name with
                    | None -> sprintf "A%d" i
                    | Some(s) -> s;
                 uri = a.uri;
                 linkerHint = a.linkerHint;
                 pragmas = a.pragmas;
                 designParams = a.designParams;
                 docStrings = a.docStrings})
        |> List.ofSeq

    if opts.verbose then
        printf "log: dnaParts dump\n"
        for a in assemblies do
            printf "log: dnaPart: %s\n" a.name
            for p in a.dnaParts do
                printf "log:      %s\n" p.description
                printf "%s\n" (format60 p.dna)

    // Check for reused pieces and number them accordingly
    // Make a list of all parts, determine unique set and assign ids
    let partIDs =
        seq {for a in assemblies do
                for p in a.dnaParts do
                    yield p.dna} // Base this on DNA for now, but could involve linkers down the road}
        |> Set.ofSeq
        |> Seq.mapi (fun i dna -> (dna,i))
        |> Map.ofSeq // TODO - could be faster with checksums
    // Relabel the pieces with IDs  - tedious, we have to reconstruct the tree
    List.map
        (fun a ->
            {a with dnaParts = List.map
                (fun (p:DNASlice) -> { p with id = Some(partIDs.[p.dna]) })
                a.dnaParts})
        assemblies

/// Just perform lexing on input file.  For testing/debugging purposes.
let lexTest opts inputFile =
    use inputF = File.OpenText inputFile
    let lexbuf = LexBuffer<_>.FromTextReader inputF

    let rec nextToken() =
        let t = gslTokenizer opts.verbose lexbuf

        if t = EOF then ()
        else
            printf "%A\n" t
            nextToken()
    nextToken()

let lexAndParse verbose gslText =
    // If the tree contains any parts that require expansion before final
    // translation to DNA then we should do that rather than emit final output
    let inBuffer = LexBuffer<_>.FromString gslText

    // Clean LET / alias tracking before parsing
    // TODO: this is global mutable state, inside the parser, leaking out.
    // Turn this inside-out.
    aliases.Clear()

    try
        if verbose then printfn "Starting tree parse...\n"
        let t = GSLParser.start (gslTokenizer verbose) inBuffer
        if verbose then printfn "Parsed tree!"
        t
    with e ->
        // Error hanlding from parsing/semantic errors
        let pos = inBuffer.EndPos
        if e.Message.Contains("@") then
            failwith e.Message // already has location information
        else
            //  General ERROR handling with no specific location mentioned
            // ==========================================================
            let lines = gslText.Replace("\r\n","\n").Split([| '\n' ; '\r'|])

            // Accumulate lines in an error report.
            let mutable msg = [
                sprintf "ERROR: near line %d col %d\nERROR: %s"
                    (pos.Line+1) (pos.Column+1) e.Message;
                "================================================================="]

            for line in max 0 (pos.Line-5) .. min (pos.Line+5) (lines.Length-1) do
                msg <- msg @ [sprintf "%s" lines.[line]]
                if line = pos.Line then
                    msg <- msg @ [sprintf "%s^" (pad pos.Column)]

            if verbose then printf "%s\n" e.StackTrace

            failwith (String.concat "\n" msg)

/// Promote long slices to regular rabits to avoid trying to build
/// impossibly long things with oligos.
let cleanLongSlices (a:AssemblyOut) =
    {a with
        dnaParts =
            a.dnaParts |> List.map (fun s ->
                if (s.sliceType = INLINEST &&
                    s.dna.Length > 30 &&
                    not (s.pragmas.ContainsKey("inline")))
                then
                    {s with
                        sliceType = REGULAR;
                        dnaSource =
                            match a.pragmas.TryGetOne("refgenome") with
                            | None -> "synthetic"
                            | Some(x) -> x;}
                else s)
    }

/// Once GSL is expanded as far as possible, go into output generation,
/// including reagent and auxillary file creation.
let writeOutput
        (opts:ParsedOptions)
        (ga:GlobalAssets)
        (outTree:Assembly list) =

    let ryseLinkers = opj opts.libDir "linkers.txt" |> loadRyseLinkers
    let hutchAncillary = opj opts.libDir "hutch.txt" |> loadThumperRef
    let markerSets = markers.loadMarkers opts.libDir

    let genomePrefixes =
        ga.rgs |> Seq.choose (fun kv->
            match kv.Value.Env.TryFind("prefix") with
            | None -> None
            | Some(prefix) -> Some(kv.Key,prefix))
        |> Map.ofSeq

    let linkerlessTree = expandTree ga.rgs opts ga.seqLibrary outTree opts.proxyURL |> List.ofSeq

    let linkerlessTreeClean = linkerlessTree |> List.map (cleanLongSlices)

    // Stitch
    // Link for different platforms
    let linkedTree =
        match opts.platform with
        | MegaStitch ->
            List.map (mapRyseLinkers opts hutchAncillary ryseLinkers) linkerlessTreeClean
        | NoPlatform -> linkerlessTreeClean

    /// Updated tree, with linkers now and with any external part references
    /// swapped in.
    let linkedReusedTree =
        match opts.proxyURL with
            | None -> linkedTree
            | Some(_) ->
                // take candidates for part reuse and actually insert their sequence and indentifiers into the 
                // slice tree.
                reuseThumperParts linkedTree
                
    // flag_new_gsl 8/12/15 Todo: At here, check if candidate can be resued.
    let primers, tweakedTree =
        if opts.noPrimers then
            None, linkedTree
        else
            let p, t = designPrimers opts linkedReusedTree
            Some(p), t

    match opts.primerFile with
    | None -> ()
    | Some(file) ->
        if opts.noPrimers then
            failwithf "ERROR: can't combine --noprimers and --primers options"
        else
            primerDump.simplePrimerDump file primers.Value tweakedTree

    // Create specific output formats according to command line options
    match opts.apeOut with
    | None -> ()
    | Some(path,tag) -> dumpAPE path tag tweakedTree

    match opts.cmOut with
    | None -> ()
    | Some(path,tag) -> dumpCM path tag tweakedTree primers

    match opts.jsonOut with
    | None -> ()
    | Some(prefix) -> dumpJsonAssemblies prefix tweakedTree

    match opts.thumperOut with
    | None -> ()
    | Some(proj) ->
        dumpThumper
            opts
            genomePrefixes
            markerSets
            ryseLinkers
            hutchAncillary
            proj
            tweakedTree
            primers.Value
            opts.sbolOut

    match opts.docStringFile with
    | None -> ()
    | Some(path) -> dumpDocStrings path tweakedTree

    // And stdout
    match opts.flatOut with
    | None -> ()
    | Some(outFile) -> dumpText outFile linkedTree // Dump flat file format

    if opts.verbose then printfn "ok"

type CompileResult =
    | FinishedAssemblies of Assembly list
    | ExpandedGSL of string

/// Run GSLC on string input.
/// Raises an exception on error.
let rec processGSL (opts:ParsedOptions) (plugins:pluginDefaults.Plugin list) (ga:GlobalAssets) gslText =
    let verbose = opts.verbose

    // If the tree contains any parts that require expansion before final
    // translation to DNA then we should do that rather than emit final output

    /// An exception will be raised here if parsing fails.
    let treeWithPragma = lexAndParse verbose gslText

    /// Build up all legal capabilities by going through plugins
    let legalCapas = plugins |> 
                     List.map (fun pi->pi.providesCapas) |>
                     List.concat |>
                     set

    // Pull out just the assembly lines
    let inTree = stuffPragmas verbose (legalCapas.Contains) initialDesignParams treeWithPragma
    let inDocTree = stuffDocstrings verbose inTree

    // Filter out everything besides assemblies
    let assemblyOnly =
        inDocTree |> List.choose (fun x ->
            match x with
            | ASSEMBLY(a) -> Some(a)
            | _ -> None)

    // Check that all submitted gene names resolve to known genes
    checkGeneNames ga.rgs ga.seqLibrary assemblyOnly

    let expNeeded = Expansion.highestExpansionPriority inDocTree

    if expNeeded <> EXP_NONE then
        // Rewrite the GSL including pragmas and split it back out again
        let newGSL =
            Expansion.gslExpansion
                (not opts.doParallel)
                verbose
                expNeeded
                plugins
                ga.rgs
                ga.codonTableCache
                inDocTree
        if verbose then
            printfn "old GSL:\n %s" gslText
            printfn "\nNewGSL:\n %s\n" newGSL
        if opts.iter then
            // recurse, go all the way
            processGSL opts plugins ga newGSL
        else
            // Take one step, write this back out and we are done
            ExpandedGSL(newGSL)
    else
        // No further expansion needed - PROCEED to output generation
        FinishedAssemblies(assemblyOnly)
