module Expansion

open shared // common routines used in several modules
open pragmaTypes
open commonTypes
open System
open constants
open thumper // Thumper output formats, dumping etc
open ryse    // RYSE architecture
open cloneManager
open ape
open l2expline
open parseTypes
open Amyris.Bio.utils
open PrettyPrint
open sgdrefformat
open DnaCreation
open alleleSwaps
open applySlices
open resolveExtPart
open pluginDefaults
open codoptSupport

/// Rewriting rules for serially expanding GSL input

// ========================================================================================================
// Expansion time - rewriting GSL serially to emit lower forms
// ========================================================================================================

/// What type of GSL expansion do we need to perform
type ExpansionLevel =
      EXP_NONE
    | EXP_MUT
    | EXP_HB
    | EXP_PROT
    | EXP_NAMING
    | EXP_INLINEROUGHAGE
    | EXP_FOR
    | EXP_FUNCTION
    | EXP_OPEN
    | EXP_L2EXPLINE
    | EXP_MULTIPART

let nameLegal =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789![]@$%^&*()'\":_-=+,.?/`~"
    |> Set.ofSeq

let whitespace c = match c with | ' ' -> true | '\t' -> true | _ -> false

let cleanHashName (s:string) =
    s |> Seq.choose (fun c ->
        if nameLegal.Contains(c) then Some(c)
        else if whitespace c then None
        else Some('_'))
    |> Array.ofSeq |> arr2seq

/// Ensure we generate a name for each part
let expandNames (tree: GSLLine list) =

    seq {
        for line in tree do
            match line with
            | ASSEMBLY(a) when a.name.IsNone -> // This one needs a name
                let literal = prettyPrintLine line |> cleanHashName
                let name = literal.Substring(0, min literal.Length maxNameLen).Replace("@","(@)")
                yield PRAGMA({definition = namePragmaDef; args = [name]})
                yield line
            | x -> yield x // Non assembly line
        } |> List.ofSeq

/// Given a choice, how high priority is a particular type of expansion
/// - e.g. do proteins before heterology blocks
let expansionPriority (e:ExpansionLevel) =
    match e with
    | EXP_NONE -> 0
    | EXP_HB -> 5
    | EXP_FUNCTION -> 8
    | EXP_PROT -> 10
    | EXP_MUT -> 20
    | EXP_FOR -> 200
    | EXP_OPEN -> 300
    | EXP_NAMING -> 500
    | EXP_L2EXPLINE -> 750
    | EXP_INLINEROUGHAGE -> 1000
    | EXP_MULTIPART -> 10000 // Highest priority to deal with b/c no syntax to pretty print this

let modIsMutation m =
    match m with
    | MUTATION(_) -> true
    | _ -> false

let modIsNotSlice m =
    match m with
    | SLICE(_) -> false
    | _ -> true

let chooseRefGenome (p:PragmaCollection) =
    match p.TryGetOne("refgenome") with
    | Some(rg) -> rg
    | None ->
        printf "%s" (refGenomeWarning())
        defaultRefGenome  // Warning - defualts to yeast codon usage

/// Check if an individual part
/// requires expansion and work out what type
let testPart (p:Part) =
    match p with
    | MULTIPART _ -> EXP_MULTIPART
    | GENEPART(gp) ->
        if gp.part.mods |> List.exists modIsMutation then EXP_MUT
        else EXP_NONE
    | INLINEDNA(_) -> EXP_NONE // no need for expansion
    | INLINEPROT(_) -> EXP_PROT
    | MARKERPART -> EXP_NONE // no need for expansion
    | HETBLOCK -> EXP_HB
    | EXPANDED(_) -> EXP_NONE
    | PARTID(x) ->
        // Can introduce a mutation into an existing part, so we check to see if
        // that expansion will be needed
        if x.mods |> List.exists modIsMutation then EXP_MUT
        else EXP_NONE
    | ERRORPART(_) -> EXP_NONE // Can't expand an error

/// Check all parts in one assembly
let rec testOnePartsList (p:PPP list) =
    match p with
    | [] -> EXP_NONE
    | hd::tl ->
        // ignore pragma arm for sake of exploration
        match testPart hd.part with
        | EXP_NONE -> testOnePartsList tl
        | _ as x -> x

/// Explore parse tree to look for things needing expansion
let rec checkNeedsExpansion (lines:GSLLine list) =
    lines |> Seq.map (fun gslLine ->
        match gslLine with
        | ROUGHAGESECTION(_) -> [EXP_INLINEROUGHAGE]
        | L2(_) -> [EXP_L2EXPLINE]
        | FORBLOCK(f) -> EXP_FOR::(checkNeedsExpansion f.body)
        | ASSEMBLY(a) ->
            let otherExp = a.parts |> List.map (fun ppp -> testPart ppp.part)
            match a.name with
            // missing assembly names take priority
            | None -> EXP_NAMING::otherExp
            // Has a name, check the parts for expansion needs
            | _ -> otherExp
        | CUTLine(c) -> checkNeedsExpansion [c.assembly]
        | DOCSTRINGSECTION(_) -> []
        | GSLFUNCTION(_) -> []
        | GSLFUNCTIONCALL(_) -> [EXP_FUNCTION]
        | OPENLINE(_) -> [EXP_OPEN]
        | GSLLINEEXPANSION(_) -> []
        | PRAGMA(_) -> []
        | LETLINE(_) -> []
        )
    |> List.concat

/// Recursively examine every GSL line and return the expansion step with
/// highest priority.
let highestExpansionPriority (lines:GSLLine list) =
    checkNeedsExpansion lines
    |> Seq.fold
        (fun highestExp v ->
            if (expansionPriority v) > (expansionPriority highestExp) then v
            else highestExp)
        EXP_NONE

/// Take inline protein sequences and expand them out to DNA sequences
let expandProtein
        useSerial
        verbose
        (rgs:GenomeDefs)
        (codonTableCache:plugins.CodonTableCache)
        (tree:GSLLine list) =

    let rewrite (codonOptParams:plugins.CodonOptParams) (refGenome:string) (p:PPP) =
        match p.part with
        | INLINEPROT(s) ->
            // Check amino acid sequence is legal
            for c in s do
                if not (aaLegal.Contains(c)) then
                    failwithf
                        "ERROR: protein sequence contains illegal amino acid '%c'"
                        c

            let refGenome' =
                match p.pr.TryGetOne("refgenome") with
                | Some(rg) -> rg
                | None -> refGenome

            match rgs.TryFind refGenome' with
            | None ->
                failwithf
                    "ERROR: unable to load refgenome %s to determine environment"
                    refGenome'
            | Some(x) ->
                // Check to see if there is a local #seed parameter and extract it else
                // fall back on the version in the codon opt parameters globally
                let codonOptData = codonTableCache.Get(x,refGenome')
    
                let seed =
                        match p.pr.TryGetOne("seed") with
                        | None -> (codonOptParams :> ICodonOptParam).Seed()
                        | Some(seed) ->
                            match Int32.TryParse seed with
                            | true,s -> s
                            | _ ->
                                failwithf
                                    "ERROR: #seed argument '%s' is not a valid integer"
                                    seed
                let result = plugins.doCodonOpt verbose codonOptParams codonOptData seed s
                {p with part = INLINEDNA(result)}
        | _ -> p

    let procOneGSLLine (a:GSLLine) =
        match a with
        | ASSEMBLY(a) ->
            let refGenome = chooseRefGenome (a.pragmas)
            
            (*
            
            let seed =
                match a.pragmas.TryGetOne("seed") with
                | None -> defaultCodonOptSeed
                | Some(s) ->
                    match Int32.TryParse s with
                    | true,s' -> s'
                    | _ -> failwithf "ERROR: invalid integer '%s' for #seed" s

            let fivePrimeWindow,globalRepeatAvoidUser =
                match a.pragmas.TryGetOne("codonopt") with
                | None -> codons.defaultCodonOptParams.fivePrimeWindow,codons.defaultCodonOptParams.globalRepeatCheck
                | Some(s) ->
                    let parts =
                        s.Split([|';'|],StringSplitOptions.RemoveEmptyEntries)
                        |> Array.map (fun s->
                            match (s.Split( [| '=' |],StringSplitOptions.RemoveEmptyEntries ) ) with
                            | [| name ; value |] -> (name,value)
                            | _ as x ->
                                failwithf
                                    "ERROR: bad name=value pair '%A' in codonopt params" x)
                        |> Map.ofArray

                    let fiveW =
                        match parts.TryFind "5window" with
                        | None -> codons.defaultCodonOptParams.fivePrimeWindow
                        | Some(v) -> int v
                    let grc =
                        match parts.TryFind "repeatcheck" with
                        | None -> codons.defaultCodonOptParams.globalRepeatCheck
                        | Some(v) ->
                            match v.ToLower() with
                                | "true" -> true
                                | "false" -> false
                                | _ as x -> failwithf "ERROR: bad #codonopt repeatcheck param '%s' - should be true or false" x
                    fiveW,grc
            *)
            let cop = plugins.parseCodonOptParams a.pragmas
            ASSEMBLY({a with parts = List.map (rewrite cop refGenome) a.parts})
        | x -> x

    // Foreach assembly, update the parts if there is a mutation present
    tree
    |> Array.ofList
    |> (if useSerial then (Array.map) else (Array.Parallel.map)) (procOneGSLLine)
    |> List.ofArray

/// Remove mutation definitions and replace with a lower level representation
let expandMut 
    verbose 
    (providers : alleleSwapPluginDefs.AlleleSwapProvider list)
    (codonTableCache:plugins.CodonTableCache) 
    (rgs:GenomeDefs) 
    (tree:GSLLine list) =
    /// Rewrite an individual part p in an assembly a.
    /// Assembly is passed in for context if needed
    let rewrite (a:Assembly) (aName:string)(p:PPP) =
        match p.part with
        | HETBLOCK -> p // don't expect this but just in case
        | INLINEDNA(_) -> p
        | INLINEPROT(_) -> p
        | MARKERPART -> p
        | EXPANDED(_) -> p
        | ERRORPART(_) -> p
        | MULTIPART(_) -> p
        | PARTID(part) ->
            // Does it contain a mutation modification
            match part.mods |> List.tryFind modIsMutation with
            | None -> p
            | Some(MUTATION(mutMod)) ->
                let rg' = getRG a rgs p.pr
                let asAACheck =
                    match a.pragmas.TryFind("warnoff") with
                    | Some(p) -> not (p.hasVal "asaacheck")
                    | None -> true
                if verbose then printfn "***** %A" a.pragmas

                // Leave pragmas intact
                {p with
                    part = EXPANDED(alleleSwaps.expandSimpleMut asAACheck rg' part mutMod)}
            | _ -> failwith "ERROR: internal error, expected Mutation match"

        | GENEPART(gp) ->
            // Does it contain a mutation modification?
            match gp.part.mods |> List.tryFind modIsMutation with
            | None -> p
            | Some(MUTATION(mutMod)) ->
                if (not (gp.part.gene.[0] = 'G' || gp.part.gene.[0] = 'g')) then
                    failwithf
                        "ERROR: allele swap gene must be g-type  e.g gABC1$x1234y.  '%c' is not a legal prefix for %s"
                        (gp.part.gene.[0]) gp.part.gene
                let rg' = getRG a rgs p.pr

                // Need to select a codon usage table
                let refGenome = chooseRefGenome (a.pragmas)

                let codonUsage =
                    codonTableCache.Get(rg',refGenome)
                    |> fun cod -> Amyris.Bio.IO.CodonUsage.prepCUT 0.0 100 cod.freq
                
                let endPref =
                    match p.pr.TryGetOne("swapend") with
                    | Some("5") -> alleleSwapPluginDefs.NTERM
                    | Some("3") -> alleleSwapPluginDefs.CTERM
                    | Some(_) -> failwithf "ERROR:  #swapend argument should be 5 or 3"
                    | None -> alleleSwapPluginDefs.NONETERM

                if verbose then
                    printf "Mutation endpreference: %s\n"
                        (match endPref with
                         | alleleSwapPluginDefs.NTERM -> "Nterm"
                         | alleleSwapPluginDefs.CTERM -> "Cterm"
                         | alleleSwapPluginDefs.NONETERM -> "No end preference")

                if not (rg'.IsValid(gp.part.gene.[1..])) then
                    failwithf "ERROR: undefined gene '%s' %s\n"
                        (gp.part.gene.[1..]) (fp gp.part.where)

                let asAACheck =
                    match a.pragmas.TryFind("warnoff") with
                    | Some(p) -> not (p.hasVal "asaacheck")
                    | None -> true

                // Check if there is a style pragma on the part itself.  User can designate
                // short or long allele swap styles at the part level
                let longStyle =
                    match p.pr.TryGetOne("style") with
                    | None -> true // default is long style
                    | Some("long") -> true // long style
                    | Some("short") -> false // short style
                    | Some(x) ->
                        failwithf
                            "ERROR: undefined allele swap style '%s', options are long or short" x

                {p with
                    part = EXPANDED(alleleSwaps.expandAS providers asAACheck aName verbose rg' codonUsage gp.part.gene mutMod endPref a.capabilities longStyle) } // Leave pragmas intact
            | _ -> failwith "ERROR: internal error, expected Mutation match"

    // Foreach assembly, update the parts if there is a mutation present
    tree
    |> List.map (fun a ->
        match a with
        | ASSEMBLY(a) ->
            ASSEMBLY(
                {a with
                    parts = List.map
                                (rewrite a (match a.name with
                                            | Some(n) -> n
                                            | None -> a.parts.Head.part.ToString()))
                                a.parts})
        | x -> x)

/// Just expand heterology blocks
let expandHB verbose (codonTableCache:plugins.CodonTableCache) (rgs:GenomeDefs) (tree: GSLLine list) =
    let getLenPragma (pr:PragmaCollection) =
        match pr.TryGetOne("len") with
        | None -> None
        | Some(v) ->
            match Int32.TryParse(v) with
            | true,i -> Some(i)
            | _ -> failwithf "ERROR: expected integer in #len pragma, not '%s'" v


    let rec scan
            (a:Assembly)
            (res:(Part*PragmaCollection*bool) list)
            (p: (Part*PragmaCollection*bool) list)  =
        // Need to select a codon usage table
        let refGenome = chooseRefGenome (a.pragmas)
        let rg = rgs.[refGenome]
        let codonUsage =
            codonTableCache.Get(rg,refGenome)
            |> fun cod -> Amyris.Bio.IO.CodonUsage.prepCUT 0.0 100 cod.freq

        match p with
        // Lone heterology block with no inlined sequence adjacent
        | (GENEPART(gpUp),pr1,fwd1)
          ::(HETBLOCK,pr2,fwd2)
          ::(GENEPART(gp),pr3,fwd3)
          ::tl ->
            let rg' = getRG a rgs pr3
            let rg'' = getRG a rgs pr1

            let sliceSeq = realizeSequence verbose fwd3 rg' gp // Get DNA sequence for this slice
            let sliceSeqUp = realizeSequence verbose fwd1 rg'' gpUp // Get DNA sequence for upstream slice
            let targetAALen = getLenPragma pr2 // Get optional length spec for the HB

            // Generate an alternative prefix for the GENEPART on RHS
            let alt =
                generateRightHB
                    codonUsage
                    minHBCodonUsage
                    targetAALen
                    a.designParams
                    sliceSeqUp
                    [||]
                    sliceSeq
            // tricky part - need to slightly adjust the slice range of gp,
            // but that's embedded down in the mod list

            // Assume they must be using a gene part next to a het block.  Bad?
            if (not (gp.part.gene.StartsWith("g"))) then
                failwithf "ERROR: heterology block must be adjacent to g part, %s not allowed" gp.part.gene
            let s = translateGenePrefix rg' GENE // Start with standard slice
            let startSlice = applySlices verbose gp.part.mods s // Apply modifiers
            let newSlice =
                {startSlice with
                    left =
                        {startSlice.left with
                            x = startSlice.left.x + (alt.Length*1<OneOffset>)}} // Chop it
            let newInline = alt|> arr2seq
            assert(alt <> sliceSeq.[0..alt.Length-1])
            // Assemble new mod list by getting rid of existing slice mods and
            // putting in new consolidated slice.
            let newMods =
                SLICE(newSlice)
                ::(gp.part.mods |> List.filter modIsNotSlice)

            // Prepend backwards as we will flip list at end - watch out, pragmas are reversed as well
            scan
                a
                ((GENEPART({gp with part = {gp.part with mods = newMods}}),pr3,fwd3)
                 ::(INLINEDNA(newInline), pr2.Add("inline"), fwd2)
                 ::(GENEPART(gpUp),pr1,fwd1)
                 ::res)
                tl

        | (GENEPART(gpUp),pr1,fwd1)
          ::(INLINEDNA(i),pr2,fwd2)
          ::(HETBLOCK,pr3,_(*fwd3*))
          ::(GENEPART(gp),pr4,fwd4)
          ::tl ->
            let ic = i.ToCharArray()
            let rg' = getRG a rgs pr4
            let rg'' = getRG a rgs pr1

            // Get DNA sequence for this slice
            let sliceSeq = realizeSequence verbose fwd4 rg' gp
            // Get DNA sequence for upstream slice
            let sliceSeqUp = realizeSequence verbose fwd1 rg'' gpUp
            // Get optional length spec for the HB
            let targetAALen = getLenPragma pr3
            // Generate an alternative prefix for the GENEPART on RHS
            let alt =
                generateRightHB
                    codonUsage minHBCodonUsage targetAALen a.designParams sliceSeqUp ic sliceSeq

            assert(alt <> sliceSeq.[0..alt.Length-1])
            if verbose then
                printf "// %s -> %s\n"
                    (arr2seq alt) (sliceSeq.[0..alt.Length-1] |> arr2seq)

            // tricky part - need to slightly adjust the slice range of gp,
            // but that's embedded down in the mod list

            // Assume they must be using a gene part next to a het block.  Bad?
            if not (gp.part.gene.[0] = 'g') then
                failwithf
                    "ERROR: slices adjacent to het block elements ~ must be gene slices - %s has '%c' gene part"
                    gp.part.gene gp.part.gene.[0]

            let s = translateGenePrefix rg'' GENE // Start with standard slice
            let startSlice = applySlices verbose gp.part.mods s // Apply modifiers
            let newSlice =
                {startSlice with
                    left =
                        {startSlice.left with
                            x = startSlice.left.x + (alt.Length*1<OneOffset>)}} // Chop it

            let newInline = Array.append ic alt |> arr2seq

            // Assemble new mod list by getting rid of existing slice mods and putting in new consolidated slice.
            let newMods =
                SLICE(newSlice)
                ::(gp.part.mods |> List.filter modIsNotSlice)
            // Prepend backwards as we will flip list at end
            // Note - currently destroy any pragmas attached to the heterology block itself
            scan
                a
                ((GENEPART({ gp with part = {gp.part with mods = newMods}} ),pr4,fwd4)
                 ::(INLINEDNA(newInline),pr2.Add("inline"),fwd2)
                 ::(GENEPART(gpUp),pr1,fwd1)
                 ::res)
                tl

        | (GENEPART(gp),pr1,fwd1)
          ::(HETBLOCK,pr2,_(*fwd2*))
          ::(INLINEDNA(i),pr3,fwd3)
          ::(GENEPART(gpDown),pr4,fwd4)
          ::tl ->
            let ic = i.ToCharArray()

            // get reference genomes warmed up
            let rg' = getRG a rgs pr1
            let rg'' = getRG a rgs pr4

            // get actual sequence for slices (with mods applied)
            let sliceSeq = realizeSequence verbose fwd1 rg' gp // Get DNA sequence for this slice
            let sliceSeqDown = realizeSequence verbose fwd4 rg'' gpDown // Get DNA sequence for this slice
            let targetAALen = getLenPragma pr2 // Get optional length spec for the HB

            // generate hetblock section off upstream slice
            // Generate an alternative prefix for the GENEPART on LHS
            let alt =
                generateLeftHB
                    codonUsage minHBCodonUsage targetAALen a.designParams sliceSeq ic sliceSeqDown
            // tricky part - need to slightly adjust the slice range of gp,
            // but that's embedded down in the mod list

            // Assume they must be using a gene part next to a het block.  Bad?
            assert(gp.part.gene.[0] = 'g')

            // now build up the slice again and apply from scratch to the gene
            let s = translateGenePrefix rg' GENE // Start with standard slice
            let startSlice = applySlices verbose gp.part.mods s // Apply modifiers

            // modify slice to take into account the bit we chopped off
            let newSlice =
                {startSlice with
                    right =
                        {startSlice.right with
                            x = startSlice.right.x - (alt.Length*1<OneOffset>)}} // Chop it
            let newInline = Array.append alt ic |> arr2seq
            assert(alt <> sliceSeq.[sliceSeq.Length-alt.Length..])
            if verbose then
                let fr = (arr2seq alt)
                let t = sliceSeq.[sliceSeq.Length-alt.Length..sliceSeq.Length-1]
                let sim = Array.map2 (fun a b -> if a = b then '.' else ' ') alt t
                let simProp =
                    (seq { for x in sim -> if x = '.' then 1.0 else 0.0} |> Seq.sum)
                  / float(alt.Length) * 100.0

                let capUsing (a:string) (b:string) =
                    Seq.zip a b
                    |> Seq.map (fun (a,b) ->
                        if a=b then b
                        else
                            match b with
                            | 'G' -> 'g'
                            | 'C' -> 'c'
                            | 'A' -> 'a'
                            | 'T' -> 't'
                            | _ as x ->x)
                    |> Array.ofSeq |> arr2seq
                printf
                    "// From: %s \n// To  : %s\n//     : %s %3.0f%%\n"
                    fr (arr2seq t |> capUsing fr) (arr2seq sim) simProp
            // Assemble new mod list by getting rid of existing slice mods and
            // putting in new consolidated slice.
            let newMods =
                SLICE(newSlice)
                ::(gp.part.mods |> List.filter modIsNotSlice)
            // Prepend backwards as we will flip list at end
            // Note - currently destroy pr2 pragmas associated with the hetblock
            scan
                a
                ((GENEPART(gpDown),pr4,fwd4)
                 ::(INLINEDNA(newInline),pr3.Add("inline"),fwd3)
                 ::(GENEPART({gp with part = {gp.part with mods = newMods}}),pr1,fwd1)
                 ::res)
                tl

        | (PARTID(pid1),pr1,fwd1)
          ::(HETBLOCK,pr2,_(*fwd2*))
          ::(INLINEDNA(i),pr3,fwd3)
          ::(PARTID(pid4),pr4,fwd4)
          ::tl ->
            // External part variation
            // ===============================================
            let ic = i.ToCharArray()
            //let rg' = getRG

            match fetchFullPartSequence(verbose) (Map.empty) pid1 with
            | EXT_FAIL(msg) -> failwithf "ERROR: fail fetching %s %s" pid1.id msg
            | EXT_FETCH_OK(part1) ->
                match fetchFullPartSequence verbose Map.empty pid4 with
                | EXT_FAIL(msg) -> failwithf "ERROR: fail fetching %s %s" pid1.id msg
                | EXT_FETCH_OK(part4) ->
                    let s1= getExtPartSlice verbose pid1
                    let s4= getExtPartSlice verbose pid4

                    // Build up upstream and downstream DNA slice
                    let sliceSeq1 = applySliceToExtSequence verbose part1 pr1 fwd1 pid1 s1
                    let sliceSeq4 = applySliceToExtSequence verbose part4 pr4 fwd4 pid4 s4

                    let targetAALen = getLenPragma pr2 // Get optional length spec for the HB

                    // generate hetblock sequence by cutting into upstream sequence
                    // Generate an alternative prefix for the GENEPART on LHS
                    let alt =
                        generateLeftHB
                            codonUsage minHBCodonUsage targetAALen
                            a.designParams sliceSeq1.dna ic sliceSeq4.dna

                    // Build up the slice mods for the upstream part from scratch
                    // Modify upstream slice to account for the part we chopped off
                    let newSlice:Slice =
                        {s1 with
                            right =
                                {s1.right with x = s1.right.x - (alt.Length*1<OneOffset>)}} // Chop it
                    let newInline = Array.append alt ic |> arr2seq
                    assert(alt <> sliceSeq1.dna.[sliceSeq1.dna.Length-alt.Length..])
                    if verbose then
                        let fr = (arr2seq alt)
                        let t = sliceSeq1.dna.[sliceSeq1.dna.Length-alt.Length..]
                        let sim = Array.map2 (fun a b -> if a = b then '.' else ' ') alt t
                        let simProp =
                            (seq { for x in sim -> if x = '.' then 1.0 else 0.0} |> Seq.sum )
                          / float(alt.Length) * 100.0

                        let capUsing (a:string) (b:string) =
                            Seq.zip a b
                            |> Seq.map (fun (a,b) ->
                                if a=b then b
                                else
                                    match b with
                                    | 'G' -> 'g'
                                    | 'C' -> 'c'
                                    | 'A' -> 'a'
                                    | 'T' -> 't'
                                    | _ as x ->x)
                            |> Array.ofSeq |> arr2seq
                        printf
                            "// From: %s \n// To  : %s\n//     : %s %3.0f%%\n"
                            fr (arr2seq t |> capUsing fr)  (arr2seq sim) simProp
                    // Assemble new mod list by getting rid of existing slice mods
                    // and putting in new consolidated slice.
                    let newMods =
                        SLICE(newSlice)
                        ::(pid1.mods |> List.filter modIsNotSlice)
                    // Prepend backwards as we will flip list at end
                    // Note - currently destroy pr2 pragmas associated with the hetblock
                    scan
                        a
                        ((PARTID(pid4),pr4,fwd4)
                         ::(INLINEDNA(newInline),pr3.Add("inline"),fwd3)
                         ::(PARTID({ pid1 with mods = newMods}),pr1,fwd1)
                         ::res)
                        tl
        | hd::tl -> scan a (hd::res) tl // do nothing
        | [] -> List.rev res
    /// Easier to process part/pragma pairs if they are explicit tuples
    let tuplePPP (ppp: PPP list) =
        ppp |> List.map (fun p -> p.part,p.pr,p.fwd)
    let untuplePPP ab =
        ab |> List.map (fun (a,b,c) -> { part = a ; pr = b ; fwd= c} )

    let res =
        tree |> List.map (fun a ->
            match a with
            | ASSEMBLY(a) ->
                ASSEMBLY({a with parts = scan a [] (tuplePPP a.parts) |> untuplePPP})
            | x -> x)

    if res |> List.exists(fun gslLine ->
        match gslLine with
        | ASSEMBLY(a) ->
            a.parts |> List.exists(fun ppp ->
                match ppp.part with | HETBLOCK -> true | _ -> false)
        | _ -> false)
    then
        failwithf
            "ERROR: attempt to expand heterology block in design %A left remaining hetblock"
            (prettyPrintTree res)
    res

/// Apply any simple sanity checks to roughage lines before we start rewriting them
let validateRoughageLine (r:Roughage) =
    // Rule 1:  must be able to work out the locus.  Locus can be either explicit (ho^) or
    //          implicit pSLN1>YNG1  but can't have just bidirectional promoters with no explicit locus  e.g.   ADH1<pGAL1-pGAL10>ADH2
    let hasLocus = r.locus <> "" || (r.parts.Length>0 && not r.parts.Head.bi)
    if not hasLocus then
        let line = match r.parts with | [] -> -1 | hd::_ -> (hd.left.Line+1)
        failwithf
            "ERROR: inline roughage construct on line %d has indeterminate locus: %s\n"
            line (prettyPrintRoughage r)
    ()

/// Take classic roughage as inline, replace it with equivalent GSL level 2 syntax
let expandInlineRoughage (_:GenomeDefs) (tree:GSLLine list) =
    /// This is the guts of the translation.  This takes one line of
    /// roughage and produces one line of GSL level 2

    /// FIXME Hard coded mapping of markers for now
    let markerMapping (s:string) =
        match s with
        | "mURA" -> "ura3"
        | "mKANA" -> "kan"
        | "mLEU2" -> "leu2"
        | "mTRP1" -> "trp1"
        | "mURA3" -> "ura3"
        | "mURA3LO" -> "ura3lo"
        | _ as x -> x // TODO: more generalized support not hard coded

    /// Emit Level 2 GSL for a list of inline roughage lines
    let rec translateRoughageLines (rLines:Roughage list) currMarker results =
        // Build the output as a string
        match rLines with
        // don't reverse because prepending these results into overall output list
        | [] -> currMarker, results
        | rLine::tl ->
            validateRoughageLine rLine
            let marker =
                match rLine.marker with
                | "" -> // No marker attached to the locus knockout
                    match rLine.parts |> List.tryFind (fun re -> re.marker.IsSome) with
                    | None -> None // No marker attached to a part either
                    | Some(x) -> x.marker
                | x -> Some(x) // Yes there is a marker attached to the locus knockout
            // What is the marker going forward?
            // For roughage, if no marker is specified, it defaults to ura3
            let newMarker = match marker with | None -> "ura3" | Some(x) -> markerMapping x

            // Expansion of one line inside the roughage inline block
            let newGSLElement =
                GSLLINEEXPANSION(
                    String.Join("",
                        seq {
                            if newMarker<> currMarker then
                                yield sprintf "#markerset %s\n" newMarker
                            if rLine.locus <> "" then
                                yield sprintf "%s^;" rLine.locus
                            for re in rLine.parts do
                                yield
                                    if re.bi then
                                        sprintf "%s>%s ; %s>%s ;"
                                            re.promoter1 re.target1 re.promoter2 re.target2
                                    else
                                        sprintf "%s>%s; "
                                            re.promoter1 re.target1
                        } ).Trim([|';' ; ' '|]) // remove trailing semicolon
                ) // End EXPANDED GSLLine definition
            // Recurse for remaining lines
            translateRoughageLines tl newMarker (newGSLElement::results)

    /// Recursively process GSL lines and unravel any Roughage sections
    /// Need to track the current marker and reset it later.
    let rec procGSLElements (lines :  GSLLine list) currMarker (results: GSLLine list) =
        match lines with
        | [] -> results |> List.rev // Reverse final list because we prepended lines as we processed them
        | ROUGHAGESECTION(rs)::tl ->
            /// Roughage section can result in more than one addition line
            /// Exand lines and keep them in the right order, remembering we will flip at the end
            let outMarker,outLines = translateRoughageLines rs currMarker []
            // let outLines = rs |> List.map (translateOneRoughageLine) |> List.rev
            if outMarker <> currMarker then
                // Need to swap the markerset back
                procGSLElements
                    tl
                    currMarker
                    (GSLLINEEXPANSION(sprintf "#markerset %s" currMarker)::outLines@results)
            else
                // Markers ended favorably, no need to change
                procGSLElements tl currMarker (outLines@results)
        | PRAGMA(p)::tl ->
            // Need to monitor pragmas, since someone might change the default
            // marker and we need to know which marker is active at any time
            if p.name = "markerset" then
                procGSLElements tl p.args.[0] (PRAGMA(p)::results)
            else
                procGSLElements tl currMarker (PRAGMA(p)::results)
        | hd::tl -> // Don't care about non roughage sections
            procGSLElements tl currMarker (hd::results)

    let startMarker = "ura3"
    procGSLElements tree startMarker []

///// Perform variable replacement in this block
//let varReplace (name:string) (replacement:string) (line:GSLLine) =
//    line // FIXFIX - unimplemented

/// Expand FOR statements, unrolling loops
let expandFor (tree:GSLLine list) =
    failwithf "ERROR: unimp"

    (*
    let rec procLines (lines:GSLLine list) (results:GSLLine list) =
        match lines with
            | [] -> List.rev results
            | FORBLOCK(f)::tl ->
                // Iterate over items in for loop,
                // replace block with substituted block as GSL EXPANSION lines
                match f.items with // TODO - match different types of list we could be expanding on
                    | PARTLIST(parts) -> ()
                    | INTLIST(il) ->
                    | INTRANGE(il,ir) ->
                        f.body |> List.map (varReplace f.name v) ] |> List.concat

            | hd::tl -> procLines tl (hd::results)

    procLines tree []
    *)

/// Expanded GSLLines that are functions or function calls
let expandFunction (tree:GSLLine list) =
    // Recursively go through GSLLine list until you bump into either a GSLFUNCTION or GSLFUNCTIONCALL
    let rec procLines
            (tree:GSLLine list)
            (funMap:Map<string, GSLFunction>)
            (newLines:GSLLine list) =
        match tree with
        | [] -> newLines
        | hd::tl ->
            match hd with
            | GSLFUNCTION(f) ->
                // put new function in funMap
                // put this same line back in the parse tree (newLines+hd)
                // recurse on the remainder of the parse tree (tl)
                let updatedNewLines = (List.append newLines [hd])
                // ^^ Decided to not put GSLFUNCTION back into the parse tree (causes infinite expansion loop)
                let updatedMap =
                    match funMap.TryFind(f.name) with
                    | Some(_) -> failwithf "ERROR: function %s already exists" f.name
                    | None -> funMap.Add(f.name, f)
                procLines tl updatedMap updatedNewLines

            | GSLFUNCTIONCALL(fc) ->
                // look for function being called in the funMap
                // add function body to new lines
                // recurse on remainder of parse tree
                let f =
                    match funMap.TryFind(fc.name) with
                    | Some(f) -> f
                    | None -> failwithf "ERROR: function %s not yet defined" fc.name

                // Check that the correct number of arguments have been provided
                if f.args.Length <> fc.args.Length then
                    failwithf "ERROR: Function %s requries %i arguments. Given %i."
                        f.name f.args.Length fc.args.Length

                let localVars =
                    seq {
                        for i in 0.. (f.args.Length - 1) do
                            yield sprintf "let %s = %s"
                                      (f.args.[i]) (printPPP (fc.args.[i]) )}
                    |> List.ofSeq
                    |> List.map (fun line -> GSLLINEEXPANSION(line))

                // put the new local variables together with the function body
                let newBody = List.append localVars f.body
                // Add all the lines back to the parse tree
                let updatedNewLines = List.append newLines newBody
                //printfn "New lines so far:\n %A" updatedNewLines
                procLines tl funMap updatedNewLines

            | _ -> procLines tl funMap (List.append newLines [hd])

    procLines tree Map.empty []

/// Expand all instances where multiple parts are
/// encoded in a single part, introduced via variables
/// Priority to get rid of these because we don't have a syntax to
/// pretty print them.  We remove them by expanding out the underlying list
/// of parts into the enclosing assembly, while observing orientation
/// and applying any multipart level pragmas to individual parts
let expandMultiparts (tree:GSLLine list) =
    let expandMultipart (ppp:PPP) =
        match ppp.part with
        | MULTIPART(m) ->
            let correctlyOriented =
                if ppp.fwd then
                    m
                else
                    // Before inverting part order, move any #fuse pragmas 1 part to the right. 
                    // This should attach it to the other side of the fused part so it is still 
                    // between them upon inverting the part order.
                    let rec shiftFuse (pppList:PPP list) (shiftedList:PPP list) (addFuse2Map:bool) =
                        match pppList with
                        | [] -> 
                            match addFuse2Map with
                            | false -> shiftedList // Done shifting fuses, just return the list
                            | true -> // 
                                failwithf "ERROR: there was a #fuse pragma attached to the right-most part of a list that got inverted... this is tricky to deal with and is not yet implemented :("
                            
                        | hd::tl ->
                            let curPragmaMap = hd.pr 
                            // if the current part's Pragma map contains a fuse, indicate that the next part should gain a fuse
                            let addFuse2NextPart = 
                                if curPragmaMap.ContainsKey("fuse") then
                                    true
                                else 
                                    false

                            // clear any fuses in this part
                            let removedFuseMap = curPragmaMap.Remove("fuse")

                            // But check if the PREVIOUS part indiciated that this part should also RECEIVE a fuse
                            let finalPragmaMap = 
                                match addFuse2Map with
                                | true -> removedFuseMap.Add("fuse")
                                | false -> removedFuseMap

                            // Now that we have the updated pragma Map for this part, make a part with this update map and add it to shiftedList
                            let partWithNewPragmas = {hd with pr = finalPragmaMap}
                            let newShiftedList = List.append shiftedList [partWithNewPragmas]

                            // call the function again with tail, the shiftedList we're accumulating, and the truth value of adding fuse to the next part

                            shiftFuse tl newShiftedList addFuse2NextPart

                    // call shift fuse on PPP list
                    let shiftedFuseList = shiftFuse m [] false  

                    // Invert part order, and orientation of the individual parts
                    // i.e pTDH3 ; mERG10 ; tCYC1 becomes !tCYC1 ; !mERG10 ; !pTDH3
                    let inverted =
                        List.rev shiftedFuseList
                        // invert the list
                        // flip part direction and inversion sensitive pragmas
                        |> List.map (fun thisppp ->
                            // get the pragmas that need to flip
                            let prToFlip =
                                thisppp.pr.pmap
                                |> Map.toSeq
                                |> Seq.choose (fun (_, pragma) ->
                                    match pragmaInverts pragma.definition with
                                    | Some(invertsTo) -> Some(pragma, invertsTo.name)
                                    | None -> None)
                                |> List.ofSeq
                            // remove pragmas we need to flip
                            let cleanedPrags =
                                prToFlip
                                |> List.fold
                                    (fun (pr:PragmaCollection) (p, _) -> pr.Remove(p.name))
                                    thisppp.pr

                            // add back the flipped versions
                            let flippedPragmas =
                                prToFlip
                                |> List.fold
                                    (fun (pr:PragmaCollection) (p, invertTo) -> pr.Add(buildPragma invertTo p.args))
                                    cleanedPrags
                            

                            // return the part with inverted direction and flipped pragmas
                            {thisppp with fwd = not thisppp.fwd; pr = flippedPragmas})

                    inverted
            
            // Finally distribute any multipart level pragmas over each individual part
            // Can imagine doing something more hierarchical here with naming etc to name a compound part
            correctlyOriented
            |> List.map (fun thisppp -> {thisppp with pr = thisppp.pr.MergeIn(ppp.pr)})
        | _ -> [ppp]

    let expandOneLine l =
        match l with
        | ASSEMBLY(a) ->
            ASSEMBLY({a with
                        parts =
                            a.parts |>
                            List.map expandMultipart |>
                            List.concat } )
        | _ as x -> x // Other lines don't have multiparts embedded in them

    tree |> List.map (expandOneLine)


/// Take input GSL and write it back out again, potentially rewriting parts of it
let gslExpansion
        useSerial
        verbose
        (level:ExpansionLevel)
        (plugins:Plugin list)
        (rg:GenomeDefs)
        (codonTableCache:plugins.CodonTableCache)
        (tree: GSLLine list) =
    prettyPrintTree
       (match level with
        | EXP_NONE -> failwith "ERROR: reached gslExpansion in non expansion mode :("
        | EXP_MULTIPART -> expandMultiparts tree
        | EXP_NAMING -> expandNames tree
        // deal with inline protein sequences
        | EXP_PROT -> expandProtein useSerial verbose rg codonTableCache tree
        // deal with heterology blocks
        | EXP_HB -> expandHB verbose codonTableCache rg tree
        // remove mutation definitions and rewrite tree
        | EXP_MUT -> 
            let alleleSwapAlgs = plugins |> List.choose (fun pi -> pi.alleleSwapAA)            
            expandMut verbose alleleSwapAlgs codonTableCache rg tree
        // Rewrite inline classic roughage as GSL Level 2
        | EXP_INLINEROUGHAGE -> expandInlineRoughage rg tree
        // Rewrite level 2 expression lines as level 1 GSL
        | EXP_L2EXPLINE -> 
            let l2Providers = plugins |> List.choose (fun provider -> provider.l2KOTitration)
            expandL2ExpLine l2Providers rg tree
        | EXP_FOR -> expandFor tree
        | EXP_FUNCTION -> expandFunction tree
        | EXP_OPEN -> failwithf "ERROR: EXP_OPEN unimplemented")