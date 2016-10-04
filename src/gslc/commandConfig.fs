/// Command line arguments, parsing, and command defaults.
module commandConfig
open System
open commonTypes
open pragmaTypes
open Amyris.Bio.utils
open semversion

let version = sprintf "%d.%d.%d" versionMajor versionMinor versionPatch

let libRoot = 
    match Environment.GetEnvironmentVariable("GSL_LIB") with
    | null -> "lib"
    | x -> x
    |> smashSlash


type CmdLineArgDef = {
    arg: string;
    param: string list;
    alias: string list;
    desc: string;
    // function which accepts the parameter arguments and incoming ParsedOptions
    // and returns a modified ParsedOptions
    proc: string list -> ParsedOptions -> ParsedOptions}


/// Define all GSLC command line arguments here.
/// An argument consists of its name, the names of its parameters, a description,
/// and a function that takes a list of passed parameters and an options record
/// and returns a modified options record.
let cmdLineArgs = [
    {arg = "ape"; param = ["outDir"; "prefix"]; alias = [];
     desc = "write APE output to prefix_##.ape to outDir\n(http://biologylabs.utah.edu/jorgensen/wayned/ape/)";
     proc = fun p opts -> {opts with apeOut = Some(p.[0], p.[1])} };

    {arg = "cm"; param = ["outDir"; "prefix"]; alias = [];
     desc="write clone manager output to prefix_##.cx5 to output directory outDir";
     proc = fun p opts -> {opts with cmOut = Some(p.[0], p.[1])} };

    {arg = "thumper"; param = ["proj"]; alias = [];
     desc = "write proj.rabit.csv proj.stitch.csv and proj.megastitch.csv files\nfor uploading to thumper - implies platform megastitch";
     proc = fun p opts -> {opts with thumperOut = Some(p.[0]) ; platform = MegaStitch}};

    {arg = "sbol"; param = []; alias = [];
     desc = "if writing thumper output, use SBOL/GBoM rather than RYCO-D";
     proc = fun _ opts -> {opts with sbolOut = true} };

    {arg = "thumper_proxy"; param = ["url"]; alias = [];
     desc = "perform part reuse from thumper given URL for thumper proxy";
     proc = fun p opts -> {opts with proxyURL = Some (p.[0])} }; 

    {arg = "rabitonly"; param = []; alias = [];
     desc = "generate just rabit output and omit stitches and megastitches";
     proc = fun _ opts -> {opts with rabitsOnly = true} };

    {arg = "reflist" ; param = [] ; alias = [];
     desc = "list available reference genomes";
     proc = fun _ opts -> {opts with refList = true}} ;

    {arg = "refdump" ; param = ["refname"] ; alias = [];
     desc = "dump available loci in reference genome";
     proc = fun p opts -> {opts with refDump = Some (p.[0])}} ;

    {arg = "flat"; param = ["outfile"]; alias = [];
     desc = "write a flat file format for results to outputfile";
     proc = fun p opts -> {opts with flatOut = Some(p.[0])} };

    {arg = "step"; param = []; alias = [];
     desc = "expand GSL just one round, and emit intermediate GSL";
     proc = fun _ opts -> {opts with iter = false}};

    {arg = "verbose"; param = []; alias = [];
     desc = "print debugging info";
     proc = fun _ opts -> {opts with verbose = true} };
    
    {arg = "version"; param = []; alias = [];
     desc = "print verison information";
     proc = fun _ opts ->
        printfn "GSL compiler version %s" version
        opts};

    {arg = "helpPragmas"; param = []; alias = [];
     desc = "print available pragmas";
     proc = fun _ opts ->
        pragmaUsage()
        opts};

    {arg = "quiet"; param = []; alias = [];
     desc = "suppress any non-essential output";
     proc = fun _ opts -> {opts with quiet = true} };

    {arg = "primers"; param = ["primerfile"]; alias = [];
     desc = "emit raw primer details (see also thumper output format)";
     proc = fun p opts -> {opts with primerFile = Some(p.[0])} };

    {arg = "noprimers"; param = []; alias = [];
     desc = "do not attempt to generate primers";
     proc = fun _ opts -> {opts with noPrimers = true} };

    {arg = "lib"; param = ["directory"]; alias = [];
     desc = "directory in which genome definitions reside\nDefault: GSL_LIB var, or 'lib' in current directory";
     proc = fun p opts -> {opts with libDir = smashSlash p.[0]}};

    {arg = "platform"; param = ["platform"]; alias = [];
     desc = "megastitch | none";
     proc = fun p opts ->
        match p.[0] with
        | "megastitch" -> {opts with platform = MegaStitch}
        | "none" -> {opts with platform = NoPlatform}
        | _ -> 
            printf "ERROR: unknown platform '%s', expected none or megastitch" p.[0]
            exit 1};

    {arg = "serial"; param = []; alias = [];
     desc = "don't run parallel operations, useful for debugging";
     proc = fun _ opts -> {opts with doParallel = false} };

    {arg = "lextest"; param = []; alias = ["tokentest"; "tokenize"];
     desc = "for debugging only, show stream of parsed tokens from input file";
     proc = fun _ opts -> {opts with lexOnly = true} };

    {arg = "docstring"; param = ["outfile"]; alias = ["docstrings"];
     desc = "log emitted documentation for each design to outfile";
     proc = fun p opts -> {opts with docStringFile = Some(p.[0])} };

    {arg = "name2id"; param = ["outfile"]; alias = [];
     desc = "filename/path for name2id mapping in rycody output.\nDefault is projname.name2id.txt";
     proc = fun p opts -> {opts with name2IdPath = Some(p.[0])} };
]

// TODO: add command aliases
let cmdLineArgIndex =
    seq {
        for a in cmdLineArgs do
            yield (a.arg, a)
            for alias in a.alias -> (alias, a)}
    |> Map.ofSeq

/// Format a command line argument and print it.
let printCmdLineArg a =
    let padSize = 35;
    let p = Seq.map (sprintf " <%s>") a.param |> Seq.fold (+) ""
    let larg = sprintf "       --%s%s" a.arg p
    let rPad = String.replicate (padSize - larg.Length) " "
    let descLines = a.desc.Split [|'\n'|]
    let firstLine = larg + rPad + "-" + descLines.[0]
    let otherLines =
        Array.map (fun (l:string) -> (String.replicate (padSize+1) " ") + l) (descLines.[1..])

    printfn "%s" firstLine
    for l in otherLines do printfn "%s" l

/// Print arg usage and help text.
let printUsage() =
    printfn "Usage:  gscl [args] input.gsl"
    for a in cmdLineArgs do printCmdLineArg a

///
let defaultOpts:ParsedOptions =
   {apeOut = None;
    quiet = false;
    libDir = libRoot;
    refStrain = "cenpk";
    cmOut = None;
    platform = NoPlatform;
    flatOut = None;
    thumperOut = None;
    sbolOut = false;
    iter = true;
    doParallel = true;
    verbose = false;
    rabitsOnly = false;
    primerFile = None;
    noPrimers = false;
    docStringFile = None;
    proxyURL = None;
    name2IdPath = None;
    lexOnly = false
    refList = false
    refDump = None
    }

/// Parse a command line arguments.  Return the parsed options and the list of
/// input files.  Raises an exception for ill-formed or unrecognized args.
let parseCommandLineArgs (argList:string list) =

    /// Parse one item and return the rest of the arg list, the modified options
    /// record, and any file names we've accumulated.
    let parseOneItem (h:string) (argList:string list) opts files =
        // Should only be a command or an input file.
        if h.StartsWith("--") then
            let arg = h.[2..]
            if arg = "help" then
                printUsage()
                (argList, opts, files)
            else
                match cmdLineArgIndex.TryFind arg with
                | None ->
                    failwithf "Unrecognied command line argument: %s" arg
                | Some(a) ->
                    // valid command line arg; parse its parameters
                    let rec getParams n argList ps =
                        // if we got all our parameters, done
                        if n = 0 then ((List.rev ps), argList)
                        else
                            match argList with
                            | [] ->
                                failwithf "Insufficient params for %s; got %A, needed %d more"
                                    h (List.rev ps) n
                            | [p] -> getParams (n-1) [] (p::ps)
                            | p::tl -> getParams (n-1) tl (p::ps)

                    let ps, argList = getParams a.param.Length argList []
                    // Make sure none of the parameters are another command
                    let badParams =
                        List.filter (fun (i:string) -> i.StartsWith("--")) ps
                    if not badParams.IsEmpty then
                        failwithf "%s got bad parameters '%A'" h badParams
                    // Call the parameter's config function and return the result
                    (argList, a.proc ps opts, files)
        // Otherwise this must be a input file (or bad input...)
        else
            (argList, opts, h::files)
    
    /// Recursively parse each input command.
    let rec parseCmdRec (argList:string list) opts files = 
        match argList with
        // nothing left to parse, we're done
        | [] -> (opts, List.rev files)
        // parse one item
        | h::tl -> parseCmdRec <||| (parseOneItem h tl opts files)

    parseCmdRec argList defaultOpts []