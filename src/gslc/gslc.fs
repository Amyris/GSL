open System.IO
open commandConfig // Command line arguments, defaults, etc
open gslcProcess // Top-level compiler operations

// Helper libs for oligo design, sequence parsing all in Amyris.Bio.dll
// These imports are only needed for the temporary primer test function below.
open Amyris.Bio
open primercore
open utils
open constants

/// Test bed for investigating primer misadventure
let testPrimer() =
    let template = "GCCAGCGATAGGAGTCCTTGGTTTAG".ToCharArray()

    for i in {15..26} do
        printfn "%d %A" i (Amyris.Bio.primercore.temp defaultParams template i)

    let fwd = true
    let pen = { primercore.defaultParams with maxLength = 30 ; tmPenalty = 3.0(* template.Length *)}
    let task : OligoTask = { tag = if fwd then "PF" else "PR" ;
                             temp = template ;
                             align = ANCHOR.LEFT ;
                             strand = STRAND.TOP ; offset =0 ; targetTemp = ryseLinkerTargetDefault;
                             sequencePenalties  = None }

    let res = oligoDesign true pen task
    printf "pen=%A \n %A %s %d\n" pen (res.Value.temp) (arr2seq res.Value.oligo) res.Value.oligo.Length

    ()

let Main() =
    // Get the input args and handle a few special cases
    let args =
        match (argv() |> List.ofArray |> List.tail) with
        // print usage info
        | [] | "--help"::_ ->
            printUsage()
            exit 1
        // temporary entry point for primer test
        | "--test"::_ ->
            testPrimer() |> ignore
            exit 1
        // parse args and run gslc
        | x -> x

    // Try to parse the arguments
    let opts, inputFiles =
        try
            parseCommandLineArgs args
        with err ->
            printfn "ERROR: %s" err.Message
            exit 1

    if not opts.quiet then printf "// GSL compiler version %s\n" version

    // For the moment we only support one file argument.
    let inputFile =
        match inputFiles with
        | [] ->
            printfn "ERROR: no input files specified"
            exit 1
        | [x] -> x
        | n ->
            printfn "ERROR: GSLC only supports one input file at a time. Got %A" n
            exit 1

    if not (File.Exists inputFile) then
        printf "ERROR: can't find file '%s'\n" inputFile
        exit 2

    // If selected, perform lexing and quit.
    if opts.lexOnly then
        lexTest opts inputFile
        exit 0

    try
        // Load static assets and initialize caches.
        let ga = loadGlobalAssets opts

        // Start with input from the users supplied file
        let input = File.ReadAllText inputFile

        /// Combined default plugins plus any extensions in plugins.fs
        let allPlugins = plugins.plugins@pluginDefaults.defaultPlugins
        // Run GSLC on the input file.
        let compilerOutput = processGSL opts allPlugins ga input
        match compilerOutput with
        | FinishedAssemblies(assems) -> writeOutput opts ga assems
        | ExpandedGSL(newGSL) -> printf "%s" newGSL
        exit 0
    with e ->
        let prependERROR (s:string) =
            if s.StartsWith("ERROR") then s
            else sprintf "ERROR: %s" s

        printfn "%s" (prependERROR e.Message)
        if e.InnerException <> null then
            printfn "%s" (prependERROR e.InnerException.Message)
            if opts.verbose then printfn "ERROR: location: \n%s" e.InnerException.StackTrace
        else
            if opts.verbose then printf "ERROR: location:\n%s" e.StackTrace

        exit 1

Main()
printf "Done!\n"