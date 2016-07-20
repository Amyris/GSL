module l2expline

///
/// Implementation of GSL Level 2 Expression Lines
/// Modelled roughly on roughage syntax  e.g. gHO^ ; a> b ; c>d etc
///
open pragmaTypes
open parseTypes
open commonTypes
open Amyris
open constants
open sgdrefformat
open System.Text
open System
open l2PluginDefs

/// Take a list of expression elements and organize them in a balanced
/// way - e.g. splitting between two halves of a megastitch
let balance (elems: L2Elem list) =
    let roundUp (x:float) =
        if x - (x|> int |> float) > Double.Epsilon then ((int x) + 1) else int x
    // Dumb implementation, doesn't handle bipromoters
    // or any clever avoidance of repeats etc, part reuse
    let countA = float(elems.Length) / 2.0 |> roundUp
    let partsA = Seq.take countA elems |> List.ofSeq
    let partsB = Seq.skip countA elems |> List.ofSeq
    partsA,partsB

/// Base implementation for level 2 knock out / promoter titration
/// give it lowest score in case someone has a preferred implementation
let l2JobScorer _ = Some 0.0<PluginScore>

/// Takes a level-2 line regarding explicit locus and returns a list.
///
/// E.g. transforms a level-2 line, gHO^ ; pA > gB ; pC > gD into
/// {"uHO"; "pA" ; "gB" ; "###" ; "!gD ; !pA" ; "dHO"}

let generateOutputsExplicitLocus (locus:L2Id) (args:l2PluginDefs.L2DesignParams) =

    let locusWithPrefix = locus.id.i
    assert locus.prefix.IsNone

    let locusWithoutPrefix = locusWithPrefix.Substring(1)
    if not (locusWithPrefix.ToUpper().StartsWith("G")) then
        failwithf "ERROR: knockout target gene %s must start with g tag (e.g. gADH1)." locusWithPrefix
    let out = seq {
                    let partsA,partsB = balance args.line.parts
                    // Emit upstream flanking region
                    yield sprintf "#name u%s__d%s\n" locusWithoutPrefix locusWithoutPrefix
                    yield sprintf "u%s" locusWithoutPrefix
                    // First half of the parts before the marker
                    for expItem in partsA do
                        yield expItem.promoter.String
                        yield sprintf "%s" (expItem.target.String)
                    if args.megastitch then yield "###" // Marker
                    // Second half of the parts after the marker
                    for expItem in partsB do
                        yield (sprintf "!%s;!%s" expItem.target.String expItem.promoter.String)
                    // Emit downstream flanking region
                    yield sprintf "d%s" locusWithoutPrefix
                } |> List.ofSeq

    // results come back as a list of strings but we need to treat first name as a separate line and ; concat remainder
    match out with
        | name::rest ->
            let partsList = String.Join(";" , rest)
            [ GSLLINEEXPANSION(name.Trim()) ; GSLLINEEXPANSION partsList ]
        | _ -> failwithf "ERROR: L2 parsing failed"


/// Takes a level-2 line regarding promoter titrations and returns a list.
///
/// E.g. transforms a level-2 line, pA>gB ; pc>gD into
/// {"uB"; "pC" ; "gD" ; "###" ; "pA" ; "gB[1:~500]"}
let generateOutputsTitrations (args:l2PluginDefs.L2DesignParams) =

    // separates the expression pGene>gGene from the rest of the line
    let locusExp,otherExp =
        match args.line.parts with
            | [] -> failwithf "ERROR: unexpected empty L2 expression construct with no locus or parts\n"
            | hd::tl -> hd,tl
    /// the titrated gene
    let locusGene = locusExp.target.id.i.Substring(1) 
    if not (locusExp.target.id.i.ToUpper().StartsWith("G")) then
        failwithf "ERROR: titrating expression target %s must start with g tag (e.g. gADH1). Variables not supported for titrations." locusExp.target.String
    let partsA,partsB = balance otherExp
    /// the flank length
    let flank = args.rgs.[args.refGenome].getFlank()
    let out = seq{
                    // Yield upstream flnaking region. 
                    yield sprintf "#name u%s_%s_d%s\n" locusGene locusExp.promoter.String locusGene
                    yield (  sprintf "u%s" locusGene) // regular locus flanking seq
                    // First half of the parts before the marker
                    for expItem in partsA do
                        yield expItem.promoter.String
                        yield expItem.target.String
                    if args.megastitch then yield "###" 
                    // Second half of the parts after the marker
                    for expItem in partsB do
                        yield (sprintf "!%s;!%s" expItem.target.String expItem.promoter.String)
                    // Finally the titrating promoter
                    yield locusExp.promoter.String
                    // Emit downstream flanking region
                    yield sprintf "%s[1:~%A]" locusExp.target.String flank
                }  |> List.ofSeq
    match out with
        | name::rest -> [GSLLINEEXPANSION(name + String.Join(";" , rest))]
        | _ -> failwithf "ERROR: L2 parsing failed"



/// Core expansion of a single L2 expression line
let expandOneL2ExpLine  
    (providers:l2PluginDefs.L2Provider list) 
    (pragmas:PragmaCollection)
    (capabilities:Capabilities)
    (rgs:GenomeDefs)
    (line : L2ExpLine) =

    // TODO TODO:   ensuring key parts must be valid GSL parts e.g. can't titrate a non g part,  can't delete a non genome part
       
    //
    // Different styles of expansion
    //
    // Explicit locus
    ////////////////////////////////////
    // No package, just deletion
    // HO^                    =====>   uHO ; ### ; dHO
    //
    // Explicit deletion locus plus expressoin
    // HO^ ; pA > gB           ====>   uHO ; pA ; gB ; ### ; dHO
    //
    // Explicit locus with two or more genes
    // HO^ ; pA > gB ; pC > gD ====>   uHO ; pA ; gB ; ### ; !gD ; !pA ; dHO
    //
    // Titrations:
    ////////////////////////////////////
    // Titration of native gene
    // pA>gB                   =====>  uB ; ### ; pA ; ~oB (DS_G_CDS in thumper parlance)
    //
    //
    // Titration of native gene with additional expression constructs
    //
    // pA>gB ; pc>gD          ======> uB ; pC ; gD ; ### ; pA ; gD[1:~500]

   
    /// Stitch or megastitch
    let megastitch = not (pragmas.ContainsKey("stitch"))
    
    /// Which reference genome are we using
    let refGenome' =
        match pragmas.TryGetOne("refgenome") with
        // specifying a different reference genome implies a non standard
        // DNA source, so we can use that too (they can override with dnasrc)
        | Some(rg) -> rg.Trim([| ' ' ; '\t' |]) 
        | None -> constants.defaultRefGenome
       
    /// Parameters to pass to specific L2 algorithm implementation  
    let designParams = { megastitch = megastitch
                         rgs = rgs
                         refGenome = refGenome' 
                         line = line
                        }

    /// List including lower-level GSL strings
    let out = 
        match line.l2Locus with
            | Some(locusWithPrefix) ->      
                // Explicit locus case.  E.g.  gADH1^::pTDH3>mERG10
                // Choose provider
                let fn = providers |> 
                             List.choose (fun provider -> 
                                                    match provider.jobScorer capabilities with
                                                        | None -> None
                                                        | Some(score) -> Some (score,provider.explicitLocusProvider)
                                            ) |>
                             List.sortWith(fun (a,_) (b,_) -> compare b a)  |> // sort largest first
                             List.head |>
                             snd

                fn locusWithPrefix designParams // returns string list
            | None -> 
                // Implicit locus case.  E.g.  gADH1^::pTDH3>mERG10
                // Choose provider
                let fn = providers |> 
                             List.choose (fun provider -> 
                                                    match provider.jobScorer capabilities with
                                                        | None -> None
                                                        | Some(score) -> Some (score,provider.implicitLocusProvider)
                                            ) |>
                             List.sortWith(fun (a,_) (b,_) -> compare b a)  |> // largest first
                             List.head |>
                             snd

                fn designParams // returns string list
    out
    
/// Look for L2 Expression Lines in GSL tree and expand those to L1 GSL, leaving other lines alone
let expandL2ExpLine (providers:l2PluginDefs.L2Provider list) (rgs:GenomeDefs) (tree:GSLLine list) =
    
    seq { for line in tree do
            match line with
                | L2(l2Line) ->
                    match l2Line.l2Design with
                        | L2EXPLINE(l2ExpLine) ->
                                // Process a single single of level 2 GSL, emiting one or more new GSL lines (likely just expanded string representation)
                                yield! expandOneL2ExpLine providers l2Line.pragmas l2Line.capabilities rgs l2ExpLine
                | _ as y -> yield y // Leave non L2 lines alone
    } |> List.ofSeq
