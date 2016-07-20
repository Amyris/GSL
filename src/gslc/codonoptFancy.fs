module codonoptFancy
open System
open System.Text.RegularExpressions
open Amyris
open Amyris.utils
open MathNet.Numerics.Random
open System.Collections.Generic
open System.IO

open System.Threading
open codoptSupport
open sgdrefformat
open pragmaTypes

(*

UUU 26.1(170666)  UCU 23.5(153557)  UAU 18.8(122728)  UGU  8.1( 52903)
UUC 18.4(120510)  UCC 14.2( 92923)  UAC 14.8( 96596)  UGC  4.8( 31095)
UUA 26.2(170884)  UCA 18.7(122028)  UAA  1.1(  6913)  UGA  0.7(  4447)
UUG 27.2(177573)  UCG  8.6( 55951)  UAG  0.5(  3312)  UGG 10.4( 67789)

CUU 12.3( 80076)  CCU 13.5( 88263)  CAU 13.6( 89007)  CGU  6.4( 41791)
etc

*)

// -----------------------------------------------------------------------------------------------------
// revTrans implementation of IDTHi like strategy
//
open Amyris.biolib
open Amyris.IO.CodonUsage
open Amyris.SuffixTree

type PopScore = { score : float ; pos : int}
type Chrom = {chrom : char array ; score : float ; problems : bool array}
type Result = { count : int ; best : Chrom ; prot : char []}

let defaultMinFreq = 0.1 //0.2 // 0.05
let defaultMaxRank = 5
let mutRate = 0.005
let mutRateHigh = 0.5
let popSize = 100
let generations = 100
let genNoProgress = 10
let selectionBias = 0.25
let defaultMerLen = 6
let defaultStartMargin = 100
let defaultCodonOptSeed = 170270

// Versions
// 1: legacy, initial implementation
// 2: includes 5prime optimization to match NYT usage
type CodonOptParams = { startMargin: int ; merLen : int ; badSeqs : string list ; minFreq : float ; maxRank : int ; 
                        codonAvoid : string list ; seed : int ; 
                        optFivePrime : bool ;
                        fivePrimeWindow : int;
                        algVersion : int;
                        prefixes : string list;
                        globalRepeatCheck : bool

                      } with
    interface ICodonOptParam with
        member x.Seed() = x.seed

let defaultCodonOptParamsV2 = { startMargin = defaultStartMargin ; badSeqs = defaultBadSeqs ; maxRank = defaultMaxRank ; 
                                merLen = defaultMerLen ; minFreq = defaultMinFreq ; codonAvoid = [] ; seed = defaultCodonOptSeed ; optFivePrime = true;
                                fivePrimeWindow = 27 ;
                                algVersion = 2;
                                prefixes = ["ACCTCCCGCGACCTCCAAAATCGAACTACCTTCACA";//A linker
                                            "ACGCACGCACACTCCCGACAGACAACTAGCTTGATA";//B linker
                                            "ACCCCACCCGAAGTCGCGCAACCAACTAACTTTACA" // C reversed
                                ];
                                globalRepeatCheck = false
                            }

let defaultCodonOptParamsV1 = { startMargin = defaultStartMargin ; badSeqs = defaultBadSeqs ; maxRank = defaultMaxRank ; 
                                merLen = 8 ; minFreq = defaultMinFreq ; codonAvoid = [] ; seed = defaultCodonOptSeed ; optFivePrime = false;
                                fivePrimeWindow = 0;
                                algVersion = 1;
                                prefixes = []
                                globalRepeatCheck = false

                            }

let defaultCodonOptParams = defaultCodonOptParamsV2

// ================= 5prime scoring support
let base2Idx c = match c with | 'A' -> 0 | 'T' -> 1 | 'C' -> 2 | 'G' ->3 | _ as x -> failwithf "ERROR: bad base %c" x

let load5Prime path =
    // Load raw nucleotide freqency data.  Should be A T C G header and frequency, one position per row
    let fivePrime = 
        Amyris.grid.TabGrid(
            File.ReadAllLines(path)
            |> Array.map (fun line -> 
                                line.Split(
                                    [| ',' ; ' ' ; '\t'|],
                                    StringSplitOptions.RemoveEmptyEntries
                                )
                            )
        )

    // Ensure base ordering in headers matches our internal base to column scheme
    "ATCG" |> Seq.iteri (fun i nt -> assert (base2Idx nt = i))
    Array2D.init 
                    ((Seq.length fivePrime.Rows)-1) 
                    4 
                    (fun row col -> 
                        match fivePrime.Row(row).[col] with
                            | "0.0" -> 999.0 // safer to match on string than float ?
                            | "0" -> 999.0
                            | _ as x -> x |> float |> log10 |> (*) (-1.0)
                    )
/// Apply 5prime scoring function to dna array starting at 'A' 'T' 'G' out to window length 'window'
let score5Prime window (fivePrimeScores : float[,]) (dnaORF:char[]) =
    match ([| for i in 3..(min (fivePrimeScores.GetLength(0)-1) (dnaORF.Length-1) |> min window) ->
                fivePrimeScores.[i,dnaORF.[i] |> base2Idx]
            |] ) with
    | [||] -> 0.0 // guard against empty dna sequences whose scores can't be averaged
    | _ as x -> x |> Array.average

let rank5Prime (fivePrimeScores : float[,]) (dnaORF:char[]) =
    let rank (c:char) pos =
        [|for i in 0..3 ->
            if (fivePrimeScores.[pos,base2Idx c] - fivePrimeScores.[pos,i]) > 0.02 then 1 else 0
        |] |> Array.sum

    [| for i in 0 .. min (fivePrimeScores.GetLength(0)-1) (dnaORF.Length-1) ->
        rank (dnaORF.[i]) i
    |]
    
// ============================================
type ScoreResult = { totalAffected:int ; probs:bool []; fivePrimeScore : float ; fiveRanks : int []}
            with member x.Score = (float x.totalAffected)+x.fivePrimeScore*20.0

/// Evaluate a dna string given a set of optimization parameters.  Returns # of problems and vector in amino acid space of problem areas
let scoreWithCoords (rng:MersenneTwister) (cop : CodonOptParams) (stPrefix:SuffixTree) (fivePrimeScores:float[,] option) (s:string) = 
        
        let st = new SuffixTree(s)
        st.buildFwdChain()
        let probs = Array.init s.Length (fun _ -> false)
        let fivePrimeScore = match fivePrimeScores with
                                | None -> 0.0
                                | Some(fps) -> 
                                    score5Prime cop.fivePrimeWindow fps (s.ToUpper().ToCharArray())
        let fiveRanks = match fivePrimeScores with
                            | None -> [||]
                            | Some(fps) ->
                                    rank5Prime fps (s.ToUpper().ToCharArray())
        for x,y in 
            seq { 
                for b in cop.badSeqs do 
                    for x in st.FindAll(b) do
                        yield x,x+b.Length-1
                for i in {0..min cop.startMargin (s.Length-cop.merLen)} do
                    let revMer = s.[i..i+cop.merLen-1].ToCharArray() |> revComp |> arr2seq 
                    let hits = st.FindAll(revMer)
                    let hitsInMargin = hits |> Array.filter (fun i -> i<cop.startMargin)
                    for x in hitsInMargin do
                        yield x,x+cop.merLen-1 // Mutate matching dest

                    

                    if hitsInMargin.Length<> 0 || stPrefix.FindAll(revMer).Length<>0 then
                        yield i,i+cop.merLen-1 // Mutate source as well
                // If needed, check for fwd fwd (tandem) repeats that could cause
                // assembly or stability problems
                if cop.globalRepeatCheck then
                    let repLen = 8
                    
                    for i in {0..(s.Length-repLen)} do
                        let thisMer = s.[i..i+repLen-1].ToCharArray()
                        for x in (st.FindAll(thisMer |> arr2seq ) |> Array.filter (fun j -> i<>j) ) do
                            if rng.Next()%10=0 then // mutate sparingly, we only need to disrupt a repeat in one place
                                yield x,x+repLen-1 // Mutate matching dest
                                yield i,i+repLen-1 // Mutate source as well

            }
            do
                for z in {x/3..y/3} do // convert to amino acid position
                    probs.[z] <- true // mark as problematic
        let totalAffected = (probs|> Array.fold (fun count x -> if x then (count+1) else count) 0 )
        {totalAffected = totalAffected ; probs = probs;fivePrimeScore = fivePrimeScore ; fiveRanks = fiveRanks}


/// Genome specific data needed to generate a codon optimized sequence for this implementation
type CodonOptData = {
    freq : Map<string,float>
    fivePrimeData:float [,] option
    codonAvoid : string list
}

//let loadCodonData (gd:GenomeDef) =
//    let env = gd.Env
//    let codonAvoid = getCodonAvoid gd
//    { codonAvoid = codonAvoid}

let parseCodonOptParams (pr:PragmaCollection) =
    let seed =
        match pr.TryGetOne("seed") with
        | None -> defaultCodonOptSeed
        | Some(s) ->
            match Int32.TryParse s with
            | true,s' -> s'
            | _ -> failwithf "ERROR: invalid integer '%s' for #seed" s

    let fivePrimeWindow,globalRepeatAvoidUser =
        match pr.TryGetOne("codonopt") with
        | None -> defaultCodonOptParams.fivePrimeWindow,defaultCodonOptParams.globalRepeatCheck
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
                | None -> defaultCodonOptParams.fivePrimeWindow
                | Some(v) -> int v
            let grc =
                match parts.TryFind "repeatcheck" with
                | None -> defaultCodonOptParams.globalRepeatCheck
                | Some(v) ->
                    match v.ToLower() with
                        | "true" -> true
                        | "false" -> false
                        | _ as x -> failwithf "ERROR: bad #codonopt repeatcheck param '%s' - should be true or false" x
            fiveW,grc

    {defaultCodonOptParams with 
        fivePrimeWindow = fivePrimeWindow ; 
        globalRepeatCheck = globalRepeatAvoidUser;
        seed = seed
    }

//let codonOptFancy verbose (cop:CodonOptParams) (codonTable':Map<string,float>) (fivePrimeScores:float[,] option) (localSeed:int) (protSeq:string) =
let doCodonOpt verbose (cop:CodonOptParams) (codoptData:CodonOptData) (localSeed:int) (protSeq:string) =

    let avoid = cop.codonAvoid |> Set.ofList

    let codonTable = codoptData.freq |> Map.filter (fun k _ -> avoid.Contains(k) |> not )
    let rng = new MersenneTwister(localSeed)

    /// Adjust preferred codon use to be even
    let renormalizeCutByFivePrime (cut:CodonLookup) =
        let byAA = cut.byAA |> Seq.map (fun kv -> 
                                            let freq = 1.0/(List.length kv.Value |> float)
                                            kv.Key,(kv.Value |> List.map (fun x -> { x with relFreq2 = freq}))
                                ) |> Map.ofSeq

        { byAA = byAA ; byCodon = cut.byCodon}

    /// Restricted codon usage table focused on top codonbs
    let cut = prepCUT cop.minFreq cop.maxRank codonTable

    /// This is a more full codon usage table for use at the 5prime end
    let cutFullPre = prepCUT 0.01 99 codonTable

    let cutFull = renormalizeCutByFivePrime cutFullPre
    /// Proportion of identical residues
    //let similarity (a:string) (b:string) = (Seq.zip a b |> Seq.fold (fun total (a,b) -> if a=b then total+1 else total) 0 |> float) / (float a.Length)
        
    // Prepare input sequence
    let p = protSeq.ToCharArray()

    /// Build an initial guess codon usgae proportionally to preference
    let randProt() = p |> Array.mapi (fun i aa -> 
                                            if i*3<cop.fivePrimeWindow then
                                                cutFull.Choose(rng.NextDouble(),aa)
                                            else
                                                cut.Choose(rng.NextDouble(),aa)
                                      ) |> Array.concat

    let mutRateFivePrime = if cop.fivePrimeWindow = 0 then 0.0 else 2.0/(float cop.fivePrimeWindow)
    /// Targeted mutation of an existing sequence  
    let mutate (mutProfile:bool array) (a : char[]) =
        let a' = Array.copy a // avoid mutating argument
        let l = a.Length/3
        for i in {0..l-1} do
            if i*3 < cop.fivePrimeWindow then
                // Five prime end case
                if rng.NextDouble() < mutRateFivePrime then
                    // Use full codon table for this part of the protein
                    a'.[i*3..i*3+2] <- cutFull.Choose(rng.NextDouble(),p.[i])
            else
                if rng.NextDouble() < (if mutProfile.[i] then mutRateHigh else mutRate) then
                    a'.[i*3..i*3+2] <- cut.Choose(rng.NextDouble(),p.[i])

        a'

    /// Take two chromosomes and matching mutation hotspot preferences and
    /// generate mutated/ crossed over versions
    let crossover (a:char[]) (mutProfileA:bool [])  (b:char[]) (mutProfileB:bool []) =
        // pick crossover point
        let i = rng.Next(a.Length/3)*3

        // Build mutated versions of chromosomes
        let a' = mutate mutProfileA a
        let b' = mutate mutProfileB b

        // Reconstruct mutated chromosome
        if i = 0 then
            (b', a')
        else
            (Array.concat [a'.[..i-1] ; b'.[i..] ]  , Array.concat [b'.[..i-1] ; a'.[i..] ] )

    
    let stPrefix = new SuffixTree(String.Join("N",cop.prefixes))
    stPrefix.buildFwdChain()
    /// Create initial random population
    let startPop = Array.init popSize (fun _ -> 
                                            let c = randProt()
                                            let r = scoreWithCoords rng cop stPrefix codoptData.fivePrimeData (arr2seq c)
                                            //let score,probs = r.totalAffected,r.probs
                                            {chrom = c ; score = (if cop.algVersion = 1 then float r.totalAffected else r.Score) ; problems = r.probs })

    /// Take through one round of mating/evolution
    let cycle (popIn : Chrom array) =
        // Large score is bad, sort so worst first
        let pop = popIn |> Array.sortWith (fun a b -> compare b.score a.score)

        // Assume minimization, worst is first after sorting, best is last
        let worst = pop.[0].score

        // Not sure why this was in there unless we want to do some scaling
        // let avg = pop |> Seq.map (fun x ->x.score) |> Seq.average

        // Invert scores so you get a bigger score for best chromosomes
        // which would previously have had low scores
        let modScore = pop |> Array.map (fun c -> c,(worst-c.score))

        // Get total modified score
        let total = modScore |> Array.fold (fun total (_,s) -> total + s) 0.0

        let choose() =
                // Chromosomes are sorted worst to best.  We pick a score and serially
                // subtract chromosomes till we find the range matching the random
                // score.  We do it smallest to largest to avoid rounding errors at the
                // end.  We bias the picking towards the top part of the population to
                // create more fitness difference between the best and the worst
                let r = (rng.NextDouble()*selectionBias + (1.0-selectionBias)) *float(total)

                let rec find i f = 
                        let chr,score = modScore.[i]
                        if score > f || i = modScore.Length-1 then chr else find (i+1) (f-score)
                let x = find 0 r
                //printf "c%d " x
                x
        let newPop = seq {
                            for _ in {0..pop.Length/2-1} do
                                let c1 = choose()
                                let c2 = choose()
                                let c3,c4 = crossover c1.chrom c1.problems c2.chrom c2.problems
                                let r3 = scoreWithCoords rng cop stPrefix codoptData.fivePrimeData (arr2seq c3)
                                //let s3,p3 = r3.totalAffected,r3.probs
                                let r4 = scoreWithCoords rng cop stPrefix codoptData.fivePrimeData (arr2seq c4)
                                //let s4,p4 = r4.totalAffected,r4.probs
                                yield {chrom=c3 ; score = (if cop.algVersion = 1 then float r3.totalAffected else r3.Score) ; problems = r3.probs}; 
                                yield {chrom = c4; score = (if cop.algVersion = 1 then float r4.totalAffected else r4.Score)  ; problems = r4.probs}
                            } |> Array.ofSeq
        newPop.[0] <- pop.[pop.Length-1] // Preserve best population member outright
        newPop

    let rec run verbose count lastScore sameFor (p : Chrom [])=
        let minScore = p |> Seq.map (fun c -> c.score) |> Seq.min
        if verbose then printf "cycle=%d min=%f\n" count minScore
        let sameFor',lastScore' = 
            match lastScore with
                | Some(ls) when abs(ls-p.[0].score)<0.000001 -> sameFor+1,Some(p.[0].score)
                | _ -> 0,Some(p.[0].score)

        if count = generations || p.[0].score <  0.000001 || sameFor' >= genNoProgress then count,p
        else 
            let newPop = cycle p
                
            run verbose (count+1) lastScore' sameFor' newPop

    let count,finalPop = run verbose 0 None 0 startPop
    let best = finalPop.[0]

    //let s = best.chrom |> arr2seq

    let trans = Amyris.biolib.translate best.chrom
    for i in 0..3..best.chrom.Length-3 do
        let codon = best.chrom.[i..i+2] |> arr2seq
        assert (avoid.Contains(codon) |> not && cutFullPre.byCodon.ContainsKey(codon))
    assert (trans = protSeq.ToCharArray())

    //printf "     %s\n" s
    //for a in args.avoid do
    //    printf "%-5.5f %s\n" (similarity s a.sequence) (a.sequence)
    if verbose then
        printfn "codopt: %s" (arr2seq best.chrom)
        printfn "        %s" (
                                String.Join (
                                   "",
                                    [| for i in 0..best.chrom.Length-1 -> string (i%3) |]
                                )
                             )
    let final = {count = count ; best = best ; prot = p}
    final.best.chrom |> arr2seq // convert to final DNA sequence as a string



/// Given a library root, load (and cache) genomes on demand
type CodonTableCache(libRoot:string) = class
    let cache = new Dictionary<string,CodonOptData>()
    let semaphore = "sequential access!"
    do
        ()

    member x.Get(gd:GenomeDef,genome:string) =
        lock semaphore (fun () ->
                        if cache.ContainsKey(genome) then cache.[genome]
                        else
                            let table = Path.Combine([| libRoot ;genome ;"codon.txt" |]) |> loadCodonTable
                            let fivePath = Path.Combine ( [| libRoot ; genome ; "5prime.csv" |])
                            let fivePrimeData = if File.Exists fivePath then Some(load5Prime fivePath) else None
                            let data = { freq = table ; fivePrimeData = fivePrimeData ; codonAvoid = getCodonAvoid gd }
                            cache.[genome] <- data
                            data
                        )
end


