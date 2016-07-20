module codoptSupport

open Amyris.Bio.utils
open Amyris.Bio
open sgdrefformat
               
/// Given an amino acid sequence and codon frequency, pick a dumb optimal codon sequence using just the most frequent
/// amino acids.  This will only really work ok for very short sequences and doesn't check for problematic sequences.                        
let codonOptDumb (codonTable:Map<string,float>) (prot:string) =
    let findBest (aa:char) =
        codonTable |> Seq.choose (fun kv -> if biolib.codon2aa (kv.Key.ToCharArray()) = aa then Some(kv.Value,kv.Key) else None)
            |> List.ofSeq |> List.sortWith (fun (v1,_) (v2,_) -> compare v2 v1) |> List.head |> fun (_,codon) -> codon
            
    prot |> Seq.map (findBest) |> Seq.concat |> Array.ofSeq |> arr2seq
       
/// typical sequences we want to exclude from synthesized DNA       
let defaultBadSeqs = [ "GCTCTTC" (*sapI*); "GAAGAGC" (* sapI RC *); "GGGGGGG" (* Trp Gly TGG GGx *); 
                    "GTTTAAAC" (*PmeI*); "GGTCTC" (* BsaI *) ; "GAGACC" (* BsaI RC *) ; 
                    "CACCTGC" (* AarI *) ; "GCAGGTG" (* Aar1 RTC *) ;
                    "AAAAAAA" ;  "TTTTTTT" ; "CCCCCCC"  ] 

/// Check if there is a codonavoid definition in the reference genome and if so,
/// return the list of codons not to touch
let getCodonAvoid (gd:GenomeDef) =
    let env = gd.Env
    match env.TryFind("codonavoid") with
    | None -> []
    | Some(x) ->
                x.Split([| ' ' ; '\t' |])
                |> List.ofArray
                |> List.map (fun s -> s.Replace('U','T'))


/// Codon opt parameters are required to implement
/// this interface
type ICodonOptParam =
    abstract member Seed:unit->int


/// Codon opt data structures for genome specific parameters
/// need to implement this interface
type ICodonOptData =
    abstract member CodonAvoid:unit->string list

/// Implementation of a codon optimizer
/// 'CodopParams describe user set parameters for the algorithm
/// 'CodoptData describe genome specific parameters for the algorithm

type CodonOptimizer<'CodoptParams,'CodoptData  when 'CodoptParams :> ICodonOptParam and 'CodoptData :> ICodonOptData> = { 
    doCodopt : bool->'CodoptParams->'CodoptData->int->string->string
    parseCodoptParams: (string*string) list->'CodoptParams
    loadCodoptData: GenomeDef->'CodoptData
}