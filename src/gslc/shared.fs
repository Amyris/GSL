module shared
/// Functions needed across several modules

open parseTypes
open commonTypes
open constants
open Amyris.utils

/// Print an integer id that might not be assigned yet 
let ambId (i :int option) = match i with None -> "?" | Some(i) -> string(i)

/// Recalculate the offsets of pieces in a list of pieces after new pieces are added in    
let recalcOffset (pieces: DNASlice list) =
    let lengths =
        pieces |> List.map (fun p -> p.dna.Length*1<ZeroOffset>)
    let _, offsets' =
        lengths |> List.fold (fun (o,r) l -> (o+l,o::r)) (0*1<ZeroOffset>,[])
    let offsets = List.rev offsets'

    List.zip pieces offsets
    |> List.map (fun (p,o) ->
        {p with destFr = o; destTo = o+(p.dna.Length-1)*1<ZeroOffset> } )

/// Produce a padding string of spaces n long
let pad n = Array.init n (fun _ -> ' ') |> arr2seq

let limitTo n (s:string) = if s.Length > n then s.Substring(0,n) else s

let extractLinker (s:string ) =
    if s.StartsWith("Linker_") then s.[7..]
    else failwithf "ERROR: unable to parse linker name '%s'" s

// List of approved linker abbreviations
let legalLinkers =
    [ '0' .. '9' ] @ [ 'A'..'E']
    |> List.map (fun c -> sprintf "%c" c) |> Set.ofSeq

let checkLinker (l:Linker) =
    if not (legalLinkers.Contains(l.l1)) then
        failwithf "ERROR: linker %s not a legal linker" l.l1
    if not (legalLinkers.Contains(l.l2)) then
        failwithf "ERROR: linker %s not a legal linker" l.l2

/// Pretty print a ParseRange line#/col# combo for error reporting
let fp (w:ParseRange) =
    sprintf "@%d,%d-%d,%d"
        (w.sp.Line+1) (w.sp.Column+1) (w.ep.Line+1) (w.ep.Column+1)