module pcrParamsParse
open System
open Amyris.Bio.primercore
open System.Text.RegularExpressions

/// Handle parsing of PCR parameters

// e.g.   ERROR: unknown pragma #pcrparams mon=50mM div=m5mM dNTP=0uM template=0uM primer=0.25uM

let re = new Regex("([^= ]*)=(\d*.?\d*)(um|mm|nm)")
let reFull = new Regex("^( *([^= ]*)=(\d*.?\d*)(um|mm|nm) *)+$")

// TODO: refactor this function to accept the pre-parsed list of strings
let parse (p:PrimerParams) (line:string) =
    if not (reFull.IsMatch(line.ToLower())) then
        failwithf "ERROR: #pcrparams should match tag=valueunit pattern (no spaces) where:
  tag   is one of  mon, div, dntp, template or primer,  
  value is a decimal number and 
  unit  is one of uM mM or nM
  e.g #pcrparams mon=50mM div=1.5mM dNTP=200uM template=0.01uM primer=0.25uM"
    else
        let rec procOne (pp : PrimerParams) (m:Match list) =
            match m with
                | [] -> pp
                | hd::tl ->
                    let v = match Single.TryParse hd.Groups.[2].Value with
                                | true,vv-> float vv
                                | false,_ ->
                                    failwithf "ERROR: parsing floating point value '%s' in #pcrparams" (hd.Groups.[2].Value)
                    let v' =  // Converted to Molar
                        match hd.Groups.[3].Value with
                            | "um" -> v*1.0<uM> |> Amyris.Bio.primercore.uM2M
                            | "mm" -> v*1.0<mM> |> Amyris.Bio.primercore.mM2M
                            | "nm" -> v*1.0<nM> |> Amyris.Bio.primercore.nM2M
                            | _ -> failwithf "inconceivable"

                    match hd.Groups.[1].Value with
                        | "mon" -> procOne { pp with monovalentConc = v' } tl
                        | "div" -> procOne { pp with divalentConc = v' } tl
                        | "dntp" -> procOne { pp with dNTPConc = v' } tl
                        | "template" -> procOne { pp with templateConc = v' } tl
                        | "primer" -> procOne { pp with primerConc = v' } tl
                        | _ as x-> failwithf "ERROR: unknown pcr parameter '%s', should be one of mon, div, dntp, template or primer" x

        [ for m in re.Matches(line.ToLower()) -> m ] |> procOne p
    

