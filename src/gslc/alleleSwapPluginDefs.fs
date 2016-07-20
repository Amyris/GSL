module alleleSwapPluginDefs
open Amyris.Bio
//open System
open constants
open parseTypes
//open commonTypes
open sgdrefformat
open IO.CodonUsage 
//open utils
//open biolib
//open primercore
//open ryse // for getrabit
type EndPref = NTERM | CTERM | NONETERM

/// Amount of extra dna adjacent to the ORF to include
let orfPlusMargin = 100

type AlleleSwapJobAccept = Capabilities->float<PluginScore> option
type AlleleSwapDesignParams = {
    verbose:bool
    longStyle:bool
    endPref:EndPref
    codonLookup:CodonLookup
    gene:string
    name:string
    rg:GenomeDef
    f:sgd.Feature
    m:Mutation
    len:int<ZeroOffset>
    mutOff:int<ZeroOffset>
    orf:char[]
    orfPlus:char[]
}

type AlleleSwapProvider = { jobScorer:AlleleSwapJobAccept ; provider:AlleleSwapDesignParams->string }
