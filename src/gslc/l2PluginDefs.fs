module l2PluginDefs
/// Type definitions for implementations providing level 2 functionality
/// i.e. knockout / gene titrations
open constants
open parseTypes
open sgdrefformat

type L2JobAccept = Capabilities->float<PluginScore> option

type L2DesignParams = {
        rgs:GenomeDefs
        megastitch : bool
        refGenome: string
        line : L2ExpLine
}
type L2Provider = { 
        jobScorer:L2JobAccept
        explicitLocusProvider:L2Id->L2DesignParams->GSLLine list
        implicitLocusProvider:L2DesignParams->GSLLine list
}

