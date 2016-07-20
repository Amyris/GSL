/// Signatures for functions implementing plugin features
module pluginDefaults

/// Customize this data structure to include extensions to compiler
type Plugin = { 
                name : string ; 
                alleleSwapAA : alleleSwapPluginDefs.AlleleSwapProvider option ;
                l2KOTitration : l2PluginDefs.L2Provider option;
                providesGlobalPragmas:string list
                providesLocalPragmas:string list
                providesCapas: string list
              }


let defaultPlugins = [
        { name = "allele swap" ; // Marker based allele swap
          alleleSwapAA = 
            Some { 
                    jobScorer = alleleSwaps.jobScorerClassicAAMut ; 
                    provider = alleleSwaps.classicAAMut
            }
          l2KOTitration = None
          providesGlobalPragmas = []
          providesLocalPragmas = []
          providesCapas = []
        } ;
        { 
            name = "level 2 ko titration"
            alleleSwapAA = None
            l2KOTitration = 
                Some {   
                    jobScorer = l2expline.l2JobScorer ;
                    explicitLocusProvider = l2expline.generateOutputsExplicitLocus
                    implicitLocusProvider = l2expline.generateOutputsTitrations
                }
            providesGlobalPragmas = []
            providesLocalPragmas = []
            providesCapas = []
        }
]




