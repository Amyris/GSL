/// Signatures for functions implementing plugin features
module pluginDefaults

/// Customize this data structure to include extensions to compiler
type Plugin = { 
                name : string ; 
                alleleSwapAA : alleleSwapPluginDefs.AlleleSwapProvider option ;
                l2KOTitration : l2PluginDefs.L2Provider option;
                providesPragmas:pragmaTypes.PragmaDef list
                providesCapas: string list
              }
/// Default empty plugin that can be selectively altered to create relevant fields
let defaultPlugin = {
            name = "example plugin" 
            alleleSwapAA = None
            l2KOTitration = None
            providesCapas = []
            providesPragmas = []
        }

let defaultPlugins = [
        { name = "allele swap" ; // Marker based allele swap
          alleleSwapAA = 
            Some { 
                    jobScorer = alleleSwaps.jobScorerClassicAAMut ; 
                    provider = alleleSwaps.classicAAMut
            }
          l2KOTitration = None
          providesPragmas = []
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
            providesPragmas = []
            providesCapas = []
        }
]




