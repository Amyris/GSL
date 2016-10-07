module jsonAssembly

open System.IO
open System.Text
open System
open commonTypes
open constants
open Amyris.Bio.utils
open parseTypes
open shared
open Newtonsoft.Json
open System.Collections.Generic

/// Represents one piece of DNA for assembly, capturing its origins and relevant details
type DNASliceJson =
   {id: string; 
    extId: string; 
    dna: string;  
    sourceChr: string; 
    sourceFr: string; 
    sourceTo: string; 
    sourceFwd: bool;
    destFr: string; 
    destTo: string;
    destFwd: bool; 
    amplified: bool; 
    sliceName: string;
    sliceType: string;
    breed: string;
    description: string ; 
}

type AssemblyOutJson = 
    { id: string;
    name: string;
    dnaSlices: DNASliceJson list;
}

///  Write out a JSON file representing the output assembly list to a given path. 
let dumpJsonAssemblies (outFile:string) (assemblies : AssemblyOut list) =

    use outF = new StreamWriter(outFile)
    let assemblyHash =
        assemblies 
            |> List.map(fun a ->
            { id = a.id.Value.ToString(); name = a.name.ToString();
            dnaSlices = (a.dnaParts |> List.map(fun d ->
            { id = (match d.id with | None -> "0" | _ -> d.id.Value.ToString());
            extId = "";
            dna = String.Join("", d.dna);
            sourceChr = d.sourceChr;
            sourceFr = d.sourceFr.ToString();
            sourceTo = d.sourceTo.ToString();
            sourceFwd = d.sourceFwd;
            destFr = d.destFr.ToString();
            destTo = d.destTo.ToString();
            destFwd = d.destFwd;
            amplified = d.amplified;
            sliceName = d.sliceName.ToString();
            sliceType = 
                (match d.sliceType with
                    | REGULAR -> "REGULAR"
                    | MARKER -> "MARKER"
                    | LINKER -> "LINKER"
                    | INLINEST -> "INLINE"
                    | FUSIONST ->"FUSION");
            breed = 
                (match d.breed with 
                    | B_PROMOTER -> "B_PROMOTER"
                    | B_TERMINATOR -> "B_TERMINATOR"
                    | B_MARKER -> "B_MARKER"
                    | B_FUSABLEORF -> "B_FUSABLEORF"
                    | B_UPSTREAM -> "B_UPSTREAM"
                    | B_DOWNSTREAM -> "B_DOWNSTREAM"
                    | B_GST -> "B_GST"
                    | G_M -> "G_M"
                    | G_STOP -> "G_STOP"
                    | B_GS -> "B_GS"
                    | B_INLINE -> "B_INLINE"
                    | B_X -> "B_X"
                    | B_VIRTUAL -> "B_VIRTUAL"
                    | B_LINKER -> "B_LINKER" );
            description = d.description.ToString()}));
            })
    
    let newstring:string = Newtonsoft.Json.JsonConvert.SerializeObject(assemblyHash, Formatting.Indented)
    outF.WriteLine(newstring)

