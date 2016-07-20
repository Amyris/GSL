/// IO routines for loading the reference file format
module sgdrefformat
open Amyris.Bio

open sgd
open utils
open parseTypes
open Amyris.Bio.SuffixTree
open constants

(*
$ ls c\:/Amyris/data/REF/cenpk/cenpk
cenpk.fsa           cenpk_features.tab

*)
open System.IO
open System

let loadEnv (p:string) =
    if File.Exists(p) then
        eachLineIn p |> Seq.map (fun x -> x.Split([|"="|],StringSplitOptions.None))
                     |> Seq.choose (fun cols -> match cols with 
                                                    | [| a; b |] -> Some(a.Trim(),b.Trim())
                                                    | _ -> printf "WARNING: bad config entry '%A'" cols ; None
                                    )
                    |> Map.ofSeq
    else
        Map.empty

/// Init genome definition, genes etc
type GenomeDef(p:string) as this = class
    let mutable fasta = None
    let mutable feats = None
    let mutable featIndex = None
    let mutable suffixTreePath : string option = None
    let mutable suffixTree : SuffixTreeDisk option = None

    let mutable env : Map<string,string> = Map.empty
    let mutable envLoaded = false

    let ensureConfigLoaded() =
        if not envLoaded then
            env <- opj p "config.txt" |> loadEnv
            envLoaded <- true

    let ensureLoaded() =
        match fasta with
            | None -> this.Load(p)
            | _ -> ()

    do
        //if p <> "" then this.Load(p)
        ()

    new() = new GenomeDef("")
    
    member x.Env = 
        ensureConfigLoaded()
        env        
    member x.Load(refDir:string) =
        ensureConfigLoaded()
        let projName = baseName refDir
        let featsPath = opj refDir (sprintf "%s_features.tab" projName)
        let fastaPath = opj refDir (sprintf "%s.fsa" projName)
        suffixTreePath <- Some(opj refDir "suffixTree.st")
        
        fasta <- Some(Amyris.Bio.biolib.readReference fastaPath)
        feats <-  Some(sgd.loadFeatures featsPath)
        
        let i1 = feats.Value |> Array.mapi (fun i f -> f.sysName,i)
        let i2 = feats.Value |> Array.mapi (fun i f -> f.gene,i)
        featIndex <- Array.concat [ i1 ; i2 ] |> Seq.filter (fun (x,_) -> x <> "") |> Map.ofSeq |> Some 
    
    member x.get(g:string) =
        ensureLoaded()
        match featIndex with
            | None -> failwith "Access to unloaded GenomeDef"
            | Some(fi) -> feats.Value.[fi.[g]]
    
    member x.Dna(errorContext:string,chr:string,l':int<ZeroOffset>,r':int<ZeroOffset>) =
            ensureLoaded()
            let l,r = l'/1<ZeroOffset> , r'/1<ZeroOffset>
            if not (r>=l) then
                failwithf "ERROR: For %s Attempt to retrieve  DNA slice with reversed coordinates %s:%d-%d\n" errorContext chr l r
                
            match fasta with
                | None -> failwithf "For %s Access to unloaded GenomeDef Fasta" errorContext
                | Some(f) ->
                    if not (f.ContainsKey(chr)  ) then failwithf "ERROR: For %s unknown chromsome '%s'" errorContext chr
                    if l < 0 || l >= f.[chr].Length then failwithf "ERROR: For %s coordinate '%d' outside chromsome %s length = %d" errorContext l chr (f.[chr].Length)
                    if r < 0 || r >= f.[chr].Length then failwithf "ERROR: For %s coordinate '%d' outside chromsome %s length = %d" errorContext r chr (f.[chr].Length)
                    f.[chr].[l..r]
                    
    member x.SuffixTree =
                    ensureLoaded()
                    match suffixTree with 
                        | Some(st) -> st 
                        | None -> if (File.Exists(suffixTreePath.Value)) then 
                                        //FIXIFX: this next line occasionally throws an 'Cannot create a file when that file already exist error
                                        (*Cannot create a file when that file already exists.

ERROR: location:
   at System.IO.__Error.WinIOError(Int32 errorCode, String maybeFullPath)
   at System.IO.MemoryMappedFiles.MemoryMappedFile.CreateCore(SafeFileHandle fileHandle, String mapName, HandleInheritability inheritability, MemoryMappedFileSecurity memoryMappedFileSecurity, MemoryMappedFileAccess access, MemoryMappedFileOptions options, Int64 capacity)
   at System.IO.MemoryMappedFiles.MemoryMappedFile.CreateFromFile(String path, FileMode mode, String mapName, Int64 capacity, MemoryMappedFileAccess access)
   at System.IO.MemoryMappedFiles.MemoryMappedFile.CreateFromFile(String path, FileMode mode, String mapName, Int64 capacity)
   at Amyris.SuffixTree.SuffixTreeDisk.openMM(String path)
   at Amyris.SuffixTree.SuffixTreeDisk..ctor(String path)
   at sgdrefformat.GenomeDef.get_SuffixTree() in C:\seq\GSL\gslc\sgdrefformat.fs:line 99

                                        *)
                                        suffixTree <- Some(new SuffixTreeDisk(suffixTreePath.Value))  
                                        suffixTree.Value
                                  else failwithf "ERROR: loading suffix tree, unable to find file %s" suffixTreePath.Value           
    member x.IsValid(f:string) =
        ensureLoaded()
        
        match featIndex with
            | None -> failwith "Uninitialized GenomeDef instance\n"
            | Some(fi) -> fi.ContainsKey(f)   
    interface System.IDisposable with
        member x.Dispose() = match suffixTree with | None -> () | Some(std) -> (std :> System.IDisposable).Dispose()
    member x.getFlank() =
        if x.Env.ContainsKey("flanklen") then
            (x.Env.["flanklen"] |> int ) * 1<OneOffset>
        else 1<OneOffset> * flanklenDefault
                 
end

type GenomeDefs = Map<string,GenomeDef>

