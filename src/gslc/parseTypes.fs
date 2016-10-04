module parseTypes
open Microsoft.FSharp.Text.Lexing
open Amyris.Bio.primercore
open constants
open System
open uri
open pragmaTypes

type Slice = {left:RelPos; lApprox:bool; right:RelPos; rApprox:bool}
type MType = AA | NT
type Mutation = {f:char; t:char; pos:int; mType:MType}

type IDLoc = {i:string; s:Position; e:Position}

type Mod = 
    | MUTATION of Mutation
    | SLICE of Slice
    | DOTMOD of IDLoc


type IDAt = {id:string; s:Position; e:Position; mods:Mod list}

type ParseRange = {sp:Position; ep:Position}        
type GenePart = {gene:string; mods:Mod list; where:ParseRange}
type Linker = {l1:string; l2:string; orient:string}

/// One line of a docstring
type DocstringLine = IDLoc

[<Measure>] type Capability 
/// Build capabilities we might be allowed to optionally use for design            
type Capabilities = Set<string>

type ParseError = {message:string; s:Position; e:Position}
type GenePartWithLinker = {part:GenePart; linker:Linker option}
type Part =
    | GENEPART of GenePartWithLinker 
    | MARKERPART
    | INLINEDNA of string 
    | INLINEPROT of string 
    | HETBLOCK 
    | EXPANDED of string 
    | PARTID of IDAt
    | ERRORPART of ParseError
    | MULTIPART of PPP list

   

/// Part plus a Pragma
and PPP = { part : Part ; pr : PragmaCollection ; fwd: bool}

type DesignParams =
   {targetTm: float<C>;
    seamlessTm: float<C>; 
    seamlessOverlapTm: float<C>;
    pp: PrimerParams;
    overlapParams: PrimerParams; 
    overlapMinLen: int}
type Assembly =
   {parts: PPP list; 
    name: string option;
    uri: Uri option;
    linkerHint: string; 
    pragmas: PragmaCollection; 
    designParams: DesignParams;
    capabilities: Capabilities; 
    docStrings: DocstringLine list list}

/// All possible types of a GSL Variable
type GSLVar = 
    | GSLV_PPP of PPP
    | GSLV_ASSEM of PPP list
    | GSLV_INT of int
    | GSLV_FLOAT of float
    | GSLV_STR of IDLoc
    | GSLV_LIST of GSLVar list


// Roughage definitions
// ================================================
/// Roughage element
type RoughageElem =
   {promoter1: string;
    target1: string; 
    promoter2: string; 
    target2: string; 
    bi: bool; 
    marker: string option; 
    left: Position; 
    right: Position}

/// One classic roughage construct
type Roughage = {locus:string; marker:string; parts:RoughageElem list}
// ================================================

// Level 2 Definitions
// ================================================

type L2Id = {prefix:IDLoc option; id:IDLoc} with 
    member x.String =
        match x.prefix with
        | None -> x.id.i 
        | Some(prefix) -> sprintf "%s.%s" prefix.i x.id.i

/// Element of a level 2 line  e.g.  pABC1>gDEF2
type L2Elem = {promoter:L2Id; target:L2Id; s:Position; e:Position}

/// L2 Top level container for the expression line  z^ ; a>b ; c > d etc
type L2ExpLine = {l2Locus:L2Id option; parts:L2Elem List}

type L2Instance =
    | L2EXPLINE of L2ExpLine
    // Allow room for expansion of L2 language e.g. for loops


/// Module name to import.  Simple for now but could be dot separated series in future
type OpenPath = {modulePath:string; pragmas:PragmaCollection}

/// L2 Top level container
type L2Line =
   {l2Design: L2Instance; 
    name: string option;
    uri: Uri option; 
    pragmas: PragmaCollection; 
    capabilities: Capabilities}

/// Instruction to generate a cut in a piece of DNA which is a series of parts
type CutTarget = {assembly:GSLLine; pragmas:PragmaCollection}

and IntListMember =
    | INTMEMBER of int
    | INTRANGE of int*int

and ForItemList =
    | INTLIST of IntListMember list
    | PARTLIST of PPP list // part plus pragma list really

and GSLFunction =
   {name: string;
    args: string list; 
    body: GSLLine list; 
    pragmas: PragmaCollection}

and GSLFunctionCall = {name:string; args:PPP list}

and GSLFor =
   {varName: string; 
    items: ForItemList; 
    body: GSLLine list; 
    pragmas: PragmaCollection}

and LetLine = {varName:string; letDefinition:GSLVar}

and GSLLine = 
    | ASSEMBLY of Assembly 
    | PRAGMA of Pragma
    | GSLLINEEXPANSION of string // When a GSLLine is wholesale rewritten, the text is emitted so
    | ROUGHAGESECTION of Roughage list
    | L2 of L2Line // Section of level 2 syntax
    | DOCSTRINGSECTION of DocstringLine list
    | GSLFUNCTION of GSLFunction
    | GSLFUNCTIONCALL of GSLFunctionCall
    | CUTLine of CutTarget
    | FORBLOCK of GSLFor
    | OPENLINE of OpenPath
    | LETLINE of LetLine


let stuffModList i s e modList (aliases:Collections.Generic.Dictionary<string,GSLVar>) = 
    let isAlias, res = aliases.TryGetValue i
    // If this alias already exists, stuff the ModList into correct field of Part inside a PPP
    if isAlias then 
         match res with
             | GSLV_ASSEM([r]) 
             | GSLV_PPP(r) -> 
                 let stuffedPart = 
                     match r.part with
                     | GENEPART(gp) -> GENEPART({gp with part = ({gp.part with mods = gp.part.mods@ List.rev modList})}) // is this necessary? 
                     | PARTID(gp) -> PARTID({gp with mods = gp.mods@ List.rev modList})
                     | _ -> 
                        failwithf "ERROR: don't yet know how to implement an 'AT ID ModList' Part other than GENEPART adn PARTID... Should this Part have a ModList?"
 
                 {r with part = stuffedPart}
            
             | _ -> failwithf "ERROR: why are we trying to stuff a ModList in a var that's not a PPP?"       
    
    // Otherwise, create a new PPP
    else 
        {part = PARTID({id = i; s = s; e = e; mods = List.rev modList});
         pr = EmptyPragmas; 
         fwd=true}   // @r1234(modification)