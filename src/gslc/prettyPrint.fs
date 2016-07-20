module PrettyPrint
open shared // common routines used in several modules
open System
open System.Text
open thumper // Thumper output formats, dumping etc
open ryse    // RYSE architecture
open cloneManager
open ape
open l2expline
open parseTypes
open Amyris.Bio.utils
open constants

/// Pretty print a RelPos
let printRP (l:RelPos) = sprintf "%A/%s" l.x (match l.relTo with |FIVEPRIME -> "S" | THREEPRIME -> "E")
let printSlice (s:Slice) = 
    sprintf "[%s%A%s:%s%A%s]" 
        (if s.lApprox then "~" else "") s.left.x (match s.left.relTo with | FIVEPRIME -> "S" | THREEPRIME -> "E") 
        (if s.rApprox then "~" else "") s.right.x (match s.right.relTo with | FIVEPRIME -> "S" | THREEPRIME -> "E")    

/// Pretty print line of Roughage
let prettyPrintRoughage (rLine:Roughage) =
    let header = 
        if rLine.locus <> "" then
            if rLine.marker <> "" then
                [sprintf "%s^[%s]" rLine.locus rLine.marker]
            else [sprintf "%s^" rLine.locus ]
        else [] // No leading section

    // Now individual elements
    let tail = rLine.parts |> List.map (fun re ->
                                            let elementHead =  
                                                if re.bi then
                                                    sprintf "%s<%s-%s>%s" re.target1 re.promoter1 re.promoter2 re.target2 
                                                else sprintf "%s>%s" re.promoter1 re.target1
                                            let elementTail = match re.marker with Some(x) -> sprintf "[%s]" x | None -> ""
                                            elementHead + elementTail
                                        )


    String.Join("::",header@tail)

/////Pretty print a param 
//let prettyPrintParam (p:Param)  = 
//    match p with
//    | PARTPARAM(s) -> s
//    | ATPARAM(s) -> "@"+s

let expandMods (ml:Mod list) =
        seq {
                for m in ml do
                    match m with
                        | MUTATION(m) -> yield sprintf "%c%c%d%c" (match m.mType with 
                                                                        | AA -> '$' 
                                                                        | NT -> '*') m.f m.pos m.t
                        | SLICE(s) -> yield printSlice s 
                        | DOTMOD(d) -> yield sprintf ".%s" d.i 
                } |> fun x -> String.Join("",x)   

let rec printPPP ppp = 
    let partOut = 
        match ppp.part with
        | HETBLOCK -> "~ " // don't do anything at this level
        | INLINEDNA(s) -> sprintf "/%s/ " s   // inline DNA sequence
        | INLINEPROT(s) -> sprintf "/$%s/ " s // inline protein sequence
        | MARKERPART -> "### "
        | PARTID(p) -> sprintf "@%s" p.id + (expandMods p.mods)
        | EXPANDED(s) -> s // Part that was already expanded into a string
        | ERRORPART(_) -> "" // Don't include anything for an error
        | MULTIPART(m) -> sprintf "(%s)" (String.Join(";",(m |> List.map printPPP)))
        | GENEPART(gp) ->
            let lOut =
                match gp.linker with
                    | None -> ""
                    | Some(l) ->
                        sprintf "%s-%s-%s-" l.l1 l.l2 l.orient // Emit linker
            let p = gp.part
                            
            let gOut = p.gene
            let modOut = expandMods p.mods
                                
            lOut + gOut + modOut // String.Join("",Array.ofSeq modOut)   
    // Now add in any inline pragma part with braces, ; separated etc
    let prOut = 
        if ppp.pr.pmap.Count=0 then "" 
        else 
            ppp.pr.pmap
            |> Seq.map (fun pv -> sprintf "#%s %s" pv.Key (pv.Value.args |> String.concat " "))
            |> fun ss -> String.Join(";",ss)
            |> sprintf "{%s}"
    (if ppp.fwd then "" else "!") + partOut + prOut
/// Pretty print line of GSL
let rec prettyPrintLine (line:GSLLine) =
    /// Pretty print parts lists
    let prettyPrintPartsList (parts:PPP list) =
        seq {
            for ppp in parts do
                yield printPPP ppp
                } |> fun x -> String.Join(";",Array.ofSeq x)

    /// Pretty print any sort of GSLVAR
    let rec prettyPrintGSLVar (var: GSLVar) = 
        match var with
        | GSLV_PPP(var) -> printPPP var
        | GSLV_ASSEM(var) -> prettyPrintPartsList var
        | GSLV_INT(var) -> sprintf "%i" var
        | GSLV_FLOAT(var) -> sprintf "%f" var
        | GSLV_STR(var) -> var.i    
        | GSLV_LIST(var) -> 
            seq {
                for v in var do
                    yield prettyPrintGSLVar v
            } |> fun x -> String.Join(";",Array.ofSeq x)


    // Main pretty print line start
    match line with
        | PRAGMA(p) -> sprintf "#%s %s" p.name (p.args |> String.concat " ") // pass pragma lines through unchanged
        | GSLLINEEXPANSION(text) -> text // Emit rewritten GSL Lines verbatim
        | OPENLINE(o) -> sprintf "open %s" o.modulePath
        | FORBLOCK(f) ->  // for x in [ 1 , 2, 3 ] do ....  end
            let formatItems (items:ForItemList) =
                sprintf "[%s]" (
                                
                                    match items with
                                        | INTLIST(il) -> String.Join (",",
                                                                il |> List.map (fun ii -> match ii with 
                                                                                            | INTMEMBER(i) -> string(i) 
                                                                                            | INTRANGE(s,e) -> sprintf "%d..%d" s e
                                                                               )
                                                         )
                                        | PARTLIST(pl) -> pl |> prettyPrintPartsList
                                )
            sprintf "for %s in %s do\n%s\nend" f.varName (formatItems f.items) (prettyPrintTree f.body)
        | CUTLine(cut) -> sprintf "cut %s" (prettyPrintLine cut.assembly)
        | ASSEMBLY(assembly) -> // Reemit assembly in standard form
            let docString = seq {
                                    for block in assembly.docStrings do
                                        for line in block do
                                            yield sprintf "/// %s" line.i
                                        yield ""
                                } |> fun x -> String.Join("\n",x)

            let a = prettyPrintPartsList assembly.parts
            docString + a
         | ROUGHAGESECTION(rList) ->
            // Expand out a roughage block which consists of individual roughage lines
            // Pretty print any roughage construct
            let sb = StringBuilder()
            sb.AppendLine("<@") |> ignore
            for rLine in rList do
                sb.AppendLine(prettyPrintRoughage rLine) |> ignore
            sb.AppendLine("@>") |> ignore
            sb.ToString()
        | L2(l2Line) ->
            // Level 2 GSL "line".  Different types potentially of L2 sections,
            // print them out accordingly
            let sb = StringBuilder()
            match l2Line.l2Design with
                | L2EXPLINE(l2ExpLine) ->
                    // Expression line  HO^ ; pA>gB etc
                    // Are there parts? - format them as a semicolon separated string
                    let partStringList = [ for part in l2ExpLine.parts -> 
                                            sprintf "%s>%s" part.promoter.String part.target.String ]
                    // Did they spec a locus?
                    match l2ExpLine.l2Locus with
                        | None -> ()
                        | Some(l) -> 
                                sb.Append(sprintf "%s^" l.String) |> ignore
                                // Might need a semi colon to separate locus and part expression string
                                if partStringList <> [] then sb.Append(";") |> ignore
                    let p = String.Join(";", partStringList)
                    sb.Append(p) |> ignore
            
            // Convert final L2 output back to a string
            sb.ToString()
        | DOCSTRINGSECTION(lines) ->
            // Docstring section - lines prefixed with ///, can re-emit them
            // verbatim without too much trouble
            let s  =lines |> List.fold (fun (sb:StringBuilder) line -> 
                                                    let x = sb.Append("///")
                                                    let y = x.Append("")
                                                    let y' = y.Append(line.i) // FIXFIX - this seems to defeat the purpose of reemitting
                                                    let z = y'.Append("\n")
                                                    z) (new StringBuilder())
            s.ToString()
        | GSLFUNCTION(f) ->
            sprintf "let %s(%s) =\n%s\nend" f.name (String.Join(",", (f.args |> List.map (fun p -> p)))) (prettyPrintTree f.body)
        | GSLFUNCTIONCALL(fc) -> 
            sprintf "%s(%s)" fc.name (String.Join(",", (fc.args |> List.map printPPP)))
        | LETLINE(l) -> sprintf "let %s=%s" l.varName (prettyPrintGSLVar l.letDefinition)
                                                                   
/// Pretty Print an assembly tree, emitting it verbatim
and prettyPrintTree (tree:GSLLine list) =
    tree |> List.map (prettyPrintLine) |> fun sl -> String.Join("\n",sl)+"\n"
            
