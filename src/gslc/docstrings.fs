/// Support for /// style docstring quotes.
module docstrings

open shared
open constants
open commonTypes
open parseTypes
open System.IO

/// Take a mixture of docstrings and
/// assembly lines and fold the matching document strings into the assemblies for use later
/// 
let stuffDocstrings _(*verbose*) (linesIn : GSLLine list) =
    let rec procOneLine (docStrings : DocstringLine list list) (lines:GSLLine list) (res: GSLLine list) =
        match lines with
            | [] ->  List.rev res // Done, flip LIFO list to restore line order
            | ASSEMBLY(a)::tl ->
                procOneLine [] tl (ASSEMBLY({ a with docStrings = docStrings})::res) // Reset list of current docstrings, push out assembly
            | DOCSTRINGSECTION(lines)::tl -> // Add to current set of docstrings
                procOneLine (lines::docStrings) tl res
            | hd::tl -> // No packing of docstrings for these cases..
                procOneLine docStrings tl (hd::res)
                
                
    procOneLine [] linesIn []
     

 /// Create output file with user or algorithm documentation of the designs
let dumpDocStrings (path:string) (assemblies:AssemblyOut list) =
    use outF = new StreamWriter(path)
    for a in assemblies do
        outF.WriteLine(sprintf "@name=%s" a.name)
        for docBlock in a.docStrings do
            for line in docBlock do
                outF.WriteLine(line.i)
            outF.WriteLine("")


