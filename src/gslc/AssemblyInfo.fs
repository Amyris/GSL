namespace System
open System.Reflection

[<assembly: AssemblyTitleAttribute("gslc")>]
[<assembly: AssemblyProductAttribute("gslc")>]
[<assembly: AssemblyDescriptionAttribute("Genotype Specification Language")>]
[<assembly: AssemblyVersionAttribute("0.2.2")>]
[<assembly: AssemblyFileVersionAttribute("0.2.2")>]
do ()

module internal AssemblyVersionInformation =
    let [<Literal>] Version = "0.2.2"
    let [<Literal>] InformationalVersion = "0.2.2"
