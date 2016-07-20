module pragmaProcess
open pragmaTypes
open shared
open constants
open commonTypes
open parseTypes
open pcrParamsParse
open Amyris.primercore // For definition of unit C

/// Starting conditions for certain pragma states
let initialPragmas = EmptyPragmas
/// Eliminate pragmas that should only affect the next assembly, but not persist
let cleanTransientPragmas (pragmas:PragmaCollection) =
    pragmas.pmap
    |> Map.filter (fun _ p -> not p.isTransient)
    |> PragmaCollection

let getTransientPragmas (pragmas:PragmaCollection) =
    pragmas.pmap
    |> Map.filter (fun _ p -> p.isTransient)
    |> PragmaCollection

/// Take a mixture of pragma instructions and
/// assembly lines and fold the appropriate pragma instructions
/// into the assembly details
let stuffPragmas _(*verbose*) (legalCapa:string->bool) (initialDesignParams:DesignParams) (linesIn : GSLLine list) =
    /// Processes the first line in lines, either by updating the environment (pragma, design parameters, capabilities,
    /// pragma stack (for push/pop pragmas), and pushed transient pragmas (for push/pop pragmas)) or by folding the 
    /// current environment into the assembly/L2 instructions.
    ///
    /// designParams::DesignParams storing the current design parameters (e.g. pcrparams)
    /// capabilities::Capabilities storing the current capabilities available for engineering
    /// pragmas::a Map object storing the current pragmas. key is the pragma name and value is the associated value
    /// pragmaStack:: list of Map object serving as a stack for push/pop pragmas. Pragmas stored in this list does not 
    ///               include transient pragmas. 
    /// pushedTransientPragmas:: Map object storing transient pragmas that were present in a previous push instruction
    ///                          that has not been applied yet. When pushedTransientPragmas is non-empty, an assembly line
    ///                          or a L2 line should use and clear those transient pragmas, and a #pop line should make
    ///                          those transients a part of the current pragma environment.
    /// lines:: GSL lines to be processed
    /// res:: GSL lines processed
    let rec procLine 
            (legalCapa:string->bool)
            (designParams:DesignParams)
            (capabilities:Capabilities) 
            (pragmas:PragmaCollection) 
            (pragmaStack:List<PragmaCollection>) 
            (pushedTransientPragmas:PragmaCollection) 
            (lines:GSLLine list) 
            res =
        match lines with
        | [] -> List.rev res
        | PRAGMA(p)::tl ->

            // replace or add pragma details
            if p.name = "warn" then
                printf "WARNING: %s\n" (p.args |> String.concat " ")
            // Updated pragma state
            let pragmas' =
                if p.name = "warnoff" then
                    // stuff warnoff entries into one
                    let newWarnOff =
                        match pragmas.TryFind("warnoff") with
                        | Some(warnoff) -> {warnoff with args = warnoff.args@p.args}
                        | None -> p
                    pragmas.Add(newWarnOff)
                elif p.name = "push" then
                    // when #push pragma is encountered, transient pragmas are separated from other pragmas.
                    // transient pragmas will be stored in pushedTransientPragmas while the rest will go on
                    // the stack and in the current pragma environment
                    cleanTransientPragmas pragmas
                elif p.name = "warn" then
                    pragmas
                elif p.name = "pop" then
                    match pragmaStack with
                    // We udate pragmaStack later, so ignore rest of stack for now
                    | mostRecentPragmaState::_ ->
                        // transients in the current pragma environment
                        let currTransients = getTransientPragmas pragmas
                        // if pushedTransientpragmas is empty, then the new pragma environment should be
                        // the combination of current transient pragmas and the pragma environment at the
                        // top of the stack.
                        if pushedTransientPragmas.pmap.IsEmpty then
                            mostRecentPragmaState.MergeIn(currTransients)
                        // if pushedTransientPragmas is not empty and there is no current transient pragmas,
                        // then the new pragma environment should be the combination of the pushedTransientPragmas
                        // and the pragma environment at the top of the stack
                        elif currTransients.pmap.IsEmpty then 
                            mostRecentPragmaState.MergeIn(pushedTransientPragmas)
                        // if the pushedTransientPragmas and the currTransients are both non-empty, then it is 
                        // unclear which transients should be used. So, it gives users a warning and sets the
                        // combination of the current transient pragmas and pragmas at the top of the stack as 
                        // the new pragma environment
                        else
                            printf
                                "WARNING: conflicting transient pragmas: %A, %A"
                                currTransients pushedTransientPragmas
                            mostRecentPragmaState.MergeIn(currTransients)
                    | [] ->
                        failwith
                            "ERROR: #pop operation cannot be done because there is no pushed state"
                else
                    // Normal pragma that replaces existing value
                    pragmas.Add(p)

            let pragmas'' =
                if p.name = "megastitch" then
                    PragmaCollection(pragmas'.pmap.Remove("stitch"))
                else pragmas'
            
            let pushedTransientPragmas' = 
                if p.name = "push" then
                    // pragmas in the current pragma environment
                    let currTransientPragmas = getTransientPragmas pragmas
                    // if there are already pragmas in pushedTransientPragmas and there are transient
                    // pragmas in the current environment, then it is unclear which transient pragmas
                    // should be stored (there should be only one set of transient pragmas at any given
                    // time because transients only apply to one line of assembly/L2 GSL). Takes the
                    // transient pragmas in the current environment to break the tie. 
                    if not (pushedTransientPragmas.pmap.IsEmpty || currTransientPragmas.pmap.IsEmpty) then 
                        printf
                            "WARNING: colliding transient pragmas %A and %A"
                            pushedTransientPragmas currTransientPragmas
                        currTransientPragmas
                    // if currTransientPragmas and/or pushedTransientPragmas is/are empty, then takes the
                    // combination of the two, which should be either one that is non-empty or empty if both
                    // are empty
                    else pushedTransientPragmas.MergeIn(currTransientPragmas)
                // if we are dealing with the #pop pragma line, there should no longer be pragmas in pushedTransientPragmas
                // because the #pop pragmas puts the pragmas in pushedTransientPragmas in the new pragmas environment
                elif p.name = "pop" then EmptyPragmas
                else pushedTransientPragmas

            let pragmaStack' = 
                // puts the current pragmas environment (minus the transients) onto the stack
                if p.name = "push" then pragmas''::pragmaStack
                // removes the pragmas environment on the top from the stack
                elif p.name = "pop" then
                    match pragmaStack with
                    |_::restPragmaHistory -> restPragmaHistory
                    |[] ->
                        failwith
                            "ERROR: #pop operation cannot be done because there is no pushed state"
                else pragmaStack

            // Update the pcr conditions if the user specifies a #pcrparams tag
            // or the targettm for seamless and regular
            let designParams' =
                match p.name, p.args with
                | "pcrparams", args ->
                    // TODO: refactor these first two to pass list of parsed args
                    {designParams with
                        pp = pcrParamsParse.parse designParams.pp (String.concat " " args)}
                | "pcrassemblyparams", args ->
                    {designParams with
                        overlapParams =
                            pcrParamsParse.parse
                                designParams.overlapParams (String.concat " " args) }
                | "targettm", v::_ ->
                    {designParams with targetTm = strToTempC v}
                | "minoverlaplen", v::_ ->
                    {designParams with overlapMinLen = int v}
                | "seamlesstm", v::_ ->
                    {designParams with seamlessTm = strToTempC v}
                | "seamlessoverlaptm", v::_ ->
                    {designParams with seamlessOverlapTm = strToTempC v}
                | "atpenalty", v::_ ->
                    {designParams with
                        pp = {designParams.pp with ATPenalty = float v * 1.0<C>}}
                | _ -> designParams

            let capabilities' =
                match p.name with
                | "capa" -> 
                    let capa = p.args.[0].Trim().ToLower()
                    if legalCapa capa then
                        capabilities.Add(capa)
                    else
                        failwithf "ERROR: unknown capability '%s'" capa
                | _ -> capabilities
            // lines kept for the next iteration in the compiler
            // keep all lines except the #warn line because the compiler should
            // only warn once
            let res' = 
                match p.name with
                |"warn" -> res
                |_ -> PRAGMA(p)::res

            procLine legalCapa designParams' capabilities' pragmas'' pragmaStack' pushedTransientPragmas' tl res'
        
        // No pragmas stuffing into GSL Lines that are text expanded
        | GSLLINEEXPANSION(s)::tl -> 
            procLine legalCapa designParams capabilities pragmas  pragmaStack pushedTransientPragmas tl (GSLLINEEXPANSION(s)::res)

        // No pragma stuffing into docstrings
        | DOCSTRINGSECTION(lines)::tl ->
            procLine legalCapa designParams capabilities pragmas  pragmaStack pushedTransientPragmas tl (DOCSTRINGSECTION(lines)::res)
        | ASSEMBLY(a)::tl ->
            let currTransients = getTransientPragmas pragmas
            let pragmas' = 
                if currTransients.pmap.IsEmpty then
                    // if current transients are empty, then the new pragmas environment should include pragmas
                    // in pushedTransientPragmas.
                    pragmas.MergeIn(pushedTransientPragmas) 
                else
                    // if both currTransients and pushedTransientPragmas are non-empty, then it is unclear which
                    // transient pragmas should be used. Takes the current transients to break the tie
                    if not pushedTransientPragmas.pmap.IsEmpty then 
                        printf "WARNING: colliding transient pragmas %A and %A"
                            pushedTransientPragmas currTransients
                    pragmas
            
            // Fold current pragma state and capability state into the assembly
            // and current design parameters
            
            let a' =
                {a with
                    name = pragmas'.TryGetOne("name");
                    uri = pragmas'.TryGetOne("uri");
                    linkerHint =
                        match pragmas'.TryGetValues("linkers") with
                        | Some(vals) -> (String.concat "" vals)
                        | None -> "";
                    pragmas = pragmas'.MergeIn(a.pragmas);
                    designParams = designParams;
                    capabilities = capabilities}

            // Eliminate pragmas that should only affect the net assembly, but not persist
            let pragmas'' = cleanTransientPragmas pragmas'
            procLine legalCapa designParams capabilities pragmas'' pragmaStack EmptyPragmas tl (ASSEMBLY(a')::res)

        | L2(l2Line)::tl ->
            let currTransients = getTransientPragmas pragmas
            let pragmas' = 
                if currTransients.pmap.IsEmpty then
                    // if current transients are empty, then the new pragmas environment should include pragmas
                    // in pushedTransientPragmas.
                    pragmas.MergeIn(pushedTransientPragmas)
                else
                    // if both currTransients and pushedTransientPragmas are non-empty, then it is unclear which
                    // transient pragmas should be used. Takes the current transients to break the tie
                    if not pushedTransientPragmas.pmap.IsEmpty then 
                        printf "WARNING: colliding transient pragmas %A and %A"
                            pushedTransientPragmas currTransients
                    pragmas
            
            // Fold current pragma state and capability state into the assembly
            // and current design parameters, as well as name and URI
            let l2' =
                {l2Line with
                    name = pragmas'.TryGetOne("name");
                    uri = pragmas'.TryGetOne("uri");
                    pragmas = pragmas'.MergeIn(l2Line.pragmas);
                    capabilities = capabilities}

            // Eliminate pragmas that should only affect the next assembly, but not persist
            let pragmas'' = pragmas' |> cleanTransientPragmas
            procLine legalCapa designParams capabilities pragmas'' pragmaStack EmptyPragmas tl (L2(l2')::res)

        | ROUGHAGESECTION(rs)::tl -> // Just pass this through, no interaction with pragmas or environment
            procLine legalCapa designParams capabilities pragmas pragmaStack pushedTransientPragmas tl (ROUGHAGESECTION(rs)::res)
        | CUTLine(c)::tl -> 
            procLine legalCapa designParams capabilities (cleanTransientPragmas pragmas) pragmaStack pushedTransientPragmas tl (CUTLine( { c with pragmas=pragmas} )::res)
        | GSLFUNCTION(c)::tl -> 
            procLine legalCapa designParams capabilities (cleanTransientPragmas pragmas) pragmaStack pushedTransientPragmas tl (GSLFUNCTION( { c with pragmas=pragmas} )::res)
        | GSLFUNCTIONCALL(c)::tl -> // No pragma stuffing for function call
            procLine legalCapa designParams capabilities (cleanTransientPragmas pragmas) pragmaStack pushedTransientPragmas tl (GSLFUNCTIONCALL( c )::res)
        | LETLINE(x)::tl -> // No pragma stuffing for let line ...?
            procLine legalCapa designParams capabilities (*(cleanTransientPragmas pragmas)*) pragmas pragmaStack pushedTransientPragmas tl (LETLINE(x)::res)
        | FORBLOCK(f)::tl -> 
            procLine legalCapa designParams capabilities (cleanTransientPragmas pragmas) pragmaStack pushedTransientPragmas tl (FORBLOCK( { f with pragmas=pragmas} )::res)
        | OPENLINE(o)::tl -> 
            procLine legalCapa designParams capabilities (cleanTransientPragmas pragmas) pragmaStack pushedTransientPragmas tl (OPENLINE( { o with pragmas=pragmas} )::res)
    
    procLine legalCapa initialDesignParams Set.empty initialPragmas [] EmptyPragmas linesIn []