
"""
VintageOptimisation(;Results=Dict(),
                     UserParameters=Dict(),
                     ObjectiveIntegrand_Input = (Sconc,Sdist,Sagg,Cconc,Cdist,t::Float64,s::Float64,Para::Dict) -> [],
                     AggregationIntegrand_Input = (Sconc,Sdist,Cconc,Cdist, t::Float64,s::Float64, Para::Dict) -> [],
                     StateDynamic_conc_Input = (Sconc,Sagg,Cconc, t::Float64, Para::Dict) -> [], 
                     StateDynamic_dist_Input = (Sconc,Sdist,Sagg,Cconc,Cdist,t::Float64, s::Float64, Para::Dict) -> [],
                     InitDist_Input = (Sconc,Sagg,Cconc, t::Float64, Para::Dict) -> [],
                     SalvageFunction_Input = (Sconc,Sdist,t::Float64,s::Float64,Para) -> [])

This function solves the two-stage optimal control problem defined by the functions

- ObjectiveIntegrand_Input -> integrand of the objective function of time and vintage
- AggregationIntegrand_Input -> integrand of the aggregated variable definition
- StateDynamic_conc_Input -> State dynamics of concentrated variables
- StateDynamic_dist_Input -> State dynamics of distributed variables
- InitDist_Input -> function describing the initial values of the distributed variables
- SalvageFunction_Input -> Salvage function in the first stage

"""
function VintageOptimisation(;Results=Dict(),UserParameters=Dict(),
                    ObjectiveIntegrand_Input = (Sconc,Sdist,Sagg,Cconc,Cdist,t::Float64,s::Float64,Para::Dict) -> [],
                    AggregationIntegrand_Input = (Sconc,Sdist,Cconc,Cdist, t::Float64,s::Float64, Para::Dict) -> [],
                    StateDynamic_conc_Input = (Sconc,Sagg,Cconc, t::Float64, Para::Dict) -> [], 
                    StateDynamic_dist_Input = (Sconc,Sdist,Sagg,Cconc,Cdist,t::Float64, s::Float64, Para::Dict) -> [],
                    InitDist_Input = (Sconc,Sagg,Cconc, t::Float64, Para::Dict) -> [],
                    SalvageFunction_Input = (Sconc,Sdist,t::Float64,s::Float64,Para) -> [])
    
    
    # Set up dictionary with all parameters either from the base values or user-specific values
    Para = Dict()
    ParametersBasics(Para,UserParameters)
    ParametersGrids(Para,UserParameters)
    ParametersControls(Para,UserParameters)
    ParametersAlgorithm(Para,UserParameters)
    ParametersInitStates(Para,UserParameters)
    ParametersPlots(Para,UserParameters)
    for kk in keys(UserParameters)
        if !(kk in keys(Para))
            Para[kk] = UserParameters[kk]
        end
    end

    # Initialize the control, state and costate variables
    Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg = InitVariables(Para)
    
    #-----------------------------------------------------------------------------------------------
    # Check wether all the function supplied by the user have the right dimensions
    SystemCheck = 1

    if !(ObjectiveIntegrand_Input(Stat_conc[1,1,:],Stat_dist[1,1,:],Stat_agg[1,1,:],Con_conc[1,1,:],Con_dist[1,1,:],0.0,0.0,Para) isa Number)
        global ObjectiveIntegrand = (Sconc,Sdist,Sagg,Cconc,Cdist,t::Float64,s::Float64,Para::Dict) -> 0.0
        println("Warning: Objective Function is not a number! Objective Function is set to 0")
        SystemCheck = 0
    else
        global ObjectiveIntegrand = ObjectiveIntegrand_Input
    end
    if size(AggregationIntegrand_Input(Stat_conc[1,1,:],Stat_dist[1,1,:],Con_conc[1,1,:],Con_dist[1,1,:],0.0,0.0,Para),1) != Para["nStat_agg"]
        global AggregationIntegrand = (Sconc,Sdist,Cconc,Cdist, t::Float64,s::Float64, Para::Dict) -> zeros(Para["nStat_agg"])
        println("Warning: Aggregation Function is not a number! Aggregation Function is set to 0")
        SystemCheck = 0
    else
        global AggregationIntegrand = AggregationIntegrand_Input
    end
    if size(StateDynamic_conc_Input(Stat_conc[1,1,:],Stat_agg[1,1,:],Con_conc[1,1,:], 0.0, Para),1) != Para["nStat_conc"]
        global StateDynamic_conc = (Sconc,Sagg,Cconc, t::Float64, Para::Dict) -> zeros(Para["nStat_conc"])
        println("Warning: Dimensions of State Dynamics (Stage 1) do not match! Dynamics are set to 0")
        SystemCheck = 0
    else
        global StateDynamic_conc = StateDynamic_conc_Input
    end
    if size(StateDynamic_dist_Input(Stat_conc[1,1,:],Stat_dist[1,1,:],Stat_agg[1,1,:],Con_conc[1,1,:],Con_dist[1,1,:],0.0,0.0,Para),1) != Para["nStat_dist"]
        global StateDynamic_dist = (Sconc,Sdist,Sagg,Cconc,Cdist,t::Float64, s::Float64, Para::Dict) -> zeros(Para["nStat_dist"])
        println("Warning: Dimensions of State Dynamics (Stage 2) do not match! Dynamics are set to 0")
        SystemCheck = 0
    else
        global StateDynamic_dist = StateDynamic_dist_Input
    end
    if size(InitDist_Input(Stat_conc[1,1,:],Stat_agg[1,1,:],Con_conc[1,1,:],0.0,Para),1) != Para["nStat_dist"]
        global InitDist = (Sconc,Sagg,Cconc, t::Float64, Para::Dict) -> zeros(Para["nStat_dist"])
        println("Warning: Shock Transition dimensions do not math! Transition is set to 0!")
        SystemCheck = 0
    else
        global InitDist = InitDist_Input
    end
    if !(SalvageFunction_Input(Stat_conc[1,1,:],Stat_dist[1,1,:],0.0,0.0,Para) isa Number)
        global SalvageFunction = (Sconc,Sdist,t::Float64,s::Float64,Para::Dict) -> 0.0
        println("Warning: Salvage Function (Stage 1) is not a number! Salvage Function is set to 0")
        SystemCheck = 0
    else
        global SalvageFunction = SalvageFunction_Input
    end
    
    # Return error message if some dimensions do not match
    if SystemCheck == 0
        println("Calculations might not work due to dimension mismatches!")
    end
   

    #-----------------------------------------------------------------------------------------------
    # Define the order of the integration method for the objective value and 
    # the aggregated variables
    if Para["IntegrationOrder"] == 1
        global Integral = integ1
    elseif Para["IntegrationOrder"] == 2
        global Integral = integ2
    elseif Para["IntegrationOrder"] == 4
        global Integral = integ4
    else
        global Integral = integ4
    end

    #-----------------------------------------------------------------------------------------------
    # Load (interpolate) the supplied controls in the results dictionary
    if Para["LoadInits"] == true
        if all(in.(["Con_conc", "Con_dist"],(keys(Results),)))
            Con_conc, Con_dist = LoadVariables(Para,Results)
        else
            println("Warning: Initial Values cannot be loaded. Con_conc and/or Con_dist are not elemet of the dictionary.")
        end
    end

    #------------------------------------------------------------------------
    #   Start Optimisation loops
    Para["HstepIter"] = 1
    Step = Para["InitLineStep"]

    # Stop calculations if the current time-step size is smaller or equal to the lower bound specified
    if Para["hstep"] > Para["hLowBound"]
        StopOuterIteration = 0
    else
        StopOuterIteration = 1
        ParaAdjust(Para["hstep"],Para)
        # Interpolate controls after hstep reduction
        Con_conc, Con_dist = ConInterpol(Con_conc,Con_dist,Para)
    end
    while StopOuterIteration == 0

        # Adjust parameters to the new time-step size
        ParaAdjust(Para["hstep"],Para)
        
        # Interpolate controls after hstep reduction
        Con_conc, Con_dist = ConInterpol(Con_conc,Con_dist,Para)

        # Initialize the gradients and other auxiliary variables
        Stat_conc,Stat_dist,Stat_agg,dHam_conc,dHam_dist,CoStat_conc,CoStat_dist,CoStat_agg = InitVariables(Para)
        Stat_conc_new,Stat_dist_new,Stat_agg_new,Con_conc_new,Con_dist_new = InitVariables(Para)
        Stat_conc_best,Stat_dist_best,Stat_agg_best,Con_conc_best,Con_dist_best = InitVariables(Para)

        # Calculate the state variables related to the current control profiles
        state_PDE_solver(Stat_conc,Stat_dist, Stat_agg, Con_conc, Con_dist, Para)
        # Calculate the corresponding objective value
        ObjValue = ObjectValue(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,Para)
        # Calculate the profiles of the costate variables
        costate_PDE_solver(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,Para)
        # Calculate the gradient of the Hamiltonian for the current control values
        GradHamiltonian(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,dHam_conc,dHam_dist,Para)
        # Adjust gradient according to user specifications
        NewDirection(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,dHam_conc,dHam_dist,Para)
        #Assign the current profiles to the Results dictionary
        AssignResults(Results, Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,dHam_conc,dHam_dist, Para)
        
        if Para["PlotResultsStart"]
            # Plot all profiles
            PlotResults(Results)
            # Wait for user key if specified
            if Para["PlotResultsWaitForKey"]
                display("Wait for key")
                readline()
            end
        end


        impr = 1
        maxabsgradient = Inf
        while Para["OptiIter"] <= Para["MaxOptiIter"] && Step > Para["stepLowBound"] && impr==1 && maxabsgradient > Para["GradBound"]

            ObjValue = ObjectValue(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,Para)
            costate_PDE_solver(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,Para)
            GradHamiltonian(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,dHam_conc,dHam_dist,Para)
            maxabsgradient = maxabsGradient(Con_conc,Con_dist,dHam_conc,dHam_dist,Para)
            NewDirection(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,dHam_conc,dHam_dist,Para)

            # Linesearch along the gradient for a better set of controls
            ObjValue, Step, impr = LineSearch(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,dHam_conc,dHam_dist,ObjValue,Step,Para,
                                              Stat_conc_new, Stat_dist_new, Stat_agg_new, Con_conc_new, Con_dist_new, 
                                              Stat_conc_best, Stat_dist_best, Stat_agg_best, Con_conc_best, Con_dist_best)

            Para["OptiIter"] = Para["OptiIter"] + 1
            
            # Plot the current results every Para["PlotResultsIntermediateFrequency"] iterations
            if Para["PlotResultsIntermediate"] == true && Para["OptiIter"]%Para["PlotResultsIntermediateFrequency"] == 0
                GradHamiltonian(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,dHam_conc,dHam_dist,Para)
                NewDirection(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,dHam_conc,dHam_dist,Para)
                AssignResults(Results,Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,dHam_conc,dHam_dist,Para)
                PlotResults(Results)
                if Para["PlotResultsWaitForKey"]
                    display("Wait for key")
                    readline()
                end
            end
            
            # Smooth the profiles of the control variables if specified by the user
            if Para["OptiIter"]%50 == 0 && Para["ConSmooth"]
                AssignResults(Results, Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,dHam_conc,dHam_dist,Para)
                for ii = 1:10
                    ConMapping(Con_conc,Para)
                    ConMapping_dist(Con_dist,Para)
                    ConSmooth(Con_conc,Con_dist,Para)
                end
                AssignResults(Results,Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,dHam_conc,dHam_dist,Para)
            end
        end

        AssignResults(Results, Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,dHam_conc,dHam_dist,Para)
        if Para["PlotResultsIntermediate"]
            PlotResults(Results)
        end

        # Decrease/Half the time-step size and reset the iteration counters
        if Para["hstep"]/2 >= Para["hLowBound"]
            Para["hstep"] = Para["hstep"] * 0.5
            Para["HstepIter"] = Para["HstepIter"] + 1
            Para["MaxOptiIter"] = floor(Int,Para["MaxOptiIter"] * 2)
            Para["UpperLineStep"] = Para["UpperLineStep"] * 1
            Step = max(Step,Para["InitLineStep"])
            StopOuterIteration = 0
            # if the smallest step-size is reached, use the original gradient for the linesearch without adjustment
            if Para["hstep"]/2 <= Para["hLowBound"] 
                Para["OptiType"] = "Gradient"
            end
        else
            StopOuterIteration = 1
        end
    end
    
    #----------------------------------------------------------------------------------------------------------
    # Calculate all final variables
    Stat_conc,Stat_dist,Stat_agg,dHam_conc,dHam_dist,CoStat_conc,CoStat_dist,CoStat_agg = InitVariables(Para)
    state_PDE_solver(Stat_conc,Stat_dist, Stat_agg, Con_conc, Con_dist, Para)
    
  
    ObjValue = ObjectValue(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,Para)
    costate_PDE_solver(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,Para)
    GradHamiltonian(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,dHam_conc,dHam_dist,Para)
    Para["OptiType"] = "Gradient"
    NewDirection(Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,dHam_conc,dHam_dist,Para)
    AssignResults(Results, Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg,dHam_conc,dHam_dist, Para)
    if Para["PlotResultsFinal"] == true
        PlotResults(Results;SavePlot = Para["SavePlot"])
    end
    return Results
end
