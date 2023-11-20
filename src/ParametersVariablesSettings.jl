
"""
    ParametersBasics(Para,UserPara)

Initialize the basic parameters necessary for the model

 - T -> time horizon of the model
 - ω -> maximum vintage age        
 - hstep -> time step between points
 - TimeDiscountRate -> time discount rate of the model
 - nCon_conc  -> number of concentrated controls
 - nCon_dist  -> number of distributed controls
 - nStat_conc -> number of concentrated state variables
 - nStat_dist -> number of distributed state variables
 - nStat_agg  -> number of aggregated state variables
"""
function ParametersBasics(Para,UserPara)
    #= --------------------------------------------------------------
        Parameters for structure/dimensions of the model
    -------------------------------------------------------------- =#
    
    Para2 = Dict()

    # time-horizont
    Para2["T"] = 1.0                       
    Para2["ω"] = 1.0                       
    
    # stepsize
    Para2["hstep"] = 0.05                     

    # time discount rate
    Para2["TimeDiscountRate"] = 0.0

    # number of concentrated controls
    Para2["nCon_conc"] = 0                   
    # number of distributed controls
    Para2["nCon_dist"] = 0         
    # number of concentrated state varibales
    Para2["nStat_conc"] = 0           
    # number of distributed state varibales
    Para2["nStat_dist"] = 0 
    # number of aggregated state varibales
    Para2["nStat_agg"] = 0 

    # Overwrite values if given by the User
    for kk in keys(Para2)
        if kk in keys(UserPara)
            Para[kk] = UserPara[kk]
        else
            Para[kk] = Para2[kk]
        end
    end
end

"""
    ParametersGrids(Para,UserPara)

Initialize the parameters associated with the time grids

- nTime -> Number of gridpoints along time
- nVintage  -> Number of gridpoints along the vintages
- tmesh -> Grid for time
- amesh -> Grid for vintages
"""
function ParametersGrids(Para,UserPara)
    Para2 = Dict()
    # number of gridpoints
    Para2["nTime"] = floor(Int,Para["T"]/Para["hstep"]) + 1 
    Para2["nVintage"] = floor(Int,Para["ω"]/Para["hstep"]) + 1 

    # grid over the time horizon
    Para2["tmesh"] = collect(0:Para["hstep"]:Para["T"])     
    Para2["amesh"] = collect(0:Para["hstep"]:Para["ω"]) 

    # Overwrite values if given by the User
    for kk in keys(Para2)
        if kk in keys(UserPara)
            Para[kk] = UserPara[kk]
        else
            Para[kk] = Para2[kk]
        end
    end
end

"""
    ParametersControls(Para,UserPara)

Initialize the parameters associated with the control variables

- LoadInits   -> Loading specific initial guesses for the control variables
- Reduce1     -> Indices of controls, which should be disabled (Stage 1)
- Reduce2     -> Indices of controls, which should be disabled (Stage 2)
- ConMin      -> Lower bound for the concentated control variables
- ConMax      -> Upper bound for the concentated control variables
- Con_distMin -> Lower bound for the distributed control variables
- Con_distMax -> Upper bound for the distributed control variables
"""
function ParametersControls(Para,UserPara)
    Para2 = Dict()
    # Loading specific initial guesses for the control variables
    Para2["LoadInits"] = false                   
    # Indices of concentrate controls, which should be disabled
    Para2["ReduceConc"] = []                       
    # Indices of distributed controls, which should be disabled
    Para2["ReduceDist"] = []       
    # Lower and Upper bounds for the control variables
    Para2["Con_concMin"] = -Inf*ones(Para["nCon_conc"])
    Para2["Con_concMax"] =  Inf*ones(Para["nCon_conc"])
    Para2["Con_distMin"] = -Inf*ones(Para["nCon_dist"])
    Para2["Con_distMax"] =  Inf*ones(Para["nCon_dist"])

    # Overwrite values if given by the User
    for kk in keys(Para2)
        if kk in keys(UserPara)
            Para[kk] = UserPara[kk]
        else
            Para[kk] = Para2[kk]
        end
    end
end

"""
    ParametersAlgorithm(Para,UserPara)

Initialize the parameters for the optimisation algorithm

* OptiType -> Adjustment of gradient of Hamiltonian
    - "Gradient" uses the standard gradient of the hamiltonian 
    - "ProbAdjust" scales the gradient by the inverse of the auxiliary variables
    - "Newton-Raphson" scales the gradient by the inverse of the Hessian
* IntegrationOrder -> Which order should the method of integration used should have.
    - 4 is the default using Simpson's rule
    - 2 is another option using the trapezoid rule
* ParallelComputing   -> Indicator whether calculations along vintages should be parallised. 
                        Needs to be implemented in the future.
* hLowBound           -> Lower bound for the step of the gridsize along time.
                        hstep gets halfed until the new value falls below hLowBound.
* InitLineStep        -> Initial step-size for the adjustment of control variables along the gradient.
* LineSearchStepIncrease -> Factor for increase in the linesearch stepsize.
* UpperLineStep       -> Upper bound for linesearch stepsize.
* stepLowBound        -> Lower bound for step in gradient search -> termination if undercut.
* GradBound           -> Lower bound for the gradient. 
                        Algorithm stops if maximum of the absolute value of the gradient is smaller.
* MaxOptiIter         -> Max iterations with same step-size h.
* MaxLineIter         -> Max iterations along the same gradient. 
* GradSmooth          -> Indicator if gradients are smoothed.
* ConSmooth           -> Indicator if controls are smoothed. 
* SearchDISP          -> Indicator if intermediate results of the search along the gradient shall be displayed.
* LineIter            -> Iteration counter
* OptiIter            -> Iteration counter
* HstepIter           -> Iteration counter
"""
function ParametersAlgorithm(Para,UserPara)
    #= ------------------------------------------------------------
        Parameters for the optimisation algorithm
    ------------------------------------------------------------ =#
    Para2 = Dict()    
    # Adjustment of gradient of Hamiltonian
    #   - "Gradient" uses the standard gradient of the hamiltonian 
    #   - "Newton-Raphson" scales the gradient by the inverse of the Hessian
    Para2["OptiType"] = "Newton-Raphson" 

    Para2["IntegrationOrder"] = 4
    Para2["ParallelComputing"] = false
    
    # Lower bound for the step-size
    Para2["hLowBound"] = Para["hstep"]*0.9 

    # Initial stepsize for linesearch
    Para2["InitLineStep"] = 1e-6          
    # Factor for increase in the linesearch stepsize
    Para2["LineSearchStepIncrease"] = 0.5
    # Upper bound for linesearch stepsize
    Para2["UpperLineStep"] = 1e1        
    # Upper bound for linesearch stepsize in one iteration. Is adjusted every step
    Para2["UpperLineStepTemporary"] = 1e1
    # Factor for the upper bound for linesearch stepsize in one iteration.
    Para2["UpperLineStepTemporaryFactor"] = 1e2
    # Lower bound for step in gradient search -> termination if undercut
    Para2["stepLowBound"] = 1e-8   
    # Lower bound for the gradient -> termination if undercut
    Para2["GradBound"] = 1e-6
    # Max iterations with same step-size h
    Para2["MaxOptiIter"] = 10000                   
    # Max iterations along the same gradient 
    Para2["MaxLineIter"] = 30                      

    # Indicator if gradients are smoothed
    Para2["GradSmooth"] = false              
    # Indicator if controls are smoothed  
    Para2["ConSmooth"] = false               

    # Indicator if intermediate results shall be displayed
    Para2["SearchDISP"] = true                 

    # Iteration counters
    Para2["LineIter"] = 0                      
    Para2["OptiIter"] = 0
    Para2["HstepIter"] = 0

    Para["GradientSteps"] = 1e-6
    Para["GradientStepsActive"] = true

    for kk in keys(Para2)
        if kk in keys(UserPara)
            Para[kk] = UserPara[kk]
        else
            Para[kk] = Para2[kk]
        end
    end
end


"""
    ParametersInitStates(Para,UserPara)

Initialize initial and boundary data

 - InitStat        -> Initial state values of concentated variables
 - InitStat_dist   -> Initial state values of distributed variables at t=t_0
 - BoundStat_dist  -> Initial state values of distributed variables at t
                    Will be used in future generalised version of the package
"""
function ParametersInitStates(Para,UserPara)
    #= ------------------------------------------------------------
        Initial and boundary data
    ------------------------------------------------------------ =#
    Para2 = Dict()
    Para2["InitStat_conc"] = zeros(1,1,Para["nStat_conc"])
    Para2["InitStat_dist"] = zeros(Para["nVintage"],1,Para["nStat_dist"])
    Para2["BoundStat_dist"] = zeros(1,Para["nTime"],Para["nStat_dist"])

    for kk in keys(Para2)
        if kk in keys(UserPara)
            Para[kk] = UserPara[kk]
        else
            Para[kk] = Para2[kk]
        end
    end
end

"""
    ParametersPlots(Para,UserPara)

Initialize parameters relevant for plotting

 - nVintagePlot            -> Number of vintages plotted for realtime solutions
 - PlotResultsStart        -> Indicator for plot of initial guess
 - PlotResultsIntermediate -> Indicator for plot of intermediate results
 - PlotResultsFinal        -> Indicator for plot of final results
 - PlotResultsIntermediateFrequency -> Number of algorithm iterations 
                                    between plot of Intermediate results
 - PlotResultsWaitForKey   -> Indicator whether user input of any key should be necessary 
                            after plot of intermediate results
 - SavePlot                -> Indicator whether plotted results should be saved
 - SavePlotPath            -> Path for saved plots

 - ControlConcLabelsLatex  -> concentated control labels using latex strings 
 - ControlConcLabelsSimple -> concentated control labels without latex strings
 - ControlDistLabelsLatex  -> distributed control labels using latex strings 
 - ControlDistLabelsSimple -> distributed control labels without latex strings 
 - StateConcLabelsLatex    -> concentated state labels using latex strings 
 - StateConcLabelsSimple   -> concentated state labels without latex strings 
 - StateDistLabelsLatex    -> distributed state labels using latex strings 
 - StateDistLabelsSimple   -> distributed state labels without latex strings 
 - StateAggLabelsLatex     -> aggregated state labels using latex strings 
 - StateAggLabelsSimple    -> aggregated state labels without latex strings 
 - CoStateConcLabelsLatex  -> concentated costate labels using latex strings 
 - CoStateConcLabelsSimple -> concentated costate labels without latex strings 
 - CoStateDistLabelsLatex  -> distributed costate labels using latex strings 
 - CoStateDistLabelsSimple -> distributed costate labels without latex strings
 - CoStateAggLabelsLatex   -> aggregated costate labels using latex strings 
 - CoStateAggLabelsSimple  -> aggregated costate labels without latex strings 
"""
function ParametersPlots(Para,UserPara)
    Para2 = Dict()
    #= ------------------------------------------------------------
        Parameters for Plots
    ------------------------------------------------------------ =#
    Para2["ControlConcLabelsLatex"] = [latexstring("u_",kk) for kk=1:Para["nCon_conc"]]
    Para2["ControlConcLabelsSimple"] = [string("u",kk) for kk=1:Para["nCon_conc"]]
    Para2["ControlDistLabelsLatex"] = [latexstring("v_",kk) for kk=1:Para["nCon_dist"]]
    Para2["ControlDistLabelsSimple"] = [string("v",kk) for kk=1:Para["nCon_dist"]]
    Para2["StateConcLabelsLatex"] = [latexstring("x_",kk) for kk=1:Para["nStat_conc"]]
    Para2["StateConcLabelsSimple"] = [string("x",kk) for kk=1:Para["nStat_conc"]]
    Para2["StateDistLabelsLatex"] = [latexstring("y_",kk) for kk=1:Para["nStat_dist"]]
    Para2["StateDistLabelsSimple"] = [string("y",kk) for kk=1:Para["nStat_dist"]]
    Para2["StateAggLabelsLatex"] = [latexstring("Q_",kk) for kk=1:Para["nStat_agg"]]
    Para2["StateAggLabelsSimple"] = [string("Q",kk) for kk=1:Para["nStat_agg"]]
    Para2["CoStateConcLabelsLatex"] = [latexstring("\\lambda_",kk) for kk=1:Para["nStat_conc"]]
    Para2["CoStateConcLabelsSimple"] = [string("lambda",kk) for kk=1:Para["nStat_conc"]]
    Para2["CoStateDistLabelsLatex"] = [latexstring("\\xi_",kk) for kk=1:Para["nStat_dist"]]
    Para2["CoStateDistLabelsSimple"] = [string("xi_",kk) for kk=1:Para["nStat_dist"]]
    Para2["CoStateAggLabelsLatex"] = [latexstring("\\nu_",kk) for kk=1:Para["nStat_agg"]]
    Para2["CoStateAggLabelsSimple"] = [string("nu_",kk) for kk=1:Para["nStat_agg"]]

    Para2["nVintagePlot"] = 10                      # Number of vintages plotted for realtime solutions
    Para2["PlotResultsStart"] = true
    Para2["PlotResultsIntermediate"] = true             # Indicator if intermediate results are plotted
    Para2["PlotResultsFinal"] = true
    Para2["PlotResultsIntermediateFrequency"] = 40
    Para2["PlotResultsWaitForKey"] = false
    Para2["SavePlot"] = false                           # Indicator if plotted results should be exported and saved
    Para2["SavePlotPath"] = "Results/Plots/"            # Path for saved plots

    for kk in keys(Para2)
        if kk in keys(UserPara)
            Para[kk] = UserPara[kk]
        else
            Para[kk] = Para2[kk]
        end
    end
end


"""
    InitVariables(Para::Dict)

Initialize control, state and costate variables for the parameters specified
"""
function InitVariables(Para::Dict)
    Stat_conc = zeros(1,               Para["nTime"],Para["nStat_conc"])
    Stat_dist = zeros(Para["nVintage"],Para["nTime"],Para["nStat_dist"])
    Stat_agg =  zeros(1,               Para["nTime"],Para["nStat_agg"])
    Con_conc =  zeros(1,               Para["nTime"],Para["nCon_conc"])
    Con_dist =  zeros(Para["nVintage"],Para["nTime"],Para["nCon_dist"])
    CoStat_conc = zeros(1,             Para["nTime"],Para["nStat_conc"])
    CoStat_dist = zeros(Para["nVintage"],Para["nTime"],Para["nStat_dist"])
    CoStat_agg = zeros(1,              Para["nTime"],Para["nStat_agg"])

    
    # Set initial guess for controls as the mean of the upper and lower bound
    for kk = 1:Para["nCon_conc"]
        if max(-Para["Con_concMin"][kk],Para["Con_concMax"][kk]) < Inf
            Con_conc[:,:,kk] .= (Para["Con_concMin"][kk] + Para["Con_concMax"][kk])/2
        elseif -Para["Con_concMin"][kk] < Inf
            Con_conc[:,:,kk] .= Para["Con_concMin"][kk] .+ 1.0
        elseif Para["Con_concMax"][kk] < Inf
            Con_conc[:,:,kk] .= Para["Con_concMax"][kk] .- 1.0
        end 
    end
    
    for kk = 1:Para["nCon_dist"]
        if max(-Para["Con_distMin"][kk],Para["Con_distMax"][kk]) < Inf
            Con_dist[:,:,kk] .= (Para["Con_distMin"][kk] + Para["Con_distMax"][kk])/2
        elseif -Para["Con_distMin"][kk] < Inf
            Con_dist[:,:,kk] .= Para["Con_distMin"][kk] .+ 1.0
        elseif Para["Con_distMax"][kk] < Inf
            Con_dist[:,:,kk] .= Para["Con_distMax"][kk] .- 1.0
        end
    end


    for kk in Para["ReduceConc"]
        Con_conc[:,:,kk] .= 0
    end
    for kk in Para["ReduceDist"]
        Con_dist[:,:,kk] .= 0
    end

    return Stat_conc,Stat_dist,Stat_agg,Con_conc,Con_dist,CoStat_conc,CoStat_dist,CoStat_agg
end


