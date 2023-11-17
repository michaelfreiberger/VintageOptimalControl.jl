
# ObjectiveIntegrand_Input(Sconc,Sdist,Sagg,Cconc,Cdist,t,s,Para) = Para["Utility"](Sagg[1]*1/Para["ω"]*(1-Cconc[1]))
# AggregationIntegrand_Input(Sconc,Sdist,Cconc,Cdist, t::Float64,s::Float64, Para::Dict) = [(Sdist[1]+1e-2)^Para["α"]*Para["A"](s)]
# StateDynamic_conc_Input(Sconc,Sagg,Cconc, t::Float64, Para::Dict) = [] 
# StateDynamic_dist_Input(Sconc,Sdist,Sagg,Cconc,Cdist,t::Float64, s::Float64, Para::Dict) = [-Para["δ"]*Sdist[1] + Sagg[1]*1/Para["ω"]*Cconc[1]]
# InitDist_Input(Sconc,Sagg,Cconc, t::Float64, Para::Dict) = [0.0]
# SalvageFunction_Input(Sconc,Sdist,t::Float64,s::Float64,Para) = 0.0

ObjectiveIntegrand_Input(Sconc,Sdist,Sagg,Cconc,Cdist,t,s,Para) = 1/Para["ω"] * (Para["R"](Sagg[1]) - Cconc[1] - 1/2*Cconc[1]^2) - Para["b"](s)*Cdist[1] - Para["c"](s)/2*Cdist[1]^2
AggregationIntegrand_Input(Sconc,Sdist,Cconc,Cdist, t::Float64,s::Float64, Para::Dict) = [Para["f"](t,s)*Para["v"](s)*Sdist[1]]
StateDynamic_conc_Input(Sconc,Sagg,Cconc, t::Float64, Para::Dict) = [] 
StateDynamic_dist_Input(Sconc,Sdist,Sagg,Cconc,Cdist,t::Float64, s::Float64, Para::Dict) = [-Para["delta"](s)*Sdist[1] + Cdist[1]]
InitDist_Input(Sconc,Sagg,Cconc, t::Float64, Para::Dict) = [Cconc[1]]
SalvageFunction_Input(Sconc,Sdist,t::Float64,s::Float64,Para) = 0.0

MyPara = Dict()
MyPara["OptiType"] = "Newton-Raphson"
#MyPara["OptiType"] = "Gradient"

MyPara["T"] = 100
MyPara["ω"] = 10
MyPara["hstep"] = 0.5

MyPara["nCon_conc"] = 1
MyPara["nCon_dist"] = 1
MyPara["nStat_conc"] = 0
MyPara["nStat_dist"] = 1
MyPara["nStat_agg"] = 1

MyPara["InitLineStep"] = 1e-5
#MyPara["LineSearchStepIncrease"] = 0.25
#MyPara["UpperLineStep"] = 1e-1
#MyPara["UpperLineStepTemporaryFactor"] = 1.5
#MyPara["hLowBound"] = 0.7
MyPara["PlotResultsIntermediateFrequency"] = 50
MyPara["PlotResultsWaitForKey"] = false

MyPara["InitStat_dist"] = [1.0,1.0]
MyPara["Con_concMin"] = [0.0]
MyPara["Con_concMax"] = [0.0]
MyPara["Con_distMin"] = [0.0]
MyPara["Con_distMax"] = [Inf]

MyPara["ConSmooth"] = false

# MyPara["Utility"] = (c) -> ((c+1e-6)^(1-1.1) - 1)/(1-1.1)
# MyPara["A"] = (s) -> 1.0 - 0.5*(s-MyPara["ω"]/2)^2/(MyPara["ω"]/2)^2
# MyPara["α"] = 0.3
# MyPara["δ"] = 0.1
# MyPara["ρ"] = 0.1
MyPara["TimeDiscountRate"] = 0.03


MyPara["R"] = (Q) -> Q - 8e-3/2*Q^2
MyPara["f"] = (t,a) -> 3.0*(t-a <= 100) + 4.0*(t-a > 100)#1 + (t-a+10)/6*log10(2)
MyPara["delta"] = (a) -> 2*a/10^2*log(5)
MyPara["v"] = (a) -> (2*sqrt(a/5)-a/5)*(a<5) + 1.0*(a>=5.0)
MyPara["b"] = (a) -> 0.8*(10-a)/10*4.931
MyPara["c"] = (a) -> 1.2*exp(-a/10)

MyPara["LoadInits"] = true
Results = Dict()
Results["Con_conc"] = 7.0*ones(1,20,MyPara["nCon_conc"])
Results["Con_dist"] = 0.4*ones(10,20,MyPara["nCon_dist"])


UserParameters = MyPara
 
#using BenchmarkTools, Profile, Coverage


Results = VintageOptimisation(Results = Results,UserParameters = MyPara,
                            ObjectiveIntegrand_Input = ObjectiveIntegrand_Input,
                            AggregationIntegrand_Input = AggregationIntegrand_Input,
                            #StateDynamic_conc_Input = StateDynamic_conc_Input, 
                            StateDynamic_dist_Input = StateDynamic_dist_Input,
                            InitDist_Input = InitDist_Input,
                            SalvageFunction_Input = SalvageFunction_Input)

############################################################################################################

ObjectiveIntegrand_Input(Sconc,Sdist,Sagg,Cconc,Cdist,t,s,Para) = Para["Utility"](Sconc[1]*(1-Cconc[1])*1/Para["ω"]*(1-Cdist[1]))
AggregationIntegrand_Input(Sconc,Sdist,Cconc,Cdist, t::Float64,s::Float64, Para::Dict) = [(Sdist[1]+1e-2)^Para["α"]*Para["A"](s)]
StateDynamic_conc_Input(Sconc,Sagg,Cconc, t::Float64, Para::Dict) = [] 
StateDynamic_dist_Input(Sconc,Sdist,Sagg,Cconc,Cdist,t::Float64, s::Float64, Para::Dict) = [-Para["δ"]*Sdist[1] + Sagg[1]*(1-Cconc[1])*1/Para["ω"]*Cdist[1]]
InitDist_Input(Sconc,Sagg,Cconc, t::Float64, Para::Dict) = [Sagg[1]*Cconc[1]]
SalvageFunction_Input(Sconc,Sdist,t::Float64,s::Float64,Para) = 0.0

MyPara = Dict()
#MyPara["OptiType"] = "Newton-Raphson"
MyPara["OptiType"] = "Gradient"

MyPara["T"] = 40
MyPara["ω"] = 10
MyPara["hstep"] = 0.25

MyPara["nCon_conc"] = 1
MyPara["nCon_dist"] = 1
MyPara["nStat_conc"] = 0
MyPara["nStat_dist"] = 1
MyPara["nStat_agg"] = 1

MyPara["InitLineStep"] = 1e-3
MyPara["LineSearchStepIncrease"] = 0.25
MyPara["UpperLineStep"] = 1e-1
#MyPara["UpperLineStepTemporaryFactor"] = 1.5
MyPara["hLowBound"] = 0.2
MyPara["PlotResultsIntermediateFrequency"] = 50
MyPara["PlotResultsWaitForKey"] = false

MyPara["InitStat_dist"] = [5.5,5.5]
MyPara["Con_concMin"] = [0.0]
MyPara["Con_concMax"] = [0.99]
MyPara["Con_distMin"] = [0.0]
MyPara["Con_distMax"] = [0.99]

MyPara["ConSmooth"] = false

MyPara["Utility"] = (c) -> ((c+1e-6)^(1-1.1) - 1)/(1-1.1)
MyPara["A"] = (s) -> 1.0 - 0.5*(s-MyPara["ω"]/2)^2/(MyPara["ω"]/2)^2
MyPara["α"] = 0.3
MyPara["δ"] = 0.1
MyPara["ρ"] = 0.1
MyPara["TimeDiscountRate"] = MyPara["ρ"]

MyPara["LoadInits"] = true
Results = Dict()
Results["Con_conc"] = 0.6*ones(1,20,MyPara["nCon_conc"])
Results["Con_dist"] = 0.3*ones(10,20,MyPara["nCon_dist"])


UserParameters = MyPara
 
#using BenchmarkTools, Profile, Coverage


Results = VintageOptimisation(Results = Results,UserParameters = MyPara,
                            ObjectiveIntegrand_Input = ObjectiveIntegrand_Input,
                            AggregationIntegrand_Input = AggregationIntegrand_Input,
                            #StateDynamic_conc_Input = StateDynamic_conc_Input, 
                            StateDynamic_dist_Input = StateDynamic_dist_Input,
                            InitDist_Input = InitDist_Input,
                            SalvageFunction_Input = SalvageFunction_Input)



SaveResults(Results,"Test")

#Results = LoadResults("Results/NewCalibration5.out")




Results = LoadResults("Results/NewCalibration12.out")
PlotResults(Results;SavePlot = true)
