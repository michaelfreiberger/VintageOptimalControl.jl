
ObjectiveIntegrand_Input(Sconc,Sdist,Sagg,Cconc,Cdist,t,s,Para) = exp(-Para["ρ"]*t) * Para["Utility"](Sagg[1]*(1-Cdist[1])*1/Para["ω"])
AggregationIntegrand_Input(Sconc,Sdist,Cconc,Cdist, t::Float64,s::Float64, Para::Dict) = [(Sdist[1]+1e-1)^Para["α"]*Para["A"](s)]
StateDynamic_conc_Input(Sconc,Sagg,Cconc, t::Float64, Para::Dict) = [] 
StateDynamic_dist_Input(Sconc,Sdist,Sagg,Cconc,Cdist,t::Float64, s::Float64, Para::Dict) = [-Para["δ"]*Sdist[1] + Sagg[1]*Cdist[1]*1/Para["ω"]]
InitDist_Input(Sconc,Sagg,Cconc, t::Float64, Para::Dict) = [0.0]
SalvageFunction_Input(Sconc,Sdist,t::Float64,s::Float64,Para) = 0.0

MyPara = Dict()
MyPara["OptiType"] = "Newton-Raphson"
#MyPara["OptiType"] = "Gradient"

MyPara["T"] = 40
MyPara["ω"] = 10
MyPara["hstep"] = 0.25

MyPara["nCon_conc"] = 0
MyPara["nCon_dist"] = 1
MyPara["nStat_conc"] = 0
MyPara["nStat_dist"] = 1
MyPara["nStat_agg"] = 1

MyPara["InitLineStep"] = 1e-7
#MyPara["LineSearchStepIncrease"] = 0.05
MyPara["UpperLineStep"] = 1e-1
MyPara["hLowBound"] = 0.2
MyPara["PlotResultsIntermediateFrequency"] = 50
MyPara["PlotResultsWaitForKey"] = false

MyPara["InitStat_dist"] = [0.5,0.8]
MyPara["Con_concMin"] = [0.00]
MyPara["Con_concMax"] = [0.99999]
MyPara["Con_distMin"] = [0.0]
MyPara["Con_distMax"] = [0.99999]

MyPara["ConSmooth"] = false

MyPara["Utility"] = (c) -> ((c+1e-6)^(1-1.1) - 1)/(1-1.1)
MyPara["A"] = (s) -> 1.0 - 0.5*(s-MyPara["ω"]/2)^2/(MyPara["ω"]/2)^2
MyPara["α"] = 0.3
MyPara["δ"] = 0.1
MyPara["ρ"] = 0.1
MyPara["TimeDiscountRate"] = MyPara["ρ"]

MyPara["LoadInits"] = true
Results = Dict()
Results["Con_conc"] = 0.1*ones(10,20,MyPara["nCon_conc"])
Results["Con_dist"] = 0.1*ones(10,20,MyPara["nCon_dist"])


UserParameters = MyPara
#
 
using BenchmarkTools, Profile, Coverage


Results = VintageOptimisation(Results = Results,UserParameters = MyPara,
                            ObjectiveIntegrand_Input = ObjectiveIntegrand_Input,
                            AggregationIntegrand_Input = AggregationIntegrand_Input,
                            #StateDynamic_conc_Input = StateDynamic_conc_Input, 
                            StateDynamic_dist_Input = StateDynamic_dist_Input,
                            InitDist_Input = InitDist_Input,
                            SalvageFunction_Input = SalvageFunction_Input)

MyPara["LoadInits"] = true
Results2 = Dict()
Results2["Con"] = copy(Results["Con"])
Results2["Con_dist"] = copy(Results["Con_dist"])
Results2["Con_dist"] .= Results["Con"]
Q(Con,Stat, t::Float64,s::Float64, Para::Dict) = Stat[2]*(Stat[1]^0.5 - Con[1]^2)

f1(Con,Stat, t::Float64, Para::Dict) = [Con[1] - 0.1*Stat[1],
                                        -0.5*Stat[2]]

Results2 = VintageOptimisation(Results = Results2,UserParameters = MyPara,
                ObjectiveIntegrand2 = U, 
                AggregationFunction2 = Q,
                StateDynamic_1_2 = f1,
                StateDynamic_2_2 = f2, 
                Shock2 = g,
                SalvageFunction_1_2=S1,
                SalvageFunction_2_2=S2)


SaveResults(Results,"Test")

#Results = LoadResults("Results/NewCalibration5.out")




Results = LoadResults("Results/NewCalibration12.out")
PlotResults(Results;SavePlot = true)
