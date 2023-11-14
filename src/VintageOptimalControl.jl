module VintageOptimalControl

    using BenchmarkTools
    using ColorSchemes
    using Coverage
    using Dierckx
    using ForwardDiff
    using Interpolations
    using LaTeXStrings
    using LinearAlgebra
    using Plots
    using Printf
    using Profile
    using Serialization
    using JLD2
    using FileIO


    include("AuxiliaryFunction.jl")
    include("LineSearch.jl")
    include("MainFunction.jl")
    include("ModelFunctions.jl")
    include("ParametersVariablesSettings.jl")
    include("ResultsHandling.jl")
    include("StateSolvers.jl")

    export VintageOptimisation
    export LineSearch
    export PlotResults
    export SaveResults
    export LoadResults

end
