# Dana Reulein, 2023

######################
# HydrogenExpansion
#####################


module P2H_CapacityExpansion
    using Pkg
    using Parameters
    using JuMP
    using Gurobi
    using XLSX
    using PlotlyJS
    using CPLEX
    using DataFrames
    using StatsBase
    using CSV
    using Clustering
    using Distances
    using YAML
    include("/cluster/home/danare/git/Clustering/TSClustering.jl")
    using .TSClustering
    using Dates
    #@reexport using TimeSeriesClustering

    const DIR = dirname(@__DIR__)
    include(joinpath("utils","post_processing.jl"))
    include(joinpath("utils","load_data.jl"))
    include(joinpath("src","set.jl"))
    include(joinpath("utils","datastructs.jl"))
    include(joinpath("src","opt.jl"))
    include(joinpath("results","plot_result.jl"))
end
