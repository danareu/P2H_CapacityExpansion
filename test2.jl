using Pkg
include("./P2H_CapacityExpansion.jl")
cd("/cluster/home/danare/git")
Pkg.activate(".")
using .P2H_CapacityExpansion
using DataFrames
using Parameters
using Flux
using Surrogates
using ScikitLearn
using LinearAlgebra, Random, Statistics
using JuMP
using XLSX
using MathOptAI
using PlotlyJS
using Ipopt
using Clustering
using CSV
using Dates
using StatsBase, MultivariateStats

file = "/cluster/home/danare/git/P2H_CapacityExpansion/results/500_scenarios_V3.txt"

df = P2H_CapacityExpansion.read_txt_file(file);
df = select(df, Not(:ENS))
size(df)

X_train, y_train, X_test, y_test = P2H_CapacityExpansion.partitionTrainTest(df, [:Cost,:Generation, :Emission], 0.7)

### scale the data ###
X_train_scaled, μX, σX  = P2H_CapacityExpansion.scaling(X_train)
X_test_scaled = (X_test .- μX) ./ σX
y_train_scaled, μy, σy  = P2H_CapacityExpansion.scaling(y_train)
    
# remove np.nan #
for i in eachindex(X_test_scaled)
    if isnan(X_test_scaled[i])
        X_test_scaled[i] = 0.0
    end
end

sg = P2H_CapacityExpansion.neural_network_model_flux(X_train_scaled, y_train_scaled, X_test_scaled, y_test, σy, μy; hidden_layer=500, epochs=100)


# load input
config = P2H_CapacityExpansion.read_yaml_file();
data = P2H_CapacityExpansion.load_cep_data(config=config);
ts_data = P2H_CapacityExpansion.load_timeseries_data_full(config=config);


cep = P2H_CapacityExpansion.run_opt(ts_data=ts_data, data=data, config=config, surrogate=true, solver=Ipopt.Optimizer)

@unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = P2H_CapacityExpansion.get_sets(cep=cep)

##################### cost optimization #####################

P2H_CapacityExpansion.setup_opt_costs_fix!(cep, config, data.data,vcat(cep.sets["non_dispatch"], cep.sets["dispatch"], 𝓢, String[s for s in cep.sets["discharging"]], cep.sets["conversion"]))

#JuMP.fix.(cep.model[:COST]["var", :, :], 0; force=true);
@variable(cep.model, COST_VAR[y ∈ 𝓨] ≥ 0);



techs = [f for f ∈ cep.sets["invest_tech"] if f != "ENS"]

for y ∈ 𝓨, r ∈ 𝓡
    x_vec = [cep.model[:TotalCapacityAnnual][r, g, y] for g ∈ techs]
    x_scaled = (x_vec ) #.- μX) ./ σX
    prediction, formulation = MathOptAI.add_predictor(cep.model, sg.model, x_scaled)
    y_rescaled = prediction #.* σy .+ μy

    ##### Cost approximation ###
    @constraint(cep.model, COST_VAR[y] .>= y_rescaled[1])

    ##### Generation  ###
    #TODO H2 and electricity combined as sum
    @constraint(cep.model, y_rescaled[2] .>= data.data["demand"][r,y,"electricity"] + data.data["demand"][r,y,"H2"])

    ##### Emission reduction  ###
    @constraint(cep.model, sum(data.data["budget"][r,y] for r ∈ 𝓡).>= y_rescaled[3] )
end


opex_discounted = sum(
    1 / ((1 + config["r"])^(y - 𝓨[1] - 10)) * (
        sum(cep.model[:COST]["fix", y, g] for g ∈ 𝓖) +
        COST_VAR[y]
    ) for y ∈ 𝓨
)

@objective(cep.model, Min, sum(
    1 / ((1 + config["r"])^(y - 𝓨[1] - 10)) *
    sum(cep.model[:COST]["cap", y, g] for g ∈ cep.sets["invest_all"])
    for y ∈ 𝓨 
) + opex_discounted)

result = P2H_CapacityExpansion.optimize_and_output(cep=cep, config=config, data=data, ts_data=ts_data, name="scenario_v3", short_sol=false)