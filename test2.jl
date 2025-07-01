using Pkg
include("./P2H_CapacityExpansion.jl")
cd("/cluster/home/danare/git")
Pkg.activate(".")
using .P2H_CapacityExpansion
using JuMP, MathOptAI
using DataFrames
using ScikitLearn
using JuMP
using Ipopt
using Parameters
using Dates




n0 = Dates.now()

file = "/cluster/home/danare/git/P2H_CapacityExpansion/results/aggregated_results/500_scenarios.txt"

df = P2H_CapacityExpansion.read_txt_file(file);

X_train, y_train, X_test, y_test = P2H_CapacityExpansion.partitionTrainTest(df, :Cost, 0.8)

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

sg = P2H_CapacityExpansion.neural_network_model_flux(X_train, y_train, X_test, y_test, μy, σy; hidden_layer=1000, epochs=1000)

n1 = Dates.now()
println(n1-n0)


# load input
config = P2H_CapacityExpansion.read_yaml_file();
data = P2H_CapacityExpansion.load_cep_data(config=config);
ts_data = P2H_CapacityExpansion.load_timeseries_data_full(config=config);

cep = P2H_CapacityExpansion.run_opt(ts_data=ts_data, data=data, config=config, surrogate=true, solver=Ipopt.Optimizer)

optimizer = optimizer_with_attributes(
    Ipopt.Optimizer,
    "max_iter" => 10000,  # or higher
    "tol" => 1e-5,
    "print_level" => 5,
)

set_optimizer(cep.model, optimizer)

@unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = P2H_CapacityExpansion.get_sets(cep=cep)
data = data.data;

##################### cost optimization #####################

P2H_CapacityExpansion.setup_opt_costs_fix!(cep, config, data,vcat(cep.sets["non_dispatch"], cep.sets["dispatch"], 𝓢, String[s for s in cep.sets["discharging"]], cep.sets["conversion"]))

#JuMP.fix.(cep.model[:COST]["var", :, :], 0; force=true);
@variable(cep.model, COST_VAR[y ∈ 𝓨] ≥ 0);

techs = [f for f ∈ cep.sets["invest_tech"] if f != "ENS"]

for y ∈ 𝓨, r ∈ 𝓡
    x_vec = [cep.model[:TotalCapacityAnnual][r, g, y] for g ∈ techs]
    prediction, formulation = MathOptAI.add_predictor(cep.model, sg.model, x_vec)
    @constraint(cep.model, COST_VAR[y] .>= prediction)
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

data = P2H_CapacityExpansion.load_cep_data(config=config);
ts_data = P2H_CapacityExpansion.load_timeseries_data_full(config=config);

result = P2H_CapacityExpansion.optimize_and_output(cep=cep, config=config, data=data, ts_data=ts_data, name="scenario_v3", short_sol=true)


n2 = Dates.now()
println(n2-n1)



#### print the results
result_f = "/cluster/home/danare/git/P2H_CapacityExpansion/results/_scenarios_V4.txt"

open(result_f, "a") do io
    for r ∈ 𝓡 ,g ∈ cep.sets["nodes"], y ∈ 𝓨
        str = "TotalCapacityAnnual[$r,$g,$y]"
        val = value.(cep.model[:TotalCapacityAnnual][r,g,y])
        println(io, "$str = $val")
    end 
end