using Pkg

# Set the CPLEX path BEFORE adding/building
ENV["CPLEX_STUDIO_BINARIES"] = "/cluster/home/danare/opt/ibm/ILOG/CPLEX_Studio2211/cplex/bin/x86-64_linux/"

# Optional: only needed once (comment if already installed)
Pkg.add("CPLEX")
Pkg.build("CPLEX")

# Activate your local project
cd("/cluster/home/danare/git/P2H_CapacityExpansion")
Pkg.activate(".")

# Load your package (after activating the project)
include("./P2H_CapacityExpansion.jl")
using .P2H_CapacityExpansion

# Load other dependencies
using AxisKeys
using Parameters
using JuMP
using Gurobi
using CPLEX
using XLSX
using PlotlyJS
using Surrogates
using Random
using Dates



# read in the data
config = P2H_CapacityExpansion.read_yaml_file();
data = P2H_CapacityExpansion.load_cep_data(config=config);
ts_data = P2H_CapacityExpansion.load_timeseries_data_full(config=config);
gas_gen = [gen for (gen, props) ∈ config["techs"] if haskey(props, "input") && get(props["input"], "fuel", nothing) == "R_Gas"]

# UNCERTAINTY SET WITH LHS 
h2_demand = [0.75, 1, 1.25]
capex_electrolyzer = [[420,260,100], [520,360,200], [620,460,300]]
gas_price = [0.015, 0.04, 0.060]
capture_rate = [0.56, 0.745, 0.93]
lifetime = [13, 16, 20]
capex_pv = [[530,430,330], [770, 638, 507], [1010,847,685]] #https://www.sciencedirect.com/science/article/pii/S0306261925005860
capex_on = [[940, 880,820], [1300,1247,1195], [1660, 1615,1570]] #https://www.sciencedirect.com/science/article/pii/S0306261925005860
capex_off = [[1700, 1500,1300], [2700,2400,2100], [3700,3300,2900]] #https://www.sciencedirect.com/science/article/pii/S0306261925005860


## not considered so far
grid_expansion = [1.5, 1.375, 1.25]


## lower and upper bounds
lb = [0.75, 420, 260, 100, 0.015, 0.56, 13, 530,430,330, 940, 880,82, 1700, 1500,1300]
ub = [1.25, 620, 460, 300, 0.06, 0.93, 20, 1010,847,685, 1660, 1615,1570, 3700,3300,2900]


# Number of samples
n = 700

# Latin Hypercube Sampling
Random.seed!(1)
scenarios = Surrogates.sample(n,lb,ub, Surrogates.LatinHypercubeSample())

k_inv = setdiff([key for (key, val) ∈ config["techs"] if get(val, "inv", "")  == true], [key for (key, val) ∈ config["techs"] if get(val, "tech_group", "")  == "transmission"] ) 
push!(k_inv, "Cost", "Generation", "Emission")

result = "/cluster/home/danare/git/P2H_CapacityExpansion/results/$(n)_scenarios_V3.txt"

### HYDROGEN DEMAND ###
open(result, "a") do io
  println(io, join(k_inv, ", "))
  for (j, s) in enumerate(scenarios)
    @info "$j scenario ..."
      
    k = 1
    for (i, y) ∈ enumerate(config["year"])
      data.data["demand"][:, y, "H2"] *= s[1]
      data.data["c_CAPEX"]["X_Electrolysis", y] = s[1+k]
      data.data["c_CAPEX"]["P_PV_Utility_Avg", y] = s[7+k]
      data.data["c_CAPEX"]["P_Wind_Onshore_Avg", y] = s[10+k]
      data.data["c_CAPEX"]["P_Wind_Offshore_Transitional", y] = s[13+k]
      k += 1
    end

    ### GAS PRICE ###
    for g ∈ gas_gen
      data.data["c_var"][g, :] = data.data["c_var"][g, :] .+ s[5]
    end
        
    ### CAPTURE RATE ###
    data.data["emission"]["X_ATR_CCS"] = data.data["emission"]["X_ATR_CCS"] * (1 - s[6])

    ### LIFE TIME ###
    data.data["lifetime"]["X_Electrolysis"] = s[7]

    cep = P2H_CapacityExpansion.run_opt(ts_data=ts_data, data=data, config=config, surrogate=false,solver=Gurobi.Optimizer)
    result = P2H_CapacityExpansion.optimize_and_output(cep=cep, config=config, data=data, ts_data=ts_data, name="$(j)_scenario_v4", short_sol=true)

    #### print the results

    for y ∈ config["year"]
      tmp_c = 0
      tmp_list = []
      for g ∈ cep.sets["invest_tech"]
        tmp_c += value.(cep.model[:COST]["var",y,g])
        tmp_capa = 0
        for r ∈ config["countries"]
          tmp_capa += value.(cep.model[:TotalCapacityAnnual][r,g,y])
        end
        push!(tmp_list, tmp_capa)
      end
      push!(tmp_list, tmp_c)

      # add the summed GENERATION
      tmp_gen = 0
      for g ∈ setdiff(cep.sets["techs"], cep.sets["storage_techs"]), r ∈ config["countries"], t ∈ axes(ts_data.ts)[3], c ∈ cep.sets["carrier"][g]
          tmp_gen += value.(cep.model[:gen][r,g,y,c,t])
      end
      push!(tmp_list, tmp_gen)

      ## CO2 Emission
      push!(tmp_list, value.(cep.model[:em][y]))

      println(io, join(tmp_list, ", "))
    end
  end
end