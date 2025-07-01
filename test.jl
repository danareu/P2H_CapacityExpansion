using Pkg
Pkg.add("AxisKeys")
using AxisKeys
include("./P2H_CapacityExpansion.jl")
cd("/cluster/home/danare/git")
Pkg.activate(".")
using .P2H_CapacityExpansion
using CPLEX
using Parameters
using JuMP
using XLSX
using PlotlyJS
using Surrogates
using Gurobi
using Dates
ENV["CPLEX_STUDIO_BINARIES"] = "/cluster/home/danare/opt/ibm/ILOG/CPLEX_Studio2211/cplex/bin/x86-64_linux/"
Pkg.add("CPLEX")



##########################################################################################################
################################################ Sampling ################################################
##########################################################################################################




# read in the data
config = P2H_CapacityExpansion.read_yaml_file();
data = P2H_CapacityExpansion.load_cep_data(config=config);
ts_data = P2H_CapacityExpansion.load_timeseries_data_full(config=config);
gas_gen = [gen for (gen, props) ∈ config["techs"] if haskey(props, "input") && get(props["input"], "fuel", nothing) == "R_Gas"]


# UNCERTAINTY SET WITH LHS 
h2_demand = [0.75, 1, 1.25]
capex = [[420,260,100], [520,360,200], [620,460,300]]
gas_price = [0.015, 0.04, 0.060]
capture_rate = [0.56, 0.745, 0.93]
lifetime = [13, 16, 20]
## not considered so far
grid_expansion = [1.5, 1.375, 1.25]


## lower and upper bounds
lb = [0.75, 420, 260, 100, 0.015, 0.56, 13]
ub = [1.25, 620, 460, 300, 0.06, 0.93, 20]


# Number of samples
n_samples = [1]


for n in n_samples
  # Latin Hypercube Sampling
  scenarios = Surrogates.sample(n,lb,ub, Surrogates.LatinHypercubeSample())

  k_inv = setdiff([key for (key, val) ∈ config["techs"] if get(val, "inv", "")  == true], [key for (key, val) ∈ config["techs"] if get(val, "tech_group", "")  == "transmission"] ) 
  push!(k_inv, "Cost")
  result = "/cluster/home/danare/git/P2H_CapacityExpansion/results/$(n)_scenarios_v3.txt"

  ### HYDROGEN DEMAND ###
  open(result, "a") do io
    println(io, join(k_inv, ", "))
    for (j, s) in enumerate(scenarios)
      @info "$j scenario ..."
      
      # check if already there
      k = 1
      for (i, y) ∈ enumerate(config["year"])
        data.data["demand"][:, y, "H2"] *= s[1]
        data.data["c_CAPEX"]["X_Electrolysis", y] = s[1+k]
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

      cep = P2H_CapacityExpansion.run_opt(ts_data=ts_data, data=data, config=config)

      result = P2H_CapacityExpansion.optimize_and_output(cep=cep, config=config, data=data, ts_data=ts_data, name="$(j)_scenario", short_sol=true)
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
        println(io, join(tmp_list, ", "))
      end
    end
  end
end