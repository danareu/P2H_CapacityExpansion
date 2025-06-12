using Pkg
Pkg.add("AxisKeys")
using AxisKeys
include("./P2H_CapacityExpansion.jl")
cd("/cluster/home/danare/git")
Pkg.activate(".")
using .P2H_CapacityExpansion
using CPLEX
using JuMP
using XLSX
using PlotlyJS
using Gurobi
using Dates
ENV["CPLEX_STUDIO_BINARIES"] = "/cluster/home/danare/opt/ibm/ILOG/CPLEX_Studio2211/cplex/bin/x86-64_linux/"
Pkg.add("CPLEX")


# read in the data
config = P2H_CapacityExpansion.read_yaml_file();
data = P2H_CapacityExpansion.load_cep_data(config=config);
ts_data = P2H_CapacityExpansion.load_timeseries_data_full(config=config);


#### sensitivities TODO remove that later 
for r in config["countries"], y in 2020:10:2050
  data.data["demand"][r, y, "electricity"] = data.data["demand"][r,y, "electricity"] *  0.7
end





if config["dispatch"]
  data.data["demand"][:, config["year"], "electricity"] = data.data["demand"][:, config["year"], "electricity"] *  0.5
  data.data["demand"]["DE", config["year"], "H2"] = 1
  data.data["cap_init"]["DE", "D_Battery_Li-Ion_in", config["year"]] = 2 
  data.data["cap_init"]["DE", "X_Electrolysis", config["year"]] = 20 
  data.data["cap_init"]["DE", "D_Battery_Li-Ion_out", config["year"]] = 2 
  data.data["cap_init"]["DE", "S_Battery_Li-Ion", config["year"]] = 6 

  data.data["cap_init"]["DE", "D_Gas_H2_in", config["year"]] = 3 
  data.data["cap_init"]["DE", "P_H2_OCGT", config["year"]] = 15 
  data.data["cap_init"]["DE", "D_Gas_H2_out", config["year"]] = 3 
  data.data["cap_init"]["DE", "S_Gas_H2", config["year"]] = 168*3 
end

# for t in ["S_Gas_H2", "D_Gas_H2_in", "D_Gas_H2_out", "P_H2_OCGT"]
#     for v in ["c_CAPEX", "c_var", "c_fix"]
#       println(v, data.data[v][t,2020])
#     end
#     println("cap", data.data["cap"]["DE", t, 2020])
#     println("cap_init", data.data["cap_init"]["DE",t,2020])
#     println("eta", data.data["eta"][t,2020])
#     println("lifetime",data.data["lifetime"][t])
# end 
  

# run the optimization model
model = P2H_CapacityExpansion.run_opt(ts_data=ts_data, data=data, config=config);

result = P2H_CapacityExpansion.optimize_and_output(cep=model, config=config, data=data, ts_data=ts_data)