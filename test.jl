using Pkg
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

# run the optimization model
model = P2H_CapacityExpansion.run_opt(ts_data=ts_data, data=data, config=config);

result = P2H_CapacityExpansion.optimize_and_output(cep=model)