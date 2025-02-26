using Pkg
include("./P2H_CapacityExpansion.jl")
cd("/cluster/home/danare/git")
Pkg.activate(".")
using .P2H_CapacityExpansion
using CPLEX
using JuMP
using PlotlyJS
using Dates
ENV["CPLEX_STUDIO_BINARIES"] = "/cluster/home/danare/opt/ibm/ILOG/CPLEX_Studio2211/cplex/bin/x86-64_linux/"
Pkg.add("CPLEX")


for ls in [1]

    SCEN_NUM = ls # number of wheater years to consider

    ## LOAD DATA ##
    config = P2H_CapacityExpansion.read_yaml_file()
    data = P2H_CapacityExpansion.load_data()
    ts_data = P2H_CapacityExpansion.load_time_series(s_num=SCEN_NUM, config=config, timesteps=1:8760)

    ## OPTIMIZATION ##
    start_time_org = now()
    model = P2H_CapacityExpansion.run_opt(ts_data=ts_data, data=data, config=config, benders=false, master=false)
    result = P2H_CapacityExpansion.optimize_and_output(cep=model)
    probabilities = result.probabilities

end

