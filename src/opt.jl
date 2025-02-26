"""
run_opt(ts_data::ClustData,opt_data::OptDataCEP,config::Dict{String,Any},optimizer::DataType)
Organizing the actual setup and run of the CEP-Problem.
Required elements are:
- `ts_data`: The time-series data.
- `opt_data`: In this case the OptDataCEP that contains information on costs, nodes, techs and for transmission also on lines.
- `config`: This includes all the settings for the design optimization problem formulation.
"""


function run_opt(; ts_data::JuMP.Containers.DenseAxisArray,
    data::Dict{Any, Any},
    config::Dict{Any, Any},
    benders::Bool=false,
    master::Bool=false,
    kwargs...
    )

    if :scenarios in keys(kwargs)
        cep=setup_opt_basic(ts_data=ts_data, scenarios=Dict(kwargs)[:scenarios])
    else
        cep=setup_opt_basic(ts_data=ts_data, scenarios=collect(1:size(ts_data)[3]))
    end

    # solve with benders benders_decomposition
    if benders && master
        set_up_equations_master(cep=cep, ts_data=ts_data, data=data, config=config, probabilities=Dict(kwargs)[:probabilities])
    elseif benders && !master
        set_up_equations_sub_problem(cep=cep, ts_data=ts_data, data=data, config=config, scenarios=Dict(kwargs)[:scenarios])
    elseif !benders && !master 
        set_up_equations(cep=cep, ts_data=ts_data, data=data, config=config)
    end
    return cep
end



function setup_opt_basic(;ts_data::JuMP.Containers.DenseAxisArray, scenarios::Vector{Int64})
    
    # Initialize model
    model =  JuMP.Model()
    set_optimizer(model, CPLEX.Optimizer)
    # Initialize solver
    optimizer = optimizer_with_attributes(
        () -> Gurobi.Optimizer(Gurobi.Env()),
        "LogToConsole" => 0,
        "OutputFlag"=> 0,
        "Method" => 2,
        "BarHomogeneous" =>  1,
        "ResultFile" =>  "Solution_julia.sol"
    )
    #set_optimizer(model, optimizer)
    # Setup set
    sets=setup_opt_set(ts_data=ts_data, scenarios=scenarios)
    # Setup Model CEP
    return OptModelCEP(model,sets,Dict())
end

function setup_opt_set(; ts_data::JuMP.Containers.DenseAxisArray, scenarios::Vector{Int64})
    sets=Dict()
    # scenarios
    sets["scenario"] = scenarios
    # generators
    sets["generators"] = ["P_Nuclear","P_Coal_Hardcoal","P_Gas","RES_Wind_Onshore_Avg","RES_Wind_Offshore_Transitional","RES_PV_Utility_Avg","X_Electrolysis", "Pumped_Hydro"];       
    # investment years
    sets["years"] = 2020:10:2050
    # time steps
    sets["timesteps"] = axes(ts_data)[4]
    return sets
end



function set_up_equations(; cep::OptModelCEP, ts_data::JuMP.Containers.DenseAxisArray, data,config, kwargs...)

    # set up sets
    𝓢 = cep.sets["scenario"]
    𝓖 = cep.sets["generators"]     
    𝓨 = cep.sets["years"]
    𝓣 = cep.sets["timesteps"]


    #set up params
    c_ll = config["cll"];                       # cost for lost load
    c_fix = data["c_fix"];                      # fixed costs
    c_var = data["c_var"];                       # variable costs
    c_capex = data["c_CAPEX"];                  # investment costs
    cap_max = data["cap"];                      # max capacity investment per year
    r = config["r"];                                   # interest rate
    η = data["eta"];                            # efficiency
    #π = Dict(s => 1/length(𝓢) for s in 𝓢)                         # probability
    if :probabilities in keys(kwargs)
        c = Dict(kwargs)
        π = c[:probabilities]
    elseif length(𝓢) == 1
        π = Dict(s => 1/length(𝓢) for s in 𝓢)
    else
        random_numbers = abs.(rand(length(𝓢)))
        π = Dict(s => random_numbers[s]/ sum(random_numbers) for s in 𝓢) 
    end
    cep.params["π"] = π
    lifetime = data["lifetime"]                 # lifetime of generators
    cap_res = data["cap_init"]                  # residual capacity in 2020
    e = data["emission"]                        # emission content per generators
    e_pen = data["emission_penalty"]            # emission penalty for fossil fuel generators
    d_power = data["d_power"]                  # power demand level annual
    d_h2 = data["d_h2"]                      # hydrogen demand https://ehb.eu/files/downloads/EHB-Analysing-the-future-demand-supply-and-transport-of-hydrogen-June-2021-v3.pdf
    timesteplength = config["season_length"]

    # set up variables
    @variable(cep.model, gen[s in 𝓢, y in 𝓨,g in 𝓖,t in 𝓣] ≥ 0)         # planned generation for generators
    @variable(cep.model, opex[s in 𝓢,y in 𝓨] ≥ 0)                      # operation costs for generators
    @variable(cep.model, capex[y in 𝓨] ≥ 0)                     # capital investment costs for generators
    @variable(cep.model, cap_new[y in 𝓨,g in 𝓖] ≥ 0)           # new capacity investments for generators
    @variable(cep.model, ll[s in 𝓢,y in 𝓨, t in 𝓣] ≥ 0)               # lost load / ENS
    @variable(cep.model, cap[y in 𝓨,g in 𝓖] ≥ 0)                # old and new capacity for generators
    @variable(cep.model, cap_acc[y in 𝓨,g in 𝓖] ≥ 0)            # accumulated capacity according to lifetime for generators
    @variable(cep.model, em[s in 𝓢,y in 𝓨] ≥ 0)                        # emission CO2 per year
    @variable(cep.model, cll[s in 𝓢,y in 𝓨] ≥ 0)  # costs for lost load yearly
    @variable(cep.model, ll_h2[s in 𝓢,y in 𝓨, t in 𝓣] ≥ 0)

    # defining constraints
    @constraint(cep.model, EnergyBalance[s in 𝓢,y in 𝓨, t in 𝓣], sum(gen[s,y, g, t] for g in config["power"]) + ll[s,y, t] - (ts_data[y, "Demand",s,t] * d_power[y])  - (sum(gen[s,y,g,t] for g in config["h2"])) == 0); # energy balance equation
    @constraint(cep.model, CAPEX[y in 𝓨], sum(cap_new[y,g]*c_capex[g,y] for g in 𝓖)==capex[y]);       # new capacity investments costs
    @constraint(cep.model, OPEX[s in 𝓢,y in 𝓨], (sum(cap[y,g]*c_fix[g,y]+sum(gen[s,y,g,t]*c_var[g,y] for t in 𝓣) for g in 𝓖)+(em[s,y]*e_pen[y]*100))==opex[s,y]);   # operational costs
    @constraint(cep.model, GenCapDisp[s in 𝓢,y in 𝓨, g in config["dispatch"], t in 𝓣], gen[s,y,g,t] ≤ cap[y,g]* η[g, y]);    # limit max generation dispatchable
    @constraint(cep.model, GenCapH2[s in 𝓢,y in 𝓨, g in config["h2"], t in 𝓣], gen[s,y,g,t] ≤ cap[y,g]);    # limit max generation hydrogen
    @constraint(cep.model,GenCapNonDisp[s in 𝓢,y in 𝓨, g in config["nondispatch"], t in 𝓣], gen[s,y,g,t] ≤ cap[y,g]*ts_data[y, g,s,t]); # max generation non dispatchable
    @constraint(cep.model, DemandH2[s in 𝓢,y in 𝓨,t in 𝓣], sum((gen[s,y,g,t]*η[g, y]+ll_h2[s,y,t]) for g in config["h2"]) == d_h2[y]/8760); # annual energy balance hydrogen
    @constraint(cep.model, InitialCap[s in 𝓢,y in 𝓨[1], g in 𝓖], cap[y,g] == (haskey(cap_res, (g,y)) ? cap_res[g,y] : 0));# initialize capacity in 2020

    for g in 𝓖
        JuMP.fix(cep.model[:cap_new][𝓨[1],g], 0; force=true);
        JuMP.fix(cep.model[:cap_acc][𝓨[1],g], 0; force=true);
    end
    @constraint(cep.model, NewCap[y in 𝓨[2:end], g in 𝓖], cap[y,g] == cap_acc[y,g]+cap_new[y,g] + (haskey(cap_res, (g,y)) ? cap_res[g,y] : 0)); # new capacity investments 
    @constraint(cep.model, AccCap[y in 𝓨[2:end], g in 𝓖], cap_acc[y,g] == sum(cap_new[hat_y,g] for hat_y in 𝓨[1]:10:y if y - 𝓨[1] ≤ lifetime[g])); # accumulated capacity  
    @constraint(cep.model, MaxCap[y in 𝓨[2:end], g in 𝓖], cap[y,g] ≤ cap_max[g,y]); # max potential capacity constraint
    @constraint(cep.model, EM[s in 𝓢,y in 𝓨], em[s,y] == sum(sum(gen[s,y,g,t]*e[g] for t in 𝓣) for g in keys(e))); # emission accounting
    @constraint(cep.model, CLL[s in 𝓢,y in 𝓨], cll[s,y]== sum((ll[s,y,t]*c_ll+ll_h2[s,y,t]*c_ll) for t in 𝓣)); # cost for lost load yearly and per scenario


    # fix fist stage decision variable
    if :evp in keys(kwargs)
        for y in axes(evp["cap_new"])[1] for g in axes(evp["cap_new"])[2]
            fix.(cap_new[y,g], evp[y,g], force=true)
            println("Are they fixed? ", is_fixed.(cap_new[y,g]))
        end end
    end

    # define the objective 
    @objective(cep.model, Min, sum(1/((1+r)^(y-𝓨[1]))*capex[y] for y in 𝓨)+ sum(π[s]*sum(1/((1+r)^(y-𝓨[1]))*(opex[s,y]+cll[s,y]) for y in 𝓨) for s in 𝓢))
end





function optimize_and_output(; cep::OptModelCEP)

    optimize!(cep.model)

    status=Symbol(termination_status(cep.model))
    println(status)
    if termination_status(cep.model) == MOI.INFEASIBLE_OR_UNBOUNDED || termination_status(cep.model) == MOI.INFEASIBLE
        JuMP.compute_conflict!(cep.model)
        list_of_conflicting_constraints = ConstraintRef[]
        for (F, S) in list_of_constraint_types(cep.model)
            for con in all_constraints(cep.model, F, S)
                if get_attribute(con, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
                    push!(list_of_conflicting_constraints, con)
                            end
            end
        end
        
        open(normpath(joinpath(dirname(@__FILE__),"..", "iis.txt")), "w") do file
            for r in list_of_conflicting_constraints
                write(file, string(r)*"\n")
            end
        end
    else
        objective = objective_value(cep.model)
        println("\n\nObjective value is $objective")
        variables = Dict()
        𝓨 = 2020:10:2050
        𝓢 = 1:16
        𝓣 = axes(value.(cep.model[:gen]))[3]
        𝓖 = ["P_Nuclear","P_Coal_Hardcoal","P_Gas","RES_Wind_Onshore_Avg","RES_Wind_Offshore_Transitional","RES_PV_Utility_Avg","X_Electrolysis", "Pumped_Hydro"];       
        variables["cap"] = value.(cep.model[:cap])
        return OptResult(cep.model, status, objective, variables, cep.params["π"])
    end  
end





function set_up_equations_master(; cep::OptModelCEP, ts_data::JuMP.Containers.DenseAxisArray, data, config, kwargs...)
    # set up sets
    𝓢 = cep.sets["scenario"]
    𝓖 = cep.sets["generators"]     
    𝓨 = cep.sets["years"]
    𝓣 = cep.sets["timesteps"]

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
    #set up parameter
    if :probabilities in keys(kwargs)
        c = Dict(kwargs)
        π = c[:probabilities]
    elseif length(𝓢) == 1
        π = Dict(s => 1/length(𝓢) for s in 𝓢)
    else
        random_numbers = abs.(rand(length(𝓢)))
        π = Dict(s => random_numbers[s]/ sum(random_numbers) for s in 𝓢) 
    end
    cep.params["π"] = π 
    c_capex = data["c_CAPEX"];                  # investment costs
    cap_max = data["cap"];                      # max capacity investment per year
    r = config["r"];                                   # interest rate
    lifetime = data["lifetime"]                 # lifetime of generators
    cap_res = data["cap_init"]                  # residual capacity in 2020
    M = -1000;


    # set up variables
    @variable(cep.model, capex[y in 𝓨] ≥ 0)                     # capital investment costs for generators
    @variable(cep.model, cap_new[y in 𝓨,g in 𝓖] ≥ 0)           # new capacity investments for generators
    @variable(cep.model, cap[y in 𝓨,g in 𝓖] ≥ 0)                # old and new capacity for generators
    @variable(cep.model, cap_acc[y in 𝓨,g in 𝓖] ≥ 0)            # accumulated capacity according to lifetime for generators
    @variable(cep.model, θ[s in 𝓢] >= M)

        # defining constraints
    @constraint(cep.model, CAPEX[y in 𝓨], sum(cap_new[y,g]*c_capex[g,y] for g in 𝓖)==capex[y]);       # new capacity investments costs
    @constraint(cep.model, InitialCap[s in 𝓢,y in 𝓨[1], g in 𝓖], cap[y,g] == (haskey(cap_res, (g,y)) ? cap_res[g,y] : 0));# initialize capacity in 2020
    @constraint(cep.model, NewCap[y in 𝓨[2:end], g in 𝓖], cap[y,g] == cap_acc[y,g]+cap_new[y,g] + (haskey(cap_res, (g,y)) ? cap_res[g,y] : 0)); # new capacity investments 
    @constraint(cep.model, AccCap[y in 𝓨[2:end], g in 𝓖], cap_acc[y,g] == sum(cap_new[hat_y,g] for hat_y in 𝓨[1]:10:y if y - 𝓨[1] ≤ lifetime[g])); # accumulated capacity  
    @constraint(cep.model, MaxCap[y in 𝓨[2:end], g in 𝓖], cap[y,g] ≤ cap_max[g,y]); # max potential capacity constraint

    for g in 𝓖
        JuMP.fix(cep.model[:cap_new][𝓨[1],g], 0; force=true);
        JuMP.fix(cep.model[:cap_acc][𝓨[1],g], 0; force=true);
    end

        # define the objective 
    @objective(cep.model, Min, sum(1/((1+r)^(y-𝓨[1]))*(capex[y]) for y in 𝓨) + sum(π[s]*θ[s] for s in 𝓢));
end




function set_up_equations_sub_problem(;cep::OptModelCEP, ts_data::JuMP.Containers.DenseAxisArray, scenarios, data, config, kwargs...)

   # set up sets
   𝓖 = cep.sets["generators"]     
   𝓨 = cep.sets["years"]
   𝓣 = cep.sets["timesteps"]

   #set up params
   c_ll = config["cll"];                       # cost for lost load
   c_fix = data["c_fix"];                      # fixed costs
   c_var = data["c_var"];                       # variable costs
   r = config["r"];                                   # interest rate
   η = data["eta"];                            # efficiency
   lifetime = data["lifetime"]                 # lifetime of generators
   cap_res = data["cap_init"]                  # residual capacity in 2020
   e = data["emission"]                        # emission content per generators
   e_pen = data["emission_penalty"]            # emission penalty for fossil fuel generators
   d_power = data["d_power"]                  # power demand level annual
   d_h2 = data["d_h2"]                      # hydrogen demand https://ehb.eu/files/downloads/EHB-Analysing-the-future-demand-supply-and-transport-of-hydrogen-June-2021-v3.pdf
   timesteplength = config["season_length"]

   # set up variables
   @variable(cep.model, gen[y in 𝓨,g in 𝓖,t in 𝓣] ≥ 0)         # planned generation for generators
   @variable(cep.model, opex[y in 𝓨] ≥ 0)                      # operation costs for generators
   @variable(cep.model, ll[y in 𝓨, t in 𝓣] ≥ 0)               # lost load / ENS
   @variable(cep.model, em[y in 𝓨] ≥ 0)                        # emission CO2 per year
   @variable(cep.model, cll[y in 𝓨] ≥ 0)  # costs for lost load yearly
   @variable(cep.model, cap[y in 𝓨,g in 𝓖] ≥ 0)
   @variable(cep.model, ll_h2[y in 𝓨, t in 𝓣] ≥ 0)  # lost load for hydrogen

   # defining constraints
   @constraint(cep.model, EnergyBalance[y in 𝓨, t in 𝓣], sum(gen[y, g, t] for g in config["power"]) + ll[y, t] - (ts_data[y, "Demand",scenarios[1],t] * d_power[y])  - (sum(gen[y,g,t] for g in config["h2"])) == 0); # energy balance equation
   @constraint(cep.model, OPEX[y in 𝓨], (sum(cap[y,g]*c_fix[g,y]+sum(gen[y,g,t]*c_var[g,y] for t in 𝓣) for g in 𝓖)+(em[y]*e_pen[y]*100))==opex[y]);   # operational costs
   @constraint(cep.model, GenCapDisp[y in 𝓨, g in config["dispatch"], t in 𝓣], gen[y,g,t] ≤ cap[y,g]* η[g, y]);    # limit max generation dispatchable
   @constraint(cep.model, GenCapH2[y in 𝓨, g in config["h2"], t in 𝓣], gen[y,g,t] ≤ cap[y,g]);    # limit max generation hydrogen
   @constraint(cep.model, GenCapNonDisp[y in 𝓨, g in config["nondispatch"], t in 𝓣], gen[y,g,t] ≤ cap[y,g]*ts_data[y,g,scenarios[1],t]); # max generation non dispatchable
   @constraint(cep.model, DemandH2[y in 𝓨, t in 𝓣], sum((gen[y,g,t]*η[g, y]+ll_h2[y, t]) for g in config["h2"]) == d_h2[y]/8760); # annual energy balance hydrogen
   @constraint(cep.model, EM[y in 𝓨], em[y] == sum(sum(gen[y,g,t]*e[g] for t in 𝓣) for g in keys(e))); # emission accounting
   @constraint(cep.model, CLL[y in 𝓨], cll[y]== sum((ll[y,t]*c_ll+ll_h2[y,t]*c_ll) for t in 𝓣)); # cost for lost load yearly and per scenario

   # define the objective 
   @objective(cep.model, Min, sum(1/((1+r)^(y-𝓨[1]))*((opex[y]+cll[y])) for y in 𝓨));
end



function solve_master_problem(; cep::OptModelCEP)

   # Solve the master problem
   optimize!(cep.model)
   status=Symbol(termination_status(cep.model))

   # Extract master problem solution
   if  termination_status(cep.model) == MOI.INFEASIBLE_OR_UNBOUNDED ||  termination_status(cep.model) == MOI.INFEASIBLE
       JuMP.compute_conflict!(cep.model)
       list_of_conflicting_constraints = ConstraintRef[]
       for (F, S) in list_of_constraint_types(cep.model)
           for con in all_constraints(cep.model, F, S)
               if get_attribute(con, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
                   push!(list_of_conflicting_constraints, con)
                           end
           end
       end
       
       open(normpath(joinpath(dirname(@__FILE__),"..", "iis.txt")), "w") do file
           for r in list_of_conflicting_constraints
               println(r)
               write(file, string(r)*"\n")
           end
       end
   else
       objective = objective_value(cep.model)
       println("\n\nObjective value is $objective")
       variables = Dict()
       variables["capex"] = value.(cep.model[:capex])
       variables["cap"] = value.(cep.model[:cap])
       return OptResult(cep.model, status, objective, variables, Dict())
   end  

end



function solve_sub_problem(; cep::OptModelCEP)
 
   optimize!(cep.model)
   status=Symbol(termination_status(cep.model))
   println(status)

   # Extract master problem solution
   if termination_status(cep.model) == MOI.INFEASIBLE_OR_UNBOUNDED || termination_status(cep.model) == MOI.INFEASIBLE
       JuMP.compute_conflict!(cep.model)
       list_of_conflicting_constraints = ConstraintRef[]
       for (F, S) in list_of_constraint_types(cep.model)
           for con in all_constraints(cep.model, F, S)
               if get_attribute(con, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
                   push!(list_of_conflicting_constraints, con)
                           end
           end
       end
       
       open(normpath(joinpath(dirname(@__FILE__),"..", "iis.txt")), "w") do file
           for r in list_of_conflicting_constraints
               println(r)
               write(file, string(r)*"\n")
           end
       end
   else
       objective = objective_value(cep.model)
       println("\n\nObjective value is $objective")
       variables = Dict()
       variables["λ"] = reduced_cost.(cep.model[:cap])
       return OptResult(cep.model, status, objective, variables, Dict())
   end 
end


function benders_decomposition(; ts_data::JuMP.Containers.DenseAxisArray, data, config, probabilities)

    MAXIMUM_ITERATIONS = 1000
    ABSOLUTE_OPTIMALITY_GAP = 1e-6

    # define sets
    𝓨 = 2020:10:2050
    𝓖 = config["generators"]
    𝓢 = axes(ts_data)[3]
    𝓣 =  axes(ts_data)[4]


    # set up master problem
    master_problem = run_opt(ts_data=ts_data, data=data, config=config, benders=true, master=true, probabilities=probabilities)

    # for each scenario store the subproblem in a dictionary 
    sub_problem = Dict(s => run_opt(ts_data=ts_data, data=data, scenarios=[s], config=config, benders=true, master=false) for s in 𝓢)

    for k in 1:MAXIMUM_ITERATIONS
        # solve master problem
        master_problem_result = solve_master_problem(cep=master_problem)
        lower_bound = master_problem_result.objective
        x_k = master_problem_result.variables["cap"]
        capex = master_problem_result.variables["capex"]

        λ = JuMP.Containers.DenseAxisArray(zeros(length(𝓨), length(𝓖), length(𝓢)),𝓨, 𝓖,1:length(𝓢))
        obj_sub = JuMP.Containers.DenseAxisArray(zeros(length(𝓢)), 1:length(𝓢))

        # for every scenario solve subproblem
        for (key, value) in sub_problem
            for y in 𝓨 for g in 𝓖
                fix.(value.model[:cap][y,g], x_k[y,g], force=true)
            end end
            sub_prob_result = solve_sub_problem(cep=value)
            for y in 𝓨 for g in 𝓖
                λ[y,g,key] = sub_prob_result.variables["λ"][y,g]
            end end 
            obj_sub[key] += sub_prob_result.objective           
        end 

        # determine upper bound and gap
        upper_bound =  (sum(1/((1+config["r"])^(y-𝓨[1]))*(capex[y]) for y in 𝓨)) + sum(obj_sub[s] *probabilities[s] for s in 𝓢)
        gap = (upper_bound - lower_bound) / upper_bound

        if gap < ABSOLUTE_OPTIMALITY_GAP
            println("Terminating with the optimal solution")
            return upper_bound
            break
        end

        # objective function depends on y and s  discounting!!!
        for s in 𝓢
            cut = @constraint(master_problem.model, master_problem.model[:θ][s] >= obj_sub[s] + sum(λ[y,g,s]*(master_problem.model[:cap][y,g]-x_k[y,g]) for g in 𝓖 for y in 𝓨))
            @info "Adding the cut $(cut)"
        end 
   end
end


 
