"""
run_opt(ts_data::ClustData,opt_data::OptDataCEP,config::Dict{String,Any},optimizer::DataType)
Organizing the actual setup and run of the CEP-Problem.
Required elements are:
- `ts_data`: The time-series data.
- `opt_data`: In this case the OptDataCEP that contains information on costs, nodes, techs and for transmission also on lines.
- `config`: This includes all the settings for the design optimization problem formulation.
"""


function run_opt(; ts_data,
    data::OptDataCEP,
    config::Dict{Any, Any},
    kwargs...
    )

    @info "Reading the data ..."
    cep=setup_opt_basic(ts_data=ts_data, config=config, data=data)

    @info "Setting up the optimization variables ..."
    setup_opt_basic_variables(cep=cep, config=config) 

    @info "Setting up the optimization equations ..."
    set_up_equations(cep=cep, ts_data=ts_data, data=data, config=config)
    setup_opt_storage!(cep, config, data)
    setup_opt_conversion!(cep, config, data)    
    set_opt_transmission!(cep, config, data)
    setup_opt_objective!(cep, config)

    return cep
end



function setup_opt_basic(;ts_data::JuMP.Containers.DenseAxisArray, 
    config::Dict{Any, Any},
    data::OptDataCEP)
    
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
    sets=Dict()
    setup_opt_set!(sets, ts_data, config, data)
    setup_opt_set_carrier!(sets, config)
    
    # Setup Model CEP
    return OptModelCEP(model,sets,Dict())
end




function setup_opt_basic_variables(; cep::OptModelCEP,  config::Dict{Any, Any})

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)

    if !config["dispatch"]
        @info "Investment mode is on."
        # capacity variables
        @variable(cep.model, TotalCapacityAnnual[r ∈ 𝓡 ,g ∈ sets["nodes"], y ∈ 𝓨] ≥ 0) # old and new capacity for generators
        @variable(cep.model, AccumulatedNewCapacity[r ∈ 𝓡 ,g ∈ sets["nodes"], y ∈ 𝓨] ≥ 0) # accumulated capacity according to lifetime for generators
        @variable(cep.model, NewCapacity[r ∈ 𝓡 ,g ∈ sets["nodes"], y ∈ 𝓨] ≥ 0)           # new capacity investments for generators      
        @variable(cep.model, capex[y ∈ 𝓨, g ∈ sets["nodes"]] ≥ 0)  # capital investment costs for generators
    end

    # generation variables
    @variable(cep.model, gen[r ∈ 𝓡 ,g ∈ 𝓖, y ∈ 𝓨, c ∈ cep.sets["carrier"][g], t ∈ 𝓣])  # planned generation for generators
    @variable(cep.model, ll[r ∈ 𝓡, y ∈ 𝓨, t ∈ 𝓣, c ∈ ["electricity", "H2"]] ≥ 0)   # lost load / ENS
    #@variable(cep.model, ll_h2[r ∈ 𝓡, y ∈ 𝓨] ≥ 0)   # lost load / ENS

    @variable(cep.model, em[y ∈ 𝓨] ≥ 0)      # emission CO2 per year ##curtail,emt??

    # cost variables
    @variable(cep.model, cll[y ∈ 𝓨] ≥ 0)  # costs for lost load yearly
    @variable(cep.model, opex[y ∈ 𝓨, g ∈ 𝓖] ≥ 0) 

    return cep
end



function set_up_equations(; cep::OptModelCEP, 
    ts_data::JuMP.Containers.DenseAxisArray, 
    data::OptDataCEP, 
    config::Dict{Any, Any}, 
    kwargs...)

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)
    data = data.data

    ## how to handle different fuels
    emitting_fuels = [g for g ∈ 𝓖 if data["emission"][g] > 0]

    # energy balance equation for each energy carrier
    @constraint(cep.model, EnergyBalance[r ∈ 𝓡, y ∈ 𝓨, t ∈ 𝓣, c ∈ 𝓒], 
    sum(cep.model[:gen][r,g,y,c,t] for g ∈ cep.sets[c]) 
    + (c ∈ ["H2", "electricity"] ? cep.model[:ll][r,y,t,c] : 0)
    - (c == "H2" ? (data["demand"][r,y,"H2"]/8760) : 0)
    - (c == "electricity" ? (ts_data[r,"Demand",t] * data["demand"][r,y,"electricity"]) : 0)
    == 0)
       
    # emission accounting
    @constraint(cep.model, EM[y ∈ 𝓨],cep.model[:em][y] == sum(cep.model[:gen][r,g,y,c,t] * data["emission"][g] for r ∈ 𝓡, g ∈ emitting_fuels, c ∈ cep.sets["carrier"][g], t ∈ 𝓣))

    # cost for lost load yearly 
    @constraint(cep.model, CLL[y ∈ 𝓨], cep.model[:cll][y] == config["cll"] * (sum(cep.model[:ll][r,y,t,c] for r ∈ 𝓡, t ∈ 𝓣, c ∈ ["H2", "electricity"])))

    # limit max and min generation dispatchable and non dispatchable
    @constraint(cep.model, GenCapDisp[r ∈ 𝓡, y ∈ 𝓨, g ∈ cep.sets["dispatch"], c ∈ cep.sets["carrier"][g], t ∈ 𝓣], cep.model[:gen][r,g,y,c,t] ≤ (config["dispatch"] ? data["cap_init"][r,g,y] : cep.model[:TotalCapacityAnnual][r,g,y]) * data["eta"][g,y])   
    @constraint(cep.model, GenCapNonDisp[r ∈ 𝓡, y ∈ 𝓨, g ∈ cep.sets["non_dispatch"], c ∈ cep.sets["carrier"][g], t ∈ 𝓣], cep.model[:gen][r,g,y,c,t] ≤ (config["dispatch"] ? data["cap_init"][r,g,y] : cep.model[:TotalCapacityAnnual][r,g,y]) * data["eta"][g,y] * ts_data[r,g,t])
    @constraint(cep.model, [r ∈ 𝓡, y ∈ 𝓨, g ∈ vcat(cep.sets["non_dispatch"],cep.sets["dispatch"]), c ∈ cep.sets["carrier"][g], t ∈ 𝓣],  0 ≤ cep.model[:gen][r,g,y,c,t])

    setup_opt_opex!(cep, config, data, vcat(cep.sets["non_dispatch"], cep.sets["dispatch"]), 1)

    if !config["dispatch"]
        # no investments in 2020
        JuMP.fix.(cep.model[:AccumulatedNewCapacity][:, :, 𝓨[1]], 0; force=true)
        JuMP.fix.(cep.model[:NewCapacity][:, :, 𝓨[1]], 0; force=true)

        setup_opt_capex!(cep, config, vcat(cep.sets["non_dispatch"], cep.sets["dispatch"]))

        # new capacity investments 
        @constraint(cep.model, NewCap[r ∈ 𝓡, g ∈ 𝓖, y ∈ 𝓨], cep.model[:TotalCapacityAnnual][r,g,y] == cep.model[:AccumulatedNewCapacity][r,g,y] + cep.model[:NewCapacity][r,g,y] + data["cap_init"][r,g,y])    
        # accumulated capacity
        @constraint(cep.model, AccCap[r ∈ 𝓡, g ∈ 𝓖, y in 𝓨[2:end]], cep.model[:AccumulatedNewCapacity][r,g,y] == sum(cep.model[:NewCapacity][r,g,hat_y] for hat_y in 𝓨[1]:10:y if y - 𝓨[1] ≤ data["lifetime"][g]))
        # max potential capacity constraint
        for r ∈ 𝓡, g ∈ 𝓖, y ∈ 𝓨[2:end]
            if data["cap"][r, g, y] > 0
                @constraint(cep.model, MaxCap[r, g, y], cep.model[:TotalCapacityAnnual][r, g, y] ≤ data["cap"][r, g, y])
            end
        end
        
        @constraint(cep.model, EM_zero[𝓨[end]], cep.model[:em][𝓨[end]] == 0)
    else
        # emission budget for each country individually
        @constraint(cep.model, EM_budget[y in 𝓨], cep.model[:em][y] ≤ sum(data["budget"][r,y] for r ∈ 𝓡))
    end
end


function setup_opt_storage!(cep::OptModelCEP, 
    config::Dict{Any, Any}, 
    data::OptDataCEP)

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)
    data = data.data

    @variable(cep.model, StorageLevel[r ∈ 𝓡, s ∈ 𝓢 , y ∈ 𝓨, t ∈ 𝓣] ≥ 0)

    # Set storage level at beginning and end of year equal
    @constraint(cep.model, SoC_Beginning[r ∈ 𝓡, s ∈ 𝓢 , y ∈ 𝓨], cep.model[:StorageLevel][r,s,y,𝓣[end]] == (config["dispatch"] ? data["cap_init"][r,s,y] : cep.model[:TotalCapacityAnnual][r,s,y]) * config["techs"][s]["constraints"]["SOC_Start"])

    # charging Soc according to max storage level 
    @constraint(cep.model, SoC[r ∈ 𝓡, s ∈ 𝓢 , y ∈ 𝓨, t ∈ 𝓣], cep.model[:StorageLevel][r,s,y,t] ≤ (config["dispatch"] ? data["cap_init"][r,s,y] : cep.model[:TotalCapacityAnnual][r,s,y]))
    setup_opt_opex!(cep, config, data, 𝓢, 1)
    
    if !config["dispatch"]
        # limit investments p/e ratio
            # Connect the previous storage level with the new storage level
        @constraint(cep.model, SoC_Balance[r ∈ 𝓡, s ∈ 𝓢 , y ∈ 𝓨, t ∈ 𝓣], 
        (t > 1 ? cep.model[:StorageLevel][r,s,y,t-1] : cep.model[:TotalCapacityAnnual][r,s,y] * config["techs"][s]["constraints"]["SOC_Start"])
        - cep.model[:gen][r,s,y,config["techs"][s]["input"]["carrier"],t]
        == cep.model[:StorageLevel][r,s,y,t] 
        )
        setup_opt_capex!(cep, config, 𝓢)

        @constraint(cep.model, P2E_ratio[r ∈ 𝓡, s ∈ 𝓢 , y ∈ 𝓨], cep.model[:TotalCapacityAnnual][r,"$(replace(s, "S_" => "D_"))_in",y] * config["techs"][s]["constraints"]["P2E"]  ≤ cep.model[:TotalCapacityAnnual][r,s,y])
        # charging and discharging investments are the same
        @constraint(cep.model, Discharg_Charge[r ∈ 𝓡, s ∈ 𝓢 , y ∈ 𝓨], cep.model[:TotalCapacityAnnual][r,"$(replace(s, "S_" => "D_"))_in",y] == cep.model[:TotalCapacityAnnual][r,"$(replace(s, "S_" => "D_"))_out",y])
    
    else
        @constraint(cep.model, SoC_Balance[r ∈ 𝓡, s ∈ 𝓢 , y ∈ 𝓨, t ∈ 𝓣], 
        (t > 1 ? cep.model[:StorageLevel][r,s,y,t-1] : data["cap_init"][r,s,y] * config["techs"][s]["constraints"]["SOC_Start"])
        - cep.model[:gen][r,s,y,config["techs"][s]["input"]["carrier"],t]
        == cep.model[:StorageLevel][r,s,y,t] 
        )   
    
    end
    return cep
end


"""
     setup_opt_conversion!((cep::OptModelCEP, config::Dict{Any, Any})

A conversion technology converts the input carrier to an output carrier with a certain efficiency.
"""

function setup_opt_conversion!(cep::OptModelCEP, 
    config::Dict{Any, Any},
    data::OptDataCEP)   

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)
    data = data.data

    # Calculate the input generation 
    @constraint(cep.model, InputConversion[r ∈ 𝓡, y ∈ 𝓨, g ∈ cep.sets["conversion"], t ∈ 𝓣], cep.model[:gen][r, g, y, config["techs"][g]["input"]["carrier"], t] ≥ (-1) * (config["dispatch"] ? data["cap_init"][r, g, y] : cep.model[:TotalCapacityAnnual][r, g, y]))
    @constraint(cep.model, InputConversion2[r ∈ 𝓡, y ∈ 𝓨, g ∈ cep.sets["conversion"], t ∈ 𝓣], cep.model[:gen][r, g, y, config["techs"][g]["input"]["carrier"], t] ≤ 0)
    
    # Calculate the output generation
    @constraint(cep.model, Outputconversion[r ∈ 𝓡, y ∈ 𝓨, g ∈ cep.sets["conversion"], t ∈ 𝓣], cep.model[:gen][r,g,y,config["techs"][g]["output"]["carrier"],t] ==  (-1) * cep.model[:gen][r,g,y,config["techs"][g]["input"]["carrier"],t] * data["eta"][g,y])

    # add the costs 
    setup_opt_opex!(cep, config, data, cep.sets["conversion"], -1)
    if !config["dispatch"]
        setup_opt_capex!(cep, config, cep.sets["conversion"])
    end

    return cep
end 





function set_opt_transmission!(cep::OptModelCEP, 
    config::Dict{Any, Any},
    data::OptDataCEP)

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)
    lines = data.lines
    data = data.data

    ## VARIABLE ##
    @variable(cep.model, FLOW[g ∈ cep.sets["transmission"], l ∈ 𝓛, dir ∈ ["uniform", "opposite"], y ∈ 𝓨, t ∈ 𝓣] >= 0)

    if !config["dispatch"]
        @variable(cep.model, NewTradeCapacityCosts[g ∈ cep.sets["transmission"], y ∈ 𝓨]  >= 0)
        @variable(cep.model, NewTradeCapacity[g ∈ cep.sets["transmission"], l ∈ 𝓛, y ∈ 𝓨]  >= 0)
        @variable(cep.model, TotalTradeCapacity[g ∈ cep.sets["transmission"], l ∈ 𝓛, y ∈ 𝓨]  >= 0)
    end

    ## TRANSMISSION TRANS ##
    @constraint(cep.model, FlowLimit[g ∈ cep.sets["transmission"], l ∈ 𝓛, dir ∈ ["uniform", "opposite"], y ∈ 𝓨, t ∈ 𝓣], cep.model[:FLOW][g,l,dir,y,t] ≤ (config["dispatch"] ? lines[(g, l)].power_lim : TotalTradeCapacity[g,l,y]))

    if !config["dispatch"]
        @constraint(cep.model, ExistingTransmCapa[g ∈ cep.sets["transmission"], l ∈ 𝓛], TotalTradeCapacity[g,l,𝓨[1]] == (c == "electricity" ? lines[(g, l)].power_lim : 0))  
        @constraint(cep.model, TransmissionExpansion[g ∈ cep.sets["transmission"], l ∈ 𝓛, i ∈ eachindex(𝓨)[2:end]], TotalTradeCapacity[g,l,𝓨[i]] == NewTradeCapacity[g,l,𝓨[i]] + TotalTradeCapacity[g,l,𝓨[i-1]])
        @constraint(cep.model, NewTradeCapacityCosts[g ∈ cep.sets["transmission"], y ∈ 𝓨[2:end]], NewTradeCapacityCosts[g,y] == sum(NewTradeCapacity[g,l,y] * lines[(g, l)].length * config["techs"][g]["investment_costs"] for l ∈ 𝓛))
        JuMP.fix.(cep.model[:NewTradeCapacity][:, :, 𝓨[1]], 0; force=true)
        JuMP.fix.(cep.model[:NewTradeCapacityCosts][:, 𝓨[1]], 0; force=true)
    end
    
    @constraint(cep.model, Nettrade[r ∈ 𝓡, g ∈ cep.sets["transmission"], y ∈ 𝓨, t ∈ 𝓣, c ∈ cep.sets["carrier"][g]], 
    cep.model[:gen][r,g,y,c,t] 
    == sum(cep.model[:FLOW][g, line_end,"uniform",y,t] - cep.model[:FLOW][g,line_end,"opposite",y,t]/lines[(g,line_end)].eff for line_end ∈ [l for ((t, l), v) ∈ lines if t == g && v.node_end == r]) + 
    sum(cep.model[:FLOW][g,line_start,"opposite",y,t] - cep.model[:FLOW][g,line_start,"uniform",y,t]/lines[(g,line_start)].eff for line_start ∈ [l for ((t,l), v) ∈ lines if t == g && v.node_start == r]))
    
    return cep
end


"""
     setup_opt_objective!(cep::OptModelCEP, config::Dict{Any, Any})
Calculate total system costs and set as objective
"""
function setup_opt_objective!(cep::OptModelCEP, 
    config::Dict{Any, Any})
    ## OBJECTIVE ##
    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)

    if !config["dispatch"]
        @objective(cep.model, Min, sum(1/((1+config["r"])^(y-𝓨[1]))*
        (cep.model[:capex][y,g] for g ∈ setdiff(𝓖, cep.sets["transmission"]) 
        + sum(NewTradeCapacityCosts[g,y] for g ∈ cep.sets["transmission"])) for y ∈ 𝓨) 
        + sum(1/((1+config["r"])^(y-𝓨[1])) * sum(sum(cep.model[:opex][y,g] for g ∈ 𝓖) + cep.model[:cll][y] for y ∈ 𝓨)))  
    else            
        @objective(cep.model, Min, sum(sum(cep.model[:opex][y,g] for g ∈ 𝓖) + cep.model[:cll][y] for y ∈ 𝓨))  
    end
  return cep
end



"""
     setup_opt_costs!(cep::OptModelCEP, config::Dict{Any, Any}, tech_group::String; sign_generation::Number=1)
add operational costs for the technology defined by `tech_group`
"""

function setup_opt_capex!(cep::OptModelCEP, 
    config::Dict{Any, Any}, 
    tech_group::String)

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)
    @constraint(cep.model, [y ∈ 𝓨, g ∈ tech_group], sum(cep.model[:NewCapacity][r,g,y] for r ∈ 𝓡) * data["c_CAPEX"][g,y] == cep.model[:capex][y,g]) 

    return cep
end


"""
     setup_opt_costs!(cep::OptModelCEP, config::Dict{Any, Any}, tech_group::String; sign_generation::Number=1)
add capacity costs for the technology defined by `tech_group`
"""

function setup_opt_opex!(cep::OptModelCEP, 
    config::Dict{Any, Any}, 
    data::Dict{Any, Any},
    tech_group::Vector{String}, 
    sign_generation::Int64)

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)
   
    # fixed and variable costs for operation
    @constraint(cep.model, [y ∈ 𝓨, g ∈ tech_group], cep.model[:opex][y,g] == sum(config["dispatch"] ? data["cap"][r,g,y] : cep.model[:TotalCapacityAnnual][r,g,y] for r ∈ 𝓡) * data["c_fix"][g, y] 
    + sign_generation * (sum(sum(cep.model[:gen][r,g,y,c,t] for r ∈ 𝓡, t ∈ 𝓣) * data["c_var"][g, y] for c ∈ cep.sets["carrier"][g])))

    return cep
end





function optimize_and_output(; cep::OptModelCEP,
    config::Dict{Any, Any}, 
    data::OptDataCEP, 
    ts_data)

    optimize!(cep.model)
    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)

    status=Symbol(termination_status(cep.model))
    println(status)
    if termination_status(cep.model) == MOI.INFEASIBLE_OR_UNBOUNDED || termination_status(cep.model) == MOI.INFEASIBLE
        JuMP.compute_conflict!(cep.model)
        list_of_conflicting_constraints = ConstraintRef[]
        for (F, S) ∈ list_of_constraint_types(cep.model)
            for con ∈ all_constraints(cep.model, F, S)
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
        # generation 
        #variables["gen"] = convert_jump_container_to_df(cep=cep, config=config)

        #plotgen(cep, config, 2030, data, ts_data)

        open(joinpath(pwd(),"P2H_CapacityExpansion","results", "solution_full.txt"), "w") do file
            println(file, "Objective = $objective")
            for v ∈ all_variables(cep.model)
                if value.(v) != 0
                    val = value.(v)
                    str = string(v)
                    variables[str] = val
                    println(file, "$str = $val")
                end
            end

            for r ∈ axes(ts_data)[1], t ∈ axes(ts_data)[3] 
                str = "Demand[$r,$t,$(config["year"])]"
                val = ts_data[r,"Demand",t] * data.data["demand"][r,config["year"],"electricity"]
                println(file, "$str = $val")
            end
            
            for r ∈ axes(data.data["cap_init"])[1], g ∈ axes(data.data["cap_init"])[2], y ∈ axes(data.data["cap_init"])[3] 
                val = data.data["cap_init"][r,g,y]
                println(file, "Capacity$r$g$y = $val") 
            end

            ## write the dual variables aka shadow prices 
            for r ∈ 𝓡, y ∈ 𝓨, t ∈ 𝓣, c ∈ 𝓒
                val = dual(cep.model[:EnergyBalance][r, y, t, c])
                println(file, "Price[$r,$y,$t,$c] = $val")
            end

        end
        return OptResult(cep.model, status, objective, variables)
    end  
end




