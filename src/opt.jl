"""
run_opt(ts_data::ClustData,opt_data::OptDataCEP,config::Dict{String,Any},optimizer::DataType)
Organizing the actual setup and run of the CEP-Problem.
Required elements are:
- `ts_data`: The time-series data.
- `opt_data`: In this case the OptDataCEP that contains information on costs, nodes, techs and for transmission also on lines.
- `config`: This includes all the settings for the design optimization problem formulation.
"""


function run_opt(; ts_data::ClustData,
    data::OptDataCEP,
    config::Dict{Any, Any},
    surrogate::Bool,
    solver::DataType, 
    kwargs...
    )

    @info "Reading the data ..."
    cep=setup_opt_basic(ts_data=ts_data, config=config, solver=solver, data=data)

    setup_opt_inv_variables!(cep, config) 
    setup_opt_investment!(cep, ts_data, data, config)

    if !surrogate
        setup_opt_op_variables!(cep, config) 
        setup_opt_operation!(cep, ts_data, data, config)
        setup_opt_storage!(cep, ts_data, config, data)
        setup_opt_conversion!(cep, config, ts_data, data)    
        set_opt_transmission!(cep, config, ts_data, data)
        setup_opt_objective!(cep, config)
    end

    return cep
end



function setup_opt_basic(;ts_data::ClustData, 
    config::Dict{Any, Any},
    solver::DataType,
    data::OptDataCEP)
    
    # Initialize model
    model =  JuMP.Model()
    set_optimizer(model, solver)
    # Initialize solver
    optimizer = optimizer_with_attributes(
        () -> Ipopt.Optimizer(Ipopt.Env()),
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





"""
    setup_opt_inv_variables(; cep::OptModelCEP, config::Dict{Any, Any}) -> OptModelCEP

Sets up the basic investment-related decision variables in the optimization model.

# Arguments
- `cep::OptModelCEP`: The optimization model container, which includes sets and the JuMP model.
- `config::Dict{Any, Any}`: Configuration dictionary that must include a `"dispatch"` key to determine model mode.

# Returns
- `OptModelCEP`: The updated model container with investment variables added.
"""

function setup_opt_inv_variables!(cep::OptModelCEP,  config::Dict{Any, Any})

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)

    if !config["dispatch"]
        @info "Investment mode is on."
        # capacity variables
        @variable(cep.model, TotalCapacityAnnual[r ∈ 𝓡 ,g ∈ cep.sets["nodes"], y ∈ 𝓨] ≥ 0) # old and new capacity for generators
        @variable(cep.model, AccumulatedNewCapacity[r ∈ 𝓡 ,g ∈ cep.sets["invest_tech"], y ∈ 𝓨] ≥ 0) # accumulated capacity according to lifetime for generators
        @variable(cep.model, NewCapacity[r ∈ 𝓡 ,g ∈ cep.sets["invest_tech"], y ∈ 𝓨] ≥ 0)           # new capacity investments for generators      
        @variable(cep.model, COST[z ∈ ["cap", "fix", "var"], y ∈ 𝓨, g ∈ 𝓖] ≥ 0) 
    else
        @variable(cep.model, COST[z ∈ ["fix", "var"], y ∈ 𝓨, g ∈ 𝓖] ≥ 0) 
    end
    return cep
end





"""
    setup_opt_op_variables(; cep::OptModelCEP, config::Dict{Any, Any}) -> OptModelCEP

Defines the basic operational decision variables in the optimization model.

# Arguments
- `cep::OptModelCEP`: The optimization model container, which includes sets and the JuMP model.
- `config::Dict{Any, Any}`: Configuration dictionary (not directly used here but kept for consistency).

# Returns
- `OptModelCEP`: The updated model container with operational variables added.
"""

function setup_opt_op_variables!(cep::OptModelCEP,  config::Dict{Any, Any})

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)

    # generation variables
    @variable(cep.model, gen[r ∈ 𝓡 , g ∈ setdiff(𝓖, cep.sets["storage_techs"]), y ∈ 𝓨, c ∈ cep.sets["carrier"][g], t ∈ 𝓣])  # planned generation for generators

    @variable(cep.model, em[y ∈ 𝓨] ≥ 0)      # emission CO2 per year 
    return cep
end



"""
    setup_opt_investment!(; cep::OptModelCEP, ts_data::ClustData, data::OptDataCEP, config::Dict{Any, Any}) -> OptModelCEP

Adds investment-related constraints to the capacity expansion model.

# Arguments
- `cep::OptModelCEP`: The optimization model container that includes sets, the JuMP model, and variable definitions.
- `ts_data::ClustData`: Time series clustering data (not directly used here but may be needed for downstream compatibility).
- `data::OptDataCEP`: Contains all scenario-specific input data such as initial capacity, lifetimes, and maximum potentials.
- `config::Dict{Any, Any}`: Configuration dictionary, passed to subfunctions (e.g., for cost setup).

# Returns
- `OptModelCEP`: The updated model container with investment constraints added.
"""


function setup_opt_investment!(cep::OptModelCEP, 
    ts_data::ClustData, 
    data::OptDataCEP, 
    config::Dict{Any, Any})

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)
    lines = data.lines
    data = data.data

    ### STORAGE AND GENERATOR CAPACITY EXPANSION ###
    # fix generation capacity where no investments are allowed to the base year
    @constraint(cep.model, NoInvestments[r ∈ 𝓡, y ∈ 𝓨, g ∈ setdiff(cep.sets["nodes"], cep.sets["invest_tech"])], cep.model[:TotalCapacityAnnual][r,g,y] == data["cap_init"][r,g,y])
    setup_opt_costs_cap!(cep, config, data, cep.sets["invest_tech"])

    # new capacity investments 
    @constraint(cep.model, NewCap[r ∈ 𝓡, g ∈ cep.sets["invest_tech"], y ∈ 𝓨], cep.model[:TotalCapacityAnnual][r,g,y] == cep.model[:AccumulatedNewCapacity][r,g,y] + data["cap_init"][r,g,y])    
    # accumulated capacity
    @constraint(cep.model, AccCap[r ∈ 𝓡, g ∈ cep.sets["invest_tech"], y in 𝓨], cep.model[:AccumulatedNewCapacity][r,g,y] == sum(cep.model[:NewCapacity][r,g,yy] for yy ∈ 𝓨 if (y - yy < data["lifetime"][g]) && (y-yy >= 0)))

    # max potential capacity constraint
    for r ∈ 𝓡, g ∈ cep.sets["invest_tech"], y ∈ 𝓨
        if data["cap"][r, g, y] > 0
            @constraint(cep.model, cep.model[:TotalCapacityAnnual][r, g, y] ≤ data["cap"][r, g, y])
        end
    end


    ### STORAGE E2P RATIO INVESTMENT RATIO ###
    @constraint(cep.model, P2E_ratio[r ∈ 𝓡, s ∈ intersect(𝓢,cep.sets["invest_tech"]) , y ∈ 𝓨], cep.model[:TotalCapacityAnnual][r,"$(replace(s, "S_" => "D_"))_out",y] * config["techs"][s]["constraints"]["P2E"]  ≤ cep.model[:TotalCapacityAnnual][r,s,y])


    ### TRANSMISSION LINE CAPACITY EXPANSION ####
    @variable(cep.model, NewTradeCapacity[g ∈ cep.sets["transmission"], l ∈ 𝓛, y ∈ 𝓨]  >= 0)
    @variable(cep.model, TotalTradeCapacity[g ∈ cep.sets["transmission"], l ∈ 𝓛, y ∈ 𝓨]  >= 0)

    @constraint(cep.model, ExistingTransmCapa[g ∈ cep.sets["transmission"], l ∈ 𝓛], TotalTradeCapacity[g,l,𝓨[1]] == lines[(g, l)].power_lim)  
    @constraint(cep.model, TransmissionExpansion[g ∈ cep.sets["transmission"], l ∈ 𝓛, i ∈ eachindex(𝓨)], TotalTradeCapacity[g,l,𝓨[i]] == NewTradeCapacity[g,l,𝓨[i]] + TotalTradeCapacity[g,l,𝓨[i-1]])    
    @constraint(cep.model, NewTradeCapacityCosts[g ∈ cep.sets["transmission"], y ∈ 𝓨], cep.model[:COST]["cap",y,g] == sum(NewTradeCapacity[g,l,y] * lines[(g, l)].length * config["techs"][g]["investment_costs"] for l ∈ 𝓛))

    return cep
end



"""
    setup_opt_operation!(; cep::OptModelCEP, ts_data::ClustData, data::OptDataCEP, config::Dict{Any, Any}) -> OptModelCEP

Adds operational constraints to the capacity expansion or dispatch model.

# Arguments
- `cep::OptModelCEP`: The optimization model object containing sets, model variables, and JuMP model.
- `ts_data::ClustData`: Clustered time series data, including weights and normalized profiles.
- `data::OptDataCEP`: Input data wrapper, including capacity, demand, efficiency, and emissions.
- `config::Dict{Any, Any}`: Configuration dictionary. The key `"dispatch"` toggles between investment and dispatch mode. Other keys include `"cll"` for the cost of lost load.

# Returns
- `OptModelCEP`: The updated model with all operational constraints defined.
"""


function setup_opt_operation!(; cep::OptModelCEP, 
    ts_data::ClustData, 
    data::OptDataCEP, 
    config::Dict{Any, Any}, )

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)
    data = data.data

    @constraint(cep.model,
    EnergyBalance[r ∈ 𝓡, y ∈ 𝓨, t ∈ 𝓣, c ∈ 𝓒],
    sum(cep.model[:gen][r,g,y,c,t] for g ∈ setdiff(cep.sets[c], cep.sets["storage_techs"])) * ts_data.weight[t]  
    == data["demand"][r,y,c] * (c == "H2" ? (1/8760) : ts_data.ts[r,"Demand",t]) * ts_data.weight[t]) 

    ## how to handle different fuels
    emitting_fuels = [g for g ∈ 𝓖 if data["emission"][g] > 0]

    # LOST LOAD ELECTRICITY ONLY
    @constraint(cep.model, [r ∈ 𝓡, y ∈ 𝓨, g ∈ ["ENS"], c ∈ cep.sets["carrier"][g], t ∈ 𝓣],  0 ≤ cep.model[:gen][r,g,y,c,t])
    @constraint(cep.model, [y ∈ 𝓨, g ∈ ["ENS"]], cep.model[:COST]["var",y,g] == sum(cep.model[:gen][r,g,y,c,t] * ts_data.weight[t] for r ∈ 𝓡, t ∈ 𝓣, c ∈ cep.sets["carrier"][g]) * config["cll"])


    # limit max and min generation dispatchable and non dispatchable
    @constraint(cep.model, GenCapDisp[r ∈ 𝓡, y ∈ 𝓨, g ∈ cep.sets["dispatch"], c ∈ cep.sets["carrier"][g], t ∈ 𝓣], cep.model[:gen][r,g,y,c,t] ≤ (config["dispatch"] ? data["cap_init"][r,g,y] : cep.model[:TotalCapacityAnnual][r,g,y]) * data["eta"][g,y])   
    @constraint(cep.model, GenCapNonDisp[r ∈ 𝓡, y ∈ 𝓨, g ∈ cep.sets["non_dispatch"], c ∈ cep.sets["carrier"][g], t ∈ 𝓣], cep.model[:gen][r,g,y,c,t] ≤ (config["dispatch"] ? data["cap_init"][r,g,y] : cep.model[:TotalCapacityAnnual][r,g,y]) * data["eta"][g,y] * ts_data.ts[r,g,t])
    @constraint(cep.model, [r ∈ 𝓡, y ∈ 𝓨, g ∈ vcat(cep.sets["non_dispatch"],cep.sets["dispatch"]), c ∈ cep.sets["carrier"][g], t ∈ 𝓣],  0 ≤ cep.model[:gen][r,g,y,c,t])

    setup_opt_costs_var!(cep, config, data, ts_data, vcat(cep.sets["non_dispatch"], cep.sets["dispatch"]), 1)
    setup_opt_costs_fix!(cep, config, data, vcat(cep.sets["non_dispatch"], cep.sets["dispatch"]))

    # EMISSION ACCOUNTING AND BUDGETS
    @constraint(cep.model, EM[y ∈ 𝓨], cep.model[:em][y] == sum(cep.model[:gen][r,g,y,c,t] * ts_data.weight[t] * data["emission"][g] for r ∈ 𝓡, g ∈ emitting_fuels, c ∈ cep.sets["carrier"][g], t ∈ 𝓣))

    if !config["dispatch"]
        @constraint(cep.model, EmissionBudget[y ∈ 𝓨], cep.model[:em][y] ≤ sum(data["budget"][r,y] for r ∈ 𝓡))
    end
    return cep
end


"""
    setup_opt_storage!(cep::OptModelCEP, ts_data::ClustData, config::Dict{Any, Any}, data::OptDataCEP) -> OptModelCEP

Adds storage-related variables and constraints to the energy system optimization model.

# Arguments
- `cep::OptModelCEP`: The optimization model container with sets, JuMP model, and defined variables.
- `ts_data::ClustData`: Contains clustered time-series data, such as inflows, time weights, and profiles.
- `config::Dict{Any, Any}`: Configuration dictionary with dispatch mode toggle and technology parameters.
- `data::OptDataCEP`: Model input data including capacity, efficiencies, and initial storage values.

# Returns
- `OptModelCEP`: The updated model with storage components added.
"""


function setup_opt_storage!(cep::OptModelCEP, 
    ts_data::ClustData,  
    config::Dict{Any, Any}, 
    data::OptDataCEP)

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)
    data = data.data

    ## STORAGE LEVEL VARIABLE ##
    @variable(cep.model, SOC[r ∈ 𝓡, s ∈ 𝓢 , y ∈ 𝓨, t ∈ 𝓣] ≥ 0)


    ## STORAGE LEVEL CONSTRAINTS ##
    # Set storage level at beginning and end of year equal
    @constraint(cep.model, SoC_Beginning[r ∈ 𝓡, s ∈ 𝓢 , y ∈ 𝓨], cep.model[:SOC][r,s,y,𝓣[end]] == (config["dispatch"] ? data["cap_init"][r,s,y] : cep.model[:TotalCapacityAnnual][r,s,y]) * config["techs"][s]["constraints"]["SOC_Start"])

    # Soc according to max storage level 
    ## add a spilling variable 
    @variable(cep.model, SPIL[r ∈ 𝓡, s ∈ 𝓢 , y ∈ 𝓨, t ∈ 𝓣] ≥ 0)

    @constraint(cep.model, SoC[r ∈ 𝓡, s ∈ 𝓢 , y ∈ 𝓨, t ∈ 𝓣], cep.model[:SOC][r,s,y,t] ≤ (config["dispatch"] ? data["cap_init"][r,s,y] : cep.model[:TotalCapacityAnnual][r,s,y]))

    # define fixed costs only
    setup_opt_costs_fix!(cep, config, data, 𝓢)
    
    # define charging and discharging limits
    setup_opt_storage_flows!(cep, config, ts_data, data)

    # storage filling
    for r ∈ 𝓡, s ∈ 𝓢, y ∈ 𝓨, c ∈ cep.sets["carrier"][s]
        soc_start = config["dispatch"] ? data["cap_init"][r,s,y] : cep.model[:TotalCapacityAnnual][r,s,y]
        for t ∈ 𝓣
            @constraint(cep.model, 
            (t > 1 ? cep.model[:SOC][r,s,y,t-1] : soc_start * config["techs"][s]["constraints"]["SOC_Start"]) +
            (("$(replace(s, "S_" => "D_"))_in" ∈ keys(config["techs"]) ? ((-1)* cep.model[:gen][r,"$(replace(s, "S_" => "D_"))_in",y,c,t] * ts_data.weight[t] * data["eta"]["$(replace(s, "S_" => "D_"))_in",y]) :  (ts_data.ts[r,"inflow",t] * ts_data.weight[t])) 
            - cep.model[:gen][r,"$(replace(s, "S_" => "D_"))_out",y,c,t] * ts_data.weight[t])
            - cep.model[:SPIL][r,s,y,t]
            == cep.model[:SOC][r,s,y,t],
            base_name="SoC_Balance$r,$s,$y,$t" 
            )   
        end
    end
    return cep
end



"""
    setup_opt_storage_flows!(cep::OptModelCEP, config::Dict{Any, Any}, data::OptDataCEP)

Adds constraints and cost components to model the energy flows into and out of storage technologies 

Returns the modified `cep` model with added constraints and cost terms.
"""
function setup_opt_storage_flows!(cep::OptModelCEP, 
    config::Dict{Any, Any}, 
    ts_data::ClustData,
    data::Dict{Any, Any},)

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)

    # charging
    @constraint(cep.model, MaxCharging[r ∈ 𝓡, y ∈ 𝓨, g ∈ cep.sets["charging"], c ∈ cep.sets["carrier"][g], t ∈ 𝓣], cep.model[:gen][r,g,y,c,t] ≥ (-1) * (config["dispatch"] ? data["cap_init"][r, "$(replace(g, "in" => "out"))", y] : cep.model[:TotalCapacityAnnual][r,"$(replace(g, "in" => "out"))", y]))
    @constraint(cep.model, MinCharging[r ∈ 𝓡, y ∈ 𝓨, g ∈ cep.sets["charging"], c ∈ cep.sets["carrier"][g], t ∈ 𝓣], cep.model[:gen][r, g, y, c, t] ≤ 0)
    
    # discharging tie it to charging capacity to avoid double costs
    @constraint(cep.model, MaxDischarging[r ∈ 𝓡, y ∈ 𝓨, g ∈ cep.sets["discharging"], t ∈ 𝓣, c ∈ cep.sets["carrier"][g]], cep.model[:gen][r,g,y,c,t] ≤ (config["dispatch"] ? data["cap_init"][r,g,y] : cep.model[:TotalCapacityAnnual][r, g, y]) * data["eta"][g,y])
    @constraint(cep.model, MinDischarging[r ∈ 𝓡, y ∈ 𝓨, g ∈ cep.sets["discharging"], t ∈ 𝓣, c ∈ cep.sets["carrier"][g]], cep.model[:gen][r,g,y,c,t] ≥ 0)

    ## add the costs for charging once
    setup_opt_costs_fix!(cep, config, data, String[s for s in cep.sets["discharging"]])
    setup_opt_costs_var!(cep, config, data, ts_data, String[s for s in cep.sets["discharging"]], 1)
    setup_opt_costs_var!(cep, config, data, ts_data, String[s for s in cep.sets["charging"]], -1)
    #JuMP.fix.(cep.model[:COST][:, :, g ∈ cep.sets["discharging"]], 0; force=true)

    return cep
end





"""
     setup_opt_conversion!((cep::OptModelCEP, config::Dict{Any, Any})

A conversion technology converts the input carrier to an output carrier with a certain efficiency.
"""

function setup_opt_conversion!(cep::OptModelCEP, 
    config::Dict{Any, Any},
    ts_data::ClustData,
    data::OptDataCEP)   

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)
    data = data.data

    # Calculate the input generation 
    @constraint(cep.model, InputConversion[r ∈ 𝓡, y ∈ 𝓨, g ∈ cep.sets["conversion"], t ∈ 𝓣], cep.model[:gen][r, g, y, config["techs"][g]["input"]["carrier"], t] ≥ (-1) * (config["dispatch"] ? data["cap_init"][r, g, y] : cep.model[:TotalCapacityAnnual][r, g, y]))
    @constraint(cep.model, InputConversion2[r ∈ 𝓡, y ∈ 𝓨, g ∈ cep.sets["conversion"], t ∈ 𝓣], cep.model[:gen][r, g, y, config["techs"][g]["input"]["carrier"], t] ≤ 0)
    
    # Calculate the output generation
    @constraint(cep.model, Outputconversion[r ∈ 𝓡, y ∈ 𝓨, g ∈ cep.sets["conversion"], t ∈ 𝓣], cep.model[:gen][r,g,y,config["techs"][g]["output"]["carrier"],t] ==  (-1) * cep.model[:gen][r,g,y,config["techs"][g]["input"]["carrier"],t] * data["eta"][g,y])

    # add the costs 
    setup_opt_costs_fix!(cep, config, data, cep.sets["conversion"])
    setup_opt_costs_var!(cep, config, data, ts_data, cep.sets["conversion"], -1)

    return cep
end 


"""
    set_opt_transmission!(cep::OptModelCEP, config::Dict{Any, Any}, ts_data::ClustData, data::OptDataCEP) -> OptModelCEP

Adds variables and constraints to model power and carrier transmission flows across regions.

# Arguments
- `cep::OptModelCEP`: The optimization model container with sets and the JuMP model.
- `config::Dict{Any, Any}`: Configuration dictionary with keys such as `"dispatch"` that toggle mode-specific behavior.
- `ts_data::ClustData`: Time series clustering data with time weights and temporal resolution.
- `data::OptDataCEP`: Contains all input data including transmission line parameters and capacity.

# Returns
- `OptModelCEP`: The updated model including transmission flow handling.
"""


function set_opt_transmission!(cep::OptModelCEP, 
    config::Dict{Any, Any},
    ts_data::ClustData,
    data::OptDataCEP)

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)
    lines = data.lines
    data = data.data

    ## VARIABLE ##
    @variable(cep.model, FLOW[g ∈ cep.sets["transmission"], l ∈ 𝓛, dir ∈ ["uniform", "opposite"], y ∈ 𝓨, t ∈ 𝓣] >= 0)

    @constraint(cep.model, Nettrade[r ∈ 𝓡, g ∈ cep.sets["transmission"], y ∈ 𝓨, t ∈ 𝓣, c ∈ cep.sets["carrier"][g]], 
    cep.model[:gen][r,g,y,c,t] 
    == sum(cep.model[:FLOW][g, line_end,"uniform",y,t] - cep.model[:FLOW][g,line_end,"opposite",y,t]/lines[(g,line_end)].eff for line_end ∈ [l for ((t, l), v) ∈ lines if t == g && v.node_end == r]) + 
    sum(cep.model[:FLOW][g,line_start,"opposite",y,t] - cep.model[:FLOW][g,line_start,"uniform",y,t]/lines[(g,line_start)].eff for line_start ∈ [l for ((t,l), v) ∈ lines if t == g && v.node_start == r]))

    setup_opt_costs_var!(cep, config, data, ts_data, cep.sets["transmission"], 1)
    JuMP.fix.(cep.model[:COST]["fix",:,cep.sets["transmission"]], 0; force=true)

    ## TRANSMISSION TRANS ##
    @constraint(cep.model, FlowLimit[g ∈ cep.sets["transmission"], l ∈ 𝓛, dir ∈ ["uniform", "opposite"], y ∈ 𝓨, t ∈ 𝓣], cep.model[:FLOW][g,l,dir,y,t] ≤ (config["dispatch"] ? lines[(g, l)].power_lim : TotalTradeCapacity[g,l,y]))
        
    return cep
end




"""
     setup_opt_costs_cap!(cep::OptModelCEP, config::Dict{Any, Any}, tech_group::String)
add capital costs for the technology defined by `tech_group`
"""

function setup_opt_costs_cap!(cep::OptModelCEP, 
    config::Dict{Any, Any}, 
    data::Dict{Any, Any},
    tech_group::Vector{String})

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)

    @constraint(cep.model, [y ∈ 𝓨, g ∈ tech_group], sum(cep.model[:NewCapacity][r,g,y] for r ∈ 𝓡) * data["c_CAPEX"][g,y] == cep.model[:COST]["cap",y,g]) 

    return cep
end


"""
     setup_opt_costs_fix!(cep::OptModelCEP, config::Dict{Any, Any}, tech_group::String)
add fixed costs for the technology defined by `tech_group`
"""

function setup_opt_costs_fix!(cep::OptModelCEP, 
    config::Dict{Any, Any}, 
    data::Dict{Any, Any},
    tech_group::Vector{String}, 
    )

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)
   
    # fixed costs for operation
    @constraint(cep.model, [y ∈ 𝓨, g ∈ tech_group], cep.model[:COST]["fix", y,g] == sum(config["dispatch"] ? data["cap"][r,g,y] : cep.model[:TotalCapacityAnnual][r,g,y] for r ∈ 𝓡) * data["c_fix"][g, y])

    return cep
end



"""
     setup_opt_costs_var!(cep::OptModelCEP, config::Dict{Any, Any}, tech_group::String; sign_generation::Number=1)
add variable costs for the technology defined by `tech_group`
"""

function setup_opt_costs_var!(cep::OptModelCEP, 
    config::Dict{Any, Any}, 
    data::Dict{Any, Any},
    ts_data::ClustData,
    tech_group::Vector{String}, 
    sign_generation::Int64)

    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)
   
    # variable costs for operation
    @constraint(cep.model, [y ∈ 𝓨, g ∈ tech_group], cep.model[:COST]["var",y,g] == sign_generation * (sum(cep.model[:gen][r,g,y,c,t] * ts_data.weight[t] for r ∈ 𝓡, t ∈ 𝓣, c ∈ cep.sets["carrier"][g]) * data["c_var"][g, y]))

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

    opex_discounted = sum(
    1 / ((1 + config["r"])^(y - 𝓨[1] - 10)) * (
        sum(cep.model[:COST]["fix", y, g] for g ∈ 𝓖) +
        sum(cep.model[:COST]["var", y, g] for g ∈ setdiff(𝓖, cep.sets["storage_techs"])) +
        sum(cep.model[:COST]["var", y, g] for g ∈ cep.sets["ENS"]) + 
        sum(cep.model[:SPIL][r,s,y,t] for r ∈ 𝓡, s ∈ 𝓢, t ∈ 𝓣) * config["cll"]
        #cep.model[:cll][y]
    ) for y ∈ 𝓨)

    if !config["dispatch"]
        @objective(cep.model, Min, sum(
            1 / ((1 + config["r"])^(y - 𝓨[1] - 10)) *
            sum(cep.model[:COST]["cap", y, g] for g ∈ cep.sets["invest_all"])
            for y ∈ 𝓨 
        ) + opex_discounted)
    else
        @objective(cep.model, Min, opex_discounted)
    end

  return cep
end



function optimize_and_output(; cep::OptModelCEP,
    config::Dict{Any, Any}, 
    data::OptDataCEP, 
    ts_data::ClustData,
    name::String,
    short_sol::Bool)

    optimize!(cep.model)
    @unpack 𝓖, 𝓨, 𝓣, 𝓡, 𝓢, 𝓛, 𝓒 = get_sets(cep=cep)

    status=Symbol(termination_status(cep.model))
    println(status)
    if termination_status(cep.model) == MOI.INFEASIBLE_OR_UNBOUNDED || termination_status(cep.model) == MOI.INFEASIBLE || termination_status(cep.model) == MOI.LOCALLY_INFEASIBLE
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
   

            if short_sol
                #open(joinpath(pwd(),"P2H_CapacityExpansion","results", "$name.txt"), "w") do file
                    #println(file, "Objective = $objective")
                # for r ∈ 𝓡, g ∈ cep.sets["nodes"], y ∈ 𝓨
                #     val = value.(cep.model[:TotalCapacityAnnual][r,g,y])
                #     str = "TotalCapacityAnnual[$r,$g,$y]"
                #     println(file, "$str = $val")
                # end

                # for y ∈ 𝓨
                #     tmp = 0
                #     for g ∈ 𝓖
                #         tmp += value.(cep.model[:COST]["var",y,g])
                #     end
                #     str = "COSTvar$y"
                #     println(file, "$str = $tmp")
                # end
                @info "Feasible"


            else
                open(joinpath(pwd(),"P2H_CapacityExpansion","results", "$name.txt"), "w") do file
                    println(file, "Objective = $objective")
                for v ∈ all_variables(cep.model)
                    if value.(v) != 0
                        val = value.(v)
                        str = string(v)
                        variables[str] = val
                        println(file, "$str = $val")
                    end
                end

                for r ∈ axes(ts_data.ts)[1], t ∈ axes(ts_data.ts)[3] 
                    str = "Demand[$r,$t,$(config["year"])]"
                    val = ts_data.ts[r,"Demand",t] * data.data["demand"][r,config["year"],"electricity"]
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
                println(file, "TimeSeriesWeight = $(ts_data.weight)") 
            end    

        end
        return OptResult(cep.model, status, objective, variables)
    end  
end




