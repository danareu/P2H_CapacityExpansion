function investigate_samples(; 
    config::Dict{Any,Any}, 
    data, 
    ts_data::JuMP.Containers.DenseAxisArray)

    # solve the expected value problem
    k = ceil(Int, size(ts_data)[3]*0.4) # 20% training set

    # run optimization problem
    ins = HydrogenExpansion.run_opt(scenarios=collect(1:k), ts_data=ts_data, data=data, config=config) 
    # fix fist-stage decision variable & run deterministic problem m times
    m = size(ts_data)[3]
    obj_avg = 0
    var = Dict()
    for i in k+1:m
        result_fix = HydrogenExpansion.run_opt(scenarios=[i], ts_data=ts_data, data=data, config=config, evp=ins.variables["cap_new"])
        obj_avg += result_fix.objective
        if i == (k+1)
            var = result_fix.variables
        else
            var_tmp = result_fix.variables
            for key in keys(var_tmp)
                if key in ["ll", "prices_el" ]
                    for s in axes(var_tmp[key])[1] for y in axes(var_tmp[key])[2] for t in axes(var_tmp[key])[3]
                        var[key][k+1,y,t] += var_tmp[key][s,y,t]
                        var[key][k+1,y,t] *= 1/2
                    end end end
                elseif key in ["gen"]
                    for s in axes(var_tmp[key])[1] for y in axes(var_tmp[key])[2] for t in axes(var_tmp[key])[3] for m in axes(var_tmp[key])[4]
                        var[key][k+1,y,t,m] += var_tmp[key][s,y,t,m]
                        var[key][k+1,y,t,m] *= 1/2
                    end end end end
                elseif key in ["cap_new", "cap"] 
                    for y in axes(var_tmp[key])[1] for g in axes(var_tmp[key])[2]
                        var[key][y,g] += var_tmp[key][y,g]
                        var[key][y,g] *= 1/2
                    end end
                elseif key in ["prices_h2"]
                    for y in axes(var_tmp[key])[1] for g in axes(var_tmp[key])[2]
                        var[key][k+1,g] += var_tmp[key][y,g]
                        var[key][k+1,g] *= 1/2
                    end end 
                end 

            end
        end
    end
    return Sample(obj_avg/(m-k), ins.objective, var,ins.variables)
end


function calculate_evpi(; evp, sp)
    # vvs=evp-sp
    evpi = round(sp - evp)
    println("EVPI: $(round(sp)) - $(round(evp)) = $evpi")
    return evpi
end

function calculate_vss(; eev, sp)
    # vvs=evp-sp
    vss = round(eev - sp)
    println("VSS: $(round(eev)) - $(round(sp)) = $vss")
    return vss
end

function calculate_eev(;config::Dict{Any,Any}, ts_data_org::JuMP.Containers.DenseAxisArray, data, probabilities::Dict{Any,Any})
    # average of all random variables
    ts_data = HydrogenExpansion.average_ts_data(ts_data=ts_data_org)
    # solve the stochastic problem
    result_first = HydrogenExpansion.run_opt(scenarios=[1], ts_data=ts_data, data=data, config=config) 
    # fix the first-stage variables & solve deterministic problem
    result_second = HydrogenExpansion.run_opt(scenarios=collect(1:size(ts_data_org)[3]), ts_data=ts_data_org, data=data, config=config, evp=result_first.variables["cap_new"])
    
    return EEV(result_first.objective, result_second.objective, result_first.variables, result_second.variables)
end



function calculate_expected_value_information(; config::Dict{Any,Any}, data, ts_data::JuMP.Containers.DenseAxisArray, probabilities::Dict{Any,Any})
    # solve each scenario deterministically and take mean value

    obj_avg = 0
    var = Dict()
    # iterate over scenarios
    for i in 1:size(ts_data)[3]
        det = HydrogenExpansion.run_opt(scenarios=[i], ts_data=ts_data, data=data, config=config)
        obj_avg += (det.objective*probabilities[i])
        if i == 1
            var = det.variables
            for key in keys(var)
                var[key] *= probabilities[i] 
            end
        else
            var_tmp = det.variables
            for key in keys(var_tmp)
                if key in ["ll", "prices_el" ]
                    for s in axes(var_tmp[key])[1] for y in axes(var_tmp[key])[2] for t in axes(var_tmp[key])[3]
                        var[key][1,y,t] += (var_tmp[key][s,y,t]*probabilities[i])
                    end end end
                elseif key in ["gen"]
                    for s in axes(var_tmp[key])[1] for y in axes(var_tmp[key])[2] for t in axes(var_tmp[key])[3] for m in axes(var_tmp[key])[4]
                        var[key][1,y,t,m] += (var_tmp[key][s,y,t,m] *probabilities[i])
                    end end end end
                elseif key in ["cap_new", "cap"] 
                    for y in axes(var_tmp[key])[1] for g in axes(var_tmp[key])[2]
                        var[key][y,g] += (var_tmp[key][y,g]*probabilities[i])
                    end end
                elseif key in ["prices_h2"]
                    for y in axes(var_tmp[key])[1] for g in axes(var_tmp[key])[2]
                        var[key][1,g] += (var_tmp[key][y,g]*probabilities[i])
                    end end 
                end 

            end
        end
    end
    return EVPI(obj_avg, var)

end

function aggregate_time_steps(; df::JuMP.Containers.DenseAxisArray)
    # aggregate to yearly values in TWh
    
    if length(size(df)) == 4
        df_gen = JuMP.Containers.DenseAxisArray(zeros(length(axes(df)[1]),length(axes(df)[2]),length(axes(df)[3])),axes(df)[1], axes(df)[2], axes(df)[3])
        for s in axes(df)[1] for y in axes(df)[2] for g in axes(df)[3]
            for t in axes(df)[4]
                df_gen[s, y, g] += df[s,y,g,t]
            end
            
        end end end
        df_gen = df_gen * (8760/size(df)[4])
    else
        df_gen = JuMP.Containers.DenseAxisArray(zeros(length(axes(df)[1]),length(axes(df)[2])),axes(df)[1], axes(df)[2])
        for s in axes(df)[1] for y in axes(df)[2]
            for t in axes(df)[3]
                df_gen[s, y] += df[s,y,t]
            end
        end end
        df_gen = df_gen * (8760/size(df)[3])

    end
    df_gen = df_gen/1000
    return df_gen
end


function aggregate_scenarios(; df::JuMP.Containers.DenseAxisArray)
    # aggregate the scenarios by taking the mean value
    # must be aggregated yearly beforehand
    
    if length(size(df)) == 3
        df_gen = JuMP.Containers.DenseAxisArray(zeros(length(axes(df)[2]),length(axes(df)[3])), axes(df)[2], axes(df)[3])
        for y in axes(df)[2] for g in axes(df)[3]
            for s in axes(df)[1]
                df_gen[y, g] += df[s,y,g]
            end
        end end
    else
        df_gen = JuMP.Containers.DenseAxisArray(zeros(length(axes(df)[2])), axes(df)[2])
        for y in axes(df)[2] 
            for s in axes(df)[1]
                df_gen[y] += df[s,y]
            end
        end
    end
    df_gen = df_gen/(size(df)[1])
    return df_gen
end