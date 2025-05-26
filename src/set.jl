function setup_opt_set!(sets::Dict{Any, Any}, 
    ts_data::JuMP.Containers.DenseAxisArray, 
    config::Dict{Any, Any},
    data)

    sets["techs"] = keys(config["techs"])      
    sets["years"] = config["dispatch"] ? config["year"] : 2020:10:2050
    sets["timesteps"] = axes(ts_data)[3]
    sets["regions"] = config["countries"]
    sets["storage_techs"] = [key for (key, val) ∈ config["techs"] if get(val, "tech_group", "") == "storage"]
    sets["dispatch"] = [key for (key, val) ∈ config["techs"]if get(val, "tech_group", "") == "dispatchable_generation"]
    sets["non_dispatch"] = [key for (key, val) ∈ config["techs"] if get(val, "tech_group", "") == "non_dispatchable_generation"]
    sets["conversion"] = [key for (key, val) ∈ config["techs"] if get(val, "tech_group", "")  == "conversion"]  
    sets["transmission"] = [key for (key, val) ∈ config["techs"] if get(val, "tech_group", "")  == "transmission"]  
    sets["nodes"] = setdiff(keys(config["techs"]), sets["transmission"])
    sets["lines"] = unique([i[2] for i ∈ keys(data.lines)])
    # the structure describes if the capacity of this tech is either setup on a node or along a line

    for i ∈ ["input", "output"]
        sets[i] = Dict(c => [tech for (tech, spec) in config["techs"] if haskey(spec, i) && get(spec[i], "carrier", nothing) == c] for c ∈ config["energy_carriers"])
    end

    return sets
end


function setup_opt_set_carrier!(sets::Dict{Any, Any}, 
    config::Dict{Any, Any})

    sets["carrier"] = Dict{String, Vector{String}}()

    for g ∈ sets["techs"]
        carriers = String[]

        if haskey(config["techs"][g]["input"], "carrier")
            input_carrier = config["techs"][g]["input"]["carrier"]
            push!(carriers, input_carrier)

            if !haskey(sets, input_carrier)
                sets[input_carrier] = String[]
            end
            push!(sets[input_carrier], g)
        end

        if haskey(config["techs"][g]["output"], "carrier")
            output_carrier = config["techs"][g]["output"]["carrier"] 
            cond = haskey(config["techs"][g]["input"], "carrier") && input_carrier == output_carrier
            if !cond
                push!(carriers, output_carrier) 
                if !haskey(sets, output_carrier)
                    sets[output_carrier] = String[]
                end
                push!(sets[output_carrier], g)
            end
        end
        sets["carrier"][g] = carriers
    
    end
    return sets
end



"""
    get_sets(cep::OptModelCEP) -> NamedTuple

Returns a NamedTuple of commonly used model sets from the `cep` model for convenience.

This includes standard sets such as:
- `𝓖`: Technologies (`"techs"`)
- `𝓨`: Investment years (`"years"`)
- `𝓣`: Timesteps (`"timesteps"`)
- `𝓡`: Regions (`"regions"`)
- `𝓢`: Storage technologies (`"storage_techs"`)
- `𝓛`: Transmission lines (`"lines"`, optional, returns `nothing` if not defined)
- `𝓒`: Energy carriers (`config["energy_carriers"]`, optional, returns `nothing` if not available)

Intended to simplify repeated unpacking of sets in model-building functions.
Use with `@unpack` for readability and maintainability.
"""
function get_sets(; cep)
    return (
        𝓖 = cep.sets["techs"],
        𝓨 = cep.sets["years"],
        𝓣 = cep.sets["timesteps"],
        𝓡 = cep.sets["regions"],
        𝓢 = cep.sets["storage_techs"],
        𝓛 = cep.sets["lines"],
    )
end



