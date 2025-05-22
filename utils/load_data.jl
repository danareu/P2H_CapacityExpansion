function sample_scenarios(; s_num::Integer, config::Dict{Any,Any},)
  # generate scenarios with sampling
  if s_num > length(config["avail_scenarios"])
    println("Choose a number between 1 and ")
  else
    return sample(config["avail_scenarios"], s_num, replace=false)
  end
end



function load_timeseries_data_provided(; config::Dict{Any,Any})

  # Generate the data path based on application and region
  CountryData = Dict(t => DataFrame(CSV.File(normpath(joinpath(dirname(@__FILE__),"..","data", "$t.csv")), delim=';', decimal=',')) for t ∈ config["timeseries"])
  for t ∈ config["timeseries"]
    select!(CountryData[t], config["countries"])
  end
  return CountryData
end



function run_clust(; CountryData::Dict{String, DataFrames.DataFrame}, representation::String, config::Dict{Any,Any}, n_clust::Integer)
  
  data = TSClustering.normalize_data(config=config, CountryData=CountryData)
  data_clustering = TSClustering.create_clustering_matrix(technology=config["timeseries"], CountryData=data)

  # define distance matrix r
  D = TSClustering.define_distance(w=0, data_clustering=data_clustering, fast_dtw=false)
  result = hclust(D, linkage=:ward)
  cl = cutree(result, k=n_clust)

  weights = Dict{Int64, Int64}()
  for i ∈ cl
      weights[i] = get(weights, i, 0) + 1
  end

  weights_yrl = JuMP.Containers.DenseAxisArray(zeros(n_clust, 24), 1:n_clust, 1:24)
  for c ∈ 1:n_clust, h ∈ 1:24
    weights_yrl[c,h] = weights[c] / 8760
  end


  # read in again because of Julia memory issues
  CountryData = P2H_CapacityExpansion.load_timeseries_data_provided(config=config)

  # calculate the representative only medoid so far available
  if representation == "medoid"
    cluster_dict_org = TSClustering.calculate_medoid(data_org=CountryData, cl=cl, config=config, K=n_clust, technology=config["timeseries"])
  else
    cluster_dict_org = TSClustering.calculate_representative(representative="centroid",data_clustering=data_clustering, cl=cl, weights=weights,k=n_clust)
  end

  sc = TSClustering.scaling(data_org=CountryData, scaled_clusters=cluster_dict_org, k=n_clust, weights=weights, config=config, technology=config["timeseries"])

  return load_timeseries_data(sc, weights_yrl)
end



function load_timeseries_data_full(; config::Dict{Any,Any})

  # Read om the data 
  CountryData = Dict(t => DataFrame(CSV.File(normpath(joinpath(dirname(@__FILE__),"..","data", "$t.csv")), delim=';', decimal=',')) for t ∈ config["timeseries"])

  sc = JuMP.Containers.DenseAxisArray(zeros(length(config["countries"]), length(config["timeseries"]), 8760), config["countries"], config["timeseries"], 1:8760) 

  # populate the array
  for t ∈ config["timeseries"], r ∈ config["countries"]
    data_col = CountryData[t][!, r]
    if length(data_col) != 8760
      error("Time series for $r in $t does not have 8.760 entries.")
    else
      sc[r,t,:] .= data_col
    end
  end
  return sc
end





function read_yaml_file()
  # return yaml file
  return YAML.load_file(joinpath(dirname(@__FILE__),"..","data","config.yml"))
end


function load_data(;config::Dict{Any,Any})

  xf = XLSX.readxlsx(config["data_path"])
  nam = XLSX.sheetnames(xf)
  data_dict = Dict()

  for n ∈ nam
    df = DataFrame(XLSX.readtable(config["data_path"], n))

    if "Unit" ∈ names(df)
      select!(df, Not("Unit"))
    end

    config["tech_mapping"] = Dict("Batt" => ["Batt_in", "Batt_out"])

    if "Generator" ∈ names(df)
      default_value_storage!(df)
    end

    # Convert to Jump Density Array
    if n ∈ ["c_CAPEX", "c_var", "c_fix"]
      data_dict[n] = create_array_from_df(df, keys(config["techs"]), 2020:10:2050)
    elseif n ∈ ["cap", "cap_init"]
      data_dict[n] = create_array_from_df(df, config["countries"], keys(config["techs"]), 2020:10:2050)
    elseif n ∈ ["lifetime", "emission"]
      data_dict[n] = create_array_from_df(df, keys(config["techs"]))
    elseif n ∈ ["demand"]
      data_dict[n] = create_array_from_df(df, config["countries"], 2020:10:2050, config["energy_carriers"])
    elseif n ∈ ["eta"]
      data_dict[n] = create_array_from_df(df, keys(config["techs"]), 2020:10:2050)
    end
  end
  return data_dict
end




function create_array_from_df(df::DataFrame, els...)

  A = JuMP.Containers.DenseAxisArray(
    zeros(length.(els)...), els...)

  # order columns 
  # Fill in values from Excel
  for r in eachrow(df)
      try
          A[r[1:end-1]...] = r.Value 
      catch err
          @debug err
      end
  end
  return A 
end


function average_ts_data(; ts_data::JuMP.Containers.DenseAxisArray)

  ts_data_2 = JuMP.Containers.DenseAxisArray(zeros(length(2020:10:2050), length(axes(ts_data)[2]), 1, length(axes(ts_data)[4])), 2020:10:2050, axes(ts_data)[2], 1, axes(ts_data)[4])

  for y in axes(ts_data)[1] for g in axes(ts_data)[2] for s in axes(ts_data)[3] for t in axes(ts_data)[4]
    ts_data_2[y,g,1,t] += ts_data[y,g,s,t]
  end end end end

  return ts_data_2 / length(axes(ts_data)[3])
end




function default_value_storage!(df)
  # map general tech => subtechs
  tech_expansions = Dict(
    "D_Battery_Li-Ion" => ["D_Battery_Li-Ion_in", "D_Battery_Li-Ion_out"]
  )

  for (base_tech, expanded) ∈ tech_expansions
    if base_tech ∈ df[!, "Generator"]
        rows = filter(row -> row["Generator"] == base_tech, df)
        for new_tech in expanded
            new_rows = deepcopy(rows)
            new_rows[!, "Generator"] .= new_tech
            append!(df, new_rows)
        end
        filter!(row -> row["Generator"] != base_tech, df)
    end
  end
  return df
end

function load_cep_data(; config::Dict{Any,Any})

  data = load_data(config=config)
  lines = load_cep_data_lines(config=config)

  return OptDataCEP(data, lines)
end 


function load_cep_data_lines(; config::Dict{Any,Any})

  # Read the sheet into a DataFrame
  tab = XLSX.readtable(config["data_path"], "trade") |> DataFrame

  techs = unique(tab[!, :Generator])
  lines_list = unique(tab[!, :Line])
  lines = Dict{Tuple{String, String}, OptDataCEPLine}()

  for tech ∈ techs, line ∈ lines_list
      row = filter(r -> r[:Generator] == tech && r[:Line] == line, tab)

      if nrow(row) == 1
          node_start = row[1, :Node_start]
          node_end   = row[1, :Node_end]
          power_lim  = row[1, :Power]
          length     = row[1, :Length]
          generator  = line  
          eff = (length^config["techs"][tech]["efficiency"])/length

          lines[(tech, line)] = OptDataCEPLine(generator, node_start, node_end, power_lim, length, eff)
      else
          @warn "No unique row found for Generator=$tech and Line=$line"
      end
  end

  return lines
end

