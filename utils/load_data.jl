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
    sc[r,t,1:8760] = CountryData[t][1:8760, r]
  end

  ## average the data if specified
  if config["average"] > 1
    @info "Each time-series is averaged in $(config["average"])-hourly steps"
    sc, weights = average_ts_data(ts_data=sc, config=config)
  else
    weights = Dict(d => 1 for d ∈ 1:8760)
  end

  return ClustData(sc, weights)
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

    # Convert to Jump Density Array
    if n ∈ ["c_CAPEX", "c_var", "c_fix"]
      data_dict[n] = create_array_from_df(df, keys(config["techs"]), config["year"]; default=0.0001)
    elseif n ∈ ["cap"]
      data_dict[n] = create_array_from_df(df, config["countries"], keys(config["techs"]), config["year"]; default=99999.0)
    elseif n ∈ ["cap_init"]
      data_dict[n] = create_array_from_df(df, config["countries"], keys(config["techs"]), config["year"]; default=0.0)
    elseif n ∈ ["lifetime", "emission"]
      data_dict[n] = create_array_from_df(df, keys(config["techs"]); default=0.0)
    elseif n ∈ ["demand"]
      data_dict[n] = create_array_from_df(df, config["countries"], config["year"], unique([config["techs"][t]["output"]["carrier"] for t ∈ keys(config["techs"])]); default=0.0)
    elseif n ∈ ["eta"]
      data_dict[n] = create_array_from_df(df, keys(config["techs"]), config["year"]; default=1.0)
    else #emission budget
      data_dict[n] = create_array_from_df(df, config["countries"], config["year"])
    end
  end
  return data_dict
end




function create_array_from_df(df::DataFrame, els...; default=0.0)

  A = JuMP.Containers.DenseAxisArray(
      fill(default, length.(els)...), els...)

  for r in eachrow(df)
      try
          A[r[1:end-1]...] = r.Value
      catch err
          @debug err
      end
  end

  return A
end



function average_ts_data(; ts_data::JuMP.Containers.DenseAxisArray, config::Dict{Any,Any})

  ts_len = Int(ceil(8760/config["average"]))

  ts_data_2 = JuMP.Containers.DenseAxisArray(zeros(length(config["countries"]), length(config["timeseries"]), ts_len), config["countries"], config["timeseries"], 1:ts_len) 

  for r ∈ config["countries"], ts ∈ config["timeseries"], t ∈ 1:ts_len
    if t == ts_len
      ts_data_2[r,ts,t] = mean(ts_data[r,ts,(t-1)*config["average"]+1:end])
    else
      ts_data_2[r,ts,t] = mean(ts_data[r,ts,(t-1)*config["average"]+1:t*config["average"]])
    end
  end

  weight = Dict(d => 8760/ts_len for d ∈ 1:ts_len)

  return ts_data_2, weight
end


function load_cep_data(; config::Dict{Any,Any})

  data = load_data(config=config)
  lines = load_cep_data_lines(config=config)

  return OptDataCEP(data, lines)
end 


function load_cep_data_lines(; config::Dict{Any,Any})

  
  # Read the sheet into a DataFrame
  df = XLSX.readtable(config["data_path"], "trade") |> DataFrame
  tab = filter(row -> row["Node_start"] ∈ config["countries"] && row["Node_end"] ∈ config["countries"], df)

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



clean_float(x) = (ismissing(x) || (x isa Number && isnan(x))) ? 0.0 :
         parse(Float64, replace(strip(string(x)), '\0' => ""))




function read_txt_file(path)
    df = CSV.read(path, DataFrame)

    # Strip whitespace from column names
    rename!(df, Dict(c => strip(c) for c in names(df)))

    # Apply clean_float to all cells
    for col in names(df)
        df[!, col] = clean_float.(df[!, col])
    end

    # Round everything to nearest integer
    df .= round.(df)

    return df
end
