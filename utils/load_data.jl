function sample_scenarios(; s_num::Integer, config::Dict{Any,Any},)
  # generate scenarios with sampling
  if s_num > length(config["avail_scenarios"])
    println("Choose a number between 1 and ")
  else
    return sample(config["avail_scenarios"], s_num, replace=false)
  end
end


function load_time_series(; s_num::Integer,
  config::Dict{Any,Any}, timesteps::UnitRange{Int})
 
    # return container with RES and Demand data
    number_of_scenarios = s_num ^ 4
    ts_data = JuMP.Containers.DenseAxisArray(zeros(length(2020:10:2050), length(config["timeseries"]), number_of_scenarios, config["season_length"]*4), 2020:10:2050, config["timeseries"], 1:number_of_scenarios, collect(1:config["season_length"]*4))
    
    scenarios = HydrogenExpansion.sample_scenarios(s_num=s_num, config=config)

    for g in config["timeseries"]
      data_path=normpath(joinpath(dirname(@__FILE__),"..","data", "$g.csv"))
      df = DataFrame(CSV.File(data_path))
      
      
      for y in 2020:10:2050
        k_beg = 1
        for t in [1:2189,2190:4379,4380:6569,6570:8760]
          repeating = Int(number_of_scenarios/s_num)
          for (i,n) in enumerate(repeat(scenarios, repeating))
            # sample one random week
            k = sample(t[1:end-config["season_length"]], 1, replace=false)[1]
            ts_data[y,g,i,k_beg:k_beg+config["season_length"]-1] = df[k:k+config["season_length"]-1, string(n)]
          end
          k_beg += config["season_length"]
        end 
      end 
    end
  return ts_data
end


function read_yaml_file()
  # return yaml file
  return YAML.load_file(joinpath(dirname(@__FILE__),"..","data","config.yml"))
end


function load_data()

  data_path=normpath(joinpath(dirname(@__FILE__),"..","data", "data.xlsx"))
  xf = XLSX.readxlsx(data_path)
  names = XLSX.sheetnames(xf)
  data_dict = Dict()

  for n in names
    df = DataFrame(XLSX.readtable(data_path, n))

    # Iterate through DataFrame rows and populate the dictionary
    temp = Dict()
    for row in eachrow(df)
        key = n in ["d_h2", "d_power", "emission_penalty"] ? (row.Year) :
        n in ["lifetime", "emission"] ? (row.Generator) :
        (row.Generator, row.Year)
        if n in ["d_h2"]
          temp[key] = row.Value*0.3
        elseif n in ["d_power"]
          temp[key] = row.Value*0.4
        else
          temp[key] = row.Value
        end

    end
    data_dict[n] = temp
  end
  return data_dict
end


function average_ts_data(; ts_data::JuMP.Containers.DenseAxisArray)

  ts_data_2 = JuMP.Containers.DenseAxisArray(zeros(length(2020:10:2050), length(axes(ts_data)[2]), 1, length(axes(ts_data)[4])), 2020:10:2050, axes(ts_data)[2], 1, axes(ts_data)[4])

  for y in axes(ts_data)[1] for g in axes(ts_data)[2] for s in axes(ts_data)[3] for t in axes(ts_data)[4]
    ts_data_2[y,g,1,t] += ts_data[y,g,s,t]
  end end end end

  return ts_data_2 / length(axes(ts_data)[3])
end

