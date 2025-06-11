
struct OptModelCEP
  model::JuMP.Model
  sets::Dict{Any,Any}
  params::Dict{Any,Any}
end

struct OptResult
  model::JuMP.Model
  status::Symbol
  objective::Float64
  variables::Dict{Any,Any}
 end

struct load_timeseries_data
  ts::JuMP.Containers.DenseAxisArray
  weights::Dict{Any,Any}
 end


struct OptDataCEPLine 
  generator::String
  node_start::String
  node_end::String
  power_lim::Number
  length::Number
  eff::Number
end

struct OptDataCEP 
  data::Dict
  lines::Dict{Tuple{String, String}, OptDataCEPLine}
end

struct ClustData
  ts::JuMP.Containers.DenseAxisArray
  weight::Dict{Any,Any}
 end