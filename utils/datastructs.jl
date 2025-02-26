
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
  probabilities::Dict{Any,Any}
 end

struct Sample
  oos_obj::Float64
  ins_obj::Float64
  oos_var::Dict{Any,Any}
  ins_var::Dict{Any,Any}
end

struct EVPI
  objective::Float64
  variables::Dict{Any,Any}
 end

 struct EEV
  ev_obj::Float64
  eev_obj::Float64
  ev_variables::Dict{Any,Any}
  eev_variables::Dict{Any,Any}
 end