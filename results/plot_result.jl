function plot_time_series(; df::JuMP.Containers.DenseAxisArray,
    year::Integer, 
    title::String,
    config::Dict)

    years = axes(df)[1]
    generators = axes(df)[2]

    plots = [scatter(
        x=1:8760, 
        y=df[year, g,:], 
        name=g, 
        stackgroup="one",
        marker_color=config["colour_codes"][g]) for g in generators]
    #append!(plots, [scatter(x=1:8760, y=z_d[year_ts,:]*d_power[year], name="Load", stackgroup="two", marker_color=config.colour_codes["Load"])])
    #append!([scatter(x=1:8760, y=ll[year,:], name="Lost Load", stackgroup="three", marker_color="red")], plots)
    plot(plots, Layout(title_text=title))
end

function plot_results_bar(; df::JuMP.Containers.DenseAxisArray, title::String, config::Dict)
    years = axes(df)[1]
    generators = axes(df)[2]
 
    plots = [bar(x=years, y=df[:, g], name=g, marker_color=config["colour_codes"][g]) for g in generators]
    #append!([(x=years, y=df[:, "LL"], name="LL", marker_color="red")], plots)
    plot(plots, Layout(barmode="relative", title_text=title))
end   


function plot_results_bar_sp_eev(; 
    df1::JuMP.Containers.DenseAxisArray,  
    df2::JuMP.Containers.DenseAxisArray,
    title1::String,
    title2::String, 
    title::String, 
    config::Dict,
    kwargs...)

    years = axes(df1)[1]
    generator = axes(df1)[2]      

    plots = [bar(x=(years, fill(title1, length(years))), y=df1[:, g], name=g, showlegend=true, marker_color=config["colour_codes"][g]) for g in generator]
    append!(plots, [bar(x=(years, fill(title2, length(years))), y=df2[:, g], name=g,  showlegend=false, marker_color=config["colour_codes"][g]) for g in generator])
    # append lost load if generation is plotted
    if :lost_load1 in keys(kwargs) && :lost_load2 in keys(kwargs)
        append!(plots, [bar(x=(years, fill(title1, length(years))), y=get(kwargs, :lost_load1, 0), showlegend=false, marker_color="red")])
        append!(plots, [bar(x=(years, fill(title2, length(years))), y=get(kwargs, :lost_load2, 0), name="Lost Load",  showlegend=true, marker_color="red")])
    end

    #display(plot(plots, Layout(barmode="relative", title_text=title)))
    return plot(plots, Layout(barmode="relative", title_text=title))
end 



function plot_results_bar_scenarios(; df::JuMP.Containers.DenseAxisArray, title::String, config::Dict, dimensions::Integer)
    scenario = axes(df)[1]
    years = axes(df)[2]

    if dimensions == 3
        legend = Dict(s => (findfirst(x -> x == s, scenario) == 1) for s in scenario)
        generators = axes(df)[3]
        plots = [bar(x=(years, repeat([findfirst(x -> x == s, scenario)], length(years))), y=df[s,:, g], name=g, showlegend=legend[s],  marker_color=config["colour_codes"][g]) for g in generators for s in scenario]
    else
        plots = [bar(x=(years, repeat([findfirst(x -> x == s, scenario)], length(years))), y=df[s,:], marker_color="orange") for s in scenario]
    end
    plot(plots, Layout(barmode="relative", title_text=title))
    return plot(plots, Layout(barmode="relative", title_text=title))
end    

