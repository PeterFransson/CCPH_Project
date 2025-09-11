
function get_cᵢ_model(modeloutputs::Vector{ModelOutput})
    cᵢ_model = Tuple{Real,Real}[]
    for modeloutput in modeloutputs
        append!(cᵢ_model,[(data_day.modelinstoutput₁.cᵢ,data_day.modelinstoutput₂.cᵢ) for data_day in modeloutput.output_day])        
    end
    return cᵢ_model
end

function get_cᵢ_mean_model(cᵢ::Vector{Tuple{R,R}},input::RawInputData) where {R<:Real}
    cᵢ_mean_model = Real[]
    for (cᵢ_day,data) in zip(cᵢ,input.weather_growth)
        cᵢ₁,cᵢ₂ = cᵢ_day
        day_nr = CCPH.Dates.dayofyear(data.date)
        daylength = CCPH.daylighthour(data.lat*pi/180,day_nr)*3600 #Seconds
        t₁,t₂,Δt₁,Δt₂ = CCPH.SDM2_get_time_points(daylength)    
        cᵢ_integral_model = 2*(cᵢ₁*Δt₁+cᵢ₂*Δt₂)
        push!(cᵢ_mean_model,cᵢ_integral_model/daylength)
    end
    return cᵢ_mean_model
end

#Calculate daily mean cᵢ
function get_cᵢ_mean_model(modeloutputs::Vector{ModelOutput},raw_input::RawInputData)
    cᵢ = get_cᵢ_model(modeloutputs)
    cᵢ_mean_model = get_cᵢ_mean_model(cᵢ,raw_input)
    return cᵢ_mean_model
end

function get_cᵢ_mean_model(par::ModelPar,raw_input::Vector{RawInputData};stand_type::Symbol=:Fertilized)
    Xₜ = Xₜ_fun.(raw_input,Ref(par))
    modeloutputs = run_week.(raw_input,Xₜ,Ref(par))    
    cᵢ_mean_model = get_cᵢ_mean_model.(modeloutputs,raw_input)
    return cᵢ_mean_model
end

function get_c_i()    
    fld = "crossval_20241105_shared_W_1_5_run_9"

    stand_type_F = JLD.load("output/"*fld*"/result_F.jld","stand_type")
    raw_input_F = RawInputData(;stand_type=stand_type_F)
    Ec_data_F = calc_Ec_data.(raw_input_F)
    GPP_data_F = get_GPP_data.(raw_input_F;stand_type=stand_type_F)   
    x_opt_F = JLD.load("output/"*fld*"/result_F.jld","x_opt")         
    par_F = ModelPar(x_opt_F;stand_type=stand_type_F)
    cᵢ_mean_model_F = get_cᵢ_mean_model(par_F,raw_input_F;stand_type=stand_type_F)

    stand_type_C = JLD.load("output/"*fld*"/result_C.jld","stand_type")
    raw_input_C = RawInputData(;stand_type=stand_type_C)
    Ec_data_C = calc_Ec_data.(raw_input_C)
    GPP_data_C = get_GPP_data.(raw_input_C;stand_type=stand_type_C)   
    x_opt_C = JLD.load("output/"*fld*"/result_C.jld","x_opt")         
    par_C = ModelPar(x_opt_C;stand_type=stand_type_C)
    cᵢ_mean_model_C = get_cᵢ_mean_model(par_C,raw_input_C;stand_type=stand_type_C)    

    pl1 = [plot(xlabel="2015",ylabel="cᵢ/cₐ",legends=false, ylims = (0.4,1.0),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (0.4,1.0),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (0.4,1.0),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (0.4,1.0),yaxis=false,guidefontsize=12)]

    save_output_fld = "./output/paper_results/plant_vs_weather"

    years = [2015,2016,2017,2018]

    for i = 1:4
        date = [weather.date for weather in raw_input_F[i].weather_growth]
        Cₐ = [weather.Cₐ for weather in raw_input_F[i].weather_growth]
        cᵢ_F = cᵢ_mean_model_F[i]
        cᵢ_C = cᵢ_mean_model_C[i]

        df_F = DataFrames.DataFrame([date, cᵢ_F], [:Date, :c_i])
        CSV.write(save_output_fld*"/c_i_data_Fertilized_$(years[i]).csv", df_F)

        df_C = DataFrames.DataFrame([date, cᵢ_C], [:Date, :c_i])
        CSV.write(save_output_fld*"/c_i_data_Control_$(years[i]).csv", df_C)

        start_tick = ""
        end_tick = ""

        plot!(pl1[i],date,cᵢ_F./Cₐ,linecolor=:blue)
        plot!(pl1[i],date,cᵢ_C./Cₐ,linecolor=:red)
        plot!(pl1[i],xticks=([date[1],date[end]],[start_tick,end_tick]))
    end

    save_fld = "./plots/paper_results/plant_vs_weather"
    isdir(save_fld)|| mkdir(save_fld)

    pl1_fin = plot(pl1[1],pl1[2],pl1[3],pl1[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)  
    savefig(pl1_fin,save_fld*"/ci_ca_ratio_date.svg")
end