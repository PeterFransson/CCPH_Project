
function get_wue_model(modeloutputs::Vector{ModelOutput})
    wue_model = Real[]
    cons = CCPH.Constants()
    for modeloutput in modeloutputs
        append!(wue_model,[(data_day.GPP/(cons.M_C*1000))/(data_day.Ec/(1000*cons.M_H2O/cons.ρ_H2O)) for data_day in modeloutput.output_day])  #mol C mol⁻¹ H2O 
    end
    return wue_model
end

function get_gₛ_model(modeloutputs::Vector{ModelOutput})
    gₛ_model = Tuple{Real,Real}[]
    for modeloutput in modeloutputs
        append!(gₛ_model,[(data_day.gₛ₁,data_day.gₛ₂) for data_day in modeloutput.optval_day])        
    end
    return gₛ_model
end

function get_gₛ_mean_model(gₛ::Vector{Tuple{R,R}},input::RawInputData) where {R<:Real}
    gₛ_mean_model = Real[]
    for (gₛ_day,data) in zip(gₛ,input.weather_growth)
        gₛ₁,gₛ₂ = gₛ_day
        day_nr = CCPH.Dates.dayofyear(data.date)
        daylength = CCPH.daylighthour(data.lat*pi/180,day_nr)*3600 #Seconds
        t₁,t₂,Δt₁,Δt₂ = CCPH.SDM2_get_time_points(daylength)    
        gₛ_integral_model = 2*(gₛ₁*Δt₁+gₛ₂*Δt₂)
        push!(gₛ_mean_model,gₛ_integral_model/daylength)
    end
    return gₛ_mean_model
end

#Calculate daily mean gₛ
function get_gₛ_mean_model(modeloutputs::Vector{ModelOutput},raw_input::RawInputData)
    gₛ = get_gₛ_model(modeloutputs)
    gₛ_mean_model = get_gₛ_mean_model(gₛ,raw_input)
    return gₛ_mean_model
end

function get_plant_var(par::ModelPar,raw_input::Vector{RawInputData};stand_type::Symbol=:Fertilized)
    Xₜ = Xₜ_fun.(raw_input,Ref(par))
    modeloutput = run_week.(raw_input,Xₜ,Ref(par))    
    wue_model = get_wue_model.(modeloutput)
    Nₘ_f_model = get_Nₘ_f_model.(modeloutput)
    gₛ_mean_model = get_gₛ_mean_model.(modeloutput,raw_input)
    return (gₛ_mean_model,Nₘ_f_model,wue_model)
end

function plant_var_vs_weather()
    fld = "crossval_20240306_shared_W_1_5_run_4"

    stand_type_F = JLD.load("output/"*fld*"/result_F.jld","stand_type")
    raw_input_F = RawInputData(;stand_type=stand_type_F)
    Ec_data_F = calc_Ec_data.(raw_input_F)
    GPP_data_F = get_GPP_data.(raw_input_F;stand_type=stand_type_F)   
    x_opt_F = JLD.load("output/"*fld*"/result_F.jld","x_opt")         
    par_F = ModelPar(x_opt_F;stand_type=stand_type_F)
    gₛ_F,Nₘ_f_F,wue_F = get_plant_var(par_F,raw_input_F;stand_type=stand_type_F)

    stand_type_C = JLD.load("output/"*fld*"/result_C.jld","stand_type")
    raw_input_C = RawInputData(;stand_type=stand_type_C)
    Ec_data_C = calc_Ec_data.(raw_input_C)
    GPP_data_C = get_GPP_data.(raw_input_C;stand_type=stand_type_C)   
    x_opt_C = JLD.load("output/"*fld*"/result_C.jld","x_opt")         
    par_C = ModelPar(x_opt_C;stand_type=stand_type_C)
    gₛ_C,Nₘ_f_C,wue_C = get_plant_var(par_C,raw_input_C;stand_type=stand_type_C)
    
    pl1 = [plot(xlabel="2015",ylabel="gₛ (mol s⁻¹ m⁻²)",legends=false, ylims = (0.04,0.18),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (0.04,0.18),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (0.04,0.18),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (0.04,0.18),yaxis=false,guidefontsize=12)]

    pl2 = [plot(xlabel="2015",ylabel="Nmf (%)",legends=false, ylims = (0,5),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (0,5),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (0,5),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (0,5),yaxis=false,guidefontsize=12)]

    pl3_1 = plot(xlabel="",ylabel="gₛ (mol s⁻¹ m⁻²)",legends=false,ylims = (0.04,0.18),xlims = (0,69),xformatter=_->"",guidefontsize=12,ytickfontsize=12)
    pl3_2 = plot(xlabel="",ylabel="",legends=false,ylims = (0.04,0.18),xlims = (0,25),xformatter=_->"",yaxis=false,guidefontsize=12,ytickfontsize=12)
    pl3_3 = plot(xlabel="",ylabel="",legends=false,ylims = (0.04,0.18),xlims = (0.0,1.5),xformatter=_->"",yaxis=false,guidefontsize=12,ytickfontsize=12)
    pl3_4 = plot(xlabel="",ylabel="",legends=false,ylims = (0.04,0.18),xlims = (6,30),xformatter=_->"",yaxis=false,guidefontsize=12,ytickfontsize=12)

    pl4_1 = plot(xlabel="I₀ (mol m⁻²)",ylabel="Nmf (%)",legends=false,ylims = (0,5),xlims = (0,69),xformatter=_->"",guidefontsize=12,ytickfontsize=12)
    pl4_2 = plot(xlabel="Tₐ (°C)",ylabel="",legends=false,ylims = (0,5),xlims = (0,25),xformatter=_->"",yaxis=false,guidefontsize=12,ytickfontsize=12)
    pl4_3 = plot(xlabel="VPD (kPA)",ylabel="",legends=false,ylims = (0,5),xlims = (0.0,1.5),xformatter=_->"",yaxis=false,guidefontsize=12,ytickfontsize=12)
    pl4_4 = plot(xlabel="θ (%)",ylabel="",legends=false,ylims = (0,5),xlims = (6,30),xformatter=_->"",yaxis=false,guidefontsize=12,ytickfontsize=12)

    pl5_1 = plot(xlabel="I₀ (mol m⁻²)",ylabel="WUE (mmol mol⁻²)",legends=false,ylims = (2,15),xlims = (0,69),guidefontsize=12,ytickfontsize=12)
    pl5_2 = plot(xlabel="Tₐ (°C)",ylabel="",legends=false,ylims = (2,15),xlims = (0,25),yaxis=false,guidefontsize=12,ytickfontsize=12)
    pl5_3 = plot(xlabel="VPD (kPA)",ylabel="",legends=false,ylims = (2,15),xlims = (0.0,1.5),yaxis=false,guidefontsize=12,ytickfontsize=12)
    pl5_4 = plot(xlabel="θ (%)",ylabel="",legends=false,ylims = (2,15),xlims = (6,30),yaxis=false,guidefontsize=12,ytickfontsize=12)

    Tmean = [[weather.Tmean for weather in raw_input_F[i].weather_growth] for i in 1:4]
    θ_F = [[weather.θₛ*100 for weather in raw_input_F[i].weather_growth] for i in 1:4]
    θ_C = [[weather.θₛ*100 for weather in raw_input_C[i].weather_growth] for i in 1:4]
    Radₜₒ = [[weather.Radₜₒ for weather in  raw_input_F[i].weather_growth] for i in 1:4] 
    I₀ = Radₜₒ*2.3*10^-6 #mol m⁻²
    VPD = [[first(get_env_from_data(weather)).VPD/1000 for weather in  raw_input_F[i].weather_growth] for i in 1:4]

    pl6 = plot(xlabel="gₛ (mol s⁻¹ m⁻²)",ylabel="Nmf (%)",legends=false, xlims = (0.04,0.18), ylims = (0,5),guidefontsize=12,ytickfontsize=12)

    pl8 = [plot(xlabel="2015",ylabel="Δgₛ (F-C) (mol s⁻¹ m⁻²)",legends=false,guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false,guidefontsize=12)]

    for i = 1:4
        date = [weather.date for weather in raw_input_F[i].weather_growth]
        
        start_tick = ""
        end_tick = ""

        plot!(pl1[i],date,gₛ_F[i],linecolor=:blue)
        plot!(pl1[i],date,gₛ_C[i],linecolor=:red)
        plot!(pl1[i],xticks=([date[1],date[end]],[start_tick,end_tick]))

        plot!(pl2[i],date,Nₘ_f_F[i]*100,linecolor=:blue)
        plot!(pl2[i],date,Nₘ_f_C[i]*100,linecolor=:red)
        plot!(pl2[i],xticks=([date[1],date[end]],[start_tick,end_tick]))       

        plot!(pl8[i],date,gₛ_F[i]-gₛ_C[i])
        plot!(pl8[i],xticks=([date[1],date[end]],[start_tick,end_tick])) 
    end   

    for i = 1:4    
        plot!(pl3_1,I₀[i],gₛ_F[i],seriestype=:scatter,markercolor=:blue)        
        plot!(pl3_2,Tmean[i],gₛ_F[i],seriestype=:scatter,markercolor=:blue)        
        plot!(pl3_3,VPD[i],gₛ_F[i],seriestype=:scatter,markercolor=:blue)        
        plot!(pl3_4,θ_F[i],gₛ_F[i],seriestype=:scatter,markercolor=:blue)
        

        plot!(pl4_1,I₀[i],Nₘ_f_F[i]*100,seriestype=:scatter,markercolor=:blue)        
        plot!(pl4_2,Tmean[i],Nₘ_f_F[i]*100,seriestype=:scatter,markercolor=:blue)        
        plot!(pl4_3,VPD[i],Nₘ_f_F[i]*100,seriestype=:scatter,markercolor=:blue)        
        plot!(pl4_4,θ_F[i],Nₘ_f_F[i]*100,seriestype=:scatter,markercolor=:blue) 
        
        plot!(pl5_1,I₀[i],wue_F[i]*1000,seriestype=:scatter,markercolor=:blue)        
        plot!(pl5_2,Tmean[i],wue_F[i]*1000,seriestype=:scatter,markercolor=:blue)        
        plot!(pl5_3,VPD[i],wue_F[i]*1000,seriestype=:scatter,markercolor=:blue)        
        plot!(pl5_4,θ_F[i],wue_F[i]*1000,seriestype=:scatter,markercolor=:blue)

        plot!(pl6,gₛ_F[i],Nₘ_f_F[i]*100,seriestype=:scatter,markercolor=:blue)
    end
    for i = 1:4
        plot!(pl3_1,I₀[i],gₛ_C[i],seriestype=:scatter,markercolor=:red)
        plot!(pl3_2,Tmean[i],gₛ_C[i],seriestype=:scatter,markercolor=:red)
        plot!(pl3_3,VPD[i],gₛ_C[i],seriestype=:scatter,markercolor=:red)
        plot!(pl3_4,θ_C[i],gₛ_C[i],seriestype=:scatter,markercolor=:red)

        plot!(pl4_1,I₀[i],Nₘ_f_C[i]*100,seriestype=:scatter,markercolor=:red)
        plot!(pl4_2,Tmean[i],Nₘ_f_C[i]*100,seriestype=:scatter,markercolor=:red)
        plot!(pl4_3,VPD[i],Nₘ_f_C[i]*100,seriestype=:scatter,markercolor=:red)
        plot!(pl4_4,θ_C[i],Nₘ_f_C[i]*100,seriestype=:scatter,markercolor=:red)

        plot!(pl5_1,I₀[i],wue_C[i]*1000,seriestype=:scatter,markercolor=:red)
        plot!(pl5_2,Tmean[i],wue_C[i]*1000,seriestype=:scatter,markercolor=:red)
        plot!(pl5_3,VPD[i],wue_C[i]*1000,seriestype=:scatter,markercolor=:red)
        plot!(pl5_4,θ_C[i],wue_C[i]*1000,seriestype=:scatter,markercolor=:red)

        plot!(pl6,gₛ_C[i],Nₘ_f_C[i]*100,seriestype=:scatter,markercolor=:red)
    end

    save_fld = "./plots/paper_results/plant_vs_weather"
    isdir(save_fld)|| mkdir(save_fld)

    #Data vs Model output
    pl1_fin = plot(pl1[1],pl1[2],pl1[3],pl1[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)    
    savefig(pl1_fin,save_fld*"/gs_date.svg")

    pl2_fin = plot(pl2[1],pl2[2],pl2[3],pl2[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl2_fin,save_fld*"/Nmf_date.svg")        
   
    plot(pl3_1,pl3_2,pl3_3,pl3_4,
    pl4_1,pl4_2,pl4_3,pl4_4,
    pl5_1,pl5_2,pl5_3,pl5_4,
    layout=(3,4),size=(900,600),margin = 4Plots.mm)
    savefig(save_fld*"/var_weather.svg")

    savefig(pl6,save_fld*"/Nmf_gs.svg")

    #Calculate cor 

    save_fld = "./output/paper_results/plant_vs_weather"
    isdir(save_fld)|| mkdir(save_fld)

    open(save_fld*"/Cor_Fertilized.txt","w") do io
        println(io,"Stand type: $(stand_type_F)")  
        println(io,"--gₛ--")  
        println(io,"I₀: $(round(calc_cor(I₀,gₛ_F),sigdigits=2))")
        println(io,"Tmean: $(round(calc_cor(Tmean,gₛ_F),sigdigits=2))")
        println(io,"VPD: $(round(calc_cor(VPD,gₛ_F),sigdigits=2))")
        println(io,"θ: $(round(calc_cor(θ_F,gₛ_F),sigdigits=2))")
        println(io,"Nₘ_f: $(round(calc_cor(Nₘ_f_F,gₛ_F),sigdigits=2))")
        println(io,"--Nₘ_f--")  
        println(io,"I₀: $(round(calc_cor(I₀,Nₘ_f_F),sigdigits=2))")
        println(io,"Tmean: $(round(calc_cor(Tmean,Nₘ_f_F),sigdigits=2))")
        println(io,"VPD: $(round(calc_cor(VPD,Nₘ_f_F),sigdigits=2))")
        println(io,"θ: $(round(calc_cor(θ_F,Nₘ_f_F),sigdigits=2))")
        println(io,"--wue--") 
        println(io,"I₀: $(round(calc_cor(I₀,wue_F),sigdigits=2))")
        println(io,"Tmean: $(round(calc_cor(Tmean,wue_F),sigdigits=2))")
        println(io,"VPD: $(round(calc_cor(VPD,wue_F),sigdigits=2))")
        println(io,"θ: $(round(calc_cor(θ_F,wue_F),sigdigits=2))")
    end

    open(save_fld*"/Cor_Control.txt","w") do io
        println(io,"Stand type: $(stand_type_C)")  
        println(io,"--gₛ--")  
        println(io,"I₀: $(round(calc_cor(I₀,gₛ_C),sigdigits=2))")
        println(io,"Tmean: $(round(calc_cor(Tmean,gₛ_C),sigdigits=2))")
        println(io,"VPD: $(round(calc_cor(VPD,gₛ_C),sigdigits=2))")
        println(io,"θ: $(round(calc_cor(θ_C,gₛ_C),sigdigits=2))")
        println(io,"Nₘ_f: $(round(calc_cor(Nₘ_f_C,gₛ_C),sigdigits=2))")
        println(io,"--Nₘ_f--") 
        println(io,"I₀: $(round(calc_cor(I₀,Nₘ_f_C),sigdigits=2))")
        println(io,"Tmean: $(round(calc_cor(Tmean,Nₘ_f_C),sigdigits=2))")
        println(io,"VPD: $(round(calc_cor(VPD,Nₘ_f_C),sigdigits=2))")
        println(io,"θ: $(round(calc_cor(θ_C,Nₘ_f_C),sigdigits=2))")
        println(io,"--wue--") 
        println(io,"I₀: $(round(calc_cor(I₀,wue_C),sigdigits=2))")
        println(io,"Tmean: $(round(calc_cor(Tmean,wue_C),sigdigits=2))")
        println(io,"VPD: $(round(calc_cor(VPD,wue_C),sigdigits=2))")
        println(io,"θ: $(round(calc_cor(θ_C,wue_C),sigdigits=2))")
    end

    plot(pl8[1],pl8[2],pl8[3],pl8[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)    
end