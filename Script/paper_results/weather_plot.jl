function create_weather_plot()
    folder_name = "weather"
    folder_path = "./plots/paper_results/"*folder_name
    isdir(folder_path)|| mkdir(folder_path)

    stand_type_F = :Fertilized
    stand_type_C = :Control
    
    raw_input_F = RawInputData(;stand_type=stand_type_F)
    Ec_data_F = calc_Ec_data.(raw_input_F)
    GPP_data_F = get_GPP_data.(raw_input_F;stand_type=stand_type_F)    

    raw_input_C = RawInputData(;stand_type=stand_type_C)
    Ec_data_C = calc_Ec_data.(raw_input_C)
    GPP_data_C = get_GPP_data.(raw_input_C;stand_type=stand_type_C)    

    #daylength = CCPH.daylighthour(data.lat*pi/180,day_nr)*3600 #Seconds
    #I₀ₜₒₜ = data.Radₜₒ*2.3*10^-6 #mol m⁻²   

    #pl1 = plot(xlabel="",ylabel="Δtg (h)",legends=false, ylims = (10,22))
    pl1 = [plot(xlabel="2015",ylabel="Δtg (h)",legends=false, ylims = (10,22),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (10,22),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (10,22),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (10,22),yaxis=false,guidefontsize=12)]
    #pl2 = plot(xlabel="",ylabel="Rₜₒₜ (MJ m⁻²)",legends=false, ylims = (0,30))
    pl2 = [plot(xlabel="2015",ylabel="Rₜₒₜ (MJ m⁻²)",legends=false, ylims = (0,30),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (0,30),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (0,30),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (0,30),yaxis=false,guidefontsize=12)]
    #pl3 = plot(xlabel="",ylabel="Tₐ (°C)",legends=false, ylims = (-6,32),aspect_ratio=10.0,xtickfontsize=16,ytickfontsize=16,yguidefontsize=16,xticks=([CCPH.Dates.Date(2016,01,01),CCPH.Dates.Date(2017,01,01),CCPH.Dates.Date(2018,01,01)],["2016","2017","2018"]))
    pl3 = [plot(xlabel="2015",ylabel="Tₐ (°C)",legends=false, ylims = (-6,32),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (-6,32),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (-6,32),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (-6,32),yaxis=false,guidefontsize=12)]
    #pl4 = plot(xlabel="",ylabel="VPD (kPa)",legends=false, ylims = (0.0,1.5))
    pl4 = [plot(xlabel="2015",ylabel="VPD (kPa)",legends=false, ylims = (0.0,1.5),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (0.0,1.5),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (0.0,1.5),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (0.0,1.5),yaxis=false,guidefontsize=12)]
    #pl5 = plot(xlabel="",ylabel="θF (%)",legends=false, ylims = (6,30),aspect_ratio=10.0,xtickfontsize=16,ytickfontsize=16,yguidefontsize=16,xticks=([CCPH.Dates.Date(2016,01,01),CCPH.Dates.Date(2017,01,01),CCPH.Dates.Date(2018,01,01)],["2016","2017","2018"]))
    pl5 = [plot(xlabel="2015",ylabel="θ (%)",legends=false, ylims = (6,30),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (6,30),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (6,30),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (6,30),yaxis=false,guidefontsize=12)]
    #pl6 = plot(xlabel="",ylabel="θC (%)",legends=false, ylims = (6,30))

    for i = 1:4
        date = [weather.date for weather in raw_input_F[i].weather_growth]
        Tmin = [weather.Tmin for weather in raw_input_F[i].weather_growth]
        Tmax = [weather.Tmax for weather in raw_input_F[i].weather_growth]
        Tmean = [weather.Tmean for weather in raw_input_F[i].weather_growth]
        θ_F = [weather.θₛ*100 for weather in raw_input_F[i].weather_growth]
        θ_C = [weather.θₛ*100 for weather in raw_input_C[i].weather_growth]
        Radₜₒ = [weather.Radₜₒ*10^-6 for weather in  raw_input_F[i].weather_growth] 
        daylength = [CCPH.daylighthour(weather.lat*pi/180,CCPH.Dates.dayofyear(weather.date)) for weather in  raw_input_F[i].weather_growth]
        VPD = [first(get_env_from_data(weather)).VPD/1000 for weather in  raw_input_F[i].weather_growth]         
        
        start_tick = ""#"$(date[1])"
        end_tick = ""#"$(date[end])"

        plot!(pl1[i],date,daylength,linecolor=:blue)
        plot!(pl1[i],xticks=([date[1],date[end]],[start_tick,end_tick]))
        plot!(pl2[i],date,Radₜₒ,linecolor=:blue)
        plot!(pl2[i],xticks=([date[1],date[end]],[start_tick,end_tick]))
        plot!(pl3[i],date,Tmin,linecolor=:green)
        plot!(pl3[i],date,Tmean,linecolor=:blue)
        plot!(pl3[i],date,Tmax,linecolor=:red)        
        plot!(pl3[i],xticks=([date[1],date[end]],[start_tick,end_tick]))
        plot!(pl4[i],date,VPD,linecolor=:blue)
        plot!(pl4[i],xticks=([date[1],date[end]],[start_tick,end_tick]))
        plot!(pl5[i],date,θ_F,linecolor=:blue)
        plot!(pl5[i],date,θ_C,linecolor=:red)
        plot!(pl5[i],xticks=([date[1],date[end]],[start_tick,end_tick]))
        #plot!(pl6,date,θ_C,linecolor=:blue)

        #=
        savefig(pl1,folder_path*"/day_len_$(i).svg")
        savefig(pl2,folder_path*"/rad_$(i).svg")
        savefig(pl3,folder_path*"/temp_$(i).svg")
        savefig(pl4,folder_path*"/VPD_$(i).svg")
        savefig(pl5,folder_path*"/swc_$(i).svg")
        =#
        #savefig(pl6,folder_path*"/swc_C_$(i).svg")
    end   

    pl1_fin = plot(pl1[1],pl1[2],pl1[3],pl1[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl1_fin,folder_path*"/daylight.svg")

    pl2_fin = plot(pl2[1],pl2[2],pl2[3],pl2[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl2_fin,folder_path*"/Rad.svg")

    pl3_fin = plot(pl3[1],pl3[2],pl3[3],pl3[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl3_fin,folder_path*"/temp.svg")  
    
    pl4_fin = plot(pl4[1],pl4[2],pl4[3],pl4[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl4_fin,folder_path*"/VPD.svg") 

    pl5_fin = plot(pl5[1],pl5[2],pl5[3],pl5[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl5_fin,folder_path*"/swc.svg") 
end