function create_weather_plot()
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

    pl1 = plot(xlabel="",ylabel="Δtg (h)",legends=false)
    pl2 = plot(xlabel="",ylabel="Rₜₒₜ (MJ m⁻²)",legends=false)
    pl3 = plot(xlabel="",ylabel="Tₐ (°C)",legends=false)
    pl4 = plot(xlabel="",ylabel="VPD (kPa)",legends=false)
    pl5 = plot(xlabel="",ylabel="θF (%)",legends=false)
    pl6 = plot(xlabel="",ylabel="θC (%)",legends=false)

    for i = 1:4
        date = [weather.date for weather in  raw_input_F[i].weather_growth]
        Tmin = [weather.Tmin for weather in  raw_input_F[i].weather_growth]
        Tmax = [weather.Tmax for weather in  raw_input_F[i].weather_growth]
        Tmean = [weather.Tmean for weather in  raw_input_F[i].weather_growth]
        θ_F = [weather.θₛ*100 for weather in  raw_input_F[i].weather_growth]
        θ_C = [weather.θₛ*100 for weather in  raw_input_C[i].weather_growth]
        Radₜₒ = [weather.Radₜₒ*10^-6 for weather in  raw_input_F[i].weather_growth] 
        daylength = [CCPH.daylighthour(weather.lat*pi/180,CCPH.Dates.dayofyear(weather.date)) for weather in  raw_input_F[i].weather_growth]
        VPD = [first(get_env_from_data(weather)).VPD/1000 for weather in  raw_input_F[i].weather_growth]   

        plot!(pl1,date,daylength,linecolor=:blue)
        plot!(pl2,date,Radₜₒ,linecolor=:blue)
        plot!(pl3,date,Tmin,linecolor=:green)
        plot!(pl3,date,Tmean,linecolor=:blue)
        plot!(pl3,date,Tmax,linecolor=:red)
        plot!(pl4,date,VPD,linecolor=:blue)
        plot!(pl5,date,θ_F,linecolor=:blue)
        plot!(pl6,date,θ_C,linecolor=:blue)
    end  

    #plot(pl1,pl2,pl3,pl4,pl5,pl6,layout=(2,3))
    plot(pl3)
end