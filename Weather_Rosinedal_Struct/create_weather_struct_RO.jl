struct RO_raw_data{T<:CSV.File}
    Weather_RO::T 
    GPP_data_RO::T
    SWC_data_RO::T 
end 

function Load_RO_weather_data(;stand_type::String="Fertilized")
    SWC_data = CSV.read("./Data/RO_data/Daily_SWC_data_$(year).csv", DataFrames.DataFrame)
    Weather_data = CSV.read("./Data/RO_data/Weather_data_$(year).csv", DataFrames.DataFrame)

    if stand_type=="Fertilized"

    elseif stand_type=="Control"
    else           
        error("Wrong input. Either \"Fertilized\" or \"Control\"") 
    end
end

function Load_RO_data(;weather_file::String="Weather_RO")
    Weather_RO = CSV.File("./Data/"*weather_file*".csv")
    GPP_data_RO = CSV.File("./Data/GPP_RO.csv")
    SWC_data_RO = CSV.File("./Data/Daily_SWC_RO.csv")  

    RO_data = RO_raw_data(Weather_RO,GPP_data_RO,SWC_data_RO)

    return RO_data
end

function create_acclimation(;X0 = -4.0,τ =  7.0,Smax = 16.0,Weather_RO_Temp)
    #Weather_RO = CSV.File("./Data/Weather_RO.csv")
    #Weather_RO_Temp = Weather_RO[8767:end]

    #photosynthesis that accounts for temperature acclimation (Mäkelä et al., 2004; Mäkelä et al., 2008). 
    #X0 = -4.0 #minium tempreture for Photosynthesis    
    #τ =  7.0 #Days
    #Smax = 16.0
    X_vec = Float64[]  
    X_growth_vec = Float64[]   
    X_growth_week_vec = Float64[]       
    Dates_growth_vec = Dates.DateTime[] 

    push!(X_vec,Weather_RO_Temp[1][2])
    if Weather_RO_Temp[1][11]==1
        push!(X_growth_vec,Weather_RO_Temp[1][2])
        push!(Dates_growth_vec,Weather_RO_Temp[1][1])
    end
    for i = 2:length(Weather_RO_Temp)
        X = X_vec[i-1]+(Weather_RO_Temp[i][2]-X_vec[i-1])/τ
        push!(X_vec,X)
        if Weather_RO_Temp[i][11]==1
            push!(X_growth_vec,X)
            push!(Dates_growth_vec,Weather_RO_Temp[i][1])
        end
    end 

    S_growth_vec = max.(X_growth_vec.-X0,Ref(0)) 
    fs_growth_vec = min.(S_growth_vec/Smax,Ref(1))
    S_growth_week_vec = Float64[]   
    fs_growth_week_vec = Float64[] 

    i = 1
    while i <=length(X_growth_vec)-6
        if Dates.year(Dates_growth_vec[i])==Dates.year(Dates_growth_vec[i+6])
            mean_val = mean(X_growth_vec[i:i+6])
            push!(X_growth_week_vec,mean_val)  
            push!(S_growth_week_vec,mean(S_growth_vec[i:i+6]))  
            push!(fs_growth_week_vec,mean(fs_growth_vec[i:i+6]))         
            i += 7
        else
            i += 1
        end        
    end

    return fs_growth_week_vec
end

function create_weather_struct_RO(RO_data::RO_raw_data;stand_type::String="Fertilized",X0 = -4.0,τ =  7.0,Smax = 16.0,
    weather_file::String="Weather_RO")
    #Creates weekly mean time series for Rosinedal 2015-2018 (Or 2015-2019 if weather_file=Weather_RO2)
    if stand_type == "Fertilized"
        GPP_ind = 2
        SWC_data_RO_ind = 3
    elseif stand_type == "Control"
        GPP_ind = 3
        SWC_data_RO_ind = 4
    else
        error("Wrong input. Either \"Fertilized\" or \"Control\"")
    end    

    Weather_RO_Temp = RO_data.Weather_RO[8767:end]
    SWC_data_RO_Temp =  RO_data.SWC_data_RO[366:end]

    Radiation_vec = Float64[]
    PAR_vec = Float64[]
    VP_vec = Float64[]
    Date_vec = Dates.DateTime[]
    GPP_vec = Float64[]
    Tmean_vec = Float64[]
    θₛ_vec = Float64[]   
    daylight_vec = Float64[]

    for i in 1:length(Weather_RO_Temp)
        row = Weather_RO_Temp[i]
        if row[11]==1
            push!(Date_vec,row[1])
            push!(Tmean_vec,row[2])
            push!(GPP_vec,RO_data.GPP_data_RO[i][GPP_ind])
            push!(θₛ_vec,SWC_data_RO_Temp[i][SWC_data_RO_ind]/100)

            val = tryparse(Float64,row[7])
            if val !== nothing
                push!(Radiation_vec,val)  
                d = Dates.dayofyear(row[1])           
                h = CCPH.daylighthour(64*pi/180,d)*3600          
                push!(PAR_vec,val*2.6/h)
                push!(daylight_vec,h)
            else
                println("Radiation data is missing at date: $(row[1]), approximating using  linear regression model")
                β = [0.0003977611889578929; -10.463723925716433]
                d = Dates.dayofyear(row[1])           
                h = CCPH.daylighthour(64*pi/180,d)*3600  
                Radiation_hat =  β[1]*h+β[2]
                push!(Radiation_vec,Radiation_hat)                         
                push!(PAR_vec,Radiation_hat*2.3/h)
                push!(daylight_vec,h)
            end     

            val = row[8]
            push!(VP_vec,val*100)
            if isnan(val)||ismissing(val)
                println("Vapor deficit data is missing (either NaN oc Missing) at date: $(row[1])")
            end          
        end
    end    

    VPD_vec = CCPH.SVPₜ.(Tmean_vec).-VP_vec

    Date_week_vec = Dates.DateTime[]
    GPP_week_vec = Float64[]
    PAR_week_vec = Float64[]    
    Tmean_week_vec = Float64[]
    VPD_week_vec = Float64[]
    θₛ_week_vec = Float64[]
    daylight_vec = Float64[]

    i = 1
    while i<=length(Date_vec)-6 
        if Dates.year(Date_vec[i])==Dates.year(Date_vec[i+6])
            push!(Date_week_vec,Date_vec[i])
            push!(Tmean_week_vec,mean(Tmean_vec[i:i+6]))
            push!(PAR_week_vec,mean(PAR_vec[i:i+6]))
            push!(VPD_week_vec,mean(VPD_vec[i:i+6]))
            push!(θₛ_week_vec,mean(θₛ_vec[i:i+6]))
            push!(GPP_week_vec,mean(GPP_vec[i:i+6]))
            d_vec = Dates.dayofyear.(Date_vec[i:i+6])
            push!(daylight_vec,sum(CCPH.daylighthour.(Ref(64*pi/180),d_vec))*3600)
            i += 7
        else
            i += 1    
        end
    end

    acclimation_fac = create_acclimation(;X0 = X0,τ =  τ,Smax = Smax,Weather_RO_Temp)    

    # Find the end of each growth period   
    if weather_file=="Weather_RO_2" 
        t_end_vec = [0,0,0,0,0]
        for j in 1:length(Date_week_vec)  
            if Dates.year(Date_week_vec[j]) == 2015 
                t_end_vec[1] += 1  
                t_end_vec[2] += 1  
                t_end_vec[3] += 1 
                t_end_vec[4] += 1   
                t_end_vec[5] += 1 
            elseif Dates.year(Date_week_vec[j]) == 2016 
                t_end_vec[2] += 1  
                t_end_vec[3] += 1 
                t_end_vec[4] += 1   
                t_end_vec[5] += 1 
            elseif Dates.year(Date_week_vec[j]) == 2017             
                t_end_vec[3] += 1 
                t_end_vec[4] += 1 
                t_end_vec[5] += 1   
            elseif Dates.year(Date_week_vec[j]) == 2018   
                t_end_vec[4] += 1
                t_end_vec[5] += 1   
            else 
                t_end_vec[5] += 1 
            end        
        end

        # calc total daylight for each year
        tot_daylight_hour = [0.0,0.0,0.0,0.0,0.0]
        tot_daylight_hour[1] = sum(daylight_vec[1:t_end_vec[1]])
        tot_daylight_hour[2] = sum(daylight_vec[t_end_vec[1] + 1:t_end_vec[2]])
        tot_daylight_hour[3] = sum(daylight_vec[t_end_vec[2] + 1:t_end_vec[3]])
        tot_daylight_hour[4] = sum(daylight_vec[t_end_vec[3] + 1:t_end_vec[4]])
        tot_daylight_hour[5] = sum(daylight_vec[t_end_vec[4] + 1:t_end_vec[5]])
            
        # Create total daylight hour time series
        tot_daylight_hour_ts = zeros(length(daylight_vec))
        tot_daylight_hour_ts[1:t_end_vec[1]] .= tot_daylight_hour[1]
        tot_daylight_hour_ts[t_end_vec[1] + 1:t_end_vec[2]] .= tot_daylight_hour[2]
        tot_daylight_hour_ts[t_end_vec[2] + 1:t_end_vec[3]] .= tot_daylight_hour[3]
        tot_daylight_hour_ts[t_end_vec[3] + 1:t_end_vec[4]] .= tot_daylight_hour[4]
        tot_daylight_hour_ts[t_end_vec[4] + 1:t_end_vec[5]] .= tot_daylight_hour[5]
    else
        t_end_vec = [0,0,0,0]
        for j in 1:length(Date_week_vec)  
            if Dates.year(Date_week_vec[j]) == 2015 
                t_end_vec[1] += 1  
                t_end_vec[2] += 1  
                t_end_vec[3] += 1 
                t_end_vec[4] += 1   
            elseif Dates.year(Date_week_vec[j]) == 2016 
                t_end_vec[2] += 1  
                t_end_vec[3] += 1 
                t_end_vec[4] += 1   
            elseif Dates.year(Date_week_vec[j]) == 2017             
                t_end_vec[3] += 1 
                t_end_vec[4] += 1   
            else
                t_end_vec[4] += 1   
            end        
        end   
    
        # calc total daylight for each year
        tot_daylight_hour = [0.0,0.0,0.0,0.0]
        tot_daylight_hour[1] = sum(daylight_vec[1:t_end_vec[1]])
        tot_daylight_hour[2] = sum(daylight_vec[t_end_vec[1] + 1:t_end_vec[2]])
        tot_daylight_hour[3] = sum(daylight_vec[t_end_vec[2] + 1:t_end_vec[3]])
        tot_daylight_hour[4] = sum(daylight_vec[t_end_vec[3] + 1:t_end_vec[4]])
            
        # Create total daylight hour time series
        tot_daylight_hour_ts = zeros(length(daylight_vec))
        tot_daylight_hour_ts[1:t_end_vec[1]] .= tot_daylight_hour[1]
        tot_daylight_hour_ts[t_end_vec[1] + 1:t_end_vec[2]] .= tot_daylight_hour[2]
        tot_daylight_hour_ts[t_end_vec[2] + 1:t_end_vec[3]] .= tot_daylight_hour[3]
        tot_daylight_hour_ts[t_end_vec[3] + 1:t_end_vec[4]] .= tot_daylight_hour[4]
    end

    weatherts = WeatherTS(Date_week_vec,daylight_vec,tot_daylight_hour_ts,PAR_week_vec,
    Tmean_week_vec,VPD_week_vec,θₛ_week_vec,acclimation_fac)

    return weatherts, GPP_week_vec
end

function plot_weather_struct_RO(weatherts::WeatherTS,GPP_week_vec::Array{Float64,1})
    pl1 = plot(weatherts.date,weatherts.daylight/3600,xlabel="Date",ylabel="Daylight (h)",seriestype=:scatter)
    pl2 = plot(weatherts.date,weatherts.PAR*10^6,xlabel="Date",ylabel="PAR (μmol s⁻¹ m⁻²)",seriestype=:scatter)
    pl3 = plot(weatherts.date,weatherts.temp,xlabel="Date",ylabel="Temp (ᵒC)",seriestype=:scatter)
    pl4 = plot(weatherts.date,weatherts.VPD/1000,xlabel="Date",ylabel="VPD (kPa)",seriestype=:scatter)
    pl5 = plot(weatherts.date,weatherts.θₛ*100,xlabel="Date",ylabel="SWC (%)",seriestype=:scatter)
    pl6 = plot(weatherts.date,CCPH.θₛ2ψₛ.(weatherts.θₛ),xlabel="Date",ylabel="SWP (MPa)",seriestype=:scatter)
    pl7 = plot(weatherts.date,GPP_week_vec,xlabel="Date",ylabel="GPP (g C m⁻² d⁻¹)",seriestype=:scatter)
    pl8 = plot(weatherts.date,weatherts.acclimation_fac,xlabel="Date",ylabel="Acclimation",seriestype=:scatter)

    plot(pl1,pl2,pl3,pl4,pl5,pl6,pl7,pl8,layout=8,legends=false)   
end
