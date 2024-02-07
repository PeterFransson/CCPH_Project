mutable struct RawInputData
    weather_raw::Vector{CCPH.WeatherDataStruct}
    weather_growth::Vector{CCPH.WeatherDataStruct}
    growth_indices::Tuple{Integer,Integer}
    growth_indices_weekly::Vector{Tuple{Integer,Integer}}
    treesize::TreeSize
    ζ::Real
end

function CCPH.WeatherDataStruct(data::DataFrames.DataFrame,
    data_idx::Integer;
    lat::Real=64,
    Cₐ::Real=400.0/10.0,
    P::Real=1.0*10^5)
    
    d = data.Date[data_idx]
    data_day = CCPH.WeatherDataStruct(d,
    lat,
    data.airTmean[data_idx],
    data.airTmin[data_idx],
    data.airTmax[data_idx]
    ,data.VP[data_idx]*100,
    data.Radiation[data_idx]*10^6,
    data.SWC[data_idx]/100,
    Cₐ,
    P)
    return data_day
end

function Load_RO_weather_data(year::Integer;stand_type::Symbol=:Fertilized)
    SWC_data = CSV.read("./Data/RO_data/Daily_SWC_data_$(year).csv", DataFrames.DataFrame)
    Weather_data = CSV.read("./Data/RO_data/Weather_data_$(year).csv", DataFrames.DataFrame)

    if stand_type==:Fertilized
        Weather_data[!,:SWC] = SWC_data.SWRos2         
    elseif stand_type==:Control
        Weather_data[!,:SWC] = SWC_data.SWRos3
    else           
        error("Wrong input. Either :Fertilized or :Control") 
    end

    n_rows =  DataFrames.nrow(Weather_data)
    raw_data = [CCPH.WeatherDataStruct(Weather_data,i) for i in 1:n_rows]

    return raw_data
end 
function Load_RO_weather_data(;stand_type::Symbol=:Fertilized)
    raw_data = Vector{Tuple{CCPH.WeatherDataStruct,Integer}}[]

    for year = 2015:2018
        data_struct_year = Load_RO_weather_data(year;stand_type=stand_type)
        push!(raw_data,data_struct_year)        
    end
    return raw_data
end

function get_weekly_indices(raw_data::Vector{CCPH.WeatherDataStruct})
    weekly_indices = Tuple{Integer,Integer}[]
    n_days = length(raw_data)
    current_day = 1
    growth_start_flag = false
    growth_start = 0
    growth_end = 0
    day_count = 0
    week_start = 0
    week_end = 0
    while current_day≤n_days
        if !isnan(raw_data[current_day].Radₜₒ)             
            if day_count==0
                week_start=current_day
                day_count+=1
                if !growth_start_flag
                    growth_start_flag = true
                    growth_start = current_day
                end
            elseif day_count==6
                week_end=current_day
                push!(weekly_indices,(week_start-growth_start+1,week_end-growth_start+1))
                day_count=0  
                if growth_start_flag
                    growth_end = current_day
                end          
            else
                day_count+=1
            end           
        end
        current_day+=1
    end
    growth_indices = (growth_start,growth_end)
    return weekly_indices,growth_indices
end

#Get stand size data and Scaling factor to upscale copy GPP to ecosystem GPP
function get_tree_size(year::Integer;stand_type::Symbol=:Fertilized,treepar::CCPH.TreePar=CCPH.TreePar())
    size_idx = Dict(2015=>1,2016=>2,2017=>3,2018=>4)

    if stand_type==:Fertilized
        
        #Stand size data Fertilized       
        H_data = [19.07,19.34,19.64,19.87] #Canopy hight (m)
        N_data = [850.0000,846.6667,846.6667,846.6667]/10000 #Stem density (# trees m⁻² ground area) 
        LAI_data = [2.45,2.42,2.38,2.44] #Above ground leaf area index (m² m⁻²) 
        K_g = 0.69 #extinction coefficient for ground vegetation layer (Tian et al. 2021)
        LAI_g = 1.0 #LAI of ground vegetation (Tian et al. 2021)       
                
    elseif stand_type==:Control

        #Stand size data Control
        H_data = [20.86,21.01,21.19,21.36] #Canopy hight (m)
        N_data = [1010.0000,1010.0000,1006.6667,1006.6667]/10000 #Stem density (# trees m⁻² ground area)
        LAI_data = [2.3,2.28,2.3,2.21] #Above ground leaf area index (m² m⁻²) 
        K_g = 0.69 #extinction coefficient for ground vegetation layer (Tian et al. 2021)
        LAI_g = 0.52 #LAI of ground vegetation (Tian et al. 2021)
        
    else
        error("Wrong input. Either :Fertilized or :Control") 
    end
    
    treesize = CCPH.TreeSize(H_data[size_idx[year]],LAI_data[size_idx[year]],N_data[size_idx[year]])

    LAI = LAI_data[size_idx[year]]

    f_par_c = 1.0.-exp(-treepar.k*LAI) #Fraction of above canopy light absorbed by canopy layer
    f_par_g = (1.0-f_par_c)*(1-exp(-K_g*LAI_g)) #Fraction of above canopy light absorbed by ground vegetation layer
    ζ = (f_par_c+f_par_g)/f_par_c #Scaling factor to upscale copy GPP to ecosystem GPP

    return treesize,ζ
end

function RawInputData(year::Integer;stand_type::Symbol=:Fertilized,treepar::CCPH.TreePar=CCPH.TreePar())
    raw_weather_data = Load_RO_weather_data(year;stand_type=stand_type)
    weekly_growth_indices,growth_indices = get_weekly_indices(raw_weather_data)
    treesize,ζ =  get_tree_size(year;stand_type=stand_type,treepar=treepar)
    return RawInputData(raw_weather_data,
    raw_weather_data[growth_indices[1]:growth_indices[2]],
    growth_indices,
    weekly_growth_indices,
    treesize,
    ζ)
end
function RawInputData(;stand_type::Symbol=:Fertilized,treepar::CCPH.TreePar=CCPH.TreePar())
    return [RawInputData(year;stand_type=stand_type,treepar=treepar) for year in 2015:2018]
end

#The delayed temperature, Sₜ, is modelled using a first order delay dynamics model
Sₜ_fun(Tₜ::Real,Sₜ₋₁::Real,τ::Real) = (1-1/τ)*Sₜ₋₁+Tₜ/τ
function Sₜ_fun(data::Vector{CCPH.WeatherDataStruct},τ::Real)
    n_data = length(data)
    S₁ = data[1].Tmean
    S_vec = [S₁]
    for i in 2:n_data
        Tₜ = data[i].Tmean   
        Sₜ₋₁ = S_vec[i-1]
        push!(S_vec,Sₜ_fun(Tₜ,Sₜ₋₁,τ))
    end
    return S_vec
end

#Photosynthetic temperature acclimation factor Xₜ (Mäkelä et al., 2004; Mäkelä et al., 2008).
function Xₜ_fun(Sₜ::Real,Smin::Real,ΔS::Real)
    if Sₜ≤Smin
        return 0
    elseif Smin<Sₜ<Smin+ΔS
        return (Sₜ-Smin)/ΔS
    elseif Sₜ≥Smin+ΔS
        return 1
    else
        error("Sₜ")
    end
end
function Xₜ_fun(raw_data::Vector{CCPH.WeatherDataStruct};Smin::Real = -4.0,ΔS::Real = 16.0,τ::Real =  7.0)
    #Calcualte photosynthetic temperature acclimation factor Xₜ (Mäkelä et al., 2004; Mäkelä et al., 2008). 
    #Smin = -4.0 #minium tempreture for Photosynthesis    
    #τ =  7.0 #Days
    #Smax = 16.0  
    
    S_vec = Sₜ_fun(raw_data,τ)
    X_vec = Xₜ_fun.(S_vec,Ref(Smin),Ref(ΔS))
    return X_vec
end
function Xₜ_fun(raw_data::RawInputData;Smin::Real = -4.0,ΔS::Real = 16.0,τ::Real =  7.0)
    Xₜ_raw = Xₜ_fun(raw_data.weather_raw;Smin=Smin,ΔS=ΔS,τ=τ)
    return Xₜ_raw[raw_data.growth_indices[1]:raw_data.growth_indices[2]]
end

meanVPD(data::CCPH.WeatherDataStruct) = CCPH.VPDₜ(data.Tmean,data.VP)
get_daylength(data::CCPH.WeatherDataStruct) = CCPH.daylighthour(data.lat*pi/180,CCPH.Dates.dayofyear(data.date))*3600

function calc_Ec_data(VPD_z::Real,Sₑ::Real;
    E_cm_ref::Real=1.812,s_VPD_z::Real=3.121,s_Sₑ::Real=18.342)
    #Tor-Ngern et al. 2017 model for estimating canopy transpiration E_c (mm/d) 
    #Default parameter values taken are form the original paper
    #for estimating E_c for Rosinedal Scots pine stand (both fertilized and control)
    #---Input---
    #VPD_z day-length normalized vapor pressure deficit (VPD*Day length/24, Pa) 
    #Sₑ effective saturation (-)
    E_cm = E_cm_ref*(1-exp(-s_VPD_z*VPD_z*10^-3))
    Ec = E_cm*(1-exp(-s_Sₑ*Sₑ))

    return Ec
end
function calc_Ec_data(data::CCPH.WeatherDataStruct;
    E_cm_ref::Real=1.812,
    s_VPD_z::Real=3.121,
    s_Sₑ::Real=18.342)

    VPD = meanVPD(data)
    daylength = get_daylength(data)    
    θₛ = data.θₛ       
    VPD_z = VPD*daylength/(24*3600)    
    Sₑ = CCPH.calcSₑ(θₛ)

    return calc_Ec_data(VPD_z,Sₑ;E_cm_ref=E_cm_ref,s_VPD_z=s_VPD_z,s_Sₑ=s_Sₑ)
end
function calc_Ec_data(data::Vector{CCPH.WeatherDataStruct};
    E_cm_ref::Real=1.812,
    s_VPD_z::Real=3.121,
    s_Sₑ::Real=18.342)

    return calc_Ec_data.(data;E_cm_ref=E_cm_ref,s_VPD_z=s_VPD_z,s_Sₑ=s_Sₑ)
end
function calc_Ec_data(data::RawInputData;
    E_cm_ref::Real=1.812,
    s_VPD_z::Real=3.121,
    s_Sₑ::Real=18.342)

    return calc_Ec_data(data.weather_growth;E_cm_ref=E_cm_ref,s_VPD_z=s_VPD_z,s_Sₑ=s_Sₑ)
end

function get_GPP_data(year::Integer;stand_type::Symbol=:Fertilized)
   GPP_data = CSV.read("./Data/RO_data/GPP_data_$(year).csv", DataFrames.DataFrame)

    if stand_type==:Fertilized
        return GPP_data.GPP_RO2       
    elseif stand_type==:Control
        return GPP_data.GPP_RO3
    else           
        error("Wrong input. Either :Fertilized or :Control") 
    end
end
function get_GPP_data(data::RawInputData;stand_type::Symbol=:Fertilized)
    year = Dates.year(data.weather_growth[1].date)
    GPP_raw = get_GPP_data(year;stand_type=stand_type)
    return GPP_raw[data.growth_indices[1]:data.growth_indices[2]]
end
