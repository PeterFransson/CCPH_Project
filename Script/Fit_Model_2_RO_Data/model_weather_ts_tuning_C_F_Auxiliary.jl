function Calctranpiration(VPD_z::T,Sₑ::T;
    E_cm_ref::T=1.812,s_VPD_z::T=3.121,s_Sₑ::T=18.342) where {T<:Float64}
    #Tor-Ngern et al. 2017 model for estimating canopy transpiration E_c (mm/d) 
    #Default parameter values taken are form the original paper
    #for estimating E_c for Rosinedal Scots pine stand
    #---Input---
    #VPD_z day-length normalized vapor pressure deficit (VPD*Day length/24, Pa) 
    #Sₑ effective saturation (-)
    E_cm = E_cm_ref*(1-exp(-s_VPD_z*VPD_z*10^-3))
    E_c = E_cm*(1-exp(-s_Sₑ*Sₑ))

    return E_c
end

function Est_gₛ(LAI::T,daylight::T;cons::CCPH.Constants=CCPH.Constants(),
    env::CCPH.EnvironmentStruct=CCPH.EnvironmentStruct(),
    E_cm_ref::T=1.812,s_VPD_z::T=3.121,s_Sₑ::T=18.342) where {T<:Float64} 
    #Estimating stomatal conductance gₛ (mol m⁻² leaf area s⁻¹) based on weather data,
    #using the Tor-Ngern et al. 2017 model for estimating canopy transpiration 
    #---Input--- 
    # LAI leaf area index (m² leaf area m⁻² ground area)    
    # daylight lenght of day (s)
    VPD = env.VPD
    θₛ = env.θₛ
    P = env.P    
    VPD_z = VPD.*daylight/(24*3600)    
    Sₑ = CCPH.CalcSₑ.(θₛ)
    E_c = Calctranpiration.(VPD_z,Sₑ;E_cm_ref=E_cm_ref,s_VPD_z=s_VPD_z,s_Sₑ=s_Sₑ) #mm/day
    E_c *= 10^-3*cons.ρ_H2O/cons.M_H2O/daylight #mol m⁻² s⁻¹
    gₛ = E_c*P/(VPD*cons.r*LAI)
    return gₛ
end

function Find_data_ind(i::Integer)
    if 1<=i<=21
        return 1
    elseif 22<=i<=44
        return 2
    elseif 45<=i<=64
        return 3
    else
        return 4
    end
end

mutable struct RoData{T<:Float64}
    H::Array{T,1}
    B::Array{T,1}
    GPP::Array{T,1}
    Wf::Array{T,1}
    N::Array{T,1}
    Hc::T
    gₛ::Array{T,1}
end

function Create_RoData_C_F(GPP_data_F::Array{T,1},GPP_data_C::Array{T,1},
    weatherts_F::CCPH.WeatherTS,weatherts_C::CCPH.WeatherTS;treepar::TreePar=TreePar()) where {T<:Float64} 
        
    H_data_F = [19.07,19.34,19.64,19.87] 
    B_data_F = [0.0452,0.0467,0.0479,0.0487]
    Wf_data_F = [8.33,8.47,8.49,8.48]
    N_data_F = [850.0000,846.6667,846.6667,846.6667]/10000
    LAI_F = Wf_data_F.*N_data_F/treepar.LMA

    Hc_F = 10.89

    H_data_C = [20.86,21.01,21.19,21.36] 
    B_data_C = [0.0348,0.0356,0.0362,0.0368]
    Wf_data_C = [5.025,5.110,5.141,5.165]
    N_data_C = [1010.0000,1010.0000,1006.6667,1006.6667]/10000
    Hc_C = 10.74
    LAI_C = Wf_data_C.*N_data_C/treepar.LMA

    gₛ_data_F = Float64[]
    gₛ_data_C = Float64[]
    for i = 1:length(weatherts_F.date)
        env_F = EnvironmentStruct(weatherts_F,i)
        env_C = EnvironmentStruct(weatherts_C,i)

        daylight_F = weatherts_F.daylight[i]/7
        daylight_C = weatherts_C.daylight[i]/7

        push!(gₛ_data_F,Est_gₛ(LAI_F[Find_data_ind(i)],daylight_F;env=env_F))
        push!(gₛ_data_C,Est_gₛ(LAI_C[Find_data_ind(i)],daylight_C;env=env_C))
    end

    data_F = RoData(H_data_F,B_data_F,GPP_data_F,Wf_data_F,N_data_F,Hc_F,gₛ_data_F)
    data_C = RoData(H_data_C,B_data_C,GPP_data_C,Wf_data_C,N_data_C,Hc_C,gₛ_data_C)

    return data_F,data_C
end

function calcR²(y_data::Array{Float64,1},f_model::Array{Float64,1})
    ymean = mean(y_data)
    SStot = sum((y_data.-ymean).^2)
    SSres = sum((y_data.-f_model).^2)
    return 1-SSres/SStot
end

function clac_GPP_R²_annual(GPP::Array{Float64,1},GPP_model::Array{Float64,1},dates::Array{Dates.DateTime,1})
    years = [Dates.year(dates[1])]    
    growth_end = Int64[]
    for i = 2:length(dates)
        if Dates.year(dates[i]) != Dates.year(dates[i-1])
            push!(growth_end,i-1)
            push!(years,Dates.year(dates[i]))
        end
    end
    push!(growth_end,length(dates))

    @show growth_end
    @show years

    println("Year: ",years[1])        
    println("R²: ",calcR²(GPP[1:growth_end[1]],GPP_model[1:growth_end[1]]))      
    println("Cor: ",cor(GPP[1:growth_end[1]],GPP_model[1:growth_end[1]])) 

    for i in 2:length(years)
        println("Year: ",years[i])        
        println("R²: ",calcR²(GPP[growth_end[i-1]+1:growth_end[i]],GPP_model[growth_end[i-1]+1:growth_end[i]]))
        println("Cor: ",cor(GPP[growth_end[i-1]+1:growth_end[i]],GPP_model[growth_end[i-1]+1:growth_end[i]]))
    end
end

function Max_k_cost(model::CCPHStruct)
    ψ₅₀,b,g,ρ_H2O,θₛ = model.hydPar.ψ₅₀,model.hydPar.b,model.cons.g,model.cons.ρ_H2O,model.env.θₛ
    #Calculate tree height
    H = model.treesize.H
    #Caluclate soil water potential
    ψₛ = CCPH.θₛ2ψₛ(θₛ)
    #Soil water potential adjusted for gravitational pressure (MPa)
    ψₛ_g = ψₛ-H*ρ_H2O*g*10^-6 

    k_cost_max = CCPH.Pfun(ψₛ_g,ψ₅₀,b)

    return k_cost_max
end

function OptimCCPH(growthlength::Float64,model::CCPH.CCPHStruct;K_cost_crit::Float64 = 0.12)
    
    gₛ_crit = CCPH.Calc_K_costᵢₙᵥ(K_cost_crit,model)

    gₛ_opt,Nₘ_f_opt = CCPH.CCPHTraitmodel(growthlength,model;gₛ_lim_hi=gₛ_crit)

    modeloutput = CCPH.CCPH_run(gₛ_opt,Nₘ_f_opt,growthlength,model)

    return modeloutput,gₛ_opt,Nₘ_f_opt
end
function OptimCCPH!(i::Integer,model::CCPH.CCPHStruct,
    weatherts::CCPH.WeatherTS,photo_kinetic::CCPH.PhotoKineticRates)
    growthlength,step_length = CCPH.Init_weather_par!(i,model,weatherts,photo_kinetic)
    modeloutput,gₛ_opt,Nₘ_f_opt = OptimCCPH(growthlength,model)

    return modeloutput,gₛ_opt,Nₘ_f_opt
end

#Calc model GPP
function CalcGPP(ccphts::CCPHTS,weatherts::WeatherTS)
    GPP = ccphts.N.*ccphts.P.*weatherts.daylight./weatherts.tot_annual_daylight*1000/7 #g C day⁻¹ m⁻² ground area

    return GPP
end
function CalcGPP(step_length::Float64,model::CCPH.CCPHStruct,modeloutput::CCPH.CCPHOutput)
    GPP = model.treesize.N*modeloutput.P*step_length*1000/7 #g C day⁻¹ m⁻² ground area

    return GPP
end
function CalcGPP!(i::Integer,model::CCPH.CCPHStruct,
    weatherts::CCPH.WeatherTS,photo_kinetic::CCPH.PhotoKineticRates)
    growthlength,step_length = CCPH.Init_weather_par!(i,model,weatherts,photo_kinetic)
    modeloutput,gₛ_opt,Nₘ_f_opt = OptimCCPH(growthlength,model)
    GPP = CalcGPP(step_length,model,modeloutput)
    
    return GPP, gₛ_opt
end

function clac_GPP_R²_annual(GPP::Array{Float64,1},GPP_model::Array{Float64,1},dates::Array{Dates.DateTime,1})
    years = [Dates.year(dates[1])]    
    growth_end = Int64[]
    for i = 2:length(dates)
        if Dates.year(dates[i]) != Dates.year(dates[i-1])
            push!(growth_end,i-1)
            push!(years,Dates.year(dates[i]))
        end
    end
    push!(growth_end,length(dates))

    @show growth_end
    @show years

    println("Year: ",years[1])        
    println("R²: ",calcR²(GPP[1:growth_end[1]],GPP_model[1:growth_end[1]]))      
    println("Cor: ",cor(GPP[1:growth_end[1]],GPP_model[1:growth_end[1]])) 

    for i in 2:length(years)
        println("Year: ",years[i])        
        println("R²: ",calcR²(GPP[growth_end[i-1]+1:growth_end[i]],GPP_model[growth_end[i-1]+1:growth_end[i]]))
        println("Cor: ",cor(GPP[growth_end[i-1]+1:growth_end[i]],GPP_model[growth_end[i-1]+1:growth_end[i]]))
    end
end

function Est_C_assimilation(GPP_Data::T,LAI::T,N::T,weatherts::CCPH.WeatherTS,
    ind::Integer;
    treepar::CCPH.TreePar=CCPH.TreePar(),
    cons::CCPH.Constants=CCPH.Constants()) where {T<:Float64}
    #Estimate per leaf area carbon assimilation rate A (mol C m⁻² leaf area s⁻¹)
    #based on GPP data (g C m⁻² ground area day⁻¹)
    Xₜ = weatherts.acclimation_fac[ind]
    daylight = weatherts.daylight[ind]

    GPP = GPP_Data*10^-3/(daylight/7) #Kg C m⁻² ground area day⁻¹
    P = GPP/N #Kg C tree⁻¹ s⁻¹
    A = P*(N*treepar.k)/(1-exp(-treepar.k*LAI))/(Xₜ*cons.M_C) #mol C m⁻² leaf area s⁻¹

    return A
end