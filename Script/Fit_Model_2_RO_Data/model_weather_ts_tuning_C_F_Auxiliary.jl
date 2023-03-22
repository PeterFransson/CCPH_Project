function Calctranpiration(VPD_z::T,Sₑ::T;
    E_cm_ref::T=1.812,s_VPD_z::T=3.121,s_Sₑ::T=18.342) where {T<:Float64}
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

function Est_Ec(daylight::T;
    env::CCPH.EnvironmentStruct=CCPH.EnvironmentStruct(),
    E_cm_ref::T=1.812,s_VPD_z::T=3.121,s_Sₑ::T=18.342) where {T<:Float64}
    #Estimating canopy transpiration E_C (mm/d) based on weather data,
    #using the Tor-Ngern et al. 2017 model for estimating canopy transpiration 
    #---Input---
    # daylight lenght of day (s)
    VPD = env.VPD
    θₛ = env.θₛ       
    VPD_z = VPD.*daylight/(24*3600)    
    Sₑ = CCPH.CalcSₑ.(θₛ)
    Ec = Calctranpiration.(VPD_z,Sₑ;E_cm_ref=E_cm_ref,s_VPD_z=s_VPD_z,s_Sₑ=s_Sₑ) #mm/day

    return Ec
end

function Est_gₛ(LAI::T,daylight::T;
    cons::CCPH.Constants=CCPH.Constants(),
    env::CCPH.EnvironmentStruct=CCPH.EnvironmentStruct(),
    treepar::CCPH.TreePar=CCPH.TreePar(),
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
    Ec = Calctranpiration.(VPD_z,Sₑ;E_cm_ref=E_cm_ref,s_VPD_z=s_VPD_z,s_Sₑ=s_Sₑ) #mm/day
    Ec *= 10^-3*cons.ρ_H2O/cons.M_H2O/daylight #mol m⁻² s⁻¹
    λ = (1-exp(-treepar.k*LAI))/treepar.k
    gₛ = Ec*P/(VPD*cons.r*λ)
    return gₛ
end

#For weekly mean GPP data, index 1-21=>2015,index 24-44=>2016, index 45-64=>2017, index 65-82=>2018
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

#Container for Rosinedal data
mutable struct RoData{T<:Float64}
    H::Array{T,1}
    B::Array{T,1}
    GPP::Array{T,1}
    Wf::Array{T,1}
    N::Array{T,1}
    Hc::T
    gₛ::Array{T,1}
    Ec::Array{T,1}
end

#Create gₛ and Ec time series for Rosinedal data (ROData)  
function Create_gₛ_Ec_Data(weatherts::CCPH.WeatherTS,LAI::Array{T,1}) where {T<:Float64}
    gₛ_data = Float64[]    
    Ec_data = Float64[]
   

    for i = 1:length(weatherts.date)
        env = EnvironmentStruct(weatherts,i)
        
        daylight = weatherts.daylight[i]/7
        
        push!(gₛ_data,Est_gₛ(LAI[Find_data_ind(i)],daylight;env=env))        
        push!(Ec_data,Est_Ec(daylight;env=env))        
    end
    return gₛ_data,Ec_data
end

function Create_RoData(GPP_data::Array{T,1},weatherts::CCPH.WeatherTS;
    treepar::TreePar=TreePar(),stand_type::String="Fertilized") where {T<:Float64}
    #Simualted size data from Hyungwoo Lim, data corresponding to mean trees size during
    #growth period: [2015,2016,2017,2018]
    if stand_type=="Fertilized"
        H_data = [19.07,19.34,19.64,19.87] 
        B_data = [0.0452,0.0467,0.0479,0.0487]
        Wf_data = [8.33,8.47,8.49,8.48]
        N_data = [850.0000,846.6667,846.6667,846.6667]/10000
        Hc = 10.89
        LAI = Wf_data.*N_data/treepar.LMA
    elseif stand_type=="Control"
        H_data = [20.86,21.01,21.19,21.36] 
        B_data = [0.0348,0.0356,0.0362,0.0368]
        Wf_data = [5.025,5.110,5.141,5.165]
        N_data = [1010.0000,1010.0000,1006.6667,1006.6667]/10000
        Hc = 10.74
        LAI = Wf_data.*N_data/treepar.LMA
    else
        error("Create_RoData: stand_type has to be either Fertilized or Control")
    end

    gₛ_data,Ec_data = Create_gₛ_Ec_Data(weatherts,LAI)

    data = RoData(H_data,B_data,GPP_data,Wf_data,N_data,Hc,gₛ_data,Ec_data)

    return data
end

#Create container for Rosinedal data for Fertilized and Control stand
function Create_RoData_C_F(GPP_data_F::Array{T,1},GPP_data_C::Array{T,1},
    weatherts_F::CCPH.WeatherTS,weatherts_C::CCPH.WeatherTS;treepar::TreePar=TreePar()) where {T<:Float64} 
    
    data_F = Create_RoData(GPP_data_F,weatherts_F;treepar=treepar,stand_type="Fertilized")
    data_C = Create_RoData(GPP_data_C,weatherts_C;treepar=treepar,stand_type="Control")  

    return data_F,data_C
end

function calcR²(y_data::Array{Float64,1},f_model::Array{Float64,1})
    ymean = mean(y_data)
    SStot = sum((y_data.-ymean).^2)
    SSres = sum((y_data-f_model).^2)
    return 1-SSres/SStot
end

function calcMAPE(y_data::Array{Float64,1},f_model::Array{Float64,1})
    return 100*mean(abs.((y_data-f_model)./y_data))
end

function calcMAPEₘₐₓ(y_data::Array{Float64,1},f_model::Array{Float64,1})
    return 100*maximum(abs.((y_data-f_model)./y_data))
end

function calcMSE(y_data::Array{Float64,1},f_model::Array{Float64,1})
    return mean((y_data-f_model).^2)
end

calcRMSE(y_data::Array{Float64,1},f_model::Array{Float64,1}) = sqrt(calcMSE(y_data,f_model))

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

#Calcualte model canopy transpiration (Ec)
function CalcEc(gₛ::Float64,growthlength::Float64,model::CCPH.CCPHStruct)
    VPD = model.env.VPD
    P = model.env.P 
    LAI = model.treesize.N*model.treesize.Wf/model.treepar.LMA
    ratio = (1-exp(-model.treepar.k*LAI))/model.treepar.k
    daylight = growthlength/7 

    Ec = model.cons.r*gₛ*VPD*ratio/P #mol m⁻² s⁻¹
    Ec *= model.cons.M_H2O*daylight*10^3/model.cons.ρ_H2O #mm day⁻¹

    return Ec 
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

function CalcModelOutput!(i::Integer,model::CCPH.CCPHStruct,
    weatherts::CCPH.WeatherTS,photo_kinetic::CCPH.PhotoKineticRates)
    growthlength,step_length = CCPH.Init_weather_par!(i,model,weatherts,photo_kinetic)    

    modeloutput,gₛ_opt,Nₘ_f_opt = OptimCCPH(growthlength,model)
    GPP = CalcGPP(step_length,model,modeloutput)

    Ec = CalcEc(gₛ_opt,growthlength*step_length,model)
    
    return GPP, Ec, modeloutput, gₛ_opt, Nₘ_f_opt
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

mutable struct ModelResult{T<:Array{Float64,1}}
    GPP::T  
    Ec::T
    gₛ::T #stomatal conductance (mol C m⁻² leaf area s⁻¹)  
    Nₘ_f::T
    A::T #leaf C assimilation time series (mol C m⁻² leaf area s⁻¹) 
    E::T #leaf transpiration time series (mol H₂O m² leaf area s⁻¹)
    cᵢ::T #intercellular carbon dioxide concentration (Pa)
end

function CreateParaDict(;αf::T=460.0,β₁::T=1.27,β₂::T=-0.27,
    Nₛ::T=0.04,rₘ::T=24.0,a_Jmax::T=0.088,b_Jmax::T=0.0,i::T=1.0,Kₓₗ₀::T=0.01,
    α_max::T=0.36,X0::T=-4.0,τ::T=7.0,Smax::T=17.4,
    a_GPP::T=0.2,b_GPP::T=0.08,a_Ec::T=0.06,b_Ec::T=0.15,
    μ_Nₘ_f::T=0.018,b_Nₘ_f::T=0.0018) where {T<:Float64}
    ParaDict = Dict(:αf=>αf,
    :β₁=>β₁,
    :β₂=>β₂,
    :Nₛ=>Nₛ,
    :rₘ=>rₘ,
    :a_Jmax=>a_Jmax,
    :b_Jmax=>b_Jmax,
    :i=>i,
    :Kₓₗ₀=>Kₓₗ₀,
    :α_max=>α_max,
    :X0=>X0,
    :τ=>τ,
    :Smax=>Smax,
    :a_GPP=>a_GPP,
    :b_GPP=>b_GPP,
    :a_Ec=>a_Ec,
    :b_Ec=>b_Ec,
    :μ_Nₘ_f=>μ_Nₘ_f,
    :b_Nₘ_f=>b_Nₘ_f)
    return ParaDict
end

function CreateParaDict(ParaDictInit::Dict{Symbol,Float64})
    ParaDict = CreateParaDict()
    for (key,value) in ParaDictInit
        ParaDict[key] = value         
    end
    return ParaDict
end

function CreateCalibParaDict(CaliParaSym::Array{Symbol,1},CaliParaVal::Array{Float64,1})
    length(CaliParaSym)==length(CaliParaVal)||error("CreateCalibParaDict: CaliParaSym and CaliParaVal must be of same length")
    CalibParaDict = Dict(zip(CaliParaSym,CaliParaVal))
    return CalibParaDict
end

function InitParaDict!(ParaDict::Dict{Symbol,Float64},CalibParaDict::Dict{Symbol,Float64})    
    for (key,value) in CalibParaDict
        ParaDict[key] = value              
    end
    return nothing
end

function CreateParaDict(CaliParaSym::Array{Symbol,1},
    CaliParaVal::Array{Float64,1};ParaDictInit=nothing)
    if typeof(ParaDictInit) == Dict{Symbol,Float64}
        ParaDict = CreateParaDict(ParaDictInit)
    else
        ParaDict = CreateParaDict()
    end

    CalibParaDict = CreateCalibParaDict(CaliParaSym,CaliParaVal)
    InitParaDict!(ParaDict,CalibParaDict)
    return ParaDict
end

function Initi_model_struct(H::T,Hc::T,N::T,Wf::T,B::T,αf::T,β₁::T,β₂::T,
    Nₛ::T,rₘ::T,a_Jmax::T,b_Jmax::T,i::T,Kₓₗ₀::T,α_max::T) where {T<:Float64}
    cons = Constants()
    env = EnvironmentStruct() 
    kinetic = PhotoKineticRates()
    photo = PhotoPar(kinetic,env.Tₐ)    
    hydPar = HydraulicsPar(;Kₓₗ₀=Kₓₗ₀,i=i)
    treepar = TreePar(;αf=αf,β₁=β₁,β₂=β₂,rₘ_ref=rₘ,Nₛ=Nₛ,a_Jmax=a_Jmax,b_Jmax=b_Jmax,α_max=α_max)   
       
    Hs = H-Hc   
    As = Wf/treepar.αf
    Ww = treepar.ρw*As*(treepar.β₁*H+treepar.β₂*Hs)    
    treesize = TreeSize(Wf,Ww,H,Hs,As,B,N)

    model = CCPHStruct(cons,env,treepar,treesize,photo,hydPar)

    return model,kinetic
end

function Initi_model_struct(H::T,Hc::T,N::T,Wf::T,B::T,ParaDict::Dict{Symbol,Float64}) where {T<:Float64}

    αf = ParaDict[:αf]
    β₁ = ParaDict[:β₁]
    β₂ = ParaDict[:β₂]
    Nₛ = ParaDict[:Nₛ]
    rₘ = ParaDict[:rₘ]
    a_Jmax = ParaDict[:a_Jmax]
    b_Jmax = ParaDict[:b_Jmax]
    i = ParaDict[:i]
    Kₓₗ₀ = ParaDict[:Kₓₗ₀]
    α_max = ParaDict[:α_max]

    model,kinetic = Initi_model_struct(H,Hc,N,Wf,B,αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₓₗ₀,α_max)
    return model,kinetic
end

function Initi_model_struct(j::Integer,data::RoData,ParaDict::Dict{Symbol,Float64}) where {T<:Float64}
    data_ind = Find_data_ind(j)
    H = data.H[data_ind]
    N = data.N[data_ind]
    Wf = data.Wf[data_ind]
    B = data.B[data_ind]
    Hc = data.Hc
    model,kinetic = Initi_model_struct(H,Hc,N,Wf,B,ParaDict)   

    return model,kinetic
end

function Run_RO_CCPH(ParaDict::Dict{Symbol,Float64},
    weatherts::WeatherTS,data::RoData;sim_steps::Integer=84) where {T<:Float64}       

    #Fertilized stand
    GPP_model = zeros(sim_steps)
    EC_model = zeros(sim_steps)  
    Nₘ_f_model = zeros(sim_steps)
    gₛ_model = zeros(sim_steps)
    A_model = zeros(sim_steps)
    E_model = zeros(sim_steps)
    cᵢ_model = zeros(sim_steps)
    for j = 1:sim_steps   
        model,kinetic = Initi_model_struct(j,data,ParaDict::Dict{Symbol,Float64})
        
        GPP_model[j],EC_model[j], modeloutput, gₛ_model[j], Nₘ_f_model[j] = CalcModelOutput!(j,model,weatherts,kinetic)
        A_model[j] = modeloutput.A
        E_model[j] = modeloutput.E
        cᵢ_model[j] = modeloutput.cᵢ
    end 
    
    modelresult = ModelResult(GPP_model,EC_model,gₛ_model,Nₘ_f_model,A_model,E_model,cᵢ_model)
    return modelresult
end

function Get_Result_RO_CCPH(ParaDict::Dict{Symbol,Float64},RO_data::RO_raw_data;stand_type::String="Fertilized",
    sim_steps::Integer = 84,weather_file::String = "Weather_RO")
    X0,τ,Smax = ParaDict[:X0],ParaDict[:τ],ParaDict[:Smax]

    weatherts,GPP_data = create_weather_struct_RO(RO_data;stand_type=stand_type,X0=X0,τ=τ,Smax=Smax,
    weather_file=weather_file)

    data = Create_RoData(GPP_data,weatherts;stand_type=stand_type)

    modelresult = Run_RO_CCPH(ParaDict,weatherts,data;sim_steps=sim_steps)

    return modelresult,weatherts,data
end

function calc_logP_term(ydata::T,ymodel::T,a::T,b::T) where {T<:Float64}
    e = abs(ydata-ymodel)
    σ = a+b*ymodel
    σ>0.0||error("Negative variance")
    return e/σ+log(σ)
end

function calc_logP_norm_term(ydata::T,ymodel::T,a::T,b::T) where {T<:Float64}
    e = abs(ydata-ymodel)
    σ = a+b*ymodel
    σ>0.0||error("Negative variance")
    return (e/σ)^2/2+log(σ)
end

function Calc_logP_GPP_Ec(model::ModelResult,data::RoData,ParaDict::Dict{Symbol,Float64};ind::Array{Int64,1}=collect(1:84)) 
    
    a_GPP,b_GPP,a_Ec,b_Ec = ParaDict[:a_GPP],ParaDict[:b_GPP],ParaDict[:a_Ec],ParaDict[:b_Ec] 
    μ_Nₘ_f,b_Nₘ_f = ParaDict[:μ_Nₘ_f],ParaDict[:b_Nₘ_f]

    logP = 0.0
    for i = ind #1:length(data.GPP)
        logP += calc_logP_term(data.GPP[i],model.GPP[i],a_GPP,b_GPP)
        logP += calc_logP_term(data.Ec[i],model.Ec[i],a_Ec,b_Ec)                 
    end    
     
    return -logP
end

function Calc_logP_GPP_Ec_Nm_f(model::ModelResult,data::RoData,ParaDict::Dict{Symbol,Float64};ind::Array{Int64,1}=collect(1:84)) 
    
    a_GPP,b_GPP,a_Ec,b_Ec = ParaDict[:a_GPP],ParaDict[:b_GPP],ParaDict[:a_Ec],ParaDict[:b_Ec] 
    μ_Nₘ_f,b_Nₘ_f = ParaDict[:μ_Nₘ_f],ParaDict[:b_Nₘ_f]

    logP = 0.0
    for i = ind #1:length(data.GPP)
        logP += calc_logP_term(data.GPP[i],model.GPP[i],a_GPP,b_GPP)
        logP += calc_logP_term(data.Ec[i],model.Ec[i],a_Ec,b_Ec)    
        logP += abs(model.Nₘ_f[i]-μ_Nₘ_f)/b_Nₘ_f          
    end    
     
    return -logP
end

function Calc_logP_GPP_Nm_f(model::ModelResult,data::RoData,ParaDict::Dict{Symbol,Float64};ind::Array{Int64,1}=collect(1:84)) 
    
    a_GPP,b_GPP,a_Ec,b_Ec = ParaDict[:a_GPP],ParaDict[:b_GPP],ParaDict[:a_Ec],ParaDict[:b_Ec] 
    μ_Nₘ_f,b_Nₘ_f = ParaDict[:μ_Nₘ_f],ParaDict[:b_Nₘ_f]

    logP = 0.0
    for i = ind #1:length(data.GPP)
        logP += calc_logP_term(data.GPP[i],model.GPP[i],a_GPP,b_GPP)            
        logP += abs(model.Nₘ_f[i]-μ_Nₘ_f)/b_Nₘ_f          
    end    
     
    return -logP
end

function Calc_logP_GPP(model::ModelResult,data::RoData,ParaDict::Dict{Symbol,Float64};ind::Array{Int64,1}=collect(1:84)) 
    
    a_GPP,b_GPP,a_Ec,b_Ec = ParaDict[:a_GPP],ParaDict[:b_GPP],ParaDict[:a_Ec],ParaDict[:b_Ec] 
    μ_Nₘ_f,b_Nₘ_f = ParaDict[:μ_Nₘ_f],ParaDict[:b_Nₘ_f]

    logP = 0.0
    for i = ind #1:length(data.GPP)
        logP += calc_logP_term(data.GPP[i],model.GPP[i],a_GPP,b_GPP)                 
    end    
     
    return -logP
end

function Calc_logP_GPP_Ec_Nm_f_sep(model::ModelResult,data::RoData,ParaDict::Dict{Symbol,Float64}) 
    
    a_GPP,b_GPP,a_Ec,b_Ec = ParaDict[:a_GPP],ParaDict[:b_GPP],ParaDict[:a_Ec],ParaDict[:b_Ec] 
    μ_Nₘ_f,b_Nₘ_f = ParaDict[:μ_Nₘ_f],ParaDict[:b_Nₘ_f]

    logP_GPP=logP_Ec=logP_Nm_f= 0.0

    for i = 1:length(data.GPP)
        logP_GPP += calc_logP_term(data.GPP[i],model.GPP[i],a_GPP,b_GPP)
        logP_Ec += calc_logP_term(data.Ec[i],model.Ec[i],a_Ec,b_Ec)    
        logP_Nm_f += abs(model.Nₘ_f[i]-μ_Nₘ_f)/b_Nₘ_f          
    end    
     
    return [-logP_GPP,-logP_Ec,-logP_Nm_f]
end

function Calc_logP_GPP_Ec_Nm_f_norm(model::ModelResult,data::RoData,ParaDict::Dict{Symbol,Float64}) 
    
    a_GPP,b_GPP,a_Ec,b_Ec = ParaDict[:a_GPP],ParaDict[:b_GPP],ParaDict[:a_Ec],ParaDict[:b_Ec] 
    μ_Nₘ_f,b_Nₘ_f = ParaDict[:μ_Nₘ_f],ParaDict[:b_Nₘ_f]

    logP = 0.0
    for i = 1:length(data.GPP)
        logP += calc_logP_norm_term(data.GPP[i],model.GPP[i],a_GPP,b_GPP)
        logP += calc_logP_norm_term(data.Ec[i],model.Ec[i],a_Ec,b_Ec)    
        logP += abs(model.Nₘ_f[i]-μ_Nₘ_f)/b_Nₘ_f          
    end    
     
    return -logP
end

function Calc_logP_GPP_Ec_Nm_f_weight(model::ModelResult,data::RoData,ParaDict::Dict{Symbol,Float64};ind::Array{Int64,1}=collect(1:84),λ_gpp::Float64=1.0) 
    
    a_GPP,b_GPP,a_Ec,b_Ec = ParaDict[:a_GPP],ParaDict[:b_GPP],ParaDict[:a_Ec],ParaDict[:b_Ec] 
    μ_Nₘ_f,b_Nₘ_f = ParaDict[:μ_Nₘ_f],ParaDict[:b_Nₘ_f]

    logP = 0.0
    for i = ind#1:length(data.GPP)
        logP += λ_gpp*calc_logP_term(data.GPP[i],model.GPP[i],a_GPP,b_GPP)
        logP += calc_logP_term(data.Ec[i],model.Ec[i],a_Ec,b_Ec)    
        logP += abs(model.Nₘ_f[i]-μ_Nₘ_f)/b_Nₘ_f          
    end    
     
    return -logP
end

struct CalibPara
    sym::Symbol
    rang::Tuple{Float64,Float64}
    sep::Bool   
end
CalibPara(Sym::Symbol,rang::Tuple{Float64,Float64};sep::Bool=false) = CalibPara(Sym,rang,sep)

function CalibParaVec(kwargs...)
    calibparavec = Array{CalibPara,1}()
    for x in kwargs
        if typeof(x) == Tuple{Symbol,Float64,Float64,Bool}
            push!(calibparavec,CalibPara(x[1],(x[2],x[3]);sep=x[4]))
        elseif typeof(x) == Tuple{Symbol,Float64,Float64}
            push!(calibparavec,CalibPara(x[1],(x[2],x[3])))
        else
            error("CalibParaVec: worng input type; Either Tuple{Symbol,Float64,Float64,Bool} or
            Tuple{Symbol,Float64,Float64}")
        end
        
    end 
    return calibparavec   
end

function CreateOptVar(calibparavec::Array{CalibPara,1})
    ranges = Array{Tuple{Float64,Float64},1}()
    para2ind = Dict{Symbol,Any}()
    ind = 1
    for para in calibparavec
        if para.sep             
            push!(ranges,para.rang,para.rang)
            merge!(para2ind,Dict(para.sym=>(ind,ind+1)))
            ind += 2
        else
            push!(ranges,para.rang)
            merge!(para2ind,Dict(para.sym=>ind))
            ind += 1
        end
    end
    return ranges,para2ind
end

function CreateParaDict(para2ind::Dict{Symbol,Any},
    CaliParaVal::Array{Float64,1};ParaDictInit_F=nothing,ParaDictInit_C=nothing)
    if typeof(ParaDictInit_F) == Dict{Symbol,Float64}
        ParaDict_F = CreateParaDict(ParaDictInit_F)
    else
        ParaDict_F = CreateParaDict()
    end
    if typeof(ParaDictInit_C) == Dict{Symbol,Float64}
        ParaDict_C = CreateParaDict(ParaDictInit_C)
    else
        ParaDict_C = CreateParaDict()
    end
    
    for (key,ind) in para2ind
        if isa(ind,Integer)
            ParaDict_F[key] = CaliParaVal[ind]
            ParaDict_C[key] = CaliParaVal[ind] 
        elseif isa(ind,Tuple{Integer,Integer})
            ParaDict_F[key] = CaliParaVal[ind[1]]
            ParaDict_C[key] = CaliParaVal[ind[2]]
        else
            error("CreateParaDict: Wrong para2ind type")
        end
    end
    
    return ParaDict_F,ParaDict_C 
end