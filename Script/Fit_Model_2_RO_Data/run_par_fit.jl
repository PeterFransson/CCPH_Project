mutable struct OptVal
    gₛ₁::Real
    gₛ₂::Real
    Nₘ_f::Real
end

mutable struct ModelPar
    Nₛ::Real
    α_max::Real
    a_Jmax::Real    
    Kₓₗ₀::Real
    τ::Real
    ΔS::Real
    a_GPP::Real
    b_GPP::Real
    μ_Nₘ_f::Real
    b_Nₘ_f::Real
    a_Ec::Real
    b_Ec::Real
end

mutable struct ModelOutput
    output_weekly::CCPH.CCPHOutput
    optval_weekly::OptVal
    output_day::Vector{CCPH.CCPHOutput}
    optval_day::Vector{OptVal}
end

function get_env_from_data(data::CCPH.WeatherDataStruct)
    day_nr = CCPH.Dates.dayofyear(data.date)
    daylength = CCPH.daylighthour(data.lat*pi/180,day_nr)*3600 #Seconds
    I₀ₜₒₜ = data.Radₜₒ*2.3*10^-6 #mol m⁻²    
    VP_Tmin = CCPH.SVPₜ(data.Tmin) #Pa
    VP = min(VP_Tmin,data.VP)

    ψₛ = CCPH.θₛ2ψₛ(data.θₛ) #MPa
    Cₐ = data.Cₐ #Pa
    P = data.P #Pa
    Tₐ = data.Tmean #°C
    VPD =  CCPH.VPDₜ(Tₐ,VP) #Pa
    I₀ = I₀ₜₒₜ/daylength #mol m⁻² s⁻¹

    ψₛ_fun(t) = ψₛ
    Cₐ_fun(t) = Cₐ
    P_fun(t) = P
    Tₐ_fun(t) = CCPH.Tₐ_fun(t,data.Tmin,data.Tmax,daylength)
    VPD_fun(t) = CCPH.VPDₜ(Tₐ_fun(t),VP)
    I₀_fun(t) = CCPH.I₀_fun(t,I₀ₜₒₜ,daylength)

    envfun = CCPH.EnvironmentFunStruct(I₀_fun,Cₐ_fun,P_fun,Tₐ_fun,VPD_fun,ψₛ_fun)
    env = CCPH.EnvironmentStruct(I₀,Cₐ,P,Tₐ,VPD,ψₛ)

    return env,envfun,daylength
end
function get_env_from_data(data_vec::Vector{CCPH.WeatherDataStruct})
    day_nrs = [CCPH.Dates.dayofyear(data.date) for data in data_vec]
    lats = [data.lat for data in data_vec]
    daylengths = CCPH.daylighthour.(lats*pi/180,day_nrs)*3600 #Seconds
    daylength = mean(daylengths)

    I₀ₜₒₜ = mean([data.Radₜₒ for data in data_vec])*2.3*10^-6 #mol m⁻²  
    Tmin_data = mean([data.Tmin for data in data_vec]) #°C
    Tmax_data = mean([data.Tmax for data in data_vec]) #°C

    VP_Tmin = CCPH.SVPₜ(Tmin_data) #Pa  
    VP_data = mean([data.VP for data in data_vec])
    VP = min(VP_Tmin,VP_data)

    ψₛ = mean([CCPH.θₛ2ψₛ(data.θₛ) for data in data_vec]) #MPa
    Cₐ = mean([data.Cₐ for data in data_vec]) #Pa
    P =  mean([data.P for data in data_vec]) #Pa
    Tₐ =  mean([data.Tmean for data in data_vec]) #°C

    VPD =  CCPH.VPDₜ(Tₐ,VP) #Pa
    I₀ = I₀ₜₒₜ/daylength #mol m⁻² s⁻¹

    ψₛ_fun(t) = ψₛ
    Cₐ_fun(t) = Cₐ
    P_fun(t) = P
    Tₐ_fun(t) = CCPH.Tₐ_fun(t,Tmin_data,Tmax_data,daylength)
    VPD_fun(t) = CCPH.VPDₜ(Tₐ_fun(t),VP)
    I₀_fun(t) = CCPH.I₀_fun(t,I₀ₜₒₜ,daylength)

    envfun = CCPH.EnvironmentFunStruct(I₀_fun,Cₐ_fun,P_fun,Tₐ_fun,VPD_fun,ψₛ_fun)
    env = CCPH.EnvironmentStruct(I₀,Cₐ,P,Tₐ,VPD,ψₛ)

    return env,envfun,daylength
end

function get_model_output(daylength::Real,model::CCPH.CCPHStruct,kinetic::CCPH.PhotoKineticRates,envfun::CCPH.EnvironmentFunStruct)
    gₛ₁_lim_hi,gₛ₂_lim_hi = CCPH.SDM2_get_gₛ_lim!(daylength,model,kinetic,envfun)     
    gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt = CCPH.CCPHOpt(daylength,
    kinetic,
    envfun,
    model;
    gₛ₁_lim_hi=gₛ₁_lim_hi,
    gₛ₂_lim_hi=gₛ₂_lim_hi,
    gₛ₁_guess=gₛ₁_lim_hi*0.9,
    gₛ₂_guess=gₛ₂_lim_hi*0.9)
    optval = OptVal(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt)
    output = CCPH.CCPH_run!(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt,daylength,kinetic,envfun,model)
    return optval,output
end
function get_model_output(Nₘ_f_opt::Real,
    gₛ₁_guess::Real,
    gₛ₂_guess::Real,
    daylength::Real,
    model::CCPH.CCPHStruct,
    kinetic::CCPH.PhotoKineticRates,
    envfun::CCPH.EnvironmentFunStruct)

    gₛ₁_lim_hi,gₛ₂_lim_hi = CCPH.SDM2_get_gₛ_lim!(daylength,model,kinetic,envfun)     
    gₛ₁_opt,gₛ₂_opt = CCPH.CCPHOpt(Nₘ_f_opt,
    daylength,
    kinetic,
    envfun,
    model;
    gₛ₁_lim_hi=gₛ₁_lim_hi,
    gₛ₂_lim_hi=gₛ₂_lim_hi,
    gₛ₁_guess=min(gₛ₁_guess,gₛ₁_lim_hi*0.9),
    gₛ₂_guess=min(gₛ₂_guess,gₛ₂_lim_hi*0.9))    
    optval = OptVal(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt)
    output = CCPH.CCPH_run!(gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt,daylength,kinetic,envfun,model)
    return optval,output
end

function Xₜ_fun(raw_data::RawInputData,par::ModelPar;Smin::Real = -4.0)
    return Xₜ_fun(raw_data;Smin=Smin,ΔS=par.ΔS,τ=par.τ)
end

function TreePar!(treepar::CCPH.TreePar,Xₜ::Real)
    treepar.Xₜ=Xₜ
    return nothing    
end
CCPH.TreePar(par::ModelPar) = CCPH.TreePar(Nₛ=par.Nₛ,α_max=par.α_max,a_Jmax=par.a_Jmax,b_Jmax=0.0)
CCPH.TreePar(par::ModelPar,Xₜ::Real) = CCPH.TreePar(Nₛ=par.Nₛ,α_max=par.α_max,a_Jmax=par.a_Jmax,b_Jmax=0.0,Xₜ=Xₜ)
CCPH.HydraulicsPar(par::ModelPar) = CCPH.HydraulicsPar(Kₓₗ₀=par.Kₓₗ₀)

function intitiate_model(data_week::Vector{CCPH.WeatherDataStruct},
    Xₜ_week::Vector{Real},
    treesize::CCPH.TreeSize,
    par::ModelPar)
    
    #Create structs
    cons = CCPH.Constants()   
    treepar = CCPH.TreePar(par,mean(Xₜ_week))    
    hydPar = CCPH.HydraulicsPar(par)
    kinetic = CCPH.PhotoKineticRates()    

    env,envfun,daylength = get_env_from_data(data_week)    
    photo = CCPH.PhotoPar(kinetic,env.Tₐ)
    model = CCPH.CCPHStruct(cons,env,treepar,treesize,photo,hydPar)

    return model,envfun,daylength,kinetic,cons,treepar,treesize,hydPar
end

function intitiate_model(data_day::CCPH.WeatherDataStruct,
    Xₜ::Real,
    kinetic::CCPH.PhotoKineticRates,
    cons::CCPH.Constants,
    treepar::CCPH.TreePar,
    treesize::TreeSize,
    hydPar::HydraulicsPar)
    
    env_day,envfun_day,daylength_day = get_env_from_data(data_day)  
    TreePar!(treepar,Xₜ)  
    photo_day = PhotoPar(kinetic,env_day.Tₐ)
    model_day = CCPHStruct(cons,env_day,treepar,treesize,photo_day,hydPar)

    return model_day,envfun_day,daylength_day
end

function run_week(data_week::Vector{CCPH.WeatherDataStruct},
    Xₜ_week::Vector{Real},
    treesize::CCPH.TreeSize,
    par::ModelPar)
    
    model,envfun,daylength,kinetic,cons,treepar,treesize,hydPar = intitiate_model(data_week,Xₜ_week,treesize,par)
    
    optval_weekly,output_weekly = get_model_output(daylength,model,kinetic,envfun)
    
    output_day_vec = Vector{CCPH.CCPHOutput}(undef,7)
    optval_day_vec =Vector{OptVal}(undef,7)

    for i = 1:7

        model_day,envfun_day,daylength_day = intitiate_model(data_week[i],Xₜ_week[i],kinetic,cons,treepar,treesize,hydPar)
                
        optval_day,output_day = get_model_output(optval_weekly.Nₘ_f,
        optval_weekly.gₛ₁,
        optval_weekly.gₛ₂,
        daylength_day,
        model_day,
        kinetic,
        envfun_day)  
        
        output_day_vec[i] = output_day
        optval_day_vec[i] = optval_day
    end
    
    return ModelOutput(output_weekly,optval_weekly,output_day_vec,optval_day_vec)
end

function run_week(input_data::RawInputData,
    Xₜ::Vector{Real},
    par::ModelPar)

    treesize = input_data.treesize    

    modeloutput = [run_week(input_data.weather_growth[idx[1]:idx[2]],Xₜ[idx[1]:idx[2]],treesize,par) for idx in input_data.growth_indices_weekly]

    return modeloutput
end

function get_GPP_model(modeloutputs::Vector{ModelOutput})
    GPP_model = Real[]
    for modeloutput in modeloutputs
        append!(GPP_model,[data_day.GPP for data_day in modeloutput.output_day])
    end
    return GPP_model
end
function get_Ec_model(modeloutputs::Vector{ModelOutput})
    Ec_model = Real[]
    for modeloutput in modeloutputs
        append!(Ec_model,[data_day.Ec for data_day in modeloutput.output_day])        
    end
    return Ec_model
end
function get_Nₘ_f_model(modeloutputs::Vector{ModelOutput})
    Nₘ_f_model = Real[]
    for modeloutput in modeloutputs
        append!(Nₘ_f_model,[data_day.Nₘ_f for data_day in modeloutput.optval_day])        
    end
    return Nₘ_f_model
end

function ModelPar(x::Vector{T};stand_type::Symbol=:Fertilized) where {T<:Real}
    Nₛ,α_max,a_Jmax,Kₓₗ₀,τ,ΔS,a_GPP,b_GPP = x   

    if stand_type==:Fertilized
        μ_Nₘ_f=0.02135
        b_Nₘ_f=0.002287
        a_Ec=0.0688
        b_Ec=0.146        
    elseif stand_type==:Control
        μ_Nₘ_f=0.01195
        b_Nₘ_f=0.002287
        a_Ec=0.0433
        b_Ec=0.179  
    else
        error("Wrong input. Either :Fertilized or :Control")
    end

    return ModelPar(Nₛ,α_max,a_Jmax,Kₓₗ₀,τ,ΔS,a_GPP,b_GPP,μ_Nₘ_f,b_Nₘ_f,a_Ec,b_Ec)
end

function calc_logP_norm_term(ydata::Real,ymodel::Real,a::Real,b::Real)
    e = abs(ydata-ymodel)
    σ = a+b*ymodel
    σ>0.0||error("Negative variance")
    return (e/σ)^2/2+log(σ)
end

function Calc_logP_GPP_Ec_Nm_f_norm(GPP_model::Vector{Vector{T}},
    Ec_model::Vector{Vector{T}},
    Nₘ_f_model::Vector{Vector{T}},
    GPP_data::Vector{Vector{R}},
    Ec_data::Vector{Vector{S}},
    par::ModelPar,
    raw_input::Vector{RawInputData}) where {T<:Real,R<:Real,S<:Real}
    
    a_GPP,b_GPP,a_Ec,b_Ec = par.a_GPP,par.b_GPP,par.a_Ec,par.b_Ec
    μ_Nₘ_f,b_Nₘ_f = par.μ_Nₘ_f,par.b_Nₘ_f

    logP = 0.0

    for i in 1:4        
        logP += sum(calc_logP_norm_term.(GPP_data[i],GPP_model[i]*raw_input[i].ζ,Ref(a_GPP),Ref(b_GPP)))
        logP += sum(calc_logP_norm_term.(Ec_data[i],Ec_model[i],Ref(a_Ec),Ref(b_Ec)))    
        logP += sum(abs.(Nₘ_f_model[i].-μ_Nₘ_f)/b_Nₘ_f)        
    end   
     
    return -logP
end

function run_model(par::ModelPar,raw_input::Vector{RawInputData};stand_type::Symbol=:Fertilized) where {T<:Real}
    Xₜ = Xₜ_fun.(raw_input,Ref(par))
    modeloutput = run_week.(raw_input,Xₜ,Ref(par))    
    GPP_model = get_GPP_model.(modeloutput)
    Ec_model = get_Ec_model.(modeloutput)
    Nₘ_f_model = get_Nₘ_f_model.(modeloutput)
    return (GPP_model,Ec_model,Nₘ_f_model)
end

function opt_par_obj(x::Vector{T},
    raw_input::Vector{RawInputData},
    GPP_data::Vector{Vector{R}},
    Ec_data::Vector{Vector{S}};
    stand_type::Symbol=:Fertilized) where{T<:Real,R<:Real,S<:Real}

    try
        par = ModelPar(x;stand_type=stand_type)
        GPP_model,Ec_model,Nₘ_f_model = run_model(par,raw_input;stand_type=stand_type)
        return -Calc_logP_GPP_Ec_Nm_f_norm(GPP_model,Ec_model,Nₘ_f_model,GPP_data,Ec_data,par,raw_input)
    catch err
        println("Parameters Error: ", err)

        return Inf
    end 
end