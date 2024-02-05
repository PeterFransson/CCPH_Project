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
    gₛ₁_opt,gₛ₂_opt,Nₘ_f_opt = CCPH.CCPHOpt(daylength,kinetic,envfun,model;gₛ₁_lim_hi=gₛ₁_lim_hi,gₛ₂_lim_hi=gₛ₂_lim_hi)
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
        #append!(Ec_model,[data_day.Ec for data_day in modeloutput.output_day])
        append!(Ec_model,[modeloutput.output_weekly.Ec for data_day in 1:7])
    end
    return Ec_model
end