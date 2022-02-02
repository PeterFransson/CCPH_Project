function Initi_model_struct(H::T,Hc::T,N::T,Wf::T,B::T,αf::T,β₁::T,β₂::T,
    Nₛ::T,rₘ::T,a_Jmax::T,b_Jmax::T,i::T,Kₓₗ₀::T) where {T<:Float64}
    cons = Constants()
    env = EnvironmentStruct() 
    kinetic = PhotoKineticRates()
    photo = PhotoPar(kinetic,env.Tₐ)    
    hydPar = HydraulicsPar(;Kₓₗ₀=Kₓₗ₀,i=i)
    treepar = TreePar(;αf=αf,β₁=β₁,β₂=β₂,rₘ=rₘ,Nₛ=Nₛ,a_Jmax=a_Jmax,b_Jmax=b_Jmax)   
       
    Hs = H-Hc   
    As = Wf/treepar.αf
    Ww = treepar.ρw*As*(treepar.β₁*H+treepar.β₂*Hs)    
    treesize = TreeSize(Wf,Ww,H,Hs,As,B,N)

    model = CCPHStruct(cons,env,treepar,treesize,photo,hydPar)

    return model,kinetic
end

function Initi_model_struct(j::Integer,data::RoData,
    αf::T,β₁::T,β₂::T,Nₛ::T,rₘ::T,a_Jmax::T,b_Jmax::T,i::T,Kₓₗ₀::T) where {T<:Float64}
    data_ind = Find_data_ind(j)
    H = data.H[data_ind]
    N = data.N[data_ind]
    Wf = data.Wf[data_ind]
    B = data.B[data_ind]
    Hc = data.Hc
    model,kinetic = Initi_model_struct(H,Hc,N,Wf,B,αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₓₗ₀)    

    return model,kinetic
end

function Run_RO_C_F_CCPH(par::Array{Float64,1},weatherts_F::WeatherTS,data_F::RoData;
    sim_steps::Integer=83)
    
    αf,β₁,β₂,b_Jmax = 460.0,1.27,-0.27,0.0
    Nₛ,rₘ,a_Jmax,Kₓₗ₀,i = par[1:5]    

    #Fertilized stand
    GPP_model_F = zeros(sim_steps+1)
    gₛ_model_F = zeros(sim_steps+1)  
    for j = 1:sim_steps+1           
       
        model,kinetic = Initi_model_struct(j,data_F,αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₓₗ₀)
        
        GPP_model_F[j],gₛ_model_F[j] = CalcGPP!(j,model,weatherts_F,kinetic)
    end  

    return GPP_model_F,gₛ_model_F
end

function RO_C_F_CCPH_ModelOutput(par::Array{Float64,1},weatherts_F::WeatherTS,
    data_F::RoData;sim_steps::Integer=83)

    αf,β₁,β₂,b_Jmax = 460.0,1.27,-0.27,0.0
    Nₛ,rₘ,a_Jmax,Kₓₗ₀,i = par[1:6] 

    #Fertilized stand
    ModelOutput_F = CCPHOutput[]
    k_cost_max_F = Float64[]
    gₛ_F = Float64[]
    gₛ_crit_F = Float64[]
    Nₘ_f_F = Float64[]
    for j = 1:sim_steps+1
       
        model,kinetic = Initi_model_struct(j,data_F,αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₓₗ₀)
        modeloutput,gₛ_opt,Nₘ_f_opt = OptimCCPH!(j,model,weatherts_F,kinetic)

        push!(ModelOutput_F,modeloutput)
        push!(k_cost_max_F,Max_k_cost(model))
        push!(gₛ_F,gₛ_opt)
        push!(gₛ_crit_F,CCPH.Calc_K_costᵢₙᵥ(0.12,model))
        push!(Nₘ_f_F,Nₘ_f_opt)
    end     

    return ModelOutput_F,k_cost_max_F,gₛ_F,gₛ_crit_F,Nₘ_f_F
end

function Calc_error_RO_C_F_CCPH(
    GPP_model_F::Array{Float64,1},gₛ_model_F::Array{Float64,1},data_F::RoData)

    n_GPP_model = length(GPP_model_F)
    n_gₛ_model = length(gₛ_model_F)
   
    error_F = 0.0
    error_F += sqrt(mean((data_F.GPP[1:n_GPP_model].-GPP_model_F).^2))
    error_F += 100*sqrt(mean((data_F.gₛ[1:n_gₛ_model].-gₛ_model_F).^2))        
     
    return  error_F
end

function Check_para_RO_C_F_CCPH(par::Array{Float64,1})
    #Nₛ,rₘ,a_Jmax,Kₓₗ₀,i
    lo = [0.001,10.0,0.01,0.001,0.1]

    up = [0.1,40.0,1.0,0.1,6.0] 

    return any(par.<lo)||any(par.>up)    
end

function Post_distri_RO_C_F_CCPH(par::Array{Float64,1},RO_data::RO_raw_data)
    if Check_para_RO_C_F_CCPH(par)
        return post = -Inf
    else
        try            
            weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized")
            weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control")           
            
            data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C,weatherts_F,weatherts_C)   

            GPP_model_F,gₛ_model_F = Run_RO_C_F_CCPH(par,weatherts_F,data_F;sim_steps = 20)              

            error_val = Calc_error_RO_C_F_CCPH(GPP_model_F,gₛ_model_F,data_F)
            
            return post = -error_val
        catch err
            println("Parameters Error: ", err)
            return post = -Inf
        end
    end
end

function run_simple_trait_model_C_F_tuning_RO_ts()
    #Init Ro data
    RO_data = Load_RO_data()    

    nsamples = 3000
    nchains = 6
    
    par_guess = [0.011976551772617215, 25.42466609569836, 0.5009242298864488,
    0.010452800544803813,1.0]

    x_current = [par_guess for i = 1:nchains]
    P_current = [Post_distri_RO_C_F_CCPH(par_guess,RO_data) for i = 1:nchains] 
        
    println(P_current)    
    
    q_vec = [0.04,24.0,0.033,0.04,1.0]/25.0
    metropolis_mcmc(x_current,P_current,
    x::Array{Float64,1}->Post_distri_RO_C_F_CCPH(x::Array{Float64,1},RO_data),q_vec,"RO__20220201";
    n_chains=nchains,n_samples=nsamples,burn_in=5000,sample_freq=5)        
end

function run_C_F_ts_mean()

    RO_data = Load_RO_data()  

    file_name = "RO__20220201"
        
    samples = load("./output/"*file_name*".jld","samples_container")    
    P_samples = load("./output/"*file_name*".jld","P_samples_container")
    Draw_trace_plot(samples,file_name)
    Draw_histogram(samples,file_name)

    par = Find_max_sample(samples,P_samples)  
  
    println(par)
    
    weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized")
    weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control")
    
    data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C,weatherts_F,weatherts_C)    

    GPP_model_F,gₛ_model_F = Run_RO_C_F_CCPH(par,weatherts_F,data_F;sim_steps = 20) 

    n_GPP_model = length(GPP_model_F)
    n_gₛ_model = length(gₛ_model_F)

    @show GPP_R²_F =  calcR²(data_F.GPP[1:n_GPP_model],GPP_model_F)
    @show GPP_cor_F = cor(data_F.GPP[1:n_GPP_model],GPP_model_F)
    @show gₛ_R²_F =  calcR²(data_F.gₛ[1:n_gₛ_model],gₛ_model_F)
    @show gₛ_cor_F = cor(data_F.gₛ[1:n_gₛ_model],gₛ_model_F)      

    plot(weatherts_F.date[1:n_GPP_model],data_F.GPP[1:n_GPP_model],label="Data")
    pl1 = plot!(weatherts_F.date[1:n_GPP_model],GPP_model_F,label="Model")
    
    plot(weatherts_F.date[1:n_GPP_model],data_F.gₛ[1:n_gₛ_model],label="Data")
    pl3 = plot!(weatherts_F.date[1:n_GPP_model],gₛ_model_F,label="Model")
    
    plot(pl1,pl3,layout=(1,2),legends=false)
end

#run_simple_trait_model_C_F_tuning_RO_ts()
run_C_F_ts_mean()