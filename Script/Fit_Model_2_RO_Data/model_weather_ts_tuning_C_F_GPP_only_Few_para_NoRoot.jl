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

function Run_RO_C_F_CCPH(par::Array{Float64,1},weatherts_F::WeatherTS,
    weatherts_C::WeatherTS,data_F::RoData,data_C::RoData;sim_steps::Integer=83)
    
    αf,β₁,β₂,b_Jmax,i = 460.0,1.27,-0.27,0.0,1.0
    Nₛ,rₘ,a_Jmax,Kₓₗ₀,Nₛ_C = par[1:5]    

    #Fertilized stand
    GPP_model_F = zeros(sim_steps+1)
    gₛ_model_F = zeros(sim_steps+1)  
    for j = 1:sim_steps+1           
       
        model,kinetic = Initi_model_struct(j,data_F,αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₓₗ₀)
        
        GPP_model_F[j],gₛ_model_F[j] = CalcGPP!(j,model,weatherts_F,kinetic)
    end  

    #Control stand
    GPP_model_C = zeros(sim_steps+1)
    gₛ_model_C = zeros(sim_steps+1)
    for j = 1:sim_steps+1        
        
        model,kinetic = Initi_model_struct(j,data_C,αf,β₁,β₂,Nₛ_C,rₘ,a_Jmax,b_Jmax,i,Kₓₗ₀)
        
        GPP_model_C[j],gₛ_model_C[j] = CalcGPP!(j,model,weatherts_C,kinetic)
    end  

    return GPP_model_F,GPP_model_C,gₛ_model_F,gₛ_model_C
end

function RO_C_F_CCPH_ModelOutput(par::Array{Float64,1},weatherts_F::WeatherTS,
    weatherts_C::WeatherTS,data_F::RoData,data_C::RoData;sim_steps::Integer=83)

    αf,β₁,β₂,b_Jmax,i = 460.0,1.27,-0.27,0.0,1.0
    Nₛ,rₘ,a_Jmax,Kₓₗ₀,Nₛ_C = par[1:5] 

    #Fertilized stand
    ModelOutput_F = CCPHOutput[]
    k_cost_max_F = Float64[]
    gₛ_F = Float64[]
    gₛ_crit_F = Float64[]
    Nₘ_f_F = Float64[]
    for j = 1:sim_steps+1
       
        model,kinetic = Initi_model_struct(j,data_C,αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₓₗ₀)
        modeloutput,gₛ_opt,Nₘ_f_opt = OptimCCPH!(j,model,weatherts_F,kinetic)

        push!(ModelOutput_F,modeloutput)
        push!(k_cost_max_F,Max_k_cost(model))
        push!(gₛ_F,gₛ_opt)
        push!(gₛ_crit_F,CCPH.Calc_K_costᵢₙᵥ(0.12,model))
        push!(Nₘ_f_F,Nₘ_f_opt)
    end  

    #Control stand
    ModelOutput_C = CCPHOutput[]
    k_cost_max_C = Float64[]
    gₛ_C = Float64[]
    gₛ_crit_C = Float64[]
    Nₘ_f_C = Float64[]
    for j = 1:sim_steps+1
        
        model,kinetic = Initi_model_struct(j,data_C,αf,β₁,β₂,Nₛ_C,rₘ,a_Jmax,b_Jmax,i,Kₓₗ₀)
        modeloutput,gₛ_opt,Nₘ_f_opt = OptimCCPH!(j,model,weatherts_C,kinetic)

        push!(ModelOutput_C,modeloutput)
        push!(k_cost_max_C,Max_k_cost(model))
        push!(gₛ_C,gₛ_opt)
        push!(gₛ_crit_C,CCPH.Calc_K_costᵢₙᵥ(0.12,model))
        push!(Nₘ_f_C,Nₘ_f_opt)
    end  

    return ModelOutput_F,ModelOutput_C,k_cost_max_F,k_cost_max_C,gₛ_F,gₛ_C,gₛ_crit_F,gₛ_crit_C,Nₘ_f_F,Nₘ_f_C
end

function Calc_error_RO_C_F_CCPH(
    GPP_model_F::Array{Float64,1},gₛ_model_F::Array{Float64,1},data_F::RoData,
    GPP_model_C::Array{Float64,1},gₛ_model_C::Array{Float64,1},data_C::RoData)   

    #[21,44,64,84]
    ind_start = [1,22,45,65]
    ind_end = [21,44,64,84]
    #Fertilized stand

    error_F = 0.0
    for i = 1:4
        error_F += 1-calcR²(data_F.GPP[ind_start[i]:ind_end[i]],GPP_model_F[ind_start[i]:ind_end[i]])
        error_F += 1-calcR²(data_F.gₛ[ind_start[i]:ind_end[i]],gₛ_model_F[ind_start[i]:ind_end[i]])
    end   

    #Control stand
    error_C = 0.0
    for i = 1:4
        error_C += 1-calcR²(data_C.GPP[ind_start[i]:ind_end[i]],GPP_model_C[ind_start[i]:ind_end[i]])
        error_C += 1-calcR²(data_C.gₛ[ind_start[i]:ind_end[i]],gₛ_model_C[ind_start[i]:ind_end[i]])
    end
     
    return  error_F+error_C
end

function Check_para_RO_C_F_CCPH(par::Array{Float64,1})
    #Nₛ,rₘ,a_Jmax,Kₓₗ₀,Nₛ_C,X0,τ,Smax,τ_C
    lo = [0.01,10.0,0.01,0.001,0.01,-8.0,3.0,10.0,3.0]

    up = [0.1,40.0,1.0,0.1,0.1,-1.0,15.0,20.0,15.0] 

    return any(par.<lo)||any(par.>up)    
end

function Post_distri_RO_C_F_CCPH(par::Array{Float64,1},RO_data::RO_raw_data)
    if Check_para_RO_C_F_CCPH(par)
        return post = 0.0
    else
        try
            X0,τ,Smax,τ_C = par[6:9]
            weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
            weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ_C,Smax=Smax)           
            
            data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C,weatherts_F,weatherts_C)   

            GPP_model_F,GPP_model_C,gₛ_model_F,gₛ_model_C = Run_RO_C_F_CCPH(par,weatherts_F,weatherts_C,data_F,data_C)              

            error_val = Calc_error_RO_C_F_CCPH(
            GPP_model_F,gₛ_model_F,data_F,
            GPP_model_C,gₛ_model_C,data_C)
            
            return post = -error_val
        catch err
            println("Parameters Error: ", err)
            return post = 0.0
        end
    end
end

function run_simple_trait_model_C_F_tuning_RO_ts()
    #Init Ro data
    RO_data = Load_RO_data()    

    nsamples = 10
    nchains = 12

    par_guess = [0.036975489745880365, 17.40739827323768, 0.5, 0.01,
    0.025675473753997178, -4.230533359595931, 12.328105896346644, 19.079489364590533,12.328105896346644]   

    x_current = [par_guess for i = 1:nchains]
    P_current = [Post_distri_RO_C_F_CCPH(par_guess,RO_data) for i = 1:nchains] 
        
    println(P_current)    
    
    q_vec = [0.04,24.0,0.033,0.01,0.04,-4.0,7.0,16.0,7.0]/25.0
    metropolis_mcmc(x_current,P_current,
    x::Array{Float64,1}->Post_distri_RO_C_F_CCPH(x::Array{Float64,1},RO_data),q_vec,"RO_C_F_GPP_Only_Few_Para_NoRoot_20220105";
    n_chains=nchains,n_samples=nsamples,burn_in=10,sample_freq=5) 
        
end

function run_C_F_ts_mean()

    RO_data = Load_RO_data()  
    
    #=
    samples = load("output/RO_C_F_GPP_Only_20211202.jld","samples_container")    
    P_samples = load("output/RO_C_F_GPP_Only_20211202.jld","P_samples_container")

    par = Find_max_sample(samples,P_samples)
    =#   
    par = [0.036975489745880365, 17.40739827323768, 0.1, 0.01,
    0.025675473753997178, -4.230533359595931, 12.328105896346644,
    19.079489364590533,12.328105896346644]  

    println(par)

    X0,τ,Smax = par[6:9]
    weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
    weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ,Smax=Smax)
    
    data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C,weatherts_F,weatherts_C)    

    GPP_model_F,GPP_model_C,gₛ_model_F,gₛ_model_C = Run_RO_C_F_CCPH(par,weatherts_F,weatherts_C,data_F,data_C) 

    @show GPP_R²_F =  calcR²(data_F.GPP,GPP_model_F)
    @show GPP_cor_F = cor(data_F.GPP,GPP_model_F)
    @show gₛ_R²_F =  calcR²(data_F.gₛ,gₛ_model_F)
    @show gₛ_cor_F = cor(data_F.gₛ,gₛ_model_F)

    @show GPP_R²_C =  calcR²(data_C.GPP,GPP_model_C)
    @show GPP_cor_C = cor(data_C.GPP,GPP_model_C)
    @show gₛ_R²_C =  calcR²(data_C.gₛ,gₛ_model_C)
    @show gₛ_cor_C = cor(data_C.gₛ,gₛ_model_C)

    println("Fertilized")
    clac_GPP_R²_annual(data_C.GPP,GPP_model_F,weatherts_F.date)
    println("Control")  
    clac_GPP_R²_annual(data_C.GPP,GPP_model_C,weatherts_C.date)  

    plot(weatherts_F.date,data_F.GPP,label="Data")
    pl1 = plot!(weatherts_F.date,GPP_model_F,label="Model")
    plot(weatherts_C.date,data_C.GPP,label="Data")
    pl2 = plot!(weatherts_C.date,GPP_model_C,label="Model")

    plot(weatherts_F.date,data_F.gₛ,label="Data")
    pl3 = plot!(weatherts_F.date,gₛ_model_F,label="Model")
    plot(weatherts_C.date,data_C.gₛ,label="Data")
    pl4 = plot!(weatherts_C.date,gₛ_model_C,label="Model")

    plot(pl1,pl3,pl2,pl4,layout=(2,2),legends=false)
end

#=
function validate_RO_2019()
    RO_data = Load_RO_data(;weather_file = "Weather_RO_2")  
    samples = load("output/RO_C_F_GPP_Only_20211119.jld","samples_container")    
    P_samples = load("output/RO_C_F_GPP_Only_20211119.jld","P_samples_container")

    par = Find_max_sample(samples,P_samples)

    println(par)

    X0,τ,Smax = par[11:13]
    weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax,
    weather_file = "Weather_RO_2")
    weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ,Smax=Smax,
    weather_file = "Weather_RO_2")

    data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C)  
    GPP_model_F,GPP_model_C = Run_RO_C_F_CCPH(par,weatherts_F,weatherts_C,data_F,data_C;sim_steps = 104) 

    ModelOutput_F,ModelOutput_C,k_cost_max_F,k_cost_max_C,gₛ_F,gₛ_C,gₛ_crit_F,gₛ_crit_C,Nₘ_f_F,Nₘ_f_C = 
    RO_C_F_CCPH_ModelOutput(par,weatherts_F,weatherts_C,data_F,data_C;sim_steps=104)
    
    K_cost_F = [Output.K_cost for Output in ModelOutput_F]
    K_cost_C = [Output.K_cost for Output in ModelOutput_C]

    αr_F = [Output.αr for Output in ModelOutput_F]
    αr_C = [Output.αr for Output in ModelOutput_C]
    
    @show length(GPP_model_F)

    println("Fertilized:")
    clac_GPP_R²_annual(data_F.GPP,GPP_model_F,weatherts_F.date)
    println("Control:")
    clac_GPP_R²_annual(data_C.GPP,GPP_model_C,weatherts_C.date)

    @show calcR²(data_F.GPP,GPP_model_F)
    @show  cor(data_F.GPP,GPP_model_F)
    @show calcR²(data_C.GPP,GPP_model_C)
    @show  cor(data_C.GPP,GPP_model_C)

    plot(weatherts_F.date,data_F.GPP,label="Data")
    pl1 = plot!(weatherts_F.date,GPP_model_F,label="Model")
    plot([2.0,9.0],[2.0,9.0],xlabel="Data",ylabel="Model")
    pl2 = plot!(data_F.GPP,GPP_model_F,seriestype=:scatter)
    plot(weatherts_C.date,data_C.GPP,label="Data")
    pl3 = plot!(weatherts_C.date,GPP_model_C,label="Model")
    plot([2.0,9.0],[2.0,9.0],xlabel="Data",ylabel="Model")
    pl4 = plot!(data_C.GPP,GPP_model_C,seriestype=:scatter)
    plot(pl1,pl2,pl3,pl4,layout=(2,2),legends=false)

    plot(weatherts_F.date,K_cost_F)
    pl1 = plot!(weatherts_F.date,k_cost_max_F) 
    plot(weatherts_C.date,K_cost_C) 
    pl2 = plot!(weatherts_C.date,k_cost_max_C)   
    plot(pl1,pl2,layout=(2,1),legends=false)  

    pl1 = plot(weatherts_F.date,gₛ_F)
    #pl1 = plot!(weatherts_F.date,gₛ_crit_F) 
    pl2 = plot(weatherts_C.date,gₛ_C) 
    #pl2 = plot!(weatherts_C.date,gₛ_crit_C)   
    plot(pl1,pl2,layout=(2,1),legends=false)  

    pl1 = plot(weatherts_F.date,Nₘ_f_F)
    pl2 = plot(weatherts_C.date,Nₘ_f_C) 
    plot(pl1,pl2,layout=(2,1),legends=false) 

    plot(weatherts_F.date,Nₘ_f_F)
    plot!(weatherts_C.date,Nₘ_f_C,legends=false)

    plot(weatherts_F.date,αr_F)
    plot!(weatherts_C.date,αr_C)

    plot(weatherts_F.date,weatherts_F.tot_annual_daylight)

    #
    H_vec = range(19.0,stop=25.0,length=50)
    αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₛᵣ₀,Nₛ_C = par[1:10]
    model,kinetic = Initi_model_struct(10,data_F,αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₛᵣ₀)
    npp_val = [first(NPP(H,0.3,0.014,250.0,9.2*10^6,model)) for H in H_vec] 
    Δ_leaf_val = [NPP(H,0.3,0.014,250.0,9.2*10^6,model)[2] for H in H_vec]
    LAI_val = [last(NPP(H,0.3,0.014,250.0,9.2*10^6,model)) for H in H_vec]

    pl1 = plot(LAI_val,npp_val)   
    pl2 = plot(LAI_val,Δ_leaf_val)   
    plot(pl1,pl2,layout=(2,1),legends=false)
end
=#

run_simple_trait_model_C_F_tuning_RO_ts()
#run_C_F_ts_mean()
#validate_RO_2019()
