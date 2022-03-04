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
    
    αf,β₁,β₂,b_Jmax = 460.0,1.27,-0.27,0.0
    Nₛ,rₘ,a_Jmax,Kₓₗ₀,Nₛ_C,i = par[1:6]    

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

    αf,β₁,β₂,b_Jmax = 460.0,1.27,-0.27,0.0
    Nₛ,rₘ,a_Jmax,Kₓₗ₀,Nₛ_C,i = par[1:6] 

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
     
    return  error_C#error_F+error_C
end

function Check_para_RO_C_F_CCPH(par::Array{Float64,1})
    #Nₛ,rₘ,a_Jmax,Kₓₗ₀,Nₛ_C,i,X0,τ,Smax,τ_C
    lo = [0.001,10.0,0.01,0.001,0.001,0.1,-8.0,1.0,10.0,1.0]

    up = [0.1,40.0,1.0,0.1,0.1,6.0,-1.0,15.0,25.0,15.0] 

    return any(par.<lo)||any(par.>up)    
end

function Post_distri_RO_C_F_CCPH(par::Array{Float64,1},RO_data::RO_raw_data)
    if Check_para_RO_C_F_CCPH(par)
        return post = -Inf
    else
        try
            X0,τ,Smax,τ_C = par[7:10]
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
            return post = -Inf
        end
    end
end

function run_simple_trait_model_C_F_tuning_RO_ts()
    #Init Ro data
    RO_data = Load_RO_data()    

    nsamples = 5000
    nchains = 12

    #par_guess = [0.036975489745880365, 17.40739827323768, 0.5, 0.01,
    #0.025675473753997178, -4.230533359595931, 12.328105896346644, 19.079489364590533,12.328105896346644]  
    
    par_guess = [0.011976551772617215, 25.42466609569836, 0.5009242298864488, 0.023503724920524185,
    0.010452800544803813,1.0,-1.312709470268466, 3.3784814235281933, 17.510523050796454, 9.566445194647766]

    x_current = [par_guess for i = 1:nchains]
    P_current = [Post_distri_RO_C_F_CCPH(par_guess,RO_data) for i = 1:nchains] 
        
    println(P_current)    
    
    q_vec = [0.04,24.0,0.033,0.01,0.04,1.0,-4.0,7.0,16.0,7.0]/25.0
    metropolis_mcmc(x_current,P_current,
    x::Array{Float64,1}->Post_distri_RO_C_F_CCPH(x::Array{Float64,1},RO_data),q_vec,"RO_C_F_GPP_Only_Few_Para_NoRoot_20220201";
    n_chains=nchains,n_samples=nsamples,burn_in=7000,sample_freq=5)        
end

function run_C_F_ts_mean()

    RO_data = Load_RO_data()  
    
    file_name = "RO_C_F_GPP_Only_Few_Para_NoRoot_20220201"

    samples = load("./output/"*file_name*".jld","samples_container")    
    P_samples = load("./output/"*file_name*".jld","P_samples_container")
    Draw_trace_plot(samples,file_name)
    Draw_histogram(samples,file_name)

    par = Find_max_sample(samples,P_samples)  
  
    println(par)

    X0,τ,Smax = par[7:10]
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

function test_adaptive_rwm()

    file_name = "adaptive_rwm_RO_C_GPP_Only_20220214"
    n_chains = 4
    n_iter = 40000   

    RO_data = Load_RO_data()

    par_guess = [0.011976551772617215, 25.42466609569836, 0.5009242298864488, 0.023503724920524185,
    0.010452800544803813,1.0,-1.312709470268466, 3.3784814235281933, 17.510523050796454, 9.566445194647766]
    
    out_mat_tot = Array{Matrix{Float64},1}(undef,n_chains)

    Threads.@threads for j = 1:n_chains
        # Run n_iter iterations of the Adaptive Metropolis method:
        out = AdaptiveMCMC.adaptive_rwm(par_guess, 
        x::Array{Float64,1}->Post_distri_RO_C_F_CCPH(x::Array{Float64,1},RO_data),
        n_iter; algorithm=:am) 

        out_mat_tot[j] = transpose(out.X)
    end   

    try
        save("./output/"*file_name*".jld","samples",out_mat_tot)

        @show typeof(out_mat_tot)
        @show size(out_mat_tot)        

        n_features = size(out_mat_tot[1])[2]

        histogram(out_mat_tot[1],layout=n_features,legend=false)
        savefig("./plots/"*file_name*".svg") 

        plot(;layout=n_features,legends=false)        
        for i = 1:n_features
            for j=1:n_chains
                plot!(out_mat_tot[j][:,i],subplot=i)           
            end
        end
        savefig("./plots/"*file_name*"_trace_plot.svg") 

        # Calculate '95% credible intervals':    
        para_stat = mapslices(x->"$(mean(x)) ± $(1.96std(x))", out_mat_tot[1], dims=1)
        for line in para_stat
            println(line)
        end
    catch err
        println("Error: ",err)
    end    
end

function find_max_adaptive_rwm()
    
    n_chains = 12
    file_name = "adaptive_rwm_RO_C_F_GPP_Only_20220203"
    out_mat_tot = load("./output/"*file_name*".jld","samples")
    out_mat = out_mat_tot[1]
    n_features = size(out_mat)[2]
    n_samples = size(out_mat)[1]

    RO_data = Load_RO_data()
    
    #=
    log_p_vec = Array{Float64,1}(undef,n_samples)

    Threads.@threads for j = 1:n_samples
        println("Done: $(j)/$(n_samples)")        
        para_val = out_mat[j,:]
        log_p_vec[j] = Post_distri_RO_C_F_CCPH(para_val,RO_data)
    end

    save("./output/"*file_name*"_log_p.jld","log_p",log_p_vec)
    =#
    
    log_p_vec = load("./output/"*file_name*"_log_p.jld","log_p")

    @show min_ind = argmax(log_p_vec) 
    
    @show size(out_mat)
    @show typeof(out_mat)

    @show par = reshape(mean(out_mat, dims=1),n_features)
    @show par = out_mat[min_ind,:]

    X0,τ,Smax,τ_C = par[7:10]
    weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
    weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ_C,Smax=Smax)
    
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

    plot(weatherts_F.date,data_F.GPP,label="Data",ylabel="GPP")
    pl1 = plot!(weatherts_F.date,GPP_model_F,label="Model")
    plot(weatherts_C.date,data_C.GPP,label="Data",ylabel="GPP")
    pl2 = plot!(weatherts_C.date,GPP_model_C,label="Model")

    plot(weatherts_F.date,data_F.gₛ,label="Data",ylabel="gₛ")
    pl3 = plot!(weatherts_F.date,gₛ_model_F,label="Model")
    plot(weatherts_C.date,data_C.gₛ,label="Data",ylabel="gₛ")
    pl4 = plot!(weatherts_C.date,gₛ_model_C,label="Model")

    plot(pl1,pl3,pl2,pl4,layout=(2,2),legends=false)
    savefig("./plots/"*file_name*"_result.svg")

    # Calculate '95% credible intervals':    
    para_stat = mapslices(x->"$(mean(x)) ± $(1.96std(x))", out_mat, dims=1)
    for line in para_stat
        println(line)
    end
end

#run_simple_trait_model_C_F_tuning_RO_ts()
#run_C_F_ts_mean()
#test_adaptive_rwm()
find_max_adaptive_rwm()