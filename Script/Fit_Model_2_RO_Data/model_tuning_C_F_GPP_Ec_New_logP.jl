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
    EC_model_F = zeros(sim_steps+1)  
    for j = 1:sim_steps+1           
       
        model,kinetic = Initi_model_struct(j,data_F,αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₓₗ₀)
        
        GPP_model_F[j],EC_model_F[j] = CalcModelOutput!(j,model,weatherts_F,kinetic)
    end  

    #Control stand
    GPP_model_C = zeros(sim_steps+1)
    EC_model_C = zeros(sim_steps+1)
    for j = 1:sim_steps+1        
        
        model,kinetic = Initi_model_struct(j,data_C,αf,β₁,β₂,Nₛ_C,rₘ,a_Jmax,b_Jmax,i,Kₓₗ₀)
        
        GPP_model_C[j],EC_model_C[j] = CalcModelOutput!(j,model,weatherts_C,kinetic)
    end  

    return GPP_model_F,GPP_model_C,EC_model_F,EC_model_C
end

function calc_logP_term(ydata::T,ymodel::T,a::T,b::T) where {T<:Float64}
    e = abs(ydata-ymodel)
    σ = a+b*ymodel
    σ>0.0||error("Negative variance")
    return e/σ+log(σ)
end

function Calc_logP(
    GPP_model_F::Array{T,1},Ec_model_F::Array{T,1},data_F::RoData,
    GPP_model_C::Array{T,1},Ec_model_C::Array{T,1},data_C::RoData,
    a_GPP::T,b_GPP::T,a_Ec::T,b_Ec::T) where {T<:Float64}     

    logP = 0.0
    for i = 1:length(data_F.GPP)
        logP += calc_logP_term(data_F.GPP[i],GPP_model_F[i],a_GPP,b_GPP)
        logP += calc_logP_term(data_C.GPP[i],GPP_model_C[i],a_GPP,b_GPP)

        logP += calc_logP_term(data_F.E_C[i],Ec_model_F[i],a_Ec,b_Ec)
        logP += calc_logP_term(data_C.E_C[i],Ec_model_C[i],a_Ec,b_Ec)        
    end    
     
    return  -logP
end

function Check_para_RO_C_F_CCPH(par::Array{Float64,1})
    #Nₛ,rₘ,a_Jmax,Kₓₗ₀,Nₛ_C,i,X0,τ,Smax,τ_C,a_GPP,b_GPP,a_Ec,b_Ec
    lo = [0.001,10.0,0.01,0.001,0.001,0.1,-18.0,1.0,10.0,1.0,0.0001,0.0001,0.0001,0.0001]

    up = [0.1,40.0,1.0,0.1,0.1,6.0,-0.9,15.0,25.0,15.0,5.0,3.0,5.0,3.0] 

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

            GPP_model_F,GPP_model_C,EC_model_F,EC_model_C = Run_RO_C_F_CCPH(par,weatherts_F,weatherts_C,data_F,data_C)              

            a_GPP,b_GPP,a_Ec,b_Ec = par[11:14]

            logP = Calc_logP(
            GPP_model_F,EC_model_F,data_F,
            GPP_model_C,EC_model_C,data_C,
            a_GPP,b_GPP,a_Ec,b_Ec)
            
            return post = logP
        catch err
            println("Parameters Error: ", err)
            return post = -Inf
        end
    end
end

function test_adaptive_rwm()

    file_name = "adaptive_rwm_RO_C_GPP_Ec_20220228"
    n_chains = 4
    n_iter = 30000   

    RO_data = Load_RO_data()

    #=
    par_guess = [0.011976551772617215, 25.42466609569836, 0.5009242298864488, 0.023503724920524185,
    0.010452800544803813,1.0,-1.312709470268466, 3.3784814235281933, 17.510523050796454, 9.566445194647766,
    0.237,0.0845,0.0688,0.146]
    =#

    par_guess = [0.004320402400586008, 27.297822676817912, 0.029549759752040473, 0.0693215277873528,
    0.0023441457252274013, 2.2501166086478843, -6.776410663775599, 13.19517957686665,
    20.362175060472445, 2.6466352227229257, 0.4645185292099437, 0.08043063207519931,
    0.14985616067841498, 0.1500009201257627]

    @show Post_distri_RO_C_F_CCPH(par_guess,RO_data)
    
    out_mat_tot = Array{Matrix{Float64},1}(undef,n_chains)

    Threads.@threads for j = 1:n_chains
        # Run n_iter iterations of the Adaptive Metropolis method:
        out = AdaptiveMCMC.adaptive_rwm(par_guess, 
        x::Array{Float64,1}->Post_distri_RO_C_F_CCPH(x::Array{Float64,1},RO_data),
        n_iter; algorithm=:ram) 

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
    
    file_name = "adaptive_rwm_RO_C_GPP_Ec_20220228"
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
    @show log_p_vec[min_ind]
    @show par = reshape(mean(out_mat, dims=1),n_features)
    @show par = out_mat[min_ind,:]  
    
    @show Post_distri_RO_C_F_CCPH(par,RO_data)

    X0,τ,Smax,τ_C = par[7:10]
    weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
    weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ_C,Smax=Smax)
    
    data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C,weatherts_F,weatherts_C)    

    GPP_model_F,GPP_model_C,EC_model_F,EC_model_C = Run_RO_C_F_CCPH(par,weatherts_F,weatherts_C,data_F,data_C) 

    a_GPP,b_GPP,a_Ec,b_Ec = par[11:14]   

    @show GPP_R²_F =  calcR²(data_F.GPP,GPP_model_F)
    @show GPP_cor_F = cor(data_F.GPP,GPP_model_F)
    @show gₛ_R²_F =  calcR²(data_F.E_C,EC_model_F)
    @show gₛ_cor_F = cor(data_F.E_C,EC_model_F)

    @show GPP_R²_C =  calcR²(data_C.GPP,GPP_model_C)
    @show GPP_cor_C = cor(data_C.GPP,GPP_model_C)
    @show gₛ_R²_C =  calcR²(data_C.E_C,EC_model_C)
    @show gₛ_cor_C = cor(data_C.E_C,EC_model_C)

    println("Fertilized")
    clac_GPP_R²_annual(data_C.GPP,GPP_model_F,weatherts_F.date)
    println("Control")  
    clac_GPP_R²_annual(data_C.GPP,GPP_model_C,weatherts_C.date)
    
    σ_GPP_F = a_GPP.+GPP_model_F*b_GPP
    σ_GPP_C = a_GPP.+GPP_model_C*b_GPP

    σ_EC_F = a_Ec.+b_Ec*EC_model_F
    σ_EC_C = a_Ec.+b_Ec*EC_model_C

    c = -log(0.025*2) #Use for calculating 95% credible intervals

    plot(weatherts_F.date,data_F.GPP,label="Data",ylabel="GPP")
    plot!(weatherts_F.date,GPP_model_F+c*σ_GPP_F)
    plot!(weatherts_F.date,GPP_model_F-c*σ_GPP_F)
    pl1 = plot!(weatherts_F.date,GPP_model_F,label="Model")
    plot(weatherts_C.date,data_C.GPP,label="Data",ylabel="GPP")
    plot!(weatherts_C.date,GPP_model_C+c*σ_GPP_C)
    plot!(weatherts_C.date,GPP_model_C-c*σ_GPP_C)
    pl2 = plot!(weatherts_C.date,GPP_model_C,label="Model")

    plot(weatherts_F.date,data_F.E_C,label="Data",ylabel="E_C")
    plot!(weatherts_F.date,EC_model_F+c*σ_EC_F)
    plot!(weatherts_F.date,EC_model_F-c*σ_EC_F)
    pl3 = plot!(weatherts_F.date,EC_model_F,label="Model")
    plot(weatherts_C.date,data_C.E_C,label="Data",ylabel="E_C")
    plot!(weatherts_C.date,EC_model_C+c*σ_EC_C)
    plot!(weatherts_C.date,EC_model_C-c*σ_EC_C)
    pl4 = plot!(weatherts_C.date,EC_model_C,label="Model")

    plot(pl1,pl3,pl2,pl4,layout=(2,2),legends=false)
    savefig("./plots/"*file_name*"_result.svg")

    # Calculate '95% credible intervals':    
    para_stat = mapslices(x->"$(mean(x)) ± $(1.96std(x))", out_mat, dims=1)
    for line in para_stat
        println(line)
    end
end

function run_simple_trait_model_C_F_tuning_RO_ts()
    #Init Ro data
    RO_data = Load_RO_data()    

    nsamples = 5000
    nchains = 12

    par_guess = [0.0028160038624083067, 27.33276184313149,
    0.028643648082103606, 0.09336572808750652,
     0.0019727630599060625, 2.2083338895077254,
      -6.817991596443347, 13.184687442117964,
       20.387232607829905, 2.6564297126274803,
        0.434540406001315, 0.06910534451564838,
         0.13628204854679268, 0.1383012394158192]

    x_current = [par_guess for i = 1:nchains]
    P_current = [Post_distri_RO_C_F_CCPH(par_guess,RO_data) for i = 1:nchains] 
        
    println(P_current)    
    
    q_vec = [0.04,24.0,0.033,0.01,0.04,1.0,4.0,7.0,16.0,7.0,
    0.237,0.0845,0.0688,0.146]/40.0

    metropolis_mcmc(x_current,P_current,
    x::Array{Float64,1}->Post_distri_RO_C_F_CCPH(x::Array{Float64,1},RO_data),q_vec,"adaptive_rwm_RO_C_GPP_Ec_20220301";
    n_chains=nchains,n_samples=nsamples,burn_in=7000,sample_freq=5)        
end

run_simple_trait_model_C_F_tuning_RO_ts()
#test_adaptive_rwm()
#find_max_adaptive_rwm()