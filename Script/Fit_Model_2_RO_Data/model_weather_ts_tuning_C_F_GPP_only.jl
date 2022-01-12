function Initi_model_struct(H,N,Wf,B,αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₛᵣ₀)
    cons = Constants()
    env = EnvironmentStruct() 
    kinetic = PhotoKineticRates()
    photo = PhotoPar(kinetic,env.Tₐ)    
    hydPar = HydraulicsPar(;Kₛᵣ₀=Kₛᵣ₀,i=i)
    treepar = TreePar(;αf=αf,β₁=β₁,β₂=β₂,rₘ=rₘ,Nₛ=Nₛ,a_Jmax=a_Jmax,b_Jmax=b_Jmax)
   
    Hc = 7.3   
    Hs = H-Hc   
    As = Wf/treepar.αf
    Ww = treepar.ρw*As*(treepar.β₁*H+treepar.β₂*Hs)    
    treesize = TreeSize(Wf,Ww,H,Hs,As,B,N)

    model = CCPHStruct(cons,env,treepar,treesize,photo,hydPar)

    return model,kinetic
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

function Initi_model_struct(j::Integer,data::RoData,
    αf::T,β₁::T,β₂::T,Nₛ::T,rₘ::T,a_Jmax::T,b_Jmax::T,i::T,Kₛᵣ₀::T) where {T<:Float64}
    data_ind = Find_data_ind(j)
    H = data.H[data_ind]
    N = data.N[data_ind]
    Wf = data.Wf[data_ind]
    B = data.B[data_ind]
    model,kinetic = Initi_model_struct(H,N,Wf,B,αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₛᵣ₀)    

    return model,kinetic
end

function Run_RO_C_F_CCPH(par::Array{Float64,1},weatherts_F::WeatherTS,
    weatherts_C::WeatherTS,data_F::RoData,data_C::RoData;sim_steps::Integer=83)

    αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₛᵣ₀,Nₛ_C = par[1:10]    

    #Fertilized stand
    GPP_model_F = zeros(sim_steps+1)
    for j = 1:sim_steps+1
        data_ind = Find_data_ind(j)
        H = data_F.H[data_ind]
        N = data_F.N[data_ind]
        Wf = data_F.Wf[data_ind]
        B = data_F.B[data_ind]
        model,kinetic = Initi_model_struct(H,N,Wf,B,αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₛᵣ₀)
        GPP_model_F[j],Transpiraiton = CalcGPP!(j,model,weatherts_F,kinetic)
    end  

    #Control stand
    GPP_model_C = zeros(sim_steps+1)
    for j = 1:sim_steps+1
        data_ind = Find_data_ind(j)
        H = data_C.H[data_ind]
        N = data_C.N[data_ind]
        Wf = data_C.Wf[data_ind]
        B = data_C.B[data_ind]
        model,kinetic = Initi_model_struct(H,N,Wf,B,αf,β₁,β₂,Nₛ_C,rₘ,a_Jmax,b_Jmax,i,Kₛᵣ₀)
        GPP_model_C[j],Transpiraiton = CalcGPP!(j,model,weatherts_C,kinetic)
    end  

    return GPP_model_F,GPP_model_C
end

function RO_C_F_CCPH_ModelOutput(par::Array{Float64,1},weatherts_F::WeatherTS,
    weatherts_C::WeatherTS,data_F::RoData,data_C::RoData;sim_steps::Integer=83)

    αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₛᵣ₀,Nₛ_C = par[1:10]

    #Fertilized stand
    ModelOutput_F = CCPHOutput[]
    k_cost_max_F = Float64[]
    gₛ_F = Float64[]
    gₛ_crit_F = Float64[]
    Nₘ_f_F = Float64[]
    for j = 1:sim_steps+1
        data_ind = Find_data_ind(j)
        H = data_F.H[data_ind]
        N = data_F.N[data_ind]
        Wf = data_F.Wf[data_ind]
        B = data_F.B[data_ind]
        model,kinetic = Initi_model_struct(H,N,Wf,B,αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₛᵣ₀)
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
        data_ind = Find_data_ind(j)
        H = data_C.H[data_ind]
        N = data_C.N[data_ind]
        Wf = data_C.Wf[data_ind]
        B = data_C.B[data_ind]
        model,kinetic = Initi_model_struct(H,N,Wf,B,αf,β₁,β₂,Nₛ_C,rₘ,a_Jmax,b_Jmax,i,Kₛᵣ₀)
        modeloutput,gₛ_opt,Nₘ_f_opt = OptimCCPH!(j,model,weatherts_C,kinetic)
        push!(ModelOutput_C,modeloutput)
        push!(k_cost_max_C,Max_k_cost(model))
        push!(gₛ_C,gₛ_opt)
        push!(gₛ_crit_C,CCPH.Calc_K_costᵢₙᵥ(0.12,model))
        push!(Nₘ_f_C,Nₘ_f_opt)
    end  

    return ModelOutput_F,ModelOutput_C,k_cost_max_F,k_cost_max_C,gₛ_F,gₛ_C,gₛ_crit_F,gₛ_crit_C,Nₘ_f_F,Nₘ_f_C
end

function Calc_error_RO_C_F_CCPH(GPP_model_F::Array{Float64,1},data_F::RoData,
    GPP_model_C::Array{Float64,1},data_C::RoData)   

    #[21,44,64,84]
    ind_start = [1,22,45,65]
    ind_end = [21,44,64,84]
    #Fertilized stand

    GPP_error_F = 0.0
    for i = 1:4
        GPP_error_F += 1-calcR²(data_F.GPP[ind_start[i]:ind_end[i]],GPP_model_F[ind_start[i]:ind_end[i]])
    end   

    #Control stand
    GPP_error_C = 0.0
    for i = 1:4
        GPP_error_C += 1-calcR²(data_C.GPP[ind_start[i]:ind_end[i]],GPP_model_C[ind_start[i]:ind_end[i]])
    end
     
    return  GPP_error_F+GPP_error_C
end

function Check_para_RO_C_F_CCPH(par::Array{Float64,1})
    lo = [390.0,1.1,-0.4,0.02,16.0,0.001,0.5e-5,0.5,5.0,0.02,-8.0,3.0,10.0,3.0]
    up = [750.0,1.5,-0.1,0.1,40.0,0.1,2.4e-5,10.0,10.0,0.1,-1.0,15.0,20.0,15.0] 

    return any(par.<lo)||any(par.>up)    
end

function Post_distri_RO_C_F_CCPH(par::Array{Float64,1},RO_data::RO_raw_data,RO_ET_data::RO_Raw_ET_Data)
    if Check_para_RO_C_F_CCPH(par)
        return post = 0.0
    else
        try
            X0,τ,Smax,τ_C = par[11:14]
            weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
            weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ_C,Smax=Smax)

            ET_data_F = RO_ET_Est(RO_ET_data,weatherts_F;stand_type="Fertilized")
            ET_data_C = RO_ET_Est(RO_ET_data,weatherts_C;stand_type="Control")
            
            data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C,ET_data_F,ET_data_C)   

            GPP_model_F,GPP_model_C = Run_RO_C_F_CCPH(par,weatherts_F,weatherts_C,data_F,data_C)              

            error_val = Calc_error_RO_C_F_CCPH(GPP_model_F,data_F,GPP_model_C,data_C)
            
            return post = exp(-error_val)
        catch err
            println("Parameters Error: ", err)
            return post = 0.0
        end
    end
end

function q_distri_RO_C_F_CCPH(par;sd::Float64=25.0)
    return par.+rand(Normal(),14).*[460.0,1.27,-0.27,0.04,24.0,0.033,1.21e-5,1.0,1.0,0.04,-4.0,7.0,16.0,7.0]/sd          
end

function run_simple_trait_model_C_F_tuning_RO_ts()
    #Init Ro data
    RO_data = Load_RO_data()    
    RO_ET_data = RO_Raw_ET_Data()

    nsamples = 3000
    nchains = 12

    par_guess = [470.328512441229, 1.2501809959672654, -0.2705395150988396,
    0.036975489745880365, 17.40739827323768, 0.015721074668478767/0.24321499065796873,
    3.4357734812576556e-6/0.24321499065796873, 1.0749476181390065, 8.621674950871828,
    0.025675473753997178, -4.230533359595931, 12.328105896346644, 19.079489364590533,12.328105896346644]   

    x_current = [par_guess for i = 1:nchains]
    P_current = [Post_distri_RO_C_F_CCPH(par_guess,RO_data,RO_ET_data) for i = 1:nchains] 
        
    println(P_current) 
    
    q_vec = [460.0,1.27,-0.27,0.04,24.0,0.033,1.21e-5,1.0,1.0,0.04,-4.0,7.0,16.0,7.0]/25.0
    metropolis_mcmc(x_current,P_current,
    x::Array{Float64,1}->Post_distri_RO_C_F_CCPH(x::Array{Float64,1},RO_data,RO_ET_data),q_vec,"RO_C_F_GPP_Only_20211202";
    n_chains=nchains,n_samples=nsamples,burn_in=4000,sample_freq=5)     
end

function run_C_F_ts_mean()

    RO_data = Load_RO_data()  
    RO_ET_data = RO_Raw_ET_Data()

    samples = load("output/RO_C_F_GPP_Only_20211202.jld","samples_container")    
    P_samples = load("output/RO_C_F_GPP_Only_20211202.jld","P_samples_container")

    par = Find_max_sample(samples,P_samples)

    println(par)

    X0,τ,Smax = par[11:13]
    weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
    weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ,Smax=Smax)

    ET_data_F = RO_ET_Est(RO_ET_data,weatherts_F;stand_type="Fertilized")
    ET_data_C = RO_ET_Est(RO_ET_data,weatherts_C;stand_type="Control")
    
    data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C,ET_data_F,ET_data_C)    

    GPP_model_F,GPP_model_C = Run_RO_C_F_CCPH(par,weatherts_F,weatherts_C,data_F,data_C) 
    
    @show GPP_R²_F =  calcR²(data_F.GPP[65:84],GPP_model_F[65:84])
    @show GPP_R²_C =  calcR²(data_C.GPP[65:84],GPP_model_C[65:84])

    @show GPP_R²_F =  calcR²(data_F.GPP,GPP_model_F)
    @show GPP_R²_C =  calcR²(data_C.GPP,GPP_model_C)

    plot(weatherts_F.date,data_F.GPP,label="Data")
    pl1 = plot!(weatherts_F.date,GPP_model_F,label="Model")
    plot(weatherts_C.date,data_C.GPP,label="Data")
    pl2 = plot!(weatherts_C.date,GPP_model_C,label="Model")
    plot(pl1,pl2,layout=(2,1),legends=false)
end

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

#run_simple_trait_model_C_F_tuning_RO_ts()
#error_RO_C_F_dist()
run_C_F_ts_mean()
#validate_RO_2019()
