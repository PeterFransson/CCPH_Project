function Run_RO_C_F_CCPH(par::Array{Float64,1},weatherts_F::WeatherTS,
    weatherts_C::WeatherTS;sim_steps::Integer=83)

    αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₛᵣ₀,Nₛ_C = par[1:10]

    #Fertilized stand
    cons = Constants()
    env = EnvironmentStruct() 
    kinetic = PhotoKineticRates()
    photo = PhotoPar(kinetic,env.Tₐ)    
    hydPar = HydraulicsPar(;Kₛᵣ₀=Kₛᵣ₀,i=i)
    treepar = TreePar(;αf=αf,β₁=β₁,β₂=β₂,rₘ=rₘ,Nₛ=Nₛ,a_Jmax=a_Jmax,b_Jmax=b_Jmax)

    H = 16.7
    Hc = 7.3   
    Hs = H-Hc
    N = 857.0/10000.0 #Initial stem density (#/m2)
    Wf = 4.54
    As = Wf/treepar.αf
    Ww = treepar.ρw*As*(treepar.β₁*H+treepar.β₂*Hs)
    B = 0.024
    treesize = TreeSize(Wf,Ww,H,Hs,As,B,N)

    model = CCPHStruct(cons,env,treepar,treesize,photo,hydPar)

    ccphts_F = CCPHStandGrowth!(model,weatherts_F,kinetic,sim_steps)

    #Control stand
    cons = Constants()
    env = EnvironmentStruct() 
    kinetic = PhotoKineticRates()
    photo = PhotoPar(kinetic,env.Tₐ)    
    hydPar = HydraulicsPar(;Kₛᵣ₀=Kₛᵣ₀,i=i)
    treepar = TreePar(;αf=αf,β₁=β₁,β₂=β₂,rₘ=rₘ,Nₛ=Nₛ_C,a_Jmax=a_Jmax,b_Jmax=b_Jmax)

    H = 16.7
    Hc = 7.3   
    Hs = H-Hc
    N = 1010.0/10000.0 #Initial stem density (#/m2)
    Wf = 4.18
    As = Wf/treepar.αf
    Ww = treepar.ρw*As*(treepar.β₁*H+treepar.β₂*Hs)
    B = 0.0225
    treesize = TreeSize(Wf,Ww,H,Hs,As,B,N)

    model = CCPHStruct(cons,env,treepar,treesize,photo,hydPar)

    ccphts_C = CCPHStandGrowth!(model,weatherts_C,kinetic,sim_steps)  

    return ccphts_F,ccphts_C
end

function Calc_error_RO_C_F_CCPH(ccphts_F::CCPHTS,data_F::RoData,weatherts_F::WeatherTS,
    ccphts_C::CCPHTS,data_C::RoData,weatherts_C::WeatherTS)

    #Fertilized stand
    GPP_model = CalcGPP(ccphts_F,weatherts_F)   

    H_model = ccphts_F.H[[21,44,64,84]]
    B_model = ccphts_F.B[[21,44,64,84]]
    Wf_model = ccphts_F.Wf[[21,44,64,84]]

    GPP_error = 1-calcR²(data_F.GPP,GPP_model)
    H_error = 1-calcR²(data_F.H,H_model)
    B_error = 1-calcR²(data_F.B,B_model)
    Wf_error = 1-calcR²(data_F.Wf,Wf_model)

    error_F = GPP_error+H_error+B_error+Wf_error

    #Control stand
    GPP_model = CalcGPP(ccphts_C,weatherts_C)

    H_model = ccphts_C.H[[21,44,64,84]]
    B_model = ccphts_C.B[[21,44,64,84]]
    Wf_model = ccphts_C.Wf[[21,44,64,84]]

    GPP_error = 1-calcR²(data_C.GPP,GPP_model)
    H_error = 1-calcR²(data_C.H,H_model)
    B_error = 1-calcR²(data_C.B,B_model)
    Wf_error = 1-calcR²(data_C.Wf,Wf_model)

    error_C = GPP_error+H_error+B_error+Wf_error
     
    return  error_F+error_C
end

function Check_para_RO_C_F_CCPH(par::Array{Float64,1})
    lo = [390.0,1.1,-0.4,0.02,16.0,0.001,0.5e-5,0.5,5.0,0.02,-8.0,3.0,10.0]
    up = [750.0,1.5,-0.1,0.1,40.0,0.1,2.4e-5,10.0,10.0,0.1,-1.0,15.0,20.0] 

    return any(par.<lo)||any(par.>up)    
end

function Log_Post_distri_RO_C_F_CCPH(par::Array{Float64,1},RO_data::RO_raw_data)
    if Check_para_RO_C_F_CCPH(par)
        return post = -Inf
    else
        try
            X0,τ,Smax = par[11:13]
            weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
            weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ,Smax=Smax)
            
            data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C)   

            ccphts_F,ccphts_C = Run_RO_C_F_CCPH(par,weatherts_F,weatherts_C)              

            error_val = Calc_error_RO_C_F_CCPH(ccphts_F,data_F,weatherts_F,ccphts_C,data_C,weatherts_C)
            
            return post = -error_val
        catch err
            println("Parameters Error: ", err)
            return post = -Inf
        end
    end
end

function q_distri_RO_C_F_CCPH(par;sd::Float64=25.0)
    return par.+rand(Normal(),13).*[460.0,1.27,-0.27,0.04,24.0,0.033,1.21e-5,1.0,1.0,0.04,-4.0,7.0,16.0]/sd          
end

function run_simple_trait_model_C_F_tuning_RO_ts()
    #Init Ro data
    RO_data = Load_RO_data()    

    nsamples = 1000
    nchains = 12

    par_guess = [470.328512441229, 1.2501809959672654, -0.2705395150988396,
    0.036975489745880365, 17.40739827323768, 0.015721074668478767/0.24321499065796873,
    3.4357734812576556e-6/0.24321499065796873, 1.0749476181390065, 8.621674950871828,
    0.025675473753997178, -4.230533359595931, 12.328105896346644, 19.079489364590533]

    x_current = [par_guess for i = 1:nchains]
    P_current = [Log_Post_distri_RO_C_F_CCPH(par_guess,RO_data) for i = 1:nchains] 
        
    println(P_current)  
    
    
    q_vec = [460.0,1.27,-0.27,0.04,24.0,0.033,1.21e-5,1.0,1.0,0.04,-4.0,7.0,16.0]/110.0
    metropolis_mcmc(x_current,P_current,
    x::Array{Float64,1}->Log_Post_distri_RO_C_F_CCPH(x::Array{Float64,1},RO_data),q_vec,"RO_C_F_w_Wf_20211116";
    n_chains=nchains,n_samples=nsamples,burn_in=3000,sample_freq=5)     
end

#=
function run_C_F_ts_mean()

    RO_data = Load_RO_data()  
    samples = load("output/RO_C_F_20211105.jld","samples_container")    
    P_samples = load("output/RO_C_F_20211105.jld","P_samples_container")

    par = Find_max_sample(samples,P_samples)

    println(par)

    X0,τ,Smax = par[11:13]
    weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
    weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ,Smax=Smax)

    data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C)

    ccphts_F,ccphts_C = Run_RO_C_F_CCPH(par,weatherts_F,weatherts_C)  

    #Fertilized stand
    GPP_model = CalcGPP(ccphts_F,weatherts_F)

    plot(weatherts_F.date,GPP_data_F,label="Data F",ylabel="g C m⁻² d⁻¹")
    pl1 = plot!(weatherts_F.date,GPP_model,label="Model F")
    
    GPP_R² = calcR²(GPP_data_F,GPP_model)   
    GPP_corr = cor(GPP_data_F,GPP_model)
    
    H_data = [17.004166666666666, 17.271666666666665, 17.59166666666667, 17.835] 
    B_data = [0.025076878666666667, 0.02603588483333334, 0.026922637833333332, 0.027666997]

    H_model = ccphts_F.H[[21,44,64,84]]
    B_model = ccphts_F.B[[21,44,64,84]]
    
    H_R² = calcR²(H_data,H_model)
    B_R² = calcR²(B_data,B_model)

    plot(2015:2018,H_data,xlabel="Date",ylabel="Height",label="Data F")
    pl2 = plot!(2015:2018,H_model,label="Model F")
    plot(2015:2018,B_data,xlabel="Date",ylabel="Basal Area",label="Data F")
    pl3 = plot!(2015:2018,B_model,label="Model F")

    println("Fertilized: ")
    println("GPP R²:",GPP_R²) 
    println("GPP Pearson correlation: ",GPP_corr)
    println("H R²:",H_R²) 
    println("B R²:",B_R²) 
    clac_GPP_R²_annual(GPP_data_F,GPP_model,weatherts_F.date)  

    #Control stand
    GPP_model = CalcGPP(ccphts_C,weatherts_C)

    plot(weatherts_C.date,GPP_data_C,label="Data C",ylabel="g C m⁻² d⁻¹")
    pl4 = plot!(weatherts_C.date,GPP_model,label="Model C")

    GPP_R² = calcR²(GPP_data_C,GPP_model)  
    GPP_corr =  cor(GPP_data_C,GPP_model) 

    H_data = [16.956,17.165,17.388,17.590] 
    B_data = [0.023015244833333334, 0.02349805816666667, 0.023938368499999998, 0.02432245283333333]

    H_model = ccphts_C.H[[21,44,64,84]]
    B_model = ccphts_C.B[[21,44,64,84]]

    H_R² = calcR²(H_data,H_model)
    B_R² = calcR²(B_data,B_model)

    plot(2015:2018,H_data,xlabel="Date",ylabel="Height",label="Data C")
    pl5 = plot!(2015:2018,H_model,label="Model C")
    plot(2015:2018,B_data,xlabel="Date",ylabel="Basal Area",label="Data C")
    pl6 = plot!(2015:2018,B_model,label="Model C")

    println("Control: ")
    println("GPP R²:",GPP_R²) 
    println("GPP Pearson correlation: ",GPP_corr)
    println("H R²:",H_R²) 
    println("B R²:",B_R²)   
    clac_GPP_R²_annual(GPP_data_C,GPP_model,weatherts_C.date)    

    plot(pl1,pl2,pl3,pl4,pl5,pl6,layout=6,legend=false)
    savefig("./plots/RO_C_F_ts_tuned_20211105.svg")
end

function validate_RO_2019()
    RO_data = Load_RO_data(;weather_file = "Weather_RO_2")  
    samples = load("output/RO_C_F_GPP_Only_20211111.jld","samples_container")    
    P_samples = load("output/RO_C_F_GPP_Only_20211111.jld","P_samples_container")

    par = Find_max_sample(samples,P_samples)

    println(par)

    X0,τ,Smax = par[11:13]
    weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax,
    weather_file = "Weather_RO_2")
    weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ,Smax=Smax,
    weather_file = "Weather_RO_2")

    data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C)  

    ccphts_F,ccphts_C = Run_RO_C_F_CCPH(par,weatherts_F,weatherts_C;sim_steps=104)
    
    GPP_F_model = CalcGPP(ccphts_F,weatherts_F)
    GPP_C_model = CalcGPP(ccphts_C,weatherts_C)

    println("Fertilized:")
    GPP_F_R² = calcR²(GPP_data_F,GPP_F_model)  
    GPP_F_corr =  cor(GPP_data_F,GPP_F_model) 

    println("R²: $(GPP_F_R²)")
    println("Corr: $(GPP_F_corr)")
    clac_GPP_R²_annual(GPP_data_F,GPP_F_model,weatherts_F.date) 

    println("Control:")
    GPP_F_R² = calcR²(GPP_data_F,GPP_F_model)  
    GPP_F_corr =  cor(GPP_data_F,GPP_F_model) 
    
    println("R²: $(GPP_F_R²)")
    println("Corr: $(GPP_F_corr)")
    clac_GPP_R²_annual(GPP_data_C,GPP_C_model,weatherts_C.date) 

    plot(weatherts_F.date,GPP_data_F)
    pl1 = plot!(weatherts_F.date,GPP_F_model)
    plot([2.0,9.0],[2.0,9.0])
    pl2 = plot!(GPP_data_F,GPP_F_model,seriestype=:scatter)

    plot(weatherts_C.date,GPP_data_C)
    pl3 = plot!(weatherts_C.date,GPP_C_model)
    plot([2.0,9.0],[2.0,9.0])
    pl4 = plot!(GPP_data_C,GPP_C_model,seriestype=:scatter)

    plot(pl1,pl2,pl3,pl4,layout=(2,2),legends=false)
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
=#

run_simple_trait_model_C_F_tuning_RO_ts()
#run_C_F_ts_mean()
#validate_RO_2019()