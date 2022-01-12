#=
mutable struct RoData{T<:Float64}
    H::Array{T,1}
    B::Array{T,1}
    GPP::Array{T,1}
    Wf::Array{T,1}
end

function Create_RoData_C_F(GPP_data_F::Array{T,1},GPP_data_C::Array{T,1}) where {T<:Float64}   

    H_data_F = [17.004166666666666, 17.271666666666665, 17.59166666666667, 17.835] 
    B_data_F = [0.025076878666666667, 0.02603588483333334, 0.026922637833333332, 0.027666997]
    Wf_data_F = [4.62,4.72,4.77,4.82]

    H_data_C = [16.956,17.165,17.388,17.590] 
    B_data_C = [0.023015244833333334, 0.02349805816666667, 0.023938368499999998, 0.02432245283333333]
    Wf_data_C = [4.18,4.20,4.20,4.2]

    data_F = RoData(H_data_F,B_data_F,GPP_data_F,Wf_data_F)
    data_C = RoData(H_data_C,B_data_C,GPP_data_C,Wf_data_C)

    return data_F,data_C
end

function calcR²(y_data,f_model)
    ymean = mean(y_data)
    SStot = sum((y_data.-ymean).^2)
    SSres = sum((y_data.-f_model).^2)
    return 1-SSres/SStot
end


#Calc model GPP
function CalcGPP(ccphts::CCPHTS,weatherts::WeatherTS)
    GPP = ccphts.N.*ccphts.P.*weatherts.daylight./weatherts.tot_annual_daylight*1000/7 #g C day⁻¹ m⁻² ground area

    return GPP
end
=#

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
    Wf = 8.27#4.54
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
    Wf = 5.03#4.2
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

    GPP_error = 1-calcR²(data_F.GPP,GPP_model)
    H_error = 1-calcR²(data_F.H,H_model)
    B_error = 1-calcR²(data_F.B,B_model)

    error_F = GPP_error+H_error+B_error

    #Control stand
    GPP_model = CalcGPP(ccphts_C,weatherts_C)

    H_model = ccphts_C.H[[21,44,64,84]]
    B_model = ccphts_C.B[[21,44,64,84]]

    GPP_error = 1-calcR²(data_C.GPP,GPP_model)
    H_error = 1-calcR²(data_C.H,H_model)
    B_error = 1-calcR²(data_C.B,B_model)

    error_C = GPP_error+H_error+B_error
     
    return  error_F+error_C
end

function Check_para_RO_C_F_CCPH(par::Array{Float64,1})
    lo = [390.0,1.1,-0.4,0.02,16.0,0.001,0.5e-5,0.5,5.0,0.02,-8.0,3.0,10.0]
    up = [750.0,1.5,-0.1,0.1,40.0,0.1,2.4e-5,10.0,10.0,0.1,-1.0,15.0,20.0] 

    return any(par.<lo)||any(par.>up)    
end

function Post_distri_RO_C_F_CCPH(par::Array{Float64,1},RO_data::RO_raw_data)
    if Check_para_RO_C_F_CCPH(par)
        return post = 0.0
    else
        try
            X0,τ,Smax = par[11:13]
            weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
            weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ,Smax=Smax)
            
            data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C)   

            ccphts_F,ccphts_C = Run_RO_C_F_CCPH(par,weatherts_F,weatherts_C)              

            error_val = Calc_error_RO_C_F_CCPH(ccphts_F,data_F,weatherts_F,ccphts_C,data_C,weatherts_C)
            
            return post = exp(-error_val)
        catch err
            println("Parameters Error: ", err)
            return post = 0.0
        end
    end
end

function q_distri_RO_C_F_CCPH(par;sd::Float64=25.0)
    return par.+rand(Normal(),13).*[460.0,1.27,-0.27,0.04,24.0,0.033,1.21e-5,1.0,1.0,0.04,-4.0,7.0,16.0]/sd          
end


function run_simple_trait_model_C_F_tuning_RO_ts()
    #Init Ro data
    RO_data = Load_RO_data()    

    nsamples = 4000
    nchains = 12

    par_guess = [470.328512441229, 1.2501809959672654, -0.2705395150988396,
    0.036975489745880365, 17.40739827323768, 0.015721074668478767/0.24321499065796873,
    3.4357734812576556e-6/0.24321499065796873, 1.0749476181390065, 8.621674950871828,
    0.025675473753997178, -4.230533359595931, 12.328105896346644, 19.079489364590533]

    x_current = [par_guess for i = 1:nchains]
    P_current = [Post_distri_RO_C_F_CCPH(par_guess,RO_data) for i = 1:nchains] 
        
    println(P_current)  

    q_vec = [460.0,1.27,-0.27,0.04,24.0,0.033,1.21e-5,1.0,1.0,0.04,-4.0,7.0,16.0]/25.0
    metropolis_mcmc(x_current,P_current,
    x::Array{Float64,1}->Post_distri_RO_C_F_CCPH(x::Array{Float64,1},RO_data),q_vec,"RO_C_F_20211124";
    n_chains=nchains,n_samples=nsamples,burn_in=4000,sample_freq=5)   
end


function run_C_F_ts_mean()

    RO_data = Load_RO_data()  
    samples = load("output/RO_C_F_20211124.jld","samples_container")    
    P_samples = load("output/RO_C_F_20211124.jld","P_samples_container")

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
    plot([2.0,10.5],[2.0,10.5],xlabel="data",ylabel="model")
    pl7 = plot!(GPP_data_F,GPP_model,seriestype=:scatter)
    
    GPP_R² = calcR²(GPP_data_F,GPP_model)   
    GPP_corr = cor(GPP_data_F,GPP_model)   

    H_model = ccphts_F.H[[21,44,64,84]]
    B_model = ccphts_F.B[[21,44,64,84]]
    Wf_model = ccphts_F.Wf[[21,44,64,84]]
    
    H_R² = calcR²(data_F.H,H_model)
    B_R² = calcR²(data_F.B,B_model)

    plot(2015:2018,data_F.H,xlabel="Date",ylabel="Height",label="Data F")
    pl2 = plot!(2015:2018,H_model,label="Model F")
    plot(2015:2018,data_F.B,xlabel="Date",ylabel="Basal Area",label="Data F")
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
    plot([2.0,10.5],[2.0,10.5],xlabel="data",ylabel="model")
    pl8 = plot!(GPP_data_C,GPP_model,seriestype=:scatter)

    GPP_R² = calcR²(GPP_data_C,GPP_model)  
    GPP_corr =  cor(GPP_data_C,GPP_model) 

    H_model = ccphts_C.H[[21,44,64,84]]
    B_model = ccphts_C.B[[21,44,64,84]]
    Wf_model = ccphts_C.Wf[[21,44,64,84]]

    H_R² = calcR²(data_C.H,H_model)
    B_R² = calcR²(data_C.B,B_model)

    plot(2015:2018,data_C.H,xlabel="Date",ylabel="Height",label="Data C")
    pl5 = plot!(2015:2018,H_model,label="Model C")
    plot(2015:2018,data_C.B,xlabel="Date",ylabel="Basal Area",label="Data C")
    pl6 = plot!(2015:2018,B_model,label="Model C")  

    println("Control: ")
    println("GPP R²:",GPP_R²) 
    println("GPP Pearson correlation: ",GPP_corr)
    println("H R²:",H_R²) 
    println("B R²:",B_R²)   
    clac_GPP_R²_annual(GPP_data_C,GPP_model,weatherts_C.date)    

    plot(pl1,pl2,pl3,pl7,pl4,pl5,pl6,pl8,layout=(2,4),legend=false)
    savefig("./plots/RO_C_F_ts_tuned_20211124.svg")    


    pl1 = plot(weatherts_F.date,ccphts_F.αr)
    pl1 = plot!(weatherts_C.date,ccphts_C.αr)
    pl2 = plot(weatherts_F.date,ccphts_F.K_cost)
    pl2 = plot!(weatherts_C.date,ccphts_C.K_cost)
    pl3 = plot(weatherts_F.date,ccphts_F.gₛ)
    pl3 = plot!(weatherts_C.date,ccphts_C.gₛ)
    pl4 = plot(weatherts_F.date,ccphts_F.Nₘ_f)
    pl4 = plot!(weatherts_C.date,ccphts_C.Nₘ_f)
    plot(pl1,pl2,pl3,pl4,layout=(2,2),legend=false)

    cons = CCPH.Constants()
    env = CCPH.EnvironmentStruct()
    treepar = CCPH.TreePar()
    kinetic = CCPH.PhotoKineticRates()
    LMA = 0.256
    P = 101325.0

    transpiration_lo = 1.0
    transpiration_hi = 4.0
    transp_F = CalcTranspiraiton(ccphts_F,weatherts_F,cons,LMA,P)
    transp_C = CalcTranspiraiton(ccphts_C,weatherts_C,cons,LMA,P)
    plot(weatherts_F.date,transp_F)
    plot!(weatherts_C.date,transp_C)

    RO_ET_F_est = RO_ET_Est(weatherts_F;stand_type="Fertilized") 
    RO_ET_C_est = RO_ET_Est(weatherts_C;stand_type="Control")

    plot(weatherts_F.date,RO_ET_F_est)
    plot!(weatherts_C.date,RO_ET_C_est)    

    RO_gₛ_F_est = EstgₛFromTranspiraiton(RO_ET_F_est,ccphts_F,weatherts_F,cons,LMA,P)
    RO_gₛ_C_est = EstgₛFromTranspiraiton(RO_ET_C_est,ccphts_C,weatherts_C,cons,LMA,P)

    plot(weatherts_F.date,RO_gₛ_F_est)
    plot!(weatherts_C.date,RO_gₛ_C_est)   
    
    RO_GPP_F_est = CalcGPP(RO_gₛ_F_est,ccphts_F,weatherts_F,cons,env,treepar,kinetic)
    RO_GPP_C_est = CalcGPP(RO_gₛ_C_est,ccphts_C,weatherts_C,cons,env,treepar,kinetic)

    plot(weatherts_F.date,data_F.GPP)
    pl1 = plot!(weatherts_F.date,RO_GPP_F_est)    
    plot(weatherts_C.date,data_C.GPP)
    pl2 = plot!(weatherts_C.date,RO_GPP_C_est)  
    
    plot(pl1,pl2,layout=(2,1),legend=false)
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

#run_simple_trait_model_C_F_tuning_RO_ts()
run_C_F_ts_mean()
#validate_RO_2019()