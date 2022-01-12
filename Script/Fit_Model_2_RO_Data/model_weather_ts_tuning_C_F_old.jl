mutable struct RoData{T<:Float64}
    H::Array{T,1}
    B::Array{T,1}
    GPP::Array{T,1}
end

function Create_RoData_C_F(GPP_data_F::Array{T,1},GPP_data_C::Array{T,1}) where {T<:Float64}   

    H_data_F = [17.004166666666666, 17.271666666666665, 17.59166666666667, 17.835] 
    B_data_F = [0.025076878666666667, 0.02603588483333334, 0.026922637833333332, 0.027666997]

    H_data_C = [16.956,17.165,17.388,17.590] 
    B_data_C = [0.023015244833333334, 0.02349805816666667, 0.023938368499999998, 0.02432245283333333]

    data_F = RoData(H_data_F,B_data_F,GPP_data_F)
    data_C = RoData(H_data_C,B_data_C,GPP_data_C)

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
    Wf = 4.2
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

function Post_distri_RO_C_F_CCPH(par::Array{Float64,1})
    if Check_para_RO_C_F_CCPH(par)
        return post = 0.0
    else
        try
            X0,τ,Smax = par[11:13]
            weatherts_F,GPP_data_F = create_weather_struct_RO(;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
            weatherts_C,GPP_data_C = create_weather_struct_RO(;stand_type="Control",X0=X0,τ=τ,Smax=Smax)
            
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

    nsamples = 3000
    nchains = 12



    par_guess = [470.328512441229, 1.2501809959672654, -0.2705395150988396,
    0.036975489745880365, 17.40739827323768, 0.015721074668478767/0.24321499065796873,
    3.4357734812576556e-6/0.24321499065796873, 1.0749476181390065, 8.621674950871828,
    0.025675473753997178, -4.230533359595931, 12.328105896346644, 19.079489364590533]

    x_current = [par_guess for i = 1:nchains]
    P_current = [Post_distri_RO_C_F_CCPH(par_guess) for i = 1:nchains] 
        
    println(P_current)  

    q_vec = [460.0,1.27,-0.27,0.04,24.0,0.033,1.21e-5,1.0,1.0,0.04,-4.0,7.0,16.0]/25.0
    metropolis_mcmc(x_current,P_current,Post_distri_RO_C_F_CCPH,
    q_vec,"RO_C_F_20211105";
    n_chains=nchains,n_samples=nsamples,burn_in=3000,sample_freq=5) 

    #=
    samples,accept_rate = metropolis_mcmc!(x_current,P_current,Post_distri_RO_C_F_CCPH,q_distri_RO_C_F_CCPH;
    n_chains=nchains,n_samples=nsamples,
    burn_in=300,sample_freq=5)
    =#
    
    #=
    metropolis_mcmc(x_current,P_current, 
    Post_distri_RO_C_F_CCPH,q_distri_RO_C_F_CCPH,"RO_C_F_20211031";
    burn_in=1000,n_samples = nsamples,n_chains = nchains,sample_freq = 5)  
    =#   

    #=
    @show accept_rate
    @show R_hat = potscalereduc(samples)    

    save("./output/mcmc_CCPH_RO_C_F_20211031.jld","samples",samples)

    samples_reshape = reshape_sample(samples)

    histogram(samples_reshape,layout=length(par_guess),legend=false)
    savefig("./plots/mcmc_CCPH_RO_C_F_20211031.svg")     
    =#                
end

#=
function run_C_F_ts_mean()

    #samples = load("./output/mcmc_simple_CCPH_RO_C_F_ts_20210831.jld","samples")
    samples = load("./output/mcmc_simple_CCPH_RO_C_F_ts_20211018.jld","samples")

    samples_reshape = reshape_sample(samples)

    par = [mean(samples_reshape[:,i]) for i = 1:size(samples_reshape,2)]
    
    #par = samples_reshape[9253,:] #C_F_ts_tuned_20210830 opt GPP + Size
    #par = samples_reshape[24587,:] #C_F_ts_tuned_20210831 opt GPP + Size
    par = samples_reshape[25088,:] #C_F_ts_tuned_20211018 opt GPP + Size

    println(par)

    X0,τ,Smax = par[11:13]
    weatherts_F,GPP_data_F,tot_daylight_hour_ts, t_end_vec = create_weather_struct_RO(;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
    weatherts_C,GPP_data_C,tot_daylight_hour_ts, t_end_vec = create_weather_struct_RO(;stand_type="Control",X0=X0,τ=τ,Smax=Smax)

    CCPHTS_F,CCPHTS_C = Run_RO_C_F_Simple_CCPH_ts(par,weatherts_F,weatherts_C,tot_daylight_hour_ts)  

    #Fertilized stand
    GPP_model = CCPHTS_F.N.*CCPHTS_F.P.*weatherts_F.daylight./tot_daylight_hour_ts*1000/7

    plot(weatherts_F.date,GPP_data_F,label="Data F",ylabel="g C m⁻² d⁻¹")
    pl1 = plot!(weatherts_F.date,GPP_model,label="Model F")
    
    GPP_R² = calcR²(GPP_data_F,GPP_model)   
    GPP_corr =  cor(GPP_data_F,GPP_model)
    
    H_data = [17.004166666666666, 17.271666666666665, 17.59166666666667, 17.835] 
    B_data = [0.025076878666666667, 0.02603588483333334, 0.026922637833333332, 0.027666997]

    H_model = CCPHTS_F.H[[21,44,64,84]]
    B_model = CCPHTS_F.B[[21,44,64,84]]
    
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
    GPP_model = CCPHTS_C.N.*CCPHTS_C.P.*weatherts_C.daylight./tot_daylight_hour_ts*1000/7

    plot(weatherts_C.date,GPP_data_C,label="Data C",ylabel="g C m⁻² d⁻¹")
    pl4 = plot!(weatherts_C.date,GPP_model,label="Model C")

    GPP_R² = calcR²(GPP_data_C,GPP_model)  
    GPP_corr =  cor(GPP_data_C,GPP_model) 

    H_data = [16.956,17.165,17.388,17.590] 
    B_data = [0.023015244833333334, 0.02349805816666667, 0.023938368499999998, 0.02432245283333333]

    H_model = CCPHTS_C.H[[21,44,64,84]]
    B_model = CCPHTS_C.B[[21,44,64,84]]

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
    savefig("./plots/simple_CCPH_RO_C_F_ts_tuned_20211018.svg")
end

function sim_RO_calc_error(par::Array{Float64, 1})
    X0,τ,Smax = par[11:13]
    weatherts_F,GPP_data_F,tot_daylight_hour_ts, t_end_vec = create_weather_struct_RO(;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
    weatherts_C,GPP_data_C,tot_daylight_hour_ts, t_end_vec = create_weather_struct_RO(;stand_type="Control",X0=X0,τ=τ,Smax=Smax)
    CCPHTS_F,CCPHTS_C = Run_RO_C_F_Simple_CCPH_ts(par,weatherts_F,weatherts_C,tot_daylight_hour_ts) 

    #Fertilized stand
    GPP_F_model = CCPHTS_F.N.*CCPHTS_F.P.*weatherts_F.daylight./tot_daylight_hour_ts*1000/7
    GPP_F_R² = calcR²(GPP_data_F,GPP_F_model)  

    #Control stand
    GPP_C_model = CCPHTS_C.N.*CCPHTS_C.P.*weatherts_C.daylight./tot_daylight_hour_ts*1000/7
    GPP_C_R² = calcR²(GPP_data_C,GPP_C_model) 

    return calc_error_RO_C_F_Simple_CCPH_ts(CCPHTS_F,GPP_data_F,CCPHTS_C,
    GPP_data_C,weatherts_F,weatherts_C,tot_daylight_hour_ts), min(GPP_F_R²,GPP_C_R²)
end

function error_RO_C_F_dist()
    samples = load("./output/mcmc_simple_CCPH_RO_C_F_ts_20211018.jld","samples")

    samples_reshape = reshape_sample(samples)

    n_samples = size(samples_reshape,1)

    error_vec = zeros(n_samples)
    GPP_filter = zeros(n_samples)    

    Threads.@threads for i = 1:n_samples    
        println("Sample Nr: $(i)/$(n_samples)")
        par = samples_reshape[i,:]

        size_error,GPP_R² = sim_RO_calc_error(par)
        size_R² = 1-size_error

        error_vec[i] = size_error

        GPP_filter[i] = GPP_R²
        
        #=
        if size_R²>0.8
            GPP_filter[i] = GPP_R²
        else
            GPP_filter[i] = 0.0
        end   
        =#     
    end

    @show 1-minimum(error_vec)
    @show min_ind = argmin(error_vec)

    @show 1-error_vec[min_ind]

    println("Minimizer:")
    par = samples_reshape[min_ind,:]
    println(par)

    @show maximum(GPP_filter)
    @show max_GPP_ind = argmax(GPP_filter)

    @show GPP_filter[max_GPP_ind]

    println("Maximizer GPP:")
    par = samples_reshape[max_GPP_ind,:]
    println(par) 

    R² = -error_vec.+1
    histogram(R²)
end

function validate_RO_2019()
    samples = load("./output/mcmc_simple_CCPH_RO_C_F_ts_20211018.jld","samples")

    samples_reshape = reshape_sample(samples)  
    #par = samples_reshape[24587,:] #date: 20210831
    par = samples_reshape[25088,:] #date: 20211018
    println(par)

    X0,τ,Smax = par[11:13]
    weatherts_F,GPP_data_F,tot_daylight_hour_ts, t_end_vec = create_weather_struct_RO(;stand_type="Fertilized",
    X0=X0,τ=τ,Smax=Smax,weather_file = "Weather_RO_2")
    weatherts_C,GPP_data_C,tot_daylight_hour_ts, t_end_vec = create_weather_struct_RO(;stand_type="Control",
    X0=X0,τ=τ,Smax=Smax,weather_file = "Weather_RO_2")

    @show length(GPP_data_F)

    CCPHTS_F,CCPHTS_C = Run_RO_C_F_Simple_CCPH_ts(par,weatherts_F,weatherts_C,tot_daylight_hour_ts;sim_steps=104)

    #Fertilized stand
    GPP_model = CCPHTS_F.N.*CCPHTS_F.P.*weatherts_F.daylight./tot_daylight_hour_ts*1000/7    

    plot(weatherts_F.date[85:end-1],GPP_data_F[85:end-1],color=1,markershape=:circle,label="Data F",ylabel="g C m⁻² d⁻¹")
    plot!(weatherts_F.date[1:84],GPP_data_F[1:84],color=1)
    plot!(weatherts_F.date[1:84],GPP_model[1:84],color=2)
    pl1 = plot!(weatherts_F.date[85:end-1],GPP_model[85:end-1],label="Model F",color=2,markershape=:circle)

    plot(GPP_data_F[85:end],GPP_model[85:end],color=1,seriestype=:scatter,xlabel="Data F",ylabel="Model F")
    plot!([1,10],[1,10],color=3)
    pl2 = plot!(GPP_data_F[1:84],GPP_model[1:84],color=2,seriestype=:scatter)
    
    GPP_R² = calcR²(GPP_data_F[85:end-1],GPP_model[85:end-1])   
    GPP_corr = cor(GPP_data_F[85:end-1],GPP_model[85:end-1])

    println("Fertilized: ")
    println("GPP R²:",GPP_R²) 
    println("GPP Pearson correlation: ",GPP_corr)

    #Control stand
    GPP_model = CCPHTS_C.N.*CCPHTS_C.P.*weatherts_C.daylight./tot_daylight_hour_ts*1000/7

    plot(weatherts_C.date[85:end-1],GPP_data_C[85:end-1],color=1,markershape=:circle,label="Data C",ylabel="g C m⁻² d⁻¹")
    plot!(weatherts_C.date[1:84],GPP_data_C[1:84],color=1)
    plot!(weatherts_C.date[1:84],GPP_model[1:84],color=2)
    pl3 = plot!(weatherts_C.date[85:end-1],GPP_model[85:end-1],label="Model C",color=2,markershape=:circle)

    plot(GPP_data_C[85:end],GPP_model[85:end],color=1,seriestype=:scatter,xlabel="Data C",ylabel="Model C") 
    plot!([1,10],[1,10],color=3)
    pl4 = plot!(GPP_data_C[1:84],GPP_model[1:84],color=2,seriestype=:scatter)


    GPP_R² = calcR²(GPP_data_C[85:end-1],GPP_model[85:end-1])  
    GPP_corr = cor(GPP_data_C[85:end-1],GPP_model[85:end-1]) 

    println("Control: ")
    println("GPP R²:",GPP_R²) 
    println("GPP Pearson correlation: ",GPP_corr)

    plot(pl1,pl2,pl3,pl4,layout=4,legend=false)
    savefig("./plots/simple_CCPH_RO_C_F_ts_validate_20211018.svg")

    plot(weatherts_F.date,CCPHTS_F.K_cost,seriestype=:scatter,ylabel="K_cost")
    pl1 = plot!(weatherts_C.date,CCPHTS_C.K_cost,seriestype=:scatter)
    plot(weatherts_F.date,CCPHTS_F.gₛ,seriestype=:scatter,ylabel="gₛ")
    pl2 = plot!(weatherts_C.date,CCPHTS_C.gₛ,seriestype=:scatter)
    plot(weatherts_F.date,CCPHTS_F.Nₘ_f,seriestype=:scatter,ylabel="Nₘ_f")
    pl3 = plot!(weatherts_C.date,CCPHTS_C.Nₘ_f,seriestype=:scatter)
    plot(weatherts_F.date,weatherts_F.θₛ,seriestype=:scatter,ylabel="θₛ")
    pl4 = plot!(weatherts_C.date,weatherts_C.θₛ,seriestype=:scatter)
    plot(weatherts_F.date,CCPHTS_F.Wf.*CCPHTS_F.N/0.256,seriestype=:scatter,ylabel="LAI")
    pl5 = plot!(weatherts_C.date,CCPHTS_C.Wf.*CCPHTS_C.N/0.256,seriestype=:scatter)
    plot(weatherts_F.date,CCPHTS_F.αr,seriestype=:scatter,ylabel="αᵣ")
    pl6 = plot!(weatherts_C.date,CCPHTS_C.αr,seriestype=:scatter)
    plot(pl1,pl2,pl3,pl4,pl5,pl6,layout=6,legend=false)
    savefig("./plots/simple_CCPH_RO_C_F_ts_validate_size_20211018.svg")
end

function clac_GPP_R²_annual(GPP::Array{Float64,1},GPP_model::Array{Float64,1},dates::Array{DateTime,1})
    years = [year(dates[1])]
    growth_end = Int64[]
    for i = 2:length(dates)
        if year(dates[i]) != year(dates[i-1])
            push!(growth_end,i-1)
            push!(years,year(dates[i]))
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
=#