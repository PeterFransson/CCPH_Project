function Initi_model_struct(H::T,Hc::T,N::T,Wf::T,B::T,αf::T,β₁::T,β₂::T,
    Nₛ::T,rₘ::T,a_Jmax::T,b_Jmax::T,i::T,Kₓₗ₀::T,r_α::T) where {T<:Float64}
    cons = Constants()
    env = EnvironmentStruct() 
    kinetic = PhotoKineticRates()
    photo = PhotoPar(kinetic,env.Tₐ)    
    hydPar = HydraulicsPar(;Kₓₗ₀=Kₓₗ₀,i=i)
    treepar = TreePar(;αf=αf,β₁=β₁,β₂=β₂,rₘ=rₘ,Nₛ=Nₛ,a_Jmax=a_Jmax,b_Jmax=b_Jmax,r_α=r_α)   
       
    Hs = H-Hc   
    As = Wf/treepar.αf
    Ww = treepar.ρw*As*(treepar.β₁*H+treepar.β₂*Hs)    
    treesize = TreeSize(Wf,Ww,H,Hs,As,B,N)

    model = CCPHStruct(cons,env,treepar,treesize,photo,hydPar)

    return model,kinetic
end

function Initi_model_struct(j::Integer,data::RoData,
    αf::T,β₁::T,β₂::T,Nₛ::T,rₘ::T,a_Jmax::T,b_Jmax::T,i::T,Kₓₗ₀::T,r_α::T) where {T<:Float64}
    data_ind = Find_data_ind(j)
    H = data.H[data_ind]
    N = data.N[data_ind]
    Wf = data.Wf[data_ind]
    B = data.B[data_ind]
    Hc = data.Hc
    model,kinetic = Initi_model_struct(H,Hc,N,Wf,B,αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₓₗ₀,r_α)   

    return model,kinetic
end

function Run_RO_CCPH(Nₛ::T,rₘ::T,a_Jmax::T,Kₓₗ₀::T,i::T,r_α::T,
    weatherts::WeatherTS,data::RoData;sim_steps::Integer=83) where {T<:Float64}
    αf,β₁,β₂,b_Jmax = 460.0,1.27,-0.27,0.0
     
    #Fertilized stand
    GPP_model = zeros(sim_steps+1)
    EC_model = zeros(sim_steps+1)  
    Nₘ_f_model = zeros(sim_steps+1)
    gₛ_model = zeros(sim_steps+1)
    for j = 1:sim_steps+1   
        model,kinetic = Initi_model_struct(j,data,αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₓₗ₀,r_α)
        
        GPP_model[j],EC_model[j], modeloutput, gₛ_model[j], Nₘ_f_model[j] = CalcModelOutput!(j,model,weatherts,kinetic)
    end 
    
    return GPP_model,EC_model,gₛ_model,Nₘ_f_model 
end


function Run_RO_CCPH(par::Array{Float64,1},weatherts::WeatherTS,data::RoData;sim_steps::Integer=83)    
    
    Nₛ,rₘ,a_Jmax,Kₓₗ₀,i,r_α = par[1:6]  

    GPP_model,EC_model,gₛ_model,Nₘ_f_model = Run_RO_CCPH(Nₛ,rₘ,a_Jmax,Kₓₗ₀,i,r_α,weatherts,data;sim_steps=sim_steps)

    return GPP_model,EC_model,gₛ_model,Nₘ_f_model
end

function calc_logP_term(ydata::T,ymodel::T,a::T,b::T) where {T<:Float64}
    e = abs(ydata-ymodel)
    σ = a+b*ymodel
    σ>0.0||error("Negative variance")
    return e/σ+log(σ)
end

function Calc_logP(
    GPP_model::Array{T,1},Ec_model::Array{T,1},Nₘ_f_model::Array{T,1},data::RoData,
    a_GPP::T,b_GPP::T,a_Ec::T,b_Ec::T,μ::T,b::T) where {T<:Float64}     

    logP = 0.0
    for i = 1:length(data.GPP)
        logP += calc_logP_term(data.GPP[i],GPP_model[i],a_GPP,b_GPP)
        logP += calc_logP_term(data.E_C[i],Ec_model[i],a_Ec,b_Ec)        
        #logP += abs(Nₘ_f_model[i]-μ)/b               
    end    
     
    return  -logP
end

function Get_Result_RO_CCPH(par::Array{Float64,1},RO_data::RO_raw_data;stand_type="Fertilized")
    X0,τ,Smax = par[7:9]

    weatherts,GPP_data = create_weather_struct_RO(RO_data;stand_type=stand_type,X0=X0,τ=τ,Smax=Smax)

    data = Create_RoData_F(GPP_data,weatherts)

    GPP_model,EC_model,gₛ_model,Nₘ_f_model = Run_RO_CCPH(par,weatherts,data)

    return GPP_model,EC_model,gₛ_model,Nₘ_f_model,weatherts,data
end

function Post_distri_RO_F_CCPH(par::Array{Float64,1},RO_data::RO_raw_data)
    try 
        GPP_model_F,EC_model_F,gₛ_model_F,Nₘ_f_model_F,weatherts_F,data_F = Get_Result_RO_CCPH(par,RO_data;stand_type="Fertilized")

        a_GPP,b_GPP,a_Ec,b_Ec = par[10:13]        

        logP = Calc_logP(GPP_model_F,EC_model_F,Nₘ_f_model_F,data_F,a_GPP,b_GPP,a_Ec,b_Ec,0.018,0.0018)
    catch err
        println("Parameters Error: ", err)
        return post = -Inf
    end           
end

function Post_distri_RO_C_CCPH(par::Array{Float64,1},RO_data::RO_raw_data)    
    try
        GPP_model_C,EC_model_C,gₛ_model_C,Nₘ_f_model_C,weatherts_C,data_C = Get_Result_RO_CCPH(par,RO_data;stand_type="Control")
        
        a_GPP,b_GPP,a_Ec,b_Ec = par[10:13]

        logP = Calc_logP(GPP_model_C,EC_model_C,Nₘ_f_model_C,data_C,a_GPP,b_GPP,a_Ec,b_Ec,0.0113,0.00063)
        
        return post = logP
    catch err
        println("Parameters Error: ", err)
        return post = -Inf
    end    
end

function run_opt()    
    @show file_name_F = "opt_RO_F_20220310"
    @show file_name_C = "opt_RO_C_20220310"

    RO_data = Load_RO_data()
    
    #Nₛ,rₘ,a_Jmax,Kₓₗ₀,i,r_α,X0,τ,Smax,a_GPP,b_GPP,a_Ec,b_Ec
    ranges = [(0.001,0.1),(10.0,40.0),(0.01,1.0),(0.001,0.1),
    (0.1,6.0),(1000,7000),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),
    (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]

    #--Fertilized--
    res = BlackBoxOptim.bboptimize(x::Array{Float64,1}->-Post_distri_RO_F_CCPH(x::Array{Float64,1},RO_data)
    ; SearchRange = ranges,PopulationSize = 50,MaxSteps=10000)

    @show xopt = BlackBoxOptim.best_candidate(res)

    save("./output/"*file_name_F*".jld","xopt",xopt)

    #--Control--
    res = BlackBoxOptim.bboptimize(x::Array{Float64,1}->-Post_distri_RO_C_CCPH(x::Array{Float64,1},RO_data)
    ; SearchRange = ranges,PopulationSize = 50,MaxSteps=10000)

    @show xopt = BlackBoxOptim.best_candidate(res)

    save("./output/"*file_name_C*".jld","xopt",xopt)
end

function run_opt_par()
    file_name_F = "opt_RO_F_20220310"
    file_name_C = "opt_RO_C_20220310"
    file_name = "opt_RO_F_C_20220310"

    par_F = load("./output/"*file_name_F*".jld","xopt")    
    par_C = load("./output/"*file_name_C*".jld","xopt")   
   
    println(file_name_F) 
    @show par_F
    println(file_name_C) 
    @show par_C

    RO_data = Load_RO_data()    

    @show Post_distri_RO_F_CCPH(par_F,RO_data)
    @show Post_distri_RO_C_CCPH(par_C,RO_data)    

    GPP_model_F,EC_model_F,gₛ_model_F,Nₘ_f_model_F,weatherts_F,data_F = Get_Result_RO_CCPH(par_F,RO_data;stand_type="Fertilized")
    GPP_model_C,EC_model_C,gₛ_model_C,Nₘ_f_model_C,weatherts_C,data_C = Get_Result_RO_CCPH(par_C,RO_data;stand_type="Control")

    println("Fertilized")
    @show GPP_R²_F =  calcR²(data_F.GPP,GPP_model_F)
    @show GPP_cor_F = cor(data_F.GPP,GPP_model_F)
    @show EC_R²_F =  calcR²(data_F.E_C,EC_model_F)
    @show EC_cor_F = cor(data_F.E_C,EC_model_F)
    println("Control")
    @show GPP_R²_C =  calcR²(data_C.GPP,GPP_model_C)
    @show GPP_cor_C = cor(data_C.GPP,GPP_model_C)
    @show EC_R²_C =  calcR²(data_C.E_C,EC_model_C)
    @show EC_cor_C = cor(data_C.E_C,EC_model_C)
    println("Fertilized")
    clac_GPP_R²_annual(data_C.GPP,GPP_model_F,weatherts_F.date)
    println("Control")  
    clac_GPP_R²_annual(data_C.GPP,GPP_model_C,weatherts_C.date)

    a_GPP_F,b_GPP_F,a_Ec_F,b_Ec_F = par_F[10:13]
    a_GPP_C,b_GPP_C,a_Ec_C,b_Ec_C = par_C[10:13]

    σ_GPP_F = a_GPP_F.+GPP_model_F*b_GPP_F
    σ_GPP_C = a_GPP_C.+GPP_model_C*b_GPP_C

    σ_EC_F = a_Ec_F.+b_Ec_F*EC_model_F
    σ_EC_C = a_Ec_C.+b_Ec_C*EC_model_C

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

    plot(pl1,pl2,pl3,pl4,layout=(2,2),legends=false)
    savefig("./plots/"*file_name*"_result.svg")

    plot(weatherts_F.date,data_F.gₛ,label="Data",ylabel="gₛ")
    pl1 = plot!(weatherts_F.date,gₛ_model_F,label="Model",ylabel="gₛ")
    plot(weatherts_C.date,data_C.gₛ,label="Data",ylabel="gₛ")
    pl2 = plot!(weatherts_C.date,gₛ_model_C,label="Model",ylabel="gₛ")

    plot(pl1,pl2,layout=(1,2),legends=false)
    savefig("./plots/"*file_name*"_result_gs.svg")

    pl1 = plot(weatherts_F.date,Nₘ_f_model_F,label="Model",ylabel="Nₘ_f (%)")    
    pl2 = plot(weatherts_C.date,Nₘ_f_model_C,label="Model",ylabel="Nₘ_f (%)")
    
    plot(pl1,pl2,layout=(1,2),legends=false)
    savefig("./plots/"*file_name*"_result_Nm_f.svg")
end

run_opt()
run_opt_par()