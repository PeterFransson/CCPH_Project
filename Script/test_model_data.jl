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
    return e/σ+log(2*σ)
end

function Calc_logP(
    GPP_model_F::Array{T,1},Ec_model_F::Array{T,1},GPP_data_F::Array{T,1},Ec_data_F::Array{T,1},
    GPP_model_C::Array{T,1},Ec_model_C::Array{T,1},GPP_data_C::Array{T,1},Ec_data_C::Array{T,1},
    a_GPP::T,b_GPP::T,a_Ec::T,b_Ec::T) where {T<:Float64}     

    logP = 0.0
    for i = 1:length(GPP_data_F)
        logP += calc_logP_term(GPP_data_F[i],GPP_model_F[i],a_GPP,b_GPP)
        logP += calc_logP_term(GPP_data_C[i],GPP_model_C[i],a_GPP,b_GPP)

        logP += calc_logP_term(Ec_data_F[i],Ec_model_F[i],a_Ec,b_Ec)
        logP += calc_logP_term(Ec_data_C[i],Ec_model_C[i],a_Ec,b_Ec)        
    end    
     
    return  -logP
end

function Post_distri_RO_C_F_CCPH(par::Array{Float64,1},RO_data::RO_raw_data,
    GPP_data_F_in::Array{Float64,1},Ec_data_F_in::Array{Float64,1},
    GPP_data_C_in::Array{Float64,1},Ec_data_C_in::Array{Float64,1})
    
        try
            X0,τ,Smax,τ_C = par[7:10]
            weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
            weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ_C,Smax=Smax)           
            
            data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C,weatherts_F,weatherts_C)   

            GPP_model_F,GPP_model_C,EC_model_F,EC_model_C = Run_RO_C_F_CCPH(par,weatherts_F,weatherts_C,data_F,data_C)              

            a_GPP,b_GPP,a_Ec,b_Ec = par[11:14]

            logP = Calc_logP(
            GPP_model_F,EC_model_F,GPP_data_F_in,Ec_data_F_in,
            GPP_model_C,EC_model_C,GPP_data_C_in,Ec_data_C_in,
            a_GPP,b_GPP,a_Ec,b_Ec)
            
            return post = logP
            
        catch err
            println("Parameters Error: ", err)
            return post = -Inf
        end     
end

function getmodeldata(par::Array{Float64,1},RO_data::RO_raw_data)
    X0,τ,Smax,τ_C = par[7:10]
    weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
    weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ_C,Smax=Smax)           
            
    data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C,weatherts_F,weatherts_C)   

    GPP_model_F,GPP_model_C,EC_model_F,EC_model_C = Run_RO_C_F_CCPH(par,weatherts_F,weatherts_C,data_F,data_C)  

    return GPP_model_F,GPP_model_C,EC_model_F,EC_model_C,weatherts_F,weatherts_C
end

function generatedatafromparameter(par::Array{Float64,1},RO_data::RO_raw_data)

    GPP_model_F,GPP_model_C,EC_model_F,EC_model_C,weatherts_F,weatherts_C = getmodeldata(par,RO_data)          

    a_GPP,b_GPP,a_Ec,b_Ec = par[11:14]

    error_GPP_F = []
    error_GPP_C = []
    error_Ec_F = []
    error_Ec_C = []

    for i = 1:length(GPP_model_F)
        push!(error_GPP_F,rand(Laplace(0.0, a_GPP+b_GPP*GPP_model_F[i])))
        push!(error_GPP_C,rand(Laplace(0.0, a_GPP+b_GPP*GPP_model_C[i])))
        push!(error_Ec_F,rand(Laplace(0.0, a_Ec+b_Ec*EC_model_F[i])))
        push!(error_Ec_C,rand(Laplace(0.0, a_Ec+b_Ec*EC_model_C[i])))
    end

    return GPP_model_F,GPP_model_C,EC_model_F,EC_model_C,error_GPP_F,error_GPP_C,error_Ec_F,error_Ec_C,weatherts_F,weatherts_C
end

function run_test()
    RO_data = Load_RO_data()

    par = [0.011976551772617215, 25.42466609569836, 0.5009242298864488, 0.023503724920524185,
    0.010452800544803813,1.0,-1.312709470268466, 3.3784814235281933, 17.510523050796454, 9.566445194647766,
    0.237,0.0845,0.0688,0.146]

    GPP_model_F,GPP_model_C,EC_model_F,EC_model_C,
    error_GPP_F,error_GPP_C,error_Ec_F,error_Ec_C,weatherts_F,weatherts_C = generatedatafromparameter(par,RO_data)

    GPP_data_F = GPP_model_F.+error_GPP_F
    GPP_data_C = GPP_model_C.+error_GPP_C
    EC_data_F = EC_model_F.+error_Ec_F
    EC_data_C = EC_model_C.+error_Ec_C
    
    plot(weatherts_F.date,GPP_model_F)
    pl1 = plot!(weatherts_F.date,GPP_model_F.+error_GPP_F)
    plot(weatherts_F.date,EC_model_F)
    pl2 = plot!(weatherts_F.date,EC_model_F.+error_Ec_F)

    plot(weatherts_C.date,GPP_model_C)
    pl3 = plot!(weatherts_C.date,GPP_model_C.+error_GPP_C)
    plot(weatherts_C.date,EC_model_C)
    pl4 = plot!(weatherts_C.date,EC_model_C.+error_Ec_C)

    @show GPP_R²_F =  calcR²(GPP_model_F.+error_GPP_F,GPP_model_F)
    @show GPP_cor_F = cor(GPP_model_F.+error_GPP_F,GPP_model_F)

    @show Ec_R²_F =  calcR²(EC_model_F.+error_Ec_F,EC_model_F)
    @show Ec_cor_F = cor(EC_model_F.+error_Ec_F,EC_model_F)

    @show GPP_R²_C =  calcR²(GPP_model_C.+error_GPP_C,GPP_model_C)
    @show GPP_cor_C = cor(GPP_model_C.+error_GPP_C,GPP_model_C)

    @show Ec_R²_C =  calcR²(EC_model_C.+error_Ec_C,EC_model_C)
    @show Ec_cor_C = cor(EC_model_C.+error_Ec_C,EC_model_C)

    plot(pl1,pl2,pl3,pl4,layout=(2,2),legends=false)    
    
    #Nₛ,rₘ,a_Jmax,Kₓₗ₀,Nₛ_C,i,X0,τ,Smax,τ_C,a_GPP,b_GPP,a_Ec,b_Ec
    ranges = [(0.001,0.1),(10.0,40.0),(0.01,1.0),(0.001,0.1),(0.001,0.1),
    (0.1,6.0),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),(1.0,15.0),
    (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]

    res = BlackBoxOptim.bboptimize(x::Array{Float64,1}->-Post_distri_RO_C_F_CCPH(x::Array{Float64,1},RO_data,
    GPP_data_F,EC_data_F,GPP_data_C,EC_data_C)
    ; SearchRange = ranges,PopulationSize = 50,MaxSteps=2000)

    @show xopt = BlackBoxOptim.best_candidate(res) 

    @show -Post_distri_RO_C_F_CCPH(par,RO_data,GPP_data_F,EC_data_F,GPP_data_C,EC_data_C)

    GPP_model_F,GPP_model_C,EC_model_F,EC_model_C,weatherts_F,weatherts_C = getmodeldata(xopt,RO_data)  

    @show GPP_R²_F =  calcR²(GPP_data_F,GPP_model_F)
    @show GPP_cor_F = cor(GPP_data_F,GPP_model_F)
    @show gₛ_R²_F =  calcR²(EC_data_F,EC_model_F)
    @show gₛ_cor_F = cor(EC_data_F,EC_model_F)

    @show GPP_R²_C =  calcR²(GPP_data_C,GPP_model_C)
    @show GPP_cor_C = cor(GPP_data_C,GPP_model_C)
    @show gₛ_R²_C =  calcR²(EC_data_C,EC_model_C)
    @show gₛ_cor_C = cor(EC_data_C,EC_model_C)    

    a_GPP,b_GPP,a_Ec,b_Ec = xopt[11:14]
    
    σ_GPP_F = a_GPP.+GPP_model_F*b_GPP
    σ_GPP_C = a_GPP.+GPP_model_C*b_GPP

    σ_EC_F = a_Ec.+b_Ec*EC_model_F
    σ_EC_C = a_Ec.+b_Ec*EC_model_C

    c = -log(0.025*2) #Use for calculating 95% credible intervals

    plot(weatherts_F.date,GPP_data_F,label="Data",ylabel="GPP")
    plot!(weatherts_F.date,GPP_model_F+c*σ_GPP_F)
    plot!(weatherts_F.date,GPP_model_F-c*σ_GPP_F)
    pl1 = plot!(weatherts_F.date,GPP_model_F,label="Model")
    plot(weatherts_C.date,GPP_data_C,label="Data",ylabel="GPP")
    plot!(weatherts_C.date,GPP_model_C+c*σ_GPP_C)
    plot!(weatherts_C.date,GPP_model_C-c*σ_GPP_C)
    pl2 = plot!(weatherts_C.date,GPP_model_C,label="Model")

    plot(weatherts_F.date,EC_data_F,label="Data",ylabel="E_C")
    plot!(weatherts_F.date,EC_model_F+c*σ_EC_F)
    plot!(weatherts_F.date,EC_model_F-c*σ_EC_F)
    pl3 = plot!(weatherts_F.date,EC_model_F,label="Model")
    plot(weatherts_C.date,EC_data_C,label="Data",ylabel="E_C")
    plot!(weatherts_C.date,EC_model_C+c*σ_EC_C)
    plot!(weatherts_C.date,EC_model_C-c*σ_EC_C)
    pl4 = plot!(weatherts_C.date,EC_model_C,label="Model")

    plot(pl1,pl3,pl2,pl4,layout=(2,2),legends=false)       
end

run_test()