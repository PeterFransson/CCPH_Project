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
        
        GPP_model_F[j],EC_model_F[j], modeloutput, gₛ_opt, Nₘ_f_opt = CalcModelOutput!(j,model,weatherts_F,kinetic)
    end  

    #Control stand
    GPP_model_C = zeros(sim_steps+1)
    EC_model_C = zeros(sim_steps+1)
    for j = 1:sim_steps+1        
        
        model,kinetic = Initi_model_struct(j,data_C,αf,β₁,β₂,Nₛ_C,rₘ,a_Jmax,b_Jmax,i,Kₓₗ₀)
        
        GPP_model_C[j],EC_model_C[j], modeloutput, gₛ_opt, Nₘ_f_opt = CalcModelOutput!(j,model,weatherts_C,kinetic)
    end  

    return GPP_model_F,GPP_model_C,EC_model_F,EC_model_C
end

function Get_Data_RO_C_F_CCPH(par::Array{Float64,1},weatherts_F::WeatherTS,
    weatherts_C::WeatherTS,data_F::RoData,data_C::RoData;sim_steps::Integer=83)
    
    αf,β₁,β₂,b_Jmax = 460.0,1.27,-0.27,0.0
    Nₛ,rₘ,a_Jmax,Kₓₗ₀,Nₛ_C,i = par[1:6]    

    #Fertilized stand
    GPP_model_F = zeros(sim_steps+1)
    EC_model_F = zeros(sim_steps+1)  
    gₛ_model_F = zeros(sim_steps+1) 
    for j = 1:sim_steps+1           
       
        model,kinetic = Initi_model_struct(j,data_F,αf,β₁,β₂,Nₛ,rₘ,a_Jmax,b_Jmax,i,Kₓₗ₀)        
        
        GPP_model_F[j],EC_model_F[j], modeloutput, gₛ_model_F[j], Nₘ_f_opt = CalcModelOutput!(j,model,weatherts_F,kinetic)
    end  

    #Control stand
    GPP_model_C = zeros(sim_steps+1)
    EC_model_C = zeros(sim_steps+1)
    gₛ_model_C = zeros(sim_steps+1) 
    for j = 1:sim_steps+1        
        
        model,kinetic = Initi_model_struct(j,data_C,αf,β₁,β₂,Nₛ_C,rₘ,a_Jmax,b_Jmax,i,Kₓₗ₀)
        
        GPP_model_C[j],EC_model_C[j], modeloutput, gₛ_model_C[j], Nₘ_f_opt = CalcModelOutput!(j,model,weatherts_C,kinetic)
    end  

    return GPP_model_F,GPP_model_C,EC_model_F,EC_model_C,gₛ_model_F,gₛ_model_C 
end

function calc_error_term(ydata::T,ymodel::T) where {T<:Float64}
    e = (ydata-ymodel)^2

    return e
end

function calc_logP_term(ydata::T,ymodel::T,a::T,b::T) where {T<:Float64}
    e = abs(ydata-ymodel)
    σ = a+b*ymodel
    σ>0.0||error("Negative variance")
    return e/σ+log(σ)
end

function obj_logP(x::Array{T,1},GPP_model_F::Array{T,1},Ec_model_F::Array{T,1},data_F::RoData,
    GPP_model_C::Array{T,1},Ec_model_C::Array{T,1},data_C::RoData) where {T<:Float64}
    a_GPP,b_GPP,a_Ec,b_Ec = x

    logP = 0.0
    for i = 1:length(data_F.GPP)
        logP += calc_logP_term(data_F.GPP[i],GPP_model_F[i],a_GPP,b_GPP)
        logP += calc_logP_term(data_C.GPP[i],GPP_model_C[i],a_GPP,b_GPP)

        logP += calc_logP_term(data_F.E_C[i],Ec_model_F[i],a_Ec,b_Ec)
        logP += calc_logP_term(data_C.E_C[i],Ec_model_C[i],a_Ec,b_Ec)        
    end    
     
    return logP
end

function calc_error(
    GPP_model_F::Array{T,1},Ec_model_F::Array{T,1},data_F::RoData,
    GPP_model_C::Array{T,1},Ec_model_C::Array{T,1},data_C::RoData) where {T<:Float64}     

    e_tot = 0.0

    #=
    for i = 1:length(data_F.GPP)
        e_tot += calc_error_term(data_F.GPP[i],GPP_model_F[i])
        e_tot += calc_error_term(data_C.GPP[i],GPP_model_C[i])

        e_tot += calc_error_term(data_F.E_C[i],Ec_model_F[i])
        e_tot += calc_error_term(data_C.E_C[i],Ec_model_C[i])        
    end 
    =#   
    e_tot += calcR²(data_F.GPP,GPP_model_F)
    e_tot += calcR²(data_C.GPP,GPP_model_C)

    e_tot += calcR²(data_F.E_C,Ec_model_F)
    e_tot += calcR²(data_C.E_C,Ec_model_C) 
     
    return -e_tot
end

function obj_RO_C_F_CCPH(par::Array{Float64,1},RO_data::RO_raw_data)    
    try
        X0,τ,Smax,τ_C = par[7:10]
        weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
        weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ_C,Smax=Smax)           
        
        data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C,weatherts_F,weatherts_C)   

        GPP_model_F,GPP_model_C,EC_model_F,EC_model_C = Run_RO_C_F_CCPH(par,weatherts_F,weatherts_C,data_F,data_C)         

        e_tot = calc_error(
        GPP_model_F,EC_model_F,data_F,
        GPP_model_C,EC_model_C,data_C)
        
        return e_tot
    catch err
        println("Parameters Error: ", err)
        return Inf
    end    
end

function run_opt()
    #file_name = "opt_no_error_RO_C_GPP_Ec_20220223_2"
    file_name = "opt_no_error_RO_C_GPP_Ec_20220224"

    RO_data = Load_RO_data()
    
    #Nₛ,rₘ,a_Jmax,Kₓₗ₀,Nₛ_C,i,X0,τ,Smax,τ_C,a_GPP,b_GPP,a_Ec,b_Ec
    ranges = [(0.001,0.1),(10.0,40.0),(0.01,1.0),(0.001,0.1),(0.001,0.1),
    (0.1,6.0),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),(1.0,15.0)]

    res = BlackBoxOptim.bboptimize(x::Array{Float64,1}->obj_RO_C_F_CCPH(x::Array{Float64,1},RO_data)
    ; SearchRange = ranges,PopulationSize = 200,MaxSteps=30000)

    @show xopt = BlackBoxOptim.best_candidate(res)

    save("./output/"*file_name*".jld","xopt",xopt)
end
function run_opt_par()
    file_name = "opt_no_error_RO_C_GPP_Ec_20220224"

    par = load("./output/"*file_name*".jld","xopt")   

    #=
    [0.012854764847976128, 14.865444332753684, 0.26557235210247054,
    0.03784509575809752, 0.013224669959225405, 0.36827028167982595,
    -14.259069117859587, 8.754639146359056, 14.402015724819176, 1.2640782289025818] 
    =#

    RO_data = Load_RO_data()

    X0,τ,Smax,τ_C = par[7:10]
    weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax)
    weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ_C,Smax=Smax)
    
    data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C,weatherts_F,weatherts_C)    

    GPP_model_F,GPP_model_C,EC_model_F,EC_model_C,gₛ_model_F,gₛ_model_C = Get_Data_RO_C_F_CCPH(par,weatherts_F,weatherts_C,data_F,data_C)       

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

    plot(weatherts_F.date,data_F.GPP,label="Data",ylabel="GPP")   
    pl1 = plot!(weatherts_F.date,GPP_model_F,label="Model")
    plot(weatherts_C.date,data_C.GPP,label="Data",ylabel="GPP")
    pl2 = plot!(weatherts_C.date,GPP_model_C,label="Model")

    plot(weatherts_F.date,data_F.E_C,label="Data",ylabel="E_C")
    pl3 = plot!(weatherts_F.date,EC_model_F,label="Model")
    plot(weatherts_C.date,data_C.E_C,label="Data",ylabel="E_C")
    pl4 = plot!(weatherts_C.date,EC_model_C,label="Model")

    plot(pl1,pl3,pl2,pl4,layout=(2,2),legends=false)
    savefig("./plots/"*file_name*"_result.svg")
    
    plot(weatherts_F.date,data_F.gₛ,label="Data",ylabel="gₛ")
    pl1 = plot!(weatherts_F.date,gₛ_model_F,label="Model",ylabel="gₛ")
    plot(weatherts_C.date,data_C.gₛ,label="Data",ylabel="gₛ")
    pl2 = plot!(weatherts_C.date,gₛ_model_C,label="Model",ylabel="gₛ")

    plot(pl1,pl2,layout=(1,2),legends=false)
    savefig("./plots/"*file_name*"_result_gs.svg")    

    #Test optimizaiton
    inner_optimizer = Optim.BFGS(linesearch=Optim.LineSearches.BackTracking())   
    
    lower = [0.0001,0.0001,0.0001,0.0001]
    upper = [5.0,3.0,5.0,3.0]
    x0 = [0.237,0.0845,0.0688,0.146]

    opt_trait = Optim.optimize(x->obj_logP(x,GPP_model_F,EC_model_F,data_F,
    GPP_model_C,EC_model_C,data_C),
    lower,upper,x0,Optim.Fminbox(inner_optimizer))

    opt_trait.minimizer
end

#run_opt()
run_opt_par()
