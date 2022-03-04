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

function run_validation_RO_2019()
    file_name = "adaptive_rwm_RO_C_GPP_Ec_20220228"
    #Nₛ,rₘ,a_Jmax,Kₓₗ₀,Nₛ_C,i,X0,τ,Smax,τ_C,a_GPP,b_GPP,a_Ec,b_Ec
    par =  [0.0028160038624083067, 27.33276184313149,
     0.028643648082103606, 0.09336572808750652,
      0.0019727630599060625, 2.2083338895077254,
       -6.817991596443347, 13.184687442117964,
        20.387232607829905, 2.6564297126274803,
         0.434540406001315, 0.06910534451564838,
          0.13628204854679268, 0.1383012394158192]

    weather_file="Weather_RO_2"  

    RO_data = Load_RO_data(weather_file=weather_file) 

    X0,τ,Smax,τ_C = par[7:10]
    weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized",X0=X0,τ=τ,Smax=Smax,weather_file= weather_file)
    weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control",X0=X0,τ=τ_C,Smax=Smax,weather_file= weather_file)   

    data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C,weatherts_F,weatherts_C)
    
    n_data = length(weatherts_F.date)

    GPP_model_F,GPP_model_C,EC_model_F,EC_model_C,gₛ_model_F,gₛ_model_C = Get_Data_RO_C_F_CCPH(par,weatherts_F,weatherts_C,data_F,data_C;sim_steps=n_data-1)
        
    ind = 85:n_data
    
    @show GPP_R²_F =  calcR²(data_F.GPP[ind],GPP_model_F[ind])
    @show GPP_cor_F = cor(data_F.GPP[ind],GPP_model_F[ind])
    @show gₛ_R²_F =  calcR²(data_F.E_C[ind],EC_model_F[ind])
    @show gₛ_cor_F = cor(data_F.E_C[ind],EC_model_F[ind])

    @show GPP_R²_C =  calcR²(data_C.GPP[ind],GPP_model_C[ind])
    @show GPP_cor_C = cor(data_C.GPP[ind],GPP_model_C[ind])
    @show gₛ_R²_C =  calcR²(data_C.E_C[ind],EC_model_C[ind])
    @show gₛ_cor_C = cor(data_C.E_C[ind],EC_model_C[ind])    

    plot(weatherts_F.date,data_F.GPP,label="Data",ylabel="GPP")   
    pl1 = plot!(weatherts_F.date,GPP_model_F,label="Model")
    plot(weatherts_C.date,data_C.GPP,label="Data",ylabel="GPP")
    pl2 = plot!(weatherts_C.date,GPP_model_C,label="Model")

    plot(weatherts_F.date,data_F.E_C,label="Data",ylabel="E_C")
    pl3 = plot!(weatherts_F.date,EC_model_F,label="Model")
    plot(weatherts_C.date,data_C.E_C,label="Data",ylabel="E_C")
    pl4 = plot!(weatherts_C.date,EC_model_C,label="Model")

    plot(pl1,pl3,pl2,pl4,layout=(2,2),legends=false) 

    a_GPP,b_GPP,a_Ec,b_Ec = par[11:14] 
    
    σ_GPP_F = a_GPP.+GPP_model_F*b_GPP
    σ_GPP_C = a_GPP.+GPP_model_C*b_GPP

    σ_EC_F = a_Ec.+b_Ec*EC_model_F
    σ_EC_C = a_Ec.+b_Ec*EC_model_C

    c = -log(0.025*2) #Use for calculating 95% credible intervals

    plot(weatherts_F.date[ind],data_F.GPP[ind],label="Data",ylabel="GPP")
    plot!(weatherts_F.date[ind],GPP_model_F[ind]+c*σ_GPP_F[ind])
    plot!(weatherts_F.date[ind],GPP_model_F[ind]-c*σ_GPP_F[ind])
    pl1 = plot!(weatherts_F.date[ind],GPP_model_F[ind],label="Model")
    plot(weatherts_C.date[ind],data_C.GPP[ind],label="Data",ylabel="GPP")
    plot!(weatherts_C.date[ind],GPP_model_C[ind]+c*σ_GPP_C[ind])
    plot!(weatherts_C.date[ind],GPP_model_C[ind]-c*σ_GPP_C[ind])
    pl2 = plot!(weatherts_C.date[ind],GPP_model_C[ind],label="Model")

    plot(weatherts_F.date[ind],data_F.E_C[ind],label="Data",ylabel="E_C")
    plot!(weatherts_F.date[ind],EC_model_F[ind]+c*σ_EC_F[ind])
    plot!(weatherts_F.date[ind],EC_model_F[ind]-c*σ_EC_F[ind])
    pl3 = plot!(weatherts_F.date[ind],EC_model_F[ind],label="Model")
    plot(weatherts_C.date[ind],data_C.E_C[ind],label="Data",ylabel="E_C")
    plot!(weatherts_C.date[ind],EC_model_C[ind]+c*σ_EC_C[ind])
    plot!(weatherts_C.date[ind],EC_model_C[ind]-c*σ_EC_C[ind])
    pl4 = plot!(weatherts_C.date[ind],EC_model_C[ind],label="Model")

    plot(pl1,pl2,pl3,pl4,layout=(2,2),legends=false)
    savefig("./plots/"*file_name*"_validation.svg")
end

run_validation_RO_2019()