mutable struct RoData{T<:Float64}
    H::Array{T,1}
    B::Array{T,1}
    GPP::Array{T,1}
    Wf::Array{T,1}
    N::Array{T,1}
    Hc::T
    ET::Array{T,1}
end

function Create_RoData_C_F(GPP_data_F::Array{T,1},GPP_data_C::Array{T,1},
    ET_data_F::Array{T,1},ET_data_C::Array{T,1}) where {T<:Float64} 
    
    #=
    H_data_F = [17.004166666666666, 17.271666666666665, 17.59166666666667, 17.835] 
    B_data_F = [0.025076878666666667, 0.02603588483333334, 0.026922637833333332, 0.027666997]
    Wf_data_F = [4.62,4.72,4.77,4.82]
    N_data_F = [850.0000,846.6667,846.6667,846.6667]/10000

    H_data_C = [16.956,17.165,17.388,17.590] 
    B_data_C = [0.023015244833333334, 0.02349805816666667, 0.023938368499999998, 0.02432245283333333]
    Wf_data_C = [4.18,4.20,4.20,4.2]
    N_data_C = [1010.0000,1010.0000,1006.6667,1006.6667]/10000
    =#
    
    H_data_F = [19.07,19.34,19.64,19.87] 
    B_data_F = [0.0452,0.0467,0.0479,0.0487]
    Wf_data_F = [8.33,8.47,8.49,8.48]
    N_data_F = [850.0000,846.6667,846.6667,846.6667]/10000
    Hc_F = 10.89

    H_data_C = [20.86,21.01,21.19,21.36] 
    B_data_C = [0.0348,0.0356,0.0362,0.0368]
    Wf_data_C = [5.025,5.110,5.141,5.165]
    N_data_C = [1010.0000,1010.0000,1006.6667,1006.6667]/10000
    Hc_C = 10.74

    data_F = RoData(H_data_F,B_data_F,GPP_data_F,Wf_data_F,N_data_F,Hc_F,ET_data_F)
    data_C = RoData(H_data_C,B_data_C,GPP_data_C,Wf_data_C,N_data_C,Hc_C,ET_data_C)

    return data_F,data_C
end

function calcR²(y_data::Array{Float64,1},f_model::Array{Float64,1})
    ymean = mean(y_data)
    SStot = sum((y_data.-ymean).^2)
    SSres = sum((y_data.-f_model).^2)
    return 1-SSres/SStot
end

function Max_k_cost(model::CCPHStruct)
    ψ₅₀,b,g,ρ_H2O,θₛ = model.hydPar.ψ₅₀,model.hydPar.b,model.cons.g,model.cons.ρ_H2O,model.env.θₛ
    #Calculate tree height
    H = model.treesize.H
    #Caluclate soil water potential
    ψₛ = CCPH.θₛ2ψₛ(θₛ)
    #Soil water potential adjusted for gravitational pressure (MPa)
    ψₛ_g = ψₛ-H*ρ_H2O*g*10^-6 

    k_cost_max = CCPH.Pfun(ψₛ_g,ψ₅₀,b)

    return k_cost_max
end

function OptimCCPH(growthlength::Float64,model::CCPH.CCPHStruct)
    K_cost_crit = 0.12
    gₛ_crit = CCPH.Calc_K_costᵢₙᵥ(K_cost_crit,model)

    gₛ_opt,Nₘ_f_opt = CCPH.CCPHTraitmodel(growthlength,model;gₛ_lim_hi=gₛ_crit)

    modeloutput = CCPH.CCPH_run(gₛ_opt,Nₘ_f_opt,growthlength,model)

    return modeloutput,gₛ_opt,Nₘ_f_opt
end
function OptimCCPH!(i::Integer,model::CCPH.CCPHStruct,
    weatherts::CCPH.WeatherTS,photo_kinetic::CCPH.PhotoKineticRates)
    growthlength,step_length = CCPH.Init_weather_par!(i,model,weatherts,photo_kinetic)
    modeloutput,gₛ_opt,Nₘ_f_opt = OptimCCPH(growthlength,model)

    return modeloutput,gₛ_opt,Nₘ_f_opt
end

#Calc model GPP
function CalcGPP(ccphts::CCPHTS,weatherts::WeatherTS)
    GPP = ccphts.N.*ccphts.P.*weatherts.daylight./weatherts.tot_annual_daylight*1000/7 #g C day⁻¹ m⁻² ground area

    return GPP
end
function CalcGPP(step_length::Float64,model::CCPH.CCPHStruct,modeloutput::CCPH.CCPHOutput)
    GPP = model.treesize.N*modeloutput.P*step_length*1000/7 #g C day⁻¹ m⁻² ground area

    return GPP
end
function CalcGPP!(i::Integer,model::CCPH.CCPHStruct,
    weatherts::CCPH.WeatherTS,photo_kinetic::CCPH.PhotoKineticRates)
    growthlength,step_length = CCPH.Init_weather_par!(i,model,weatherts,photo_kinetic)
    modeloutput,gₛ_opt,Nₘ_f_opt = OptimCCPH(growthlength,model)
    GPP = CalcGPP(step_length,model,modeloutput)
    Transpiraiton = CalcTranspiraiton(gₛ_opt,model,growthlength*step_length)
    return GPP, Transpiraiton 
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

function NPP(H::T,gₛ::T,Nₘ_f::T,αr::T,growthlength::T,model::CCPHStruct) where {T<:Float64}

    A = model.treesize.As*((H-model.treesize.Hs)/(model.treesize.H-model.treesize.Hs))^model.treepar.z
    Ww = model.treepar.ρw*A*(model.treepar.β₁*H+model.treepar.β₂*model.treesize.Hs) 
    Wr = αr*A
    Wf = model.treepar.αf*A
    LAI = Wf*model.treesize.N/model.treepar.LMA

    #Calculate per sapwood mass nitrogen concentration
    Nₘ_w = CCPH.Calc_Nₘ_w(Nₘ_f,model)
    #Calcualte per fine roots mass nitrogen concentration
    Nₘ_r = CCPH.Calc_Nₘ_r(Nₘ_f,model)
    #Calcualte per leaf area nitrogen concentration    
    Nₐ = CCPH.Calc_Nₐ(Nₘ_f,model)
    #Calculate total maintenance respiration 
    Rₘ = CCPH.Calc_Rₘ(Nₘ_f,Nₘ_w,Nₘ_r,Wf,Ww,Wr,model)
    #calcualte total conductance
    gₜ = CCPH.Calc_gₜ(gₛ,model)    
    #Calculate Jmax
    Jmax = CCPH.Calc_Jmax(Nₐ,model.treepar.a_Jmax,model.treepar.b_Jmax,model.photopar.b_Jmax)
    #Irradiance incident on a leaf at canopy top
    Iᵢ = CCPH.Calc_Iᵢ(model.env.I₀,model)
    #Calculate per tree carbon assimilation
    P = CCPH.GPP(gₜ,Iᵢ,Jmax,LAI,growthlength,model) 

    Δ_leaf = CCPH.Calc_Δ_leaf(gₜ,Iᵢ,LAI,growthlength,Nₘ_f,Jmax,model)

    return P-Rₘ,Δ_leaf,LAI
end

function CalcTranspiraiton(gₛ::Float64,model::CCPH.CCPHStruct,daylight::Float64)
    LAI = model.treesize.Wf.*model.treesize.N/model.treepar.LMA
    tranpiraiton = model.cons.r*model.cons.M_H2O/model.cons.ρ_H2O*
    model.env.VPD/model.env.P*gₛ*LAI.*daylight*1000/7
    return tranpiraiton
end

function CalcTranspiraiton(ccphts::CCPHTS,weatherts::WeatherTS,cons::CCPH.Constants,LMA::Float64,P::Float64)
    LAI = ccphts.Wf.*ccphts.N/LMA    
    tranpiraiton = cons.r*cons.M_H2O/cons.ρ_H2O.*weatherts.VPD/P.*ccphts.gₛ.*LAI.*weatherts.daylight*1000/7
    return tranpiraiton
end

function EstgₛFromTranspiraiton(tranpiraiton::Float64,ccphts::CCPHTS,
    weatherts::WeatherTS,cons::CCPH.Constants,LMA::Float64,P::Float64)
    LAI = ccphts.Wf.*ccphts.N/LMA 
    gₛ = tranpiraiton*cons.ρ_H2O/(cons.r*cons.M_H2O)*P./weatherts.VPD./(LAI.*weatherts.daylight)*7/1000

    return gₛ
end
function EstgₛFromTranspiraiton(tranpiraiton::Array{T,1},ccphts::CCPHTS,
    weatherts::WeatherTS,cons::CCPH.Constants,LMA::T,P::T) where {T<:Float64}
    LAI = ccphts.Wf.*ccphts.N/LMA 
    gₛ = tranpiraiton.*cons.ρ_H2O/(cons.r*cons.M_H2O)*P./weatherts.VPD./(LAI.*weatherts.daylight)*7/1000

    return gₛ
end

function CalcGPP(gₛ::Array{T,1},ccphts::CCPHTS,weatherts::WeatherTS,cons::CCPH.Constants,
    env::CCPH.EnvironmentStruct,treepar::CCPH.TreePar,kinetic::CCPH.PhotoKineticRates) where {T<:Float64}
   
    GPP_vec = Float64[]
    for i = 1:length(gₛ)

        Temp = weatherts.temp[i]
        Xₜ = weatherts.acclimation_fac[i]
        growthlength = weatherts.daylight[i]
        #Irradiance
        I₀ = weatherts.PAR[i]
        #Create Photosynthesis parameters
        photopar = PhotoPar(kinetic,Temp)
        LAI = ccphts.Wf[i].*ccphts.N[i]/treepar.LMA
        #calcualte total conductance
        gₜ = gₛ[i]*treepar.r_gₛ
        #Calc J_max
        Jmax = photopar.b_Jmax*kinetic.Jmax.Kₒₚₜ
        #Irradiance incident on a leaf at canopy top
        Iᵢ = I₀*treepar.k/(1-treepar.m)
        #Calcualte leaf C assimilation
        A_C = CCPH.Farquhar(gₜ,Iᵢ,Jmax,photopar,env)[1]*Xₜ*cons.M_C*growthlength*1000/7 

        GPP = A_C*(1-exp(-treepar.k*LAI))/(ccphts.N[i]*treepar.k) #g C day⁻¹ m⁻² ground area

        push!(GPP_vec,GPP) 
    end
    return GPP_vec
end

struct RO_Raw_ET_Data{T<:CSV.File}
    RO_ET_F_2015::T
    RO_ET_F_2016::T
    RO_ET_F_2017::T
    RO_ET_F_2018::T
    RO_ET_C_2015::T
    RO_ET_C_2016::T
    RO_ET_C_2017::T
    RO_ET_C_2018::T
end
function RO_Raw_ET_Data()
    RO_ET_F_2015 = CSV.File("./Data/RO_ET_Fertilized_2015.csv")
    RO_ET_F_2016 = CSV.File("./Data/RO_ET_Fertilized_2016.csv")
    RO_ET_F_2017 = CSV.File("./Data/RO_ET_Fertilized_2017.csv")
    RO_ET_F_2018 = CSV.File("./Data/RO_ET_Fertilized_2018.csv")
    RO_ET_C_2015 = CSV.File("./Data/RO_ET_Control_2015.csv")
    RO_ET_C_2016 = CSV.File("./Data/RO_ET_Control_2016.csv")
    RO_ET_C_2017 = CSV.File("./Data/RO_ET_Control_2017.csv")
    RO_ET_C_2018 = CSV.File("./Data/RO_ET_Control_2018.csv")

    return RO_Raw_ET_Data(RO_ET_F_2015,RO_ET_F_2016,RO_ET_F_2017,RO_ET_F_2018,
    RO_ET_C_2015,RO_ET_C_2016,RO_ET_C_2017,RO_ET_C_2018)
end

function RO_ET_Est(ET_data::RO_Raw_ET_Data,weatherts::WeatherTS;stand_type::String="Fertilized")
    if stand_type == "Fertilized"
        RO_ET_2015 = ET_data.RO_ET_F_2015
        RO_ET_2016 = ET_data.RO_ET_F_2016
        RO_ET_2017 = ET_data.RO_ET_F_2017
        RO_ET_2018 = ET_data.RO_ET_F_2017         
    elseif stand_type == "Control"
        RO_ET_2015 = ET_data.RO_ET_C_2015
        RO_ET_2016 = ET_data.RO_ET_C_2016
        RO_ET_2017 = ET_data.RO_ET_C_2017
        RO_ET_2018 = ET_data.RO_ET_C_2017 
    else
        error("Wrong input. Either \"Fertilized\" or \"Control\"")
    end 

    ET_Lin_Int_2015 = Interpolations.LinearInterpolation(RO_ET_2015.Day, RO_ET_2015.ET)
    ET_Lin_Int_2016 = Interpolations.LinearInterpolation(RO_ET_2016.Day, RO_ET_2016.ET)
    ET_Lin_Int_2017 = Interpolations.LinearInterpolation(RO_ET_2017.Day, RO_ET_2017.ET)
    ET_Lin_Int_2018 = Interpolations.LinearInterpolation(RO_ET_2018.Day, RO_ET_2018.ET)

    ET_est = Float64[]
    
    for d in weatherts.date
        y = Dates.year(d)
        d_nr = Dates.dayofyear(d)

        if y == 2015
            push!(ET_est,ET_Lin_Int_2015(d_nr))
        elseif y == 2016
            push!(ET_est,ET_Lin_Int_2016(d_nr))
        elseif y == 2017
            push!(ET_est,ET_Lin_Int_2017(d_nr))
        elseif y == 2018
            push!(ET_est,ET_Lin_Int_2018(d_nr))
        else
            error("Date must be between 2015-2018")
        end    
    end
    return ET_est
end

function RO_ET_Est(weatherts::WeatherTS;stand_type::String="Fertilized")
    if stand_type == "Fertilized"
        RO_ET_2015 = CSV.File("./Data/RO_ET_Fertilized_2015.csv")
        RO_ET_2016 = CSV.File("./Data/RO_ET_Fertilized_2016.csv")
        RO_ET_2017 = CSV.File("./Data/RO_ET_Fertilized_2017.csv")
        RO_ET_2018 = CSV.File("./Data/RO_ET_Fertilized_2018.csv")         
    elseif stand_type == "Control"
        RO_ET_2015 = CSV.File("./Data/RO_ET_Control_2015.csv")
        RO_ET_2016 = CSV.File("./Data/RO_ET_Control_2016.csv")
        RO_ET_2017 = CSV.File("./Data/RO_ET_Control_2017.csv")
        RO_ET_2018 = CSV.File("./Data/RO_ET_Control_2018.csv") 
    else
        error("Wrong input. Either \"Fertilized\" or \"Control\"")
    end 

    ET_Lin_Int_2015 = Interpolations.LinearInterpolation(RO_ET_2015.Day, RO_ET_2015.ET)
    ET_Lin_Int_2016 = Interpolations.LinearInterpolation(RO_ET_2016.Day, RO_ET_2016.ET)
    ET_Lin_Int_2017 = Interpolations.LinearInterpolation(RO_ET_2017.Day, RO_ET_2017.ET)
    ET_Lin_Int_2018 = Interpolations.LinearInterpolation(RO_ET_2018.Day, RO_ET_2018.ET)

    ET_est = Float64[]
    
    for d in weatherts.date
        y = Dates.year(d)
        d_nr = Dates.dayofyear(d)

        if y == 2015
            push!(ET_est,ET_Lin_Int_2015(d_nr))
        elseif y == 2016
            push!(ET_est,ET_Lin_Int_2016(d_nr))
        elseif y == 2017
            push!(ET_est,ET_Lin_Int_2017(d_nr))
        elseif y == 2018
            push!(ET_est,ET_Lin_Int_2018(d_nr))
        else
            error("Date must be between 2015-2018")
        end    
    end
    return ET_est
end