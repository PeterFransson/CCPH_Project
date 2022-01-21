using Pkg
Pkg.activate(".")
using Revise, CCPH, Statistics, Distributions, Plots, JLD
import Dates, CSV, Interpolations

include("./mcmc/metropolis_log_P.jl")
include("./mcmc/analysis.jl")
include("./Weather_Rosinedal_Struct/create_weather_struct_RO.jl")
include("./Script/Fit_Model_2_RO_Data/model_weather_ts_tuning_C_F_Auxiliary.jl")
include("./Transpiration/transpiration.jl")

#=
RO_data = Load_RO_data() 
weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized")
weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control")

data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C,weatherts_F,weatherts_C) 

treepar = TreePar()
cons = CCPH.Constants()

A_c_F = Float64[]
A_c_C = Float64[]
A_c_F_calc = Float64[]
GPP_F_calc = Float64[]

J_max_opt = 184.0e-6
J_max_opt = 400.0e-6
for i = 1:length(weatherts_F.date)

    photorates = CCPH.PhotoKineticRates()
    env = CCPH.EnvironmentStruct(weatherts_F,i)    
    photopar = PhotoPar(photorates,env.Tₐ)
    Jₘₐₓ = J_max_opt*photopar.b_Jmax
    Iᵢ = CCPH.Calc_Iᵢ(env.I₀,treepar)
    gₜ = CCPH.Calc_gₜ(data_F.gₛ[i],treepar)

    A_c = CCPH.Farquhar(gₜ,Iᵢ,Jₘₐₓ,photopar,env)[1]

    push!(A_c_F_calc,A_c)
   
    ind = Find_data_ind(i)

    LAI_F = data_F.Wf[ind]*data_F.N[ind]/treepar.LMA
    N_F = data_F.N[ind]
    push!(A_c_F,Est_C_assimilation(GPP_data_F[i],LAI_F,N_F,weatherts_F,i))
    
    LAI_F = data_C.Wf[ind]*data_C.N[ind]/treepar.LMA
    N_F = data_C.N[ind]
    push!(A_c_C,Est_C_assimilation(GPP_data_C[i],LAI_F,N_F,weatherts_C,i)[1])

    Xₜ = weatherts_F.acclimation_fac[i]
    daylight = weatherts_F.daylight[i]/7
    GPP = A_c*(1-exp(-treepar.k*LAI_F))/treepar.k*Xₜ*cons.M_C*daylight*10^3 #gC m⁻² ground area day⁻¹
    push!(GPP_F_calc,GPP)        
end

plot(weatherts_F.date,data_F.gₛ)
plot!(weatherts_C.date,data_C.gₛ)

@show calcR²(A_c_F,A_c_F_calc)
@show cor(A_c_F,A_c_F_calc)

plot(weatherts_F.date,A_c_F*10^6)
plot!(weatherts_C.date,A_c_C*10^6)
pl1 = plot!(weatherts_C.date,A_c_F_calc*10^6)

plot(A_c_C*10^6,A_c_F_calc*10^6,seriestype=:scatter,xlabel="Data")
pl2 = plot!([3.0,11.0],[3.0,11.0])

plot(pl1,pl2,layout=(1,2))

@show GPP_R²_F =  calcR²(data_F.GPP,GPP_F_calc)
@show GPP_cor_F = cor(data_F.GPP,GPP_F_calc)

plot(weatherts_F.date,data_F.GPP)
plot!(weatherts_C.date,data_C.GPP)
pl1 = plot!(weatherts_F.date,GPP_F_calc)
plot(data_F.GPP,GPP_F_calc,seriestype=:scatter,xlabel="Data")
pl2 = plot!([1.0,9.0],[1.0,9.0])
plot(pl1,pl2,layout=(1,2))
=#

#plot(data_F.gₛ,A_c_F*10^6,seriestype=:scatter)
#plot!(data_C.gₛ,A_c_C*10^6,seriestype=:scatter)

#plot(weatherts_F.date,weatherts_F.PAR*10^6)

include("Script/Fit_Model_2_RO_Data/model_weather_ts_tuning_C_F_GPP_only_Few_para_NoRoot.jl")