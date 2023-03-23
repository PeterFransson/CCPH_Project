stand_type="Fertilized"
time_idx = 15

parasym = [:Nₛ,
    :a_Jmax,
    :Kₓₗ₀,    
    :α_max,
    :τ]

par_F = [0.015,0.014,0.00079,0.15,14.6]

ParaDict = CreateParaDict(parasym,par_F;ParaDictInit=nothing)

RO_data = Load_RO_data()

X0,τ,Smax = ParaDict[:X0],ParaDict[:τ],ParaDict[:Smax]

weatherts,GPP_data = create_weather_struct_RO(RO_data;stand_type=stand_type,X0=X0,τ=τ,Smax=Smax,
    weather_file="Weather_RO")
data = Create_RoData(GPP_data,weatherts;stand_type=stand_type)

model,kinetic=Initi_model_struct(time_idx::Integer,data::RoData,ParaDict::Dict{Symbol,Float64})
growthlength,step_length = CCPH.Init_weather_par!(time_idx,model,weatherts,kinetic) 

@show model.env.I₀
@show model.env.VPD
@show model.env.θₛ

n_val = 100
I₀_vec = range(0.0004,stop=0.001,length=n_val)
VPD_vec = range(250,stop=1000,length=n_val)
θₛ_vec =  range(0.17,stop=0.26,length=n_val)

gₛ_vec = zeros(n_val)
Nₘ_f_vec = zeros(n_val)

for i in 1:n_val  
    #@show model.env.I₀ = I₀_vec[i]
    @show model.env.VPD = VPD_vec[i]
    #@show model.env.θₛ = θₛ_vec[i]

    modeloutput,gₛ_opt,Nₘ_f_opt = OptimCCPH(growthlength,model)

    gₛ_vec[i]=gₛ_opt
    Nₘ_f_vec[i]=Nₘ_f_opt
end
#xvec = I₀_vec*10^6 
xvec = VPD_vec*10^-3
#xvec = θₛ_vec*100

pl1 = plot(xvec,gₛ_vec,xlabel="",ylabel="gₛ")
pl2 = plot(xvec,Nₘ_f_vec,xlabel="",ylabel="Nₘ_f")
plot(pl1,pl2,layout=2,legends=false)