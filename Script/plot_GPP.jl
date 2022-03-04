RO_data = Load_RO_data() 
weatherts_F,GPP_data_F = create_weather_struct_RO(RO_data;stand_type="Fertilized")
weatherts_C,GPP_data_C = create_weather_struct_RO(RO_data;stand_type="Control")

data_F,data_C = Create_RoData_C_F(GPP_data_F,GPP_data_C,weatherts_F,weatherts_C) 

treepar = TreePar()
cons = Constants()

gₛ_alt = 0.04 #mol CO₂ m⁻² s⁻¹
J_max_opt_season = 184.0e-6

A_c_F = Float64[]
Cᵢ_F = Float64[]
GPP_F = Float64[]
C_comp = Float64[]

for i = 1:length(weatherts_F.date)

    photorates = PhotoKineticRates()
    env = EnvironmentStruct(weatherts_F,i)      
    photopar = PhotoPar(photorates,env.Tₐ)
    xₜ = weatherts_F.acclimation_fac[i]

    C_crit = (photopar.Γ+photopar.b_r)/(1-photopar.b_r)

    ind = Find_data_ind(i)
    LAI = data_F.Wf[ind]*data_F.N[ind]/treepar.LMA    
    N = data_F.N[ind]

    daylight = weatherts_F.daylight[i]/7

    Iᵢ = CCPH.Calc_Iᵢ(env.I₀,treepar)
    gₜ = CCPH.Calc_gₜ(data_F.gₛ[i],treepar)
    Jₘₐₓ = xₜ*photopar.b_Jmax*J_max_opt_season
    photopar.α = CCPH.Calc_α(Jₘₐₓ,treepar.r_α)

    A_c,Cᵢ = CCPH.Farquhar(gₜ,Iᵢ,Jₘₐₓ,photopar,env)    

    GPP = cons.M_C*(1-exp(-treepar.k*LAI))/(treepar.k)*A_c*daylight*10^3

    push!(A_c_F,A_c)
    push!(GPP_F,GPP)
    push!(Cᵢ_F,Cᵢ)
    push!(C_comp,C_crit)
end

plot(weatherts_F.date,data_F.GPP,ylabel="GPP")
plot!(weatherts_F.date,GPP_F)

@show cor(data_F.GPP,GPP_F)

plot(weatherts_F.date,A_c_F,ylabel="A_c")

plot(weatherts_F.date,data_F.gₛ,ylabel="gₛ")

plot(weatherts_F.date,Cᵢ_F,ylabel="Cᵢ")
plot!(weatherts_F.date,C_comp)
plot!([weatherts_F.date[1],weatherts_F.date[end]],[40,40])