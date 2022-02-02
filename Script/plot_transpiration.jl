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
E_c_F = Calctranpiration.(weatherts_F.VPD.*weatherts_F.daylight/(7*24*3600),CCPH.CalcSₑ.(weatherts_F.θₛ))
E_c_C = Calctranpiration.(weatherts_C.VPD.*weatherts_C.daylight/(7*24*3600),CCPH.CalcSₑ.(weatherts_C.θₛ))
 
g_c = 4 #mm s⁻¹
J_max_opt = 184.0e-6
J_max_opt = 400.0e-6

E_c_Calc = g_c*weatherts_C.VPD/10^5.0.*weatherts_F.daylight/7*(cons.ρ_vapor/cons.ρ_H2O)
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

plot(weatherts_F.date,E_c_F,xlabel="Date",ylabel="Canopy tranpiraiton (mm day⁻¹)",label="F")
plot!(weatherts_C.date,E_c_C,label="C")
plot!(weatherts_F.date,E_c_Calc,label="g_c = 4 mm s⁻¹")
savefig("plots/Ec.svg")

plot(weatherts_F.date,data_F.gₛ,xlabel="Date",
ylabel="Stomatal conductance (mol C m⁻² leaf area s⁻¹)",label="F")
plot!(weatherts_C.date,data_C.gₛ,label="C")
savefig("plots/gs.svg")

plot(weatherts_F.date,data_F.GPP,xlabel="Date",
ylabel="GPP (g C m⁻² ground area day⁻¹)",label="F data")
plot!(weatherts_C.date,data_C.GPP,label="C data")
plot!(weatherts_F.date,GPP_F_calc,label="F calc")
savefig("plots/GPP.svg")

pl1 = plot(weatherts_F.date,data_F.GPP,xlabel="Date",ylabel="GPP",label="F data")
plot!(weatherts_F.date,GPP_F_calc)
pl2 = plot(weatherts_F.date,E_c_F,ylabel="E_c",label="F")
pl3 = plot(weatherts_F.date,data_F.gₛ,ylabel="gₛ",label="F")
pl4 = plot(weatherts_F.date,weatherts_F.VPD,ylabel="VPD",label="F")
pl5 = plot(weatherts_F.date,weatherts_F.θₛ,ylabel="θₛ",label="F")
pl6 = plot(weatherts_F.date,A_c_F_calc,ylabel="A_c",label="F")
pl7 = plot(weatherts_F.date,GPP_F_calc./A_c_F_calc,ylabel="Γ")


plot(pl1,pl6,pl7,pl3,layout=(4,1),legends=false)
savefig("plots/fig1.svg")
plot(pl2,pl3,pl4,pl5,layout=(4,1),legends=false)
savefig("plots/fig2.svg")