#Create water use efficiency use plot
function create_wue_plot(fig_name::String,model_result::ModelResult,weatherts::WeatherTS)    
    wue = model_result.GPP./model_result.Ec*10 #kg C ha⁻¹ mm⁻¹
    wue_mol = model_result.A./model_result.E*1000 #mmol C mol⁻¹ H₂O
    
    pl1 = plot(weatherts.PAR*10^6,wue_mol,xlabel="I₀ (μmol s⁻¹ m⁻²)",ylabel="WUE (mmol mol⁻¹)",seriestype=:scatter)
    pl2 = plot(weatherts.temp,wue_mol,xlabel="Tₐ (ᵒC)",ylabel="WUE (mmol mol⁻¹)",seriestype=:scatter)
    pl3 = plot(weatherts.VPD/1000,wue_mol,xlabel="VPD (kPa)",ylabel="WUE (mmol mol⁻¹)",seriestype=:scatter)
    pl4 = plot(weatherts.θₛ*100,wue_mol,xlabel="θ (%)",ylabel="WUE (mmol mol⁻¹)",seriestype=:scatter)

    plot(pl1,pl2,pl3,pl4,layout=4,legends=false)
    savefig(fig_name*".svg")
end

#Create assimilation plot
function create_A_plot(fig_name::String,model_result::ModelResult,weatherts::WeatherTS)    
    A = model_result.A.*10^6 #μmol C m⁻² leaf area s⁻¹
    
    pl1 = plot(weatherts.PAR*10^6,A,xlabel="I₀ (μmol s⁻¹ m⁻²)",ylabel="A (μmol C m⁻² s⁻¹)",seriestype=:scatter)
    pl2 = plot(weatherts.temp,A,xlabel="Tₐ (ᵒC)",ylabel="A (μmol C m⁻² s⁻¹)",seriestype=:scatter)
    pl3 = plot(weatherts.VPD/1000,A,xlabel="VPD (kPa)",ylabel="A (μmol C m⁻² s⁻¹)",seriestype=:scatter)
    pl4 = plot(weatherts.θₛ*100,A,xlabel="θ (%)",ylabel="A (μmol C m⁻² s⁻¹)",seriestype=:scatter)

    plot(pl1,pl2,pl3,pl4,layout=4,legends=false)
    savefig(fig_name*".svg")
end

#Create cᵢ plot
function create_cᵢ_plot(fig_name::String,model_result::ModelResult,weatherts::WeatherTS)    
    cᵢ = model_result.cᵢ #Pa
    
    pl1 = plot(weatherts.PAR*10^6,cᵢ,xlabel="I₀ (μmol s⁻¹ m⁻²)",ylabel="cᵢ (Pa)",seriestype=:scatter)
    pl2 = plot(weatherts.temp,cᵢ,xlabel="Tₐ (ᵒC)",ylabel="cᵢ (Pa)",seriestype=:scatter)
    pl3 = plot(weatherts.VPD/1000,cᵢ,xlabel="VPD (kPa)",ylabel="cᵢ (Pa)",seriestype=:scatter)
    pl4 = plot(weatherts.θₛ*100,cᵢ,xlabel="θ (%)",ylabel="cᵢ (Pa)",seriestype=:scatter)

    plot(pl1,pl2,pl3,pl4,layout=4,legends=false)
    savefig(fig_name*".svg")
end

#Optimal traits
function create_trait_plots(fig_name::String,model_result::ModelResult,weatherts::WeatherTS)
    
    pl1 = plot(weatherts.PAR*10^6,model_result.gₛ,xlabel="I₀ (μmol s⁻¹ m⁻²)",ylabel="gₛ (mol s⁻¹ m⁻²)",seriestype=:scatter)
    pl2 = plot(weatherts.temp,model_result.gₛ,xlabel="Tₐ (ᵒC)",ylabel="gₛ (mol s⁻¹ m⁻²)",seriestype=:scatter)
    pl3 = plot(weatherts.VPD/1000,model_result.gₛ,xlabel="VPD (kPa)",ylabel="gₛ (mol s⁻¹ m⁻²)",seriestype=:scatter)
    pl4 = plot(weatherts.θₛ*100,model_result.gₛ,xlabel="θ (%)",ylabel="gₛ (mol s⁻¹ m⁻²)",seriestype=:scatter)    
    
    plot(pl1,pl2,pl3,pl4,layout=4,legends=false)
    savefig(fig_name*"_gₛ_weather.svg")


    Nₘ_f = model_result.Nₘ_f*100
    pl1 = plot(weatherts.PAR*10^6,Nₘ_f,xlabel="I₀ (μmol s⁻¹ m⁻²)",ylabel="Nₘ_f (%)",seriestype=:scatter)
    pl2 = plot(weatherts.temp,Nₘ_f,xlabel="Tₐ (ᵒC)",ylabel="Nₘ_f (%)",seriestype=:scatter)
    pl3 = plot(weatherts.VPD/1000,Nₘ_f,xlabel="VPD (kPa)",ylabel="Nₘ_f (%)",seriestype=:scatter)
    pl4 = plot(weatherts.θₛ*100,Nₘ_f,xlabel="θ (%)",ylabel="Nₘ_f (%)",seriestype=:scatter)    
    
    plot(pl1,pl2,pl3,pl4,layout=4,legends=false)
    savefig(fig_name*"_Nₘ_f_weather.svg")

    plot(model_result.gₛ,Nₘ_f,xlabel="gₛ (mol s⁻¹ m⁻²)",ylabel="Nₘ_f (%)",legends=false,seriestype=:scatter)
    savefig(fig_name*"_Nₘ_f_gₛ.svg")
end

#Optimal traits for both treatments (F and C)
function create_trait_plots(fig_name::String,
    model_result_F::ModelResult,
    weatherts_F::WeatherTS,
    model_result_C::ModelResult,
    weatherts_C::WeatherTS)
    
    pl1 = plot(weatherts_F.PAR*10^6,model_result_F.gₛ,xlabel="I₀ (μmol s⁻¹ m⁻²)",ylabel="gₛ (mol s⁻¹ m⁻²)",seriestype=:scatter)
    pl1 = plot!(weatherts_C.PAR*10^6,model_result_C.gₛ,seriestype=:scatter)
    pl2 = plot(weatherts_F.temp,model_result_F.gₛ,xlabel="Tₐ (ᵒC)",ylabel="gₛ (mol s⁻¹ m⁻²)",seriestype=:scatter)
    pl2 = plot!(weatherts_C.temp,model_result_C.gₛ,seriestype=:scatter)
    pl3 = plot(weatherts_F.VPD/1000,model_result_F.gₛ,xlabel="VPD (kPa)",ylabel="gₛ (mol s⁻¹ m⁻²)",seriestype=:scatter)
    pl3 = plot!(weatherts_C.VPD/1000,model_result_C.gₛ,seriestype=:scatter)
    pl4 = plot(weatherts_F.θₛ*100,model_result_F.gₛ,xlabel="θ (%)",ylabel="gₛ (mol s⁻¹ m⁻²)",seriestype=:scatter)    
    pl4 = plot!(weatherts_C.θₛ*100,model_result_C.gₛ,seriestype=:scatter)

    plot(pl1,pl2,pl3,pl4,layout=4,legends=false)
    savefig(fig_name*"_gₛ_weather.svg")

    Nₘ_f = model_result_F.Nₘ_f*100    
    pl1 = plot(weatherts_F.PAR*10^6,Nₘ_f,xlabel="I₀ (μmol s⁻¹ m⁻²)",ylabel="Nₘ_f (%)",seriestype=:scatter)
    pl2 = plot(weatherts_F.temp,Nₘ_f,xlabel="Tₐ (ᵒC)",ylabel="Nₘ_f (%)",seriestype=:scatter)
    pl3 = plot(weatherts_F.VPD/1000,Nₘ_f,xlabel="VPD (kPa)",ylabel="Nₘ_f (%)",seriestype=:scatter)
    pl4 = plot(weatherts_F.θₛ*100,Nₘ_f,xlabel="θ (%)",ylabel="Nₘ_f (%)",seriestype=:scatter)    
    
    Nₘ_f = model_result_C.Nₘ_f*100    
    plot!(pl1,weatherts_C.PAR*10^6,Nₘ_f,seriestype=:scatter)
    plot!(pl2,weatherts_C.temp,Nₘ_f,seriestype=:scatter)
    plot!(pl3,weatherts_C.VPD/1000,Nₘ_f,seriestype=:scatter)
    plot!(pl4,weatherts_C.θₛ*100,Nₘ_f,seriestype=:scatter)

    plot(pl1,pl2,pl3,pl4,layout=4,legends=false)
    savefig(fig_name*"_Nₘ_f_weather.svg")

    Nₘ_f_F = model_result_F.Nₘ_f*100
    Nₘ_f_C = model_result_C.Nₘ_f*100 
    plot(model_result_F.gₛ,Nₘ_f_F,xlabel="gₛ (mol s⁻¹ m⁻²)",ylabel="Nₘ_f (%)",legends=false,seriestype=:scatter)
    plot!(model_result_C.gₛ,Nₘ_f_C,seriestype=:scatter)
    savefig(fig_name*"_Nₘ_f_gₛ.svg")

    plot(weatherts_F.date,model_result_F.gₛ,ylabel="gₛ (mol s⁻¹ m⁻²)",legends=false)
    pl1 = plot!(weatherts_C.date,model_result_C.gₛ)
    plot(weatherts_F.date,Nₘ_f_F,ylabel="Nₘ_f (%)",legends=false)
    pl2 = plot!(weatherts_C.date,Nₘ_f_C)
    
    savefig(pl1,fig_name*"_gₛ_date.svg")
    savefig(pl2,fig_name*"_Nₘ_f_date.svg")
end

#Create water use efficiency use plot for both treatments (F and C)
function create_wue_plot(fig_name::String,
    model_result_F::ModelResult,
    weatherts_F::WeatherTS,
    model_result_C::ModelResult,
    weatherts_C::WeatherTS)    
    
    
    wue = model_result_F.GPP./model_result_F.Ec*10 #kg C ha⁻¹ mm⁻¹
    wue_mol = model_result_F.A./model_result_F.E*1000 #mmol C mol⁻¹ H₂O
    
    pl1 = plot(weatherts_F.PAR*10^6,wue_mol,xlabel="I₀ (μmol s⁻¹ m⁻²)",ylabel="WUE (mmol mol⁻¹)",seriestype=:scatter)
    pl2 = plot(weatherts_F.temp,wue_mol,xlabel="Tₐ (ᵒC)",ylabel="WUE (mmol mol⁻¹)",seriestype=:scatter)
    pl3 = plot(weatherts_F.VPD/1000,wue_mol,xlabel="VPD (kPa)",ylabel="WUE (mmol mol⁻¹)",seriestype=:scatter)
    pl4 = plot(weatherts_F.θₛ*100,wue_mol,xlabel="θ (%)",ylabel="WUE (mmol mol⁻¹)",seriestype=:scatter)
    pl5 = plot(model_result_F.GPP,wue_mol,xlabel="GPP (g C day⁻¹ m⁻²)",ylabel="WUE (mmol mol⁻¹)",seriestype=:scatter)
    pl6 = plot(model_result_F.Ec,wue_mol,xlabel="Ec (mm day⁻¹)",ylabel="WUE (mmol mol⁻¹)",seriestype=:scatter)

    wue = model_result_C.GPP./model_result_C.Ec*10 #kg C ha⁻¹ mm⁻¹
    wue_mol = model_result_C.A./model_result_C.E*1000 #mmol C mol⁻¹ H₂O

    plot!(pl1,weatherts_C.PAR*10^6,wue_mol,seriestype=:scatter)
    plot!(pl2,weatherts_C.temp,wue_mol,seriestype=:scatter)
    plot!(pl3,weatherts_C.VPD/1000,wue_mol,seriestype=:scatter)
    plot!(pl4,weatherts_C.θₛ*100,wue_mol,seriestype=:scatter)
    plot!(pl5,model_result_C.GPP,wue_mol,seriestype=:scatter)
    plot!(pl6,model_result_C.Ec,wue_mol,seriestype=:scatter)

    plot(pl1,pl2,pl3,pl4,pl5,pl6,layout=(3,2),legends=false)
    savefig(fig_name*".svg")

    wue_mol_F = model_result_F.A./model_result_F.E*1000 #mmol C mol⁻¹ H₂O
    wue_mol_E = model_result_C.A./model_result_C.E*1000 #mmol C mol⁻¹ H₂O
    ψₛ_F = CCPH.θₛ2ψₛ.(weatherts_F.θₛ)
    ψₛ_C = CCPH.θₛ2ψₛ.(weatherts_C.θₛ)
    pl1 = plot(ψₛ_F,wue_mol_F,xlabel="ψₛ (MPa)",ylabel="WUE (mmol mol⁻¹)",seriestype=:scatter)
    plot!(pl1,ψₛ_C,wue_mol_E,seriestype=:scatter) 

    ψₛ_vec = range(-1,stop=-0.0001,length=300)
    Pval_vec = CCPH.Pfun.(ψₛ_vec,Ref(-2.0),(2.0))
    Pval_F = CCPH.Pfun.(ψₛ_F,Ref(-2.0),(2.0))
    Pval_C = CCPH.Pfun.(ψₛ_C,Ref(-2.0),(2.0))
    

    pl2 = plot(ψₛ_F,Pval_F,ψₛ_vec,xlabel="ψₛ (MPa)",ylabel="P",seriestype=:scatter,xaxis=:flip)
    plot!(pl2,ψₛ_C,Pval_C,seriestype=:scatter)
    plot!(pl2,ψₛ_vec,Pval_vec)
    plot!(pl2,[-0.3,-0.3],[first(Pval_vec),last(Pval_vec)])
    
    plot(pl1,pl2,layout=(1,2),legend=false)
    savefig(fig_name*"_SP.svg")
end
