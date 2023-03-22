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