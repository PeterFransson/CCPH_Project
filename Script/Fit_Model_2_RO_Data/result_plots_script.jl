#Create water use efficiency use plot
function create_wue_plot(fig_name::String,model_result::ModelResult,weatherts::WeatherTS)    
    wue = model_result.GPP./model_result.Ec #g C m⁻² mm⁻¹
    plot(weatherts.date,wue,seriestype=:scatter,ylabel="WUE")
    savefig(fig_name)
end

#Optimal traits
function create_trait_plots(fig_name::String,model_result::ModelResult,weatherts::WeatherTS)
    
    pl1 = plot(weatherts.PAR*10^6,model_result.gₛ,xlabel="PAR (μmol s⁻¹ m⁻²)",ylabel="gₛ",seriestype=:scatter)
    pl2 = plot(weatherts.temp,model_result.gₛ,xlabel="Temp (ᵒC)",ylabel="gₛ",seriestype=:scatter)
    pl3 = plot(weatherts.VPD/1000,model_result.gₛ,xlabel="VPD (kPa)",ylabel="gₛ",seriestype=:scatter)
    pl4 = plot(CCPH.θₛ2ψₛ.(weatherts.θₛ),model_result.gₛ,xlabel="SWP (MPa)",ylabel="gₛ",seriestype=:scatter)    
    
    plot(pl1,pl2,pl3,pl4,layout=4,legends=false)
    savefig(fig_name*"_gₛ_weather.svg")


    Nₘ_f = model_result.Nₘ_f*100
    pl1 = plot(weatherts.PAR*10^6,Nₘ_f,xlabel="PAR (μmol s⁻¹ m⁻²)",ylabel="Nₘ_f (%)",seriestype=:scatter)
    pl2 = plot(weatherts.temp,Nₘ_f,xlabel="Temp (ᵒC)",ylabel="Nₘ_f (%)",seriestype=:scatter)
    pl3 = plot(weatherts.VPD/1000,Nₘ_f,xlabel="VPD (kPa)",ylabel="Nₘ_f (%)",seriestype=:scatter)
    pl4 = plot(CCPH.θₛ2ψₛ.(weatherts.θₛ),Nₘ_f,xlabel="SWP (MPa)",ylabel="Nₘ_f (%)",seriestype=:scatter)    
    
    plot(pl1,pl2,pl3,pl4,layout=4,legends=false)
    savefig(fig_name*"_Nₘ_f_weather.svg")
end