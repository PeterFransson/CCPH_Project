#Calculate water use efficiency statistics
function create_wue_stat(file_name::String,model_result::ModelResult,weatherts::WeatherTS)    
    wue = model_result.GPP./model_result.Ec*10 #kg C ha⁻¹ mm⁻¹
    wue_mol = model_result.A./model_result.E*1000 #mmol C mol⁻¹ H₂O     

    open(file_name*".txt","w") do io
        println(io,"--Water use efficiency (wue) statistics--") 
        println(io,"cor(PAR,wue) = $(cor(weatherts.PAR*10^6,wue_mol))")
        println(io,"cor(temp,wue) = $(cor(weatherts.temp,wue_mol))")
        println(io,"cor(VPD,wue) = $(cor(weatherts.VPD/1000,wue_mol))")
        println(io,"cor(θₛ,wue) = $(cor(weatherts.θₛ*100,wue_mol))")
    end
end

#Optimal traits stats
function create_trait_stat(file_name::String,model_result::ModelResult,weatherts::WeatherTS)
    Nₘ_f = model_result.Nₘ_f*100
    gₛ = model_result.gₛ

    open(file_name*".txt","w") do io
        println(io,"--Optimal traits  (Nₘ_f & gₛ) statistics--")
        println(io,"--Nₘ_f--")
        println(io,"cor(PAR,Nₘ_f) = $(cor(weatherts.PAR*10^6,Nₘ_f))")
        println(io,"cor(temp,Nₘ_f) = $(cor(weatherts.temp,Nₘ_f))")
        println(io,"cor(VPD,Nₘ_f) = $(cor(weatherts.VPD/1000,Nₘ_f))")
        println(io,"cor(θₛ,Nₘ_f) = $(cor(weatherts.θₛ*100,Nₘ_f))")
        println(io,"--gₛ--")
        println(io,"cor(PAR,gₛ) = $(cor(weatherts.PAR*10^6,gₛ))")
        println(io,"cor(temp,gₛ) = $(cor(weatherts.temp,gₛ))")
        println(io,"cor(VPD,gₛ) = $(cor(weatherts.VPD/1000,gₛ))")
        println(io,"cor(θₛ,gₛ) = $(cor(weatherts.θₛ*100,gₛ))")
        println(io,"--gₛ_&_Nₘ_f--")
        println(io,"cor(Nₘ_f,gₛ) = $(cor(Nₘ_f,gₛ))")
    end
end