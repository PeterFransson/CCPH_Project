function run_validation_RO(ParaDict::Dict{Symbol,Float64},file_name::String;
    stand_type::String="Fertilized",ind::Array{Int64,1}=collect(65:84))
    if stand_type == "Fertilized"
        type = "_F"
    elseif stand_type == "Control"
        type = "_C"
    else
        error("Wrong input. Either \"Fertilized\" or \"Control\"")
    end 

    weather_file="Weather_RO"  

    RO_data = Load_RO_data(weather_file=weather_file)           
    
    model,weatherts,data = Get_Result_RO_CCPH(ParaDict,RO_data;stand_type=stand_type,
    sim_steps=84,weather_file=weather_file)        
    
    open("./output/"*file_name*type*"_validation.txt","w") do io
        println(io,"--Statistics--")        
        println(io,stand_type)        
        println(io,"  GPP: R² = $(calcR²(data.GPP[ind],model.GPP[ind])), Cor = $(cor(data.GPP[ind],model.GPP[ind]))")
        println(io,"            MAPE = $(calcMAPE(data.GPP[ind],model.GPP[ind])), MAPEₘₐₓ = $(calcMAPEₘₐₓ(data.GPP[ind],model.GPP[ind])),")
        println(io,"            MSE = $(calcMSE(data.GPP[ind],model.GPP[ind])), RMSE = $(calcRMSE(data.GPP[ind],model.GPP[ind]))")
        println(io,"  Ec: R² = $(calcR²(data.Ec[ind],model.Ec[ind])), Cor = $(cor(data.Ec[ind],model.Ec[ind])),")
        println(io,"            MAPE = $(calcMAPE(data.Ec[ind],model.Ec[ind])), MAPEₘₐₓ = $(calcMAPEₘₐₓ(data.Ec[ind],model.Ec[ind])),")
        println(io,"            MSE = $(calcMSE(data.Ec[ind],model.Ec[ind])), RMSE = $(calcRMSE(data.Ec[ind],model.Ec[ind]))")
        println(io,"");println(io,"")
        println(io,"--Parameters--")
        println(io,"Fertilized")
        for (key,value) in ParaDict            
            println(io,"  $(key): $(value)")                    
        end
    end

    a_GPP,b_GPP,a_Ec,b_Ec = ParaDict[:a_GPP],ParaDict[:b_GPP],ParaDict[:a_Ec],ParaDict[:b_Ec]
    
    σ_GPP = a_GPP.+model.GPP*b_GPP 
    σ_EC = a_Ec.+b_Ec*model.Ec
   
    c = -log(0.025*2) #Use for calculating 95% credible intervals
    
    plot(weatherts.date[ind],data.GPP[ind],label="Data",ylabel="GPP")
    plot!(weatherts.date[ind],model.GPP[ind]+c*σ_GPP[ind])
    plot!(weatherts.date[ind],model.GPP[ind]-c*σ_GPP[ind])
    pl1 = plot!(weatherts.date[ind],model.GPP[ind],label="Model")

    plot(weatherts.date[ind],data.Ec[ind],label="Data",ylabel="E_C")
    plot!(weatherts.date[ind],model.Ec[ind]+c*σ_EC[ind])
    plot!(weatherts.date[ind],model.Ec[ind]-c*σ_EC[ind])
    pl2 = plot!(weatherts.date[ind],model.Ec[ind],label="Model")
    
    plot(pl1,pl2,layout=(2,1),legends=false)
    savefig("./plots/"*file_name*type*"_validation.svg")
end

function run_validation_RO(file_name::String,ranges::Array{Tuple{Float64, Float64},1},
    parasym::Array{Symbol,1};ParaDictInit_F=nothing,ParaDictInit_C=nothing)
    file_name_F = file_name*"_F"
    file_name_C = file_name*"_C"

    par_F = load("./output/"*file_name_F*".jld","xopt")     
    par_C = load("./output/"*file_name_C*".jld","xopt") 
    ind_train = load("./output/"*file_name_F*".jld","ind_train")   
    ind = setdiff(1:84, ind_train)

    ParaDict_F = CreateParaDict(parasym,par_F;ParaDictInit=ParaDictInit_F)
    run_validation_RO(ParaDict_F,file_name;stand_type="Fertilized",ind=ind)

    ParaDict_C = CreateParaDict(parasym,par_C;ParaDictInit=ParaDictInit_C)
    run_validation_RO(ParaDict_C,file_name;stand_type="Control",ind=ind)   
end

function run_validation_RO(file_name::String,ranges::Array{Tuple{Float64, Float64},1},
    para2ind::Dict{Symbol,Any};ParaDictInit_F=nothing,ParaDictInit_C=nothing)

    par = load("./output/"*file_name*".jld","xopt")
    ind_train = load("./output/"*file_name*".jld","ind_train")   
    ind = setdiff(1:84, ind_train)

    ParaDict_F,ParaDict_C = CreateParaDict(para2ind,
    par;
    ParaDictInit_F=ParaDictInit_F,
    ParaDictInit_C=ParaDictInit_C)
    
    run_validation_RO(ParaDict_F,file_name;stand_type="Fertilized",ind=ind)    
    run_validation_RO(ParaDict_C,file_name;stand_type="Control",ind=ind)  
end