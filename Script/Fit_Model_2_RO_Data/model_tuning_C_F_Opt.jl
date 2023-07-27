function log_post_distri_RO_CCPH(para::Array{Float64,1},para2ind::Dict{Symbol,Any},
    RO_data::RO_raw_data;Calc_logP::Function=Calc_logP_GPP_Ec,
    ParaDictInit_F=nothing,ParaDictInit_C=nothing,ind::Array{Int64,1}=collect(1:84))
    try 
        ParaDict_F,ParaDict_C = CreateParaDict(para2ind,para;
        ParaDictInit_F=ParaDictInit_F,ParaDictInit_C=ParaDictInit_C)

        model_F,weatherts_F,data_F = Get_Result_RO_CCPH(ParaDict_F,RO_data;stand_type = "Fertilized")
        model_C,weatherts_C,data_C = Get_Result_RO_CCPH(ParaDict_C,RO_data;stand_type = "Control")
        
        logP_F = Calc_logP(model_F,data_F,ParaDict_F,ind=ind) 
        logP_C = Calc_logP(model_C,data_C,ParaDict_C,ind=ind)

        post = logP_F+logP_C
        return post      
    catch err
        println("Parameters Error: ", err)
        post = -Inf
        return post
    end           
end

function run_opt(file_name::String,
    ranges::Array{Tuple{Float64, Float64},1},
    para2ind::Dict{Symbol,Any},
    RO_data::RO_raw_data;
    PopulationSize::Integer=50,
    MaxSteps::Integer=10000,
    Calc_logP::Function=Calc_logP_GPP_Ec_Nm_f,
    ParaDictInit_F=nothing,
    ParaDictInit_C=nothing,
    ind::Array{Int64,1}=collect(1:84))   
    
    res = BlackBoxOptim.bboptimize(x::Array{Float64,1}->
    -log_post_distri_RO_CCPH(x,para2ind,RO_data;
    Calc_logP=Calc_logP,ParaDictInit_F=ParaDictInit_F,ParaDictInit_C=ParaDictInit_C,ind=ind)
    ; SearchRange = ranges,PopulationSize = PopulationSize, MaxSteps = MaxSteps)

    
    xopt = BlackBoxOptim.best_candidate(res)
    log_likelihood = -BlackBoxOptim.best_fitness(res)
     
    log_likelihood_fun = "$(Calc_logP)"
       
    save("./output/"*file_name*".jld",
    "xopt",xopt,
    "log_likelihood",log_likelihood,
    "log_likelihood_fun",log_likelihood_fun,
    "para2ind",para2ind,
    "ParaDictInit_F",ParaDictInit_F,
    "ParaDictInit_C",ParaDictInit_C) 
     
     return (xopt,log_likelihood)
end

function save_run_opt(file_name::String,
    xopt::Array{Float64},
    log_likelihood::Float64,
    para2ind::Dict{Symbol,Any};
    Calc_logP::Function=Calc_logP_GPP_Ec_Nm_f,
    ParaDictInit_F=nothing,
    ParaDictInit_C=nothing,
    ind::Array{Int64,1}=collect(1:84))

    log_likelihood_fun = "$(Calc_logP)"

    save("./output/"*file_name*".jld",
    "xopt",xopt,
    "log_likelihood",log_likelihood,
    "log_likelihood_fun",log_likelihood_fun,
    "para2ind",para2ind,
    "ParaDictInit_F",ParaDictInit_F,
    "ParaDictInit_C",ParaDictInit_C,
    "ind_train",ind) 
    return nothing
end

function run_opt_par(file_name::String,ranges::Array{Tuple{Float64, Float64},1},
    para2ind::Dict{Symbol,Any};ParaDictInit_F=nothing,ParaDictInit_C=nothing)
    
    file_name_output = file_name*"_output"

    par = load("./output/"*file_name*".jld","xopt") 
    ind = load("./output/"*file_name*".jld","ind_train")  
    log_likelihood = load("./output/"*file_name*".jld","log_likelihood") 
    log_likelihood_fun = load("./output/"*file_name*".jld","log_likelihood_fun") 

    RO_data = Load_RO_data()

    ParaDict_F,ParaDict_C = CreateParaDict(para2ind,par;
    ParaDictInit_F=ParaDictInit_F,ParaDictInit_C=ParaDictInit_C)
    
    model_F,weatherts_F,data_F = Get_Result_RO_CCPH(ParaDict_F,RO_data;stand_type="Fertilized")
    
    model_C,weatherts_C,data_C = Get_Result_RO_CCPH(ParaDict_C,RO_data;stand_type="Control")

    open("./output/"*file_name_output*".txt","w") do io
        println(io,"--Statistics--")        
        println(io,"Fertilized")        
        println(io,"  GPP: R² = $(calcR²(data_F.GPP[ind],model_F.GPP[ind])), Cor = $(cor(data_F.GPP[ind],model_F.GPP[ind]))")
        println(io,"            MAPE = $(calcMAPE(data_F.GPP[ind],model_F.GPP[ind])), MAPEₘₐₓ = $(calcMAPEₘₐₓ(data_F.GPP[ind],model_F.GPP[ind])),")
        println(io,"            MSE = $(calcMSE(data_F.GPP[ind],model_F.GPP[ind])), RMSE = $(calcRMSE(data_F.GPP[ind],model_F.GPP[ind]))")
        println(io,"  Ec: R² = $(calcR²(data_F.Ec[ind],model_F.Ec[ind])), Cor = $(cor(data_F.Ec[ind],model_F.Ec[ind]))")
        println(io,"            MAPE = $(calcMAPE(data_F.Ec[ind],model_F.Ec[ind])), MAPEₘₐₓ = $(calcMAPEₘₐₓ(data_F.Ec[ind],model_F.Ec[ind])),")
        println(io,"            MSE = $(calcMSE(data_F.Ec[ind],model_F.Ec[ind])), RMSE = $(calcRMSE(data_F.Ec[ind],model_F.Ec[ind]))")
        println(io,"  Log likelihood: $(log_likelihood)")
        println(io,"  Log likelihood function selected: $(log_likelihood_fun)")
        println(io,"")
        println(io,"Control")        
        println(io,"  GPP: R² = $(calcR²(data_C.GPP[ind],model_C.GPP[ind])), Cor = $(cor(data_C.GPP[ind],model_C.GPP[ind]))")
        println(io,"            MAPE = $(calcMAPE(data_C.GPP[ind],model_C.GPP[ind])), MAPEₘₐₓ = $(calcMAPEₘₐₓ(data_C.GPP[ind],model_C.GPP[ind])),")
        println(io,"            MSE = $(calcMSE(data_C.GPP[ind],model_C.GPP[ind])), RMSE = $(calcRMSE(data_C.GPP[ind],model_C.GPP[ind]))")
        println(io,"  Ec: R² = $(calcR²(data_C.Ec[ind],model_C.Ec[ind])), Cor = $(cor(data_C.Ec[ind],model_C.Ec[ind]))")
        println(io,"            MAPE = $(calcMAPE(data_C.Ec[ind],model_C.Ec[ind])), MAPEₘₐₓ = $(calcMAPEₘₐₓ(data_C.Ec[ind],model_C.Ec[ind])),")
        println(io,"            MSE = $(calcMSE(data_C.Ec[ind],model_C.Ec[ind])), RMSE = $(calcRMSE(data_C.Ec[ind],model_C.Ec[ind]))")
        println(io,"  Log likelihood: $(log_likelihood)")
        println(io,"  Log likelihood function selected: $(log_likelihood_fun)")
        println(io,"");println(io,"")
        println(io,"--Parameters--")
        println(io,"Fertilized")
        for (key,value) in ParaDict_F
            if haskey(para2ind,key) 
                if isa(para2ind[key],Integer)
                    indx = para2ind[key] 
                    println(io,"  $(key): $(value) $(ranges[indx])")
                else
                    indx = para2ind[key][1]
                    println(io,"  $(key): $(value) $(ranges[indx]), Separate")
                end               
            else
                println(io,"  $(key): $(value)")
            end            
        end
        println(io,"")
        println(io,"Control") 
        for (key,value) in ParaDict_C
            if haskey(para2ind,key)
                if isa(para2ind[key],Integer)
                    indx = para2ind[key] 
                    println(io,"  $(key): $(value) $(ranges[indx])")
                else
                    indx = para2ind[key][2]
                    println(io,"  $(key): $(value) $(ranges[indx]), Separate")
                end                
            else
                println(io,"  $(key): $(value)")
            end            
        end
    end    
    
    #=
    a_GPP_F,b_GPP_F,a_Ec_F,b_Ec_F = ParaDict_F[:a_GPP],ParaDict_F[:b_GPP],ParaDict_F[:a_Ec],ParaDict_F[:b_Ec]
    a_GPP_C,b_GPP_C,a_Ec_C,b_Ec_C = ParaDict_C[:a_GPP],ParaDict_C[:b_GPP],ParaDict_C[:a_Ec],ParaDict_C[:b_Ec]

    σ_GPP_F = a_GPP_F.+model_F.GPP*b_GPP_F
    σ_GPP_C = a_GPP_C.+model_C.GPP*b_GPP_C

    σ_EC_F = a_Ec_F.+b_Ec_F*model_F.Ec
    σ_EC_C = a_Ec_C.+b_Ec_C*model_C.Ec

    c = -log(0.025*2) #Use for calculating 95% credible intervals

    plot(weatherts_F.date[ind],data_F.GPP[ind],label="Data",ylabel="GPP")
    plot!(weatherts_F.date[ind],model_F.GPP[ind]+c*σ_GPP_F[ind])
    plot!(weatherts_F.date[ind],model_F.GPP[ind]-c*σ_GPP_F[ind])
    pl1 = plot!(weatherts_F.date[ind],model_F.GPP[ind],label="Model")
    plot(weatherts_C.date[ind],data_C.GPP[ind],label="Data",ylabel="GPP")
    plot!(weatherts_C.date[ind],model_C.GPP[ind]+c*σ_GPP_C[ind])
    plot!(weatherts_C.date[ind],model_C.GPP[ind]-c*σ_GPP_C[ind])
    pl2 = plot!(weatherts_C.date[ind],model_C.GPP[ind],label="Model")

    plot(weatherts_F.date[ind],data_F.Ec[ind],label="Data",ylabel="E_C")
    plot!(weatherts_F.date[ind],model_F.Ec[ind]+c*σ_EC_F[ind])
    plot!(weatherts_F.date[ind],model_F.Ec[ind]-c*σ_EC_F[ind])
    pl3 = plot!(weatherts_F.date[ind],model_F.Ec[ind],label="Model")
    plot(weatherts_C.date[ind],data_C.Ec[ind],label="Data",ylabel="E_C")
    plot!(weatherts_C.date[ind],model_C.Ec[ind]+c*σ_EC_C[ind])
    plot!(weatherts_C.date[ind],model_C.Ec[ind]-c*σ_EC_C[ind])
    pl4 = plot!(weatherts_C.date[ind],model_C.Ec[ind],label="Model")

    plot(pl1,pl2,pl3,pl4,layout=(2,2),legends=false)
    savefig("./plots/"*file_name*"_result.svg")

    plot(weatherts_F.date[ind],data_F.gₛ[ind],label="Data",ylabel="gₛ")
    pl1 = plot!(weatherts_F.date[ind],model_F.gₛ[ind],label="Model")
    plot(weatherts_C.date[ind],data_C.gₛ[ind],label="Data",ylabel="gₛ")
    pl2 = plot!(weatherts_C.date[ind],model_C.gₛ[ind],label="Model")

    plot(pl1,pl2,layout=(1,2),legends=false)
    savefig("./plots/"*file_name*"_result_gs.svg")

    pl1 = plot(weatherts_F.date[ind],model_F.Nₘ_f[ind]*100,label="Model",ylabel="Nₘ_f (%)")    
    pl2 = plot(weatherts_C.date[ind],model_C.Nₘ_f[ind]*100,label="Model",ylabel="Nₘ_f (%)")
    
    plot(pl1,pl2,layout=(1,2),legends=false)
    savefig("./plots/"*file_name*"_result_Nm_f.svg")

    ind_val = setdiff(1:84, ind)
    plot([0.0,15.0],[0.0,15.0],xlabel="Data GPP",ylabel="Model GPP")
    plot!(data_F.GPP[ind],model_F.GPP[ind],seriestype=:scatter)
    pl1 = plot!(data_F.GPP[ind_val],model_F.GPP[ind_val],seriestype=:scatter)
    plot([0.0,15.0],[0.0,15.0],xlabel="Data GPP",ylabel="Model GPP")
    plot!(data_C.GPP[ind],model_C.GPP[ind],seriestype=:scatter)
    pl2 = plot!(data_C.GPP[ind_val],model_C.GPP[ind_val],seriestype=:scatter)

    plot([0.0,2.0],[0.0,2.0],xlabel="Data Ec",ylabel="Model Ec")
    plot!(data_F.Ec[ind],model_F.Ec[ind],seriestype=:scatter)
    pl3 = plot!(data_F.Ec[ind_val],model_F.Ec[ind_val],seriestype=:scatter)
    plot([0.0,2.0],[0.0,2.0],xlabel="Data Ec",ylabel="Model Ec")
    plot!(data_C.Ec[ind],model_C.Ec[ind],seriestype=:scatter)
    pl4 = plot!(data_C.Ec[ind_val],model_C.Ec[ind_val],seriestype=:scatter)

    plot(pl1,pl2,pl3,pl4,layout=(2,2),legends=false)    
    savefig("./plots/"*file_name*"_annual_error.svg")

    plot([0.0,15.0],[0.0,0.0],xlabel="Data GPP",ylabel="Res GPP")
    plot!(data_F.GPP[ind],model_F.GPP[ind]-data_F.GPP[ind],seriestype=:scatter)
    pl1 = plot!(data_F.GPP[ind_val],model_F.GPP[ind_val]-data_F.GPP[ind_val],seriestype=:scatter)
    plot([0.0,15.0],[0.0,0.0],xlabel="Data GPP",ylabel="Res GPP")
    plot!(data_C.GPP[ind],model_C.GPP[ind]-data_C.GPP[ind],seriestype=:scatter)
    pl2 = plot!(data_C.GPP[ind_val],model_C.GPP[ind_val]-data_C.GPP[ind_val],seriestype=:scatter)

    plot([0.0,2.0],[0.0,0.0],xlabel="Data Ec",ylabel="Res Ec")
    plot!(data_F.Ec[ind],model_F.Ec[ind]-data_F.Ec[ind],seriestype=:scatter)
    pl3 = plot!(data_F.Ec[ind_val],model_F.Ec[ind_val]-data_F.Ec[ind_val],seriestype=:scatter)
    plot([0.0,2.0],[0.0,0.0],xlabel="Data Ec",ylabel="Res Ec")
    plot!(data_C.Ec[ind],model_C.Ec[ind]-data_C.Ec[ind],seriestype=:scatter)
    pl4 = plot!(data_C.Ec[ind_val],model_C.Ec[ind_val]-data_C.Ec[ind_val],seriestype=:scatter)

    plot(pl1,pl2,pl3,pl4,layout=(2,2),legends=false)    
    savefig("./plots/"*file_name*"_residuals.svg")

    create_wue_plot("./plots/"*file_name*"_wue_F",model_F,weatherts_F)  
    create_wue_plot("./plots/"*file_name*"_wue_C",model_C,weatherts_C) 

    create_trait_plots("./plots/"*file_name*"_trait_F",model_F,weatherts_F)
    create_trait_plots("./plots/"*file_name*"_trait_C",model_C,weatherts_C)
    =#
end

function run_calibration(file_name::String,
    ranges::Array{Tuple{Float64, Float64},1},
    para2ind::Dict{Symbol,Any},
    RO_data::RO_raw_data;
    PopulationSize::Integer=50,
    MaxSteps::Integer=10000,
    Calc_logP::Function=Calc_logP_GPP_Ec,
    ParaDictInit_F=nothing,
    ParaDictInit_C=nothing,
    ind::Array{Int64,1}=collect(1:84))

    run_opt(file_name,
    ranges,
    para2ind,
    RO_data;
    PopulationSize=PopulationSize,
    MaxSteps=MaxSteps,
    Calc_logP=Calc_logP,
    ParaDictInit_F=ParaDictInit_F,
    ParaDictInit_C=ParaDictInit_C,
    ind=ind)

    run_opt_par(file_name,ranges,para2ind;ParaDictInit_F=ParaDictInit_F,ParaDictInit_C=ParaDictInit_C)
    return nothing
end