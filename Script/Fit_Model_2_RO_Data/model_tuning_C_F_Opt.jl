function log_post_distri_RO_CCPH(para::Array{Float64,1},para2ind::Dict{Symbol,Any},
    RO_data::RO_raw_data;Calc_logP::Function=Calc_logP_GPP_Ec,
    ParaDictInit_F=nothing,ParaDictInit_C=nothing)
    try 
        ParaDict_F,ParaDict_C = CreateParaDict(para2ind,para;
        ParaDictInit_F=ParaDictInit_F,ParaDictInit_C=ParaDictInit_C)

        model_F,weatherts_F,data_F = Get_Result_RO_CCPH(ParaDict_F,RO_data;stand_type = "Fertilized")
        model_C,weatherts_C,data_C = Get_Result_RO_CCPH(ParaDict_C,RO_data;stand_type = "Control")
        
        logP_F = Calc_logP(model_F,data_F,ParaDict_F) 
        logP_C = Calc_logP(model_C,data_C,ParaDict_C)

        post = logP_F+logP_C
        return post      
    catch err
        println("Parameters Error: ", err)
        post = -Inf
        return post
    end           
end

function run_opt(file_name::String,ranges::Array{Tuple{Float64, Float64},1},para2ind::Dict{Symbol,Any};
    PopulationSize::Integer=50,MaxSteps::Integer=10000,Calc_logP::Function=Calc_logP_GPP_Ec,
    ParaDictInit_F=nothing,ParaDictInit_C=nothing)    
    
    RO_data = Load_RO_data()
    
    res = BlackBoxOptim.bboptimize(x::Array{Float64,1}->
    -log_post_distri_RO_CCPH(x,para2ind,RO_data;
    Calc_logP=Calc_logP,ParaDictInit_F=ParaDictInit_F,ParaDictInit_C=ParaDictInit_C)
    ; SearchRange = ranges,PopulationSize = PopulationSize, MaxSteps = MaxSteps)

    xopt = BlackBoxOptim.best_candidate(res)
    log_likelihood = -BlackBoxOptim.best_fitness(res)  
    log_likelihood_fun = "$(Calc_logP)"
    save("./output/"*file_name*".jld","xopt",xopt,"log_likelihood",log_likelihood,
    "log_likelihood_fun",log_likelihood_fun) 
end

function run_opt_par(file_name::String,ranges::Array{Tuple{Float64, Float64},1},
    para2ind::Dict{Symbol,Any};ParaDictInit_F=nothing,ParaDictInit_C=nothing)
    
    file_name_output = file_name*"_output"

    par = load("./output/"*file_name*".jld","xopt") 
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
        println(io,"  GPP: R² = $(calcR²(data_F.GPP,model_F.GPP)), Cor = $(cor(data_F.GPP,model_F.GPP))")
        println(io,"  Ec: R² = $(calcR²(data_F.Ec,model_F.Ec)), Cor = $(cor(data_F.Ec,model_F.Ec))")
        println(io,"  Log likelihood: $(log_likelihood)")
        println(io,"  Log likelihood function selected: $(log_likelihood_fun)")
        println(io,"")
        println(io,"Control")        
        println(io,"  GPP: R² = $(calcR²(data_C.GPP,model_C.GPP)), Cor = $(cor(data_C.GPP,model_C.GPP))")
        println(io,"  Ec: R² = $(calcR²(data_C.Ec,model_C.Ec)), Cor = $(cor(data_C.Ec,model_C.Ec))")
        println(io,"  Log likelihood: $(log_likelihood)")
        println(io,"  Log likelihood function selected: $(log_likelihood_fun)")
        println(io,"");println(io,"")
        println(io,"--Parameters--")
        println(io,"Fertilized")
        for (key,value) in ParaDict_F
            if haskey(para2ind,key) 
                if isa(para2ind[key],Integer)
                    ind = para2ind[key] 
                    println(io,"  $(key): $(value) $(ranges[ind])")
                else
                    ind = para2ind[key][1]
                    println(io,"  $(key): $(value) $(ranges[ind]), Separate")
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
                    ind = para2ind[key] 
                    println(io,"  $(key): $(value) $(ranges[ind])")
                else
                    ind = para2ind[key][2]
                    println(io,"  $(key): $(value) $(ranges[ind]), Separate")
                end                
            else
                println(io,"  $(key): $(value)")
            end            
        end
    end    
    
    a_GPP_F,b_GPP_F,a_Ec_F,b_Ec_F = ParaDict_F[:a_GPP],ParaDict_F[:b_GPP],ParaDict_F[:a_Ec],ParaDict_F[:b_Ec]
    a_GPP_C,b_GPP_C,a_Ec_C,b_Ec_C = ParaDict_C[:a_GPP],ParaDict_C[:b_GPP],ParaDict_C[:a_Ec],ParaDict_C[:b_Ec]

    σ_GPP_F = a_GPP_F.+model_F.GPP*b_GPP_F
    σ_GPP_C = a_GPP_C.+model_C.GPP*b_GPP_C

    σ_EC_F = a_Ec_F.+b_Ec_F*model_F.Ec
    σ_EC_C = a_Ec_C.+b_Ec_C*model_C.Ec

    c = -log(0.025*2) #Use for calculating 95% credible intervals

    plot(weatherts_F.date,data_F.GPP,label="Data",ylabel="GPP")
    plot!(weatherts_F.date,model_F.GPP+c*σ_GPP_F)
    plot!(weatherts_F.date,model_F.GPP-c*σ_GPP_F)
    pl1 = plot!(weatherts_F.date,model_F.GPP,label="Model")
    plot(weatherts_C.date,data_C.GPP,label="Data",ylabel="GPP")
    plot!(weatherts_C.date,model_C.GPP+c*σ_GPP_C)
    plot!(weatherts_C.date,model_C.GPP-c*σ_GPP_C)
    pl2 = plot!(weatherts_C.date,model_C.GPP,label="Model")

    plot(weatherts_F.date,data_F.Ec,label="Data",ylabel="E_C")
    plot!(weatherts_F.date,model_F.Ec+c*σ_EC_F)
    plot!(weatherts_F.date,model_F.Ec-c*σ_EC_F)
    pl3 = plot!(weatherts_F.date,model_F.Ec,label="Model")
    plot(weatherts_C.date,data_C.Ec,label="Data",ylabel="E_C")
    plot!(weatherts_C.date,model_C.Ec+c*σ_EC_C)
    plot!(weatherts_C.date,model_C.Ec-c*σ_EC_C)
    pl4 = plot!(weatherts_C.date,model_C.Ec,label="Model")

    plot(pl1,pl2,pl3,pl4,layout=(2,2),legends=false)
    savefig("./plots/"*file_name*"_result.svg")

    plot(weatherts_F.date,data_F.gₛ,label="Data",ylabel="gₛ")
    pl1 = plot!(weatherts_F.date,model_F.gₛ,label="Model")
    plot(weatherts_C.date,data_C.gₛ,label="Data",ylabel="gₛ")
    pl2 = plot!(weatherts_C.date,model_C.gₛ,label="Model")

    plot(pl1,pl2,layout=(1,2),legends=false)
    savefig("./plots/"*file_name*"_result_gs.svg")

    pl1 = plot(weatherts_F.date,model_F.Nₘ_f*100,label="Model",ylabel="Nₘ_f (%)")    
    pl2 = plot(weatherts_C.date,model_C.Nₘ_f*100,label="Model",ylabel="Nₘ_f (%)")
    
    plot(pl1,pl2,layout=(1,2),legends=false)
    savefig("./plots/"*file_name*"_result_Nm_f.svg")
end

function run_calibration(file_name::String,ranges::Array{Tuple{Float64, Float64},1},para2ind::Dict{Symbol,Any};
    PopulationSize::Integer=50,MaxSteps::Integer=10000,Calc_logP::Function=Calc_logP_GPP_Ec,
    ParaDictInit_F=nothing,ParaDictInit_C=nothing)
    run_opt(file_name,ranges,para2ind;
    PopulationSize=PopulationSize,MaxSteps=MaxSteps,Calc_logP=Calc_logP,
    ParaDictInit_F=ParaDictInit_F,ParaDictInit_C=ParaDictInit_C)
    run_opt_par(file_name,ranges,para2ind;ParaDictInit_F=ParaDictInit_F,ParaDictInit_C=ParaDictInit_C)
    return nothing
end