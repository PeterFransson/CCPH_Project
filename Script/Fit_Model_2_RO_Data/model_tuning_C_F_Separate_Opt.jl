function Calc_logP_GPP_Ec(model::ModelResult,data::RoData,ParaDict::Dict{Symbol,Float64}) 
    
    a_GPP,b_GPP,a_Ec,b_Ec = ParaDict[:a_GPP],ParaDict[:b_GPP],ParaDict[:a_Ec],ParaDict[:b_Ec] 

    logP = 0.0
    for i = 1:length(data.GPP)
        logP += calc_logP_term(data.GPP[i],model.GPP[i],a_GPP,b_GPP)
        logP += calc_logP_term(data.Ec[i],model.Ec[i],a_Ec,b_Ec)                      
    end    
     
    return -logP
end

function Calc_logP_GPP_Ec_Nm_f(model::ModelResult,data::RoData,ParaDict::Dict{Symbol,Float64}) 
    
    a_GPP,b_GPP,a_Ec,b_Ec = ParaDict[:a_GPP],ParaDict[:b_GPP],ParaDict[:a_Ec],ParaDict[:b_Ec] 
    μ_Nₘ_f,b_Nₘ_f = ParaDict[:μ_Nₘ_f],ParaDict[:b_Nₘ_f]

    logP = 0.0
    for i = 1:length(data.GPP)
        logP += calc_logP_term(data.GPP[i],model.GPP[i],a_GPP,b_GPP)
        logP += calc_logP_term(data.Ec[i],model.Ec[i],a_Ec,b_Ec)    
        logP += abs(model.Nₘ_f[i]-μ_Nₘ_f)/b_Nₘ_f          
    end    
     
    return -logP
end


function log_post_distri_RO_CCPH(para::Array{Float64,1},parasym::Array{Symbol,1},
    RO_data::RO_raw_data;stand_type::String="Fertilized",Calc_logP::Function=Calc_logP_GPP_Ec,
    ParaDictInit=nothing)
    try 
        ParaDict = CreateParaDict(parasym,para;ParaDictInit=ParaDictInit)

        model,weatherts,data = Get_Result_RO_CCPH(ParaDict,RO_data;stand_type)
            
        logP = Calc_logP(model,data,ParaDict) 
        return logP       
    catch err
        println("Parameters Error: ", err)
        post = -Inf
        return post
    end           
end

function run_opt(file_name::String,ranges::Array{Tuple{Float64, Float64},1},parasym::Array{Symbol,1};
    PopulationSize::Integer=50,MaxSteps::Integer=10000,Calc_logP::Function=Calc_logP_GPP_Ec,
    ParaDictInit_F=nothing,ParaDictInit_C=nothing)    
    file_name_F = file_name*"_F"
    file_name_C = file_name*"_C"

    RO_data = Load_RO_data()  

    #--Fertilized--
    res = BlackBoxOptim.bboptimize(x::Array{Float64,1}->
    -log_post_distri_RO_CCPH(x,parasym,RO_data;stand_type="Fertilized",
    Calc_logP=Calc_logP,ParaDictInit=ParaDictInit_F)
    ; SearchRange = ranges,PopulationSize = PopulationSize, MaxSteps = MaxSteps)

    xopt = BlackBoxOptim.best_candidate(res)
    log_likelihood = -BlackBoxOptim.best_fitness(res)  
    log_likelihood_fun = "$(Calc_logP)"
    save("./output/"*file_name_F*".jld","xopt",xopt,"log_likelihood",log_likelihood,
    "log_likelihood_fun",log_likelihood_fun)

    #--Control--
    res = BlackBoxOptim.bboptimize(x::Array{Float64,1}->
    -log_post_distri_RO_CCPH(x,parasym,RO_data;stand_type="Control",
    Calc_logP=Calc_logP,ParaDictInit=ParaDictInit_C)
    ; SearchRange = ranges,PopulationSize = PopulationSize, MaxSteps = MaxSteps)

    xopt = BlackBoxOptim.best_candidate(res)
    log_likelihood = -BlackBoxOptim.best_fitness(res)  
    log_likelihood_fun = "$(Calc_logP)"

    save("./output/"*file_name_C*".jld","xopt",xopt,"log_likelihood",log_likelihood,
    "log_likelihood_fun",log_likelihood_fun)
end

function run_opt_par(file_name::String,ranges::Array{Tuple{Float64, Float64},1},
    parasym::Array{Symbol,1};ParaDictInit_F=nothing,ParaDictInit_C=nothing)
    file_name_F = file_name*"_F"
    file_name_C = file_name*"_C"
    file_name_output = file_name*"_output"

    par_F = load("./output/"*file_name_F*".jld","xopt") 
    log_likelihood_F = load("./output/"*file_name_F*".jld","log_likelihood") 
    log_likelihood_fun_F = load("./output/"*file_name_F*".jld","log_likelihood_fun") 
    par_C = load("./output/"*file_name_C*".jld","xopt") 
    log_likelihood_C = load("./output/"*file_name_C*".jld","log_likelihood") 
    log_likelihood_fun_C = load("./output/"*file_name_C*".jld","log_likelihood_fun") 

    RO_data = Load_RO_data()       

    ParaDict_F = CreateParaDict(parasym,par_F;ParaDictInit=ParaDictInit_F)
    model_F,weatherts_F,data_F = Get_Result_RO_CCPH(ParaDict_F,RO_data;stand_type="Fertilized")

    ParaDict_C = CreateParaDict(parasym,par_C;ParaDictInit=ParaDictInit_C)
    model_C,weatherts_C,data_C = Get_Result_RO_CCPH(ParaDict_C,RO_data;stand_type="Control")

    open("./output/"*file_name_output*".txt","w") do io
        println(io,"--Statistics--")        
        println(io,"Fertilized")        
        println(io,"  GPP: R² = $(calcR²(data_F.GPP,model_F.GPP)), Cor = $(cor(data_F.GPP,model_F.GPP))")
        println(io,"  Ec: R² = $(calcR²(data_F.Ec,model_F.Ec)), Cor = $(cor(data_F.Ec,model_F.Ec))")
        println(io,"  Log likelihood: $(log_likelihood_F)")
        println(io,"  Log likelihood function selected: $(log_likelihood_fun_F)")
        println(io,"")
        println(io,"Control")        
        println(io,"  GPP: R² = $(calcR²(data_C.GPP,model_C.GPP)), Cor = $(cor(data_C.GPP,model_C.GPP))")
        println(io,"  Ec: R² = $(calcR²(data_C.Ec,model_C.Ec)), Cor = $(cor(data_C.Ec,model_C.Ec))")
        println(io,"  Log likelihood: $(log_likelihood_C)")
        println(io,"  Log likelihood function selected: $(log_likelihood_fun_C)")
        println(io,"");println(io,"")
        println(io,"--Parameters--")
        println(io,"Fertilized")
        for (key,value) in ParaDict_F
            if any(key.==parasym)
                ind = findfirst(x->x==key,parasym)
                println(io,"  $(key): $(value) $(ranges[ind])")
            else
                println(io,"  $(key): $(value)")
            end            
        end
        println(io,"")
        println(io,"Control") 
        for (key,value) in ParaDict_C
            if any(key.==parasym)
                ind = findfirst(x->x==key,parasym)
                println(io,"  $(key): $(value) $(ranges[ind])")
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

function run_calibration(file_name::String,ranges::Array{Tuple{Float64, Float64},1},parasym::Array{Symbol,1};
    PopulationSize::Integer=50,MaxSteps::Integer=10000,Calc_logP::Function=Calc_logP_GPP_Ec,
    ParaDictInit_F=nothing,ParaDictInit_C=nothing)
    run_opt(file_name,ranges,parasym;
    PopulationSize=PopulationSize,MaxSteps=MaxSteps,Calc_logP=Calc_logP,
    ParaDictInit_F=ParaDictInit_F,ParaDictInit_C=ParaDictInit_C)
    run_opt_par(file_name,ranges,parasym;ParaDictInit_F=ParaDictInit_F,ParaDictInit_C=ParaDictInit_C)
    return nothing
end

function work_list()
    
    parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:X0,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
    ranges = [(0.001,0.1),(10.0,40.0),(0.01,1.0),(0.001,0.1),
    (0.1,6.0),(0.1,0.5),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),
    (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
    file_name = "RO_Opt_Separate_GPP_EC_20220311" 
    run_calibration(file_name,ranges,parasym;
    PopulationSize=50,MaxSteps=100,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C)
end

work_list()