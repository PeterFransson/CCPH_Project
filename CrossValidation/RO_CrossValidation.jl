mutable struct ModelStat{T<:Float64}
    R²::T
    PearsonCor::T
    MAPE::T
    MSE::T
    RMSE::T
    R²_tot::T
end
function (K::ModelStat)(Stat_type::String)
    if Stat_type=="R²"
        return K.R²
    elseif Stat_type=="PearsonCor"
        return K.PearsonCor
    elseif Stat_type=="MAPE"
        return K.MAPE
    elseif Stat_type=="MSE"
        return K.MSE
    elseif Stat_type=="RMSE" 
        return K.RMSE
    elseif Stat_type=="R²_tot"
        return K.R²_tot
    else
        error("Stat_type either: \"R²\", \"PearsonCor\", \"MAPE\", \"MSE\", or \"RMSE\"")           
    end
end

mutable struct ModelStatSum{T<:ModelStat}
    GPP::T
    Ec::T
end
function (K::ModelStatSum)(Data_type::String)
    if Data_type=="GPP"
        return K.GPP
    elseif Data_type=="Ec"
        return K.Ec
    else
        error("Data_type either: \"GPP\" or \"Ec\"")
    end
end

function GetStatArray(modelstatsum::Array{ModelStatSum,1},
    Data_type::String,
    Stat_type::String)

    return [x(Data_type)(Stat_type) for x in modelstatsum]
end

function Get_Model_Stat(ParaDict::Dict{Symbol,Float64};
    stand_type::String="Fertilized",
    ind_train::Array{Int64,1}=collect(1:64))
    
    ind_val= setdiff(1:84, ind_train)
    RO_data = Load_RO_data()
    model,weatherts,data = Get_Result_RO_CCPH(ParaDict,RO_data;stand_type=stand_type)

    modeltrain_GPP = ModelStat(calcR²(data.GPP[ind_train],model.GPP[ind_train]),
    cor(data.GPP[ind_train],model.GPP[ind_train]),
    calcMAPE(data.GPP[ind_train],model.GPP[ind_train]),
    calcMSE(data.GPP[ind_train],model.GPP[ind_train]),
    calcRMSE(data.GPP[ind_train],model.GPP[ind_train]),
    calcR²(data.GPP,model.GPP))
    
    modelval_GPP = ModelStat(calcR²(data.GPP[ind_val],model.GPP[ind_val]),
    cor(data.GPP[ind_val],model.GPP[ind_val]),
    calcMAPE(data.GPP[ind_val],model.GPP[ind_val]),
    calcMSE(data.GPP[ind_val],model.GPP[ind_val]),
    calcRMSE(data.GPP[ind_val],model.GPP[ind_val]),
    calcR²(data.GPP,model.GPP))

    modeltrain_Ec = ModelStat(calcR²(data.Ec[ind_train],model.Ec[ind_train]),
    cor(data.Ec[ind_train],model.Ec[ind_train]),
    calcMAPE(data.Ec[ind_train],model.Ec[ind_train]),
    calcMSE(data.Ec[ind_train],model.Ec[ind_train]),
    calcRMSE(data.Ec[ind_train],model.Ec[ind_train]),
    calcR²(data.Ec,model.Ec))
    
    modelval_Ec = ModelStat(calcR²(data.Ec[ind_val],model.Ec[ind_val]),
    cor(data.Ec[ind_val],model.Ec[ind_val]),
    calcMAPE(data.Ec[ind_val],model.Ec[ind_val]),
    calcMSE(data.Ec[ind_val],model.Ec[ind_val]),
    calcRMSE(data.Ec[ind_val],model.Ec[ind_val]),
    calcR²(data.Ec,model.Ec))

    modeltrain = ModelStatSum(modeltrain_GPP,modeltrain_Ec)
    modelval = ModelStatSum(modelval_GPP,modelval_Ec)    

    return (modeltrain,modelval)
end

function Get_Model_Stat(file_name::String,
    ranges::Array{Tuple{Float64,Float64},1},
    para2ind::Dict{Symbol,Any};
    ParaDictInit_F::Union{Nothing,Dict{Symbol,Float64}}=nothing,
    ParaDictInit_C::Union{Nothing,Dict{Symbol,Float64}}=nothing)
    
    par = load("./output/"*file_name*".jld","xopt")   
    ind_train = load("./output/"*file_name*".jld","ind_train")   

    ParaDict_F,ParaDict_C = CreateParaDict(para2ind,par;
    ParaDictInit_F=ParaDictInit_F,ParaDictInit_C=ParaDictInit_C)
      
    modeltrain_F,modelval_F = Get_Model_Stat(ParaDict_F;stand_type="Fertilized",ind_train=ind_train)       
    modeltrain_C,modelval_C = Get_Model_Stat(ParaDict_C;stand_type="Control",ind_train=ind_train)

    return (modeltrain_F,modelval_F,modeltrain_C,modelval_C)
end

function Get_Model_Stat(file_name::String,
    ranges::Array{Tuple{Float64,Float64},1},
    parasym::Array{Symbol,1};
    ParaDictInit_F::Union{Nothing,Dict{Symbol,Float64}}=nothing,
    ParaDictInit_C::Union{Nothing,Dict{Symbol,Float64}}=nothing)

    file_name_F = file_name*"_F"
    file_name_C = file_name*"_C"
    ind_train = load("./output/"*file_name_F*".jld","ind_train")
    par_F = load("./output/"*file_name_F*".jld","xopt") 
    par_C = load("./output/"*file_name_C*".jld","xopt")

    ParaDict_F = CreateParaDict(parasym,par_F;ParaDictInit=ParaDictInit_F)
    ParaDict_C = CreateParaDict(parasym,par_C;ParaDictInit=ParaDictInit_C)
      
    modeltrain_F,modelval_F = Get_Model_Stat(ParaDict_F;stand_type="Fertilized",ind_train=ind_train)       
    modeltrain_C,modelval_C = Get_Model_Stat(ParaDict_C;stand_type="Control",ind_train=ind_train)

    return (modeltrain_F,modelval_F,modeltrain_C,modelval_C)
end

function writestat(io::IO,modeltrain::Array{ModelStatSum,1},
    modelval::Array{ModelStatSum,1},
    Data_type::String,
    Stat_type::String)

    println(io,"    Training: $(round.(GetStatArray(modeltrain,Data_type,Stat_type),digits=2))")
    println(io,"  Validation: $(round.(GetStatArray(modelval,Data_type,Stat_type),digits=2))")
    return nothing
end
function writestandstat(io::IO,
    modeltrain::Array{ModelStatSum,1},
    modelval::Array{ModelStatSum,1},
    Data_type::String)

    println(io,"R²:")
    writestat(io,modeltrain,modelval,Data_type,"R²")
    println(io,"MAPE:")
    writestat(io,modeltrain,modelval,Data_type,"MAPE")
    println(io,"RMSE:")
    writestat(io,modeltrain,modelval,Data_type,"RMSE")
    println(io,"Cor:")
    writestat(io,modeltrain,modelval,Data_type,"PearsonCor")
    println(io,"")
    println(io,"R² (tot): $(round.((GetStatArray(modeltrain,Data_type,"R²_tot")),digits=2))")
    print(io,"Mean RMSE: $(round(mean(GetStatArray(modeltrain,Data_type,"RMSE")),digits=2)) ")
    println(io,"($(round(mean(GetStatArray(modelval,Data_type,"RMSE")),digits=2)))")
    print(io,"Mean MAPE: $(round(mean(GetStatArray(modeltrain,Data_type,"MAPE")),digits=2)) ")
    println(io,"($(round(mean(GetStatArray(modelval,Data_type,"MAPE")),digits=2)))")
    return nothing
end

function writestat2file(file_name::String,
    modeltrain_F::Array{ModelStatSum,1},
    modelval_F::Array{ModelStatSum,1},
    modeltrain_C::Array{ModelStatSum,1},
    modelval_C::Array{ModelStatSum,1})

    open("./output/"*file_name*".txt","w") do io
        println(io,"--Statistics--")
        println(io,"Info: validation values in parentheses")
        println(io,"");println(io,"Fertilized")
        println(io,"");println(io,"GPP:")
        writestandstat(io,modeltrain_F,modelval_F,"GPP")
        println(io,"");println(io,"Ec:")
        writestandstat(io,modeltrain_F,modelval_F,"Ec")
        println(io,"");println(io,"");println(io,"Control")
        println(io,"");println(io,"GPP:")
        writestandstat(io,modeltrain_C,modelval_C,"GPP")
        println(io,"");println(io,"Ec:")
        writestandstat(io,modeltrain_C,modelval_C,"Ec")
    end
    return nothing
end

function CreateStatPlot(file_name::String,
    n_runs::Integer,
    modeltrain::Array{ModelStatSum,1},
    modelval::Array{ModelStatSum,1})

    sx = repeat(["Training", "Validation"], inner = n_runs)    
    nam = repeat(string.(1:n_runs), outer = 2)

    R²_GPP = append!(GetStatArray(modeltrain,"GPP","R²"),GetStatArray(modelval,"GPP","R²"))
    MAPE_GPP = append!(GetStatArray(modeltrain,"GPP","MAPE"),GetStatArray(modelval,"GPP","MAPE"))
    RMSE_GPP = append!(GetStatArray(modeltrain,"GPP","RMSE"),GetStatArray(modelval,"GPP","RMSE"))
    Cor_GPP = append!(GetStatArray(modeltrain,"GPP","PearsonCor"),GetStatArray(modelval,"GPP","PearsonCor"))
    
    R²_Ec = append!(GetStatArray(modeltrain,"Ec","R²"),GetStatArray(modelval,"Ec","R²"))
    MAPE_Ec = append!(GetStatArray(modeltrain,"Ec","MAPE"),GetStatArray(modelval,"Ec","MAPE"))
    RMSE_Ec = append!(GetStatArray(modeltrain,"Ec","RMSE"),GetStatArray(modelval,"Ec","RMSE"))
    Cor_Ec = append!(GetStatArray(modeltrain,"Ec","PearsonCor"),GetStatArray(modelval,"Ec","PearsonCor"))

    p_R²_GPP = groupedbar(nam, R²_GPP, group = sx, ylabel = "R²",title="GPP")
    p_MAPE_GPP = groupedbar(nam, MAPE_GPP, group = sx, ylabel = "MAPE")
    p_RMSE_GPP = groupedbar(nam, RMSE_GPP, group = sx, ylabel = "RMSE")  
    p_cor_GPP = groupedbar(nam, Cor_GPP, group = sx, ylabel = "Cor") 

    p_R²_Ec = groupedbar(nam, R²_Ec, group = sx, ylabel = "R²",title="Ec")
    p_MAPE_Ec = groupedbar(nam, MAPE_Ec, group = sx, ylabel = "MAPE")
    p_RMSE_Ec = groupedbar(nam, RMSE_Ec, group = sx, ylabel = "RMSE")  
    p_cor_Ec = groupedbar(nam, Cor_Ec, group = sx, ylabel = "Cor") 

    statplot = plot(p_R²_GPP,p_R²_Ec,p_MAPE_GPP,p_MAPE_Ec,
    p_RMSE_GPP,p_RMSE_Ec,p_cor_GPP,p_cor_Ec,layout=(4,2),legends=false)
    savefig(statplot,"./plots/"*file_name*".svg")
    return nothing
end
function CreateStatPlot(file_name::String,
    n_runs::Integer,
    modeltrain_F::Array{ModelStatSum,1},
    modelval_F::Array{ModelStatSum,1},
    modeltrain_C::Array{ModelStatSum,1},
    modelval_C::Array{ModelStatSum,1})
    
    CreateStatPlot(file_name*"_Stat_Fertilized",
    n_runs,
    modeltrain_F,
    modelval_F)

    CreateStatPlot(file_name*"_Stat_Control",
    n_runs,
    modeltrain_C,
    modelval_C)
    return nothing
end

function CreateTrainnValInd(Ind_sort::Array{Int64,1};
    val_prop::Float64=0.2)

    n_ind = length(Ind_sort)
    n_val = round(Int,n_ind*val_prop)
    ind_val = sort(shuffle(Ind_sort)[1:n_val])
    ind_train = setdiff(Ind_sort,ind_val)
    return (ind_train,ind_val)
end

function run_CrossValidation(file_name::String,
    n_runs::Integer,
    ranges::Array{Tuple{Float64, Float64},1},
    para2ind::Dict{Symbol,Any};
    PopulationSize::Integer=50,
    MaxSteps::Integer=10000,
    Calc_logP::Function=Calc_logP_GPP_Ec_Nm_f,
    ParaDictInit_F::Union{Nothing,Dict{Symbol,Float64}}=nothing,
    ParaDictInit_C::Union{Nothing,Dict{Symbol,Float64}}=nothing)
    
    modeltrain_F = Array{ModelStatSum,1}(undef,n_runs)
    modelval_F = Array{ModelStatSum,1}(undef,n_runs)
    modeltrain_C = Array{ModelStatSum,1}(undef,n_runs)
    modelval_C = Array{ModelStatSum,1}(undef,n_runs)

    xopt = Array{Array{Float64,1},1}(undef,n_runs)
    log_likelihood = Array{Float64,1}(undef,n_runs)
    ind_train_vec = Array{Array{Int64,1},1}(undef,n_runs)

    RO_data = Load_RO_data()

    isdir("./output/"*file_name)|| mkdir("./output/"*file_name)
    isdir("./plots/"*file_name) || mkdir("./plots/"*file_name)
    
    Threads.@threads for ix = 1:n_runs
        file_name_save = file_name*"/"*file_name*"_$(ix)"
        ind_train,ind_val = CreateTrainnValInd(collect(1:84))        

        ind_train_vec[ix] = ind_train
        xopt[ix],log_likelihood[ix] = run_opt(file_name_save,
        ranges,
        para2ind,
        RO_data;
        PopulationSize=PopulationSize,
        MaxSteps=MaxSteps,
        Calc_logP=Calc_logP,
        ParaDictInit_F=ParaDictInit_F,
        ParaDictInit_C=ParaDictInit_C,
        ind=ind_train)
    end
        
    for ix = 1:n_runs
        file_name_save = file_name*"/"*file_name*"_$(ix)"

        save_run_opt(file_name_save,
        xopt[ix],
        log_likelihood[ix],        
        para2ind;
        Calc_logP=Calc_logP,
        ParaDictInit_F=ParaDictInit_F,
        ParaDictInit_C=ParaDictInit_C,
        ind=ind_train_vec[ix])
        
        run_opt_par(file_name_save,
        ranges,
        para2ind;
        ParaDictInit_F=ParaDictInit_F,
        ParaDictInit_C=ParaDictInit_C)

        run_validation_RO(file_name_save,
        ranges,
        para2ind;
        ParaDictInit_F=ParaDictInit_F,
        ParaDictInit_C=ParaDictInit_C)

        modeltrain_F[ix],modelval_F[ix],modeltrain_C[ix],modelval_C[ix] = Get_Model_Stat(file_name_save,
        ranges,
        para2ind;
        ParaDictInit_F=ParaDictInit_F,
        ParaDictInit_C=ParaDictInit_C)
    end    

    writestat2file(file_name*"/"*file_name,
    modeltrain_F,
    modelval_F,
    modeltrain_C,
    modelval_C)    
    
    CreateStatPlot(file_name*"/"*file_name,
    n_runs,
    modeltrain_F,
    modelval_F,
    modeltrain_C,
    modelval_C)
end

function run_CrossValidation(file_name::String,
    n_runs::Integer,
    ranges::Array{Tuple{Float64, Float64},1},
    parasym::Array{Symbol,1};
    PopulationSize::Integer=50,
    MaxSteps::Integer=10000,
    Calc_logP::Function=Calc_logP_GPP_Ec_Nm_f,
    ParaDictInit_F::Union{Nothing,Dict{Symbol,Float64}}=nothing,
    ParaDictInit_C::Union{Nothing,Dict{Symbol,Float64}}=nothing)
    
    modeltrain_F = Array{ModelStatSum,1}(undef,n_runs)
    modelval_F = Array{ModelStatSum,1}(undef,n_runs)
    modeltrain_C = Array{ModelStatSum,1}(undef,n_runs)
    modelval_C = Array{ModelStatSum,1}(undef,n_runs)

    xopt_F = Array{Array{Float64,1},1}(undef,n_runs)
    log_likelihood_F = Array{Float64,1}(undef,n_runs)
    xopt_C = Array{Array{Float64,1},1}(undef,n_runs)
    log_likelihood_C = Array{Float64,1}(undef,n_runs)
    ind_train_vec = Array{Array{Int64,1},1}(undef,n_runs)

    RO_data = Load_RO_data()

    isdir("./output/"*file_name)|| mkdir("./output/"*file_name)
    isdir("./plots/"*file_name) || mkdir("./plots/"*file_name)
    #=
    Threads.@threads for ix = 1:n_runs
        file_name_save = file_name*"/"*file_name*"_$(ix)"
        ind_train,ind_val = CreateTrainnValInd(collect(1:84))        

        ind_train_vec[ix] = ind_train

        xopt_F[ix],log_likelihood_F[ix],xopt_C[ix],log_likelihood_C[ix] = run_opt(file_name_save,
        ranges,
        parasym,
        RO_data;
        PopulationSize=PopulationSize,
        MaxSteps=MaxSteps,
        Calc_logP=Calc_logP,
        ParaDictInit_F=ParaDictInit_F,
        ParaDictInit_C=ParaDictInit_C,
        ind=ind_train)        
    end
    =#    
    for ix = 1:n_runs
        file_name_save = file_name*"/"*file_name*"_$(ix)"        
        #=
        save_run_opt(file_name_save,
        xopt_F[ix],
        xopt_C[ix],
        log_likelihood_F[ix],
        log_likelihood_C[ix],
        parasym;
        Calc_logP=Calc_logP,
        ParaDictInit_F=ParaDictInit_F,
        ParaDictInit_C=ParaDictInit_C,
        ind=ind_train_vec[ix])
        =#
        
        run_opt_par(file_name_save,
        ranges,
        parasym;
        ParaDictInit_F=ParaDictInit_F,
        ParaDictInit_C=ParaDictInit_C)        

        run_validation_RO(file_name_save,
        ranges,
        parasym;
        ParaDictInit_F=ParaDictInit_F,
        ParaDictInit_C=ParaDictInit_C)

        modeltrain_F[ix],modelval_F[ix],modeltrain_C[ix],modelval_C[ix] = Get_Model_Stat(file_name_save,
        ranges,
        parasym;
        ParaDictInit_F=ParaDictInit_F,
        ParaDictInit_C=ParaDictInit_C)
    end    

    writestat2file(file_name*"/"*file_name,
    modeltrain_F,
    modelval_F,
    modeltrain_C,
    modelval_C)    
    
    CreateStatPlot(file_name*"/"*file_name,
    n_runs,
    modeltrain_F,
    modelval_F,
    modeltrain_C,
    modelval_C)
end