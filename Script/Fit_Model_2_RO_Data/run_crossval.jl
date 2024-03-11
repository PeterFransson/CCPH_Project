function CreateTrainValSet(raw_input::Vector{RawInputData};
    val_prop::Float64=0.2)

    n_years = length(raw_input)
    @show n_weeks = [length(input.growth_indices_weekly) for input in raw_input]    

    train_set = [ones(Bool,week*7) for week in n_weeks]
    val_set = [zeros(Bool,week*7) for week in n_weeks]

    train_set_weekly = [ones(Bool,week) for week in n_weeks]
    val_set_weekly = [zeros(Bool,week) for week in n_weeks]

    total_set = [(year_idx,week_idx) for year_idx in 1:n_years for week_idx in 1:n_weeks[year_idx]]
  
    sum(n_weeks)==length(total_set)||error("sum(n_weeks)!=length(total_set)")

    n_val = round(Int,sum(n_weeks)*val_prop)
    val_idx = sort(shuffle(total_set)[1:n_val])    

    for idx in val_idx 

        start_idx,end_idx = raw_input[idx[1]].growth_indices_weekly[idx[2]]

        train_set_weekly[idx[1]][idx[2]] = false
        val_set_weekly[idx[1]][idx[2]] = true

        train_set[idx[1]][start_idx:end_idx] .= false
        val_set[idx[1]][start_idx:end_idx] .= true
    end
    
    return (train_set,val_set,train_set_weekly,val_set_weekly)
end

#Non-shared parameters
function opt_par_obj(x::Vector{T},
    raw_input::Vector{RawInputData},
    GPP_data::Vector{Vector{R}},
    Ec_data::Vector{Vector{S}},
    train_set::Vector{Vector{Bool}};
    stand_type::Symbol=:Fertilized,
    weight_GPP::Real=1.0) where{T<:Real,R<:Real,S<:Real}

    try
        par = ModelPar(x;stand_type=stand_type)
        GPP_model,Ec_model,Nₘ_f_model = run_model(par,raw_input;stand_type=stand_type)

        GPP_data_train = [GPP_data[i][train_set[i]] for i in 1:4] 
        Ec_data_train = [Ec_data[i][train_set[i]] for i in 1:4] 

        GPP_model_train = [GPP_model[i][train_set[i]] for i in 1:4]
        Ec_model_train = [Ec_model[i][train_set[i]] for i in 1:4]
        Nₘ_f_model_train = [Nₘ_f_model[i][train_set[i]] for i in 1:4]

        return -Calc_logP_GPP_Ec_Nm_f(GPP_model_train,
        Ec_model_train,
        Nₘ_f_model_train,
        GPP_data_train,
        Ec_data_train,
        par,
        raw_input;
        weight_GPP=weight_GPP)
    catch err
        println("Parameters Error: ", err)

        return Inf
    end 
end

#Non-shared parameters
function run_crossval_opt(folder_name::String,
    stand_type::Symbol,
    weight_GPP::Real)
    
    filename = "result"

    isdir("./output/"*folder_name)|| mkdir("./output/"*folder_name)

    #stand_type=:Control
    raw_input = RawInputData(;stand_type=stand_type)
    Ec_data = calc_Ec_data.(raw_input)
    GPP_data = get_GPP_data.(raw_input;stand_type=stand_type)    

    (train_set,val_set,train_set_weekly,val_set_weekly) = CreateTrainValSet(raw_input) 
    
    #Nₛ,α_max,a_Jmax,Kₓₗ₀,τ,ΔS,a_GPP,b_GPP = x    
    range = [(0.0001,0.1),
    (0.1,0.5),
    (0.001,1.0),
    (0.0004,0.1),
    (1.0,15.0),
    (10.0,25.0),
    (0.0001,5.0),
    (0.0001,3.0)]      

    res = BlackBoxOptim.bboptimize(x->opt_par_obj(x,
    raw_input,
    GPP_data,
    Ec_data,
    train_set;
    stand_type=stand_type,
    weight_GPP=weight_GPP); SearchRange = range,NThreads=Threads.nthreads()-1)
    x_opt = BlackBoxOptim.best_candidate(res) 
    par_opt = ModelPar(x_opt;stand_type=stand_type)
    
    JLD.save("output/"*folder_name*"/"*filename*".jld",
    "x_opt",x_opt,
    "par_opt",par_opt,
    "stand_type",stand_type,
    "weight_GPP",weight_GPP,
    "train_set",train_set,
    "val_set",val_set,
    "train_set_weekly",train_set_weekly,
    "val_set_weekly",val_set_weekly)
    return nothing
end  

#Shared parameters
function opt_par_obj(x::Vector{T},
    raw_input_F::Vector{RawInputData},
    GPP_data_F::Vector{Vector{R}},
    Ec_data_F::Vector{Vector{S}},
    train_set_F::Vector{Vector{Bool}},    
    raw_input_C::Vector{RawInputData},
    GPP_data_C::Vector{Vector{R}},
    Ec_data_C::Vector{Vector{S}},
    train_set_C::Vector{Vector{Bool}};
    weight_GPP::Real=1.0) where{T<:Real,R<:Real,S<:Real}

    x_F = x[1:8]
    x_C = x[[9,2,3,10,5,6,7,8]]

    try
        obj_F = opt_par_obj(x_F,raw_input_F,GPP_data_F,Ec_data_F,train_set_F;stand_type=:Fertilized,weight_GPP=weight_GPP)
        obj_C = opt_par_obj(x_C,raw_input_C,GPP_data_C,Ec_data_C,train_set_C;stand_type=:Control,weight_GPP=weight_GPP)

        return obj_F+obj_C
    catch err
        println("Parameters Error: ", err)

        return Inf
    end
end 


#Shared parameters
function run_crossval_opt(folder_name::String,    
    weight_GPP::Real;
    MaxFuncEvals::Integer=10000)
    
    filename = "result"

    isdir("./output/"*folder_name)|| mkdir("./output/"*folder_name)

    stand_type_F = :Fertilized
    stand_type_C = :Control
    
    raw_input_F = RawInputData(;stand_type=stand_type_F)
    Ec_data_F = calc_Ec_data.(raw_input_F)
    GPP_data_F = get_GPP_data.(raw_input_F;stand_type=stand_type_F)    

    raw_input_C = RawInputData(;stand_type=stand_type_C)
    Ec_data_C = calc_Ec_data.(raw_input_C)
    GPP_data_C = get_GPP_data.(raw_input_C;stand_type=stand_type_C)    

    (train_set,val_set,train_set_weekly,val_set_weekly) = CreateTrainValSet(raw_input_F) 
    
    #Shared paramters
    #Nₛ_F,α_max,a_Jmax,Kₓₗ₀_F,τ,ΔS,a_GPP,b_GPP,Nₛ_C,Kₓₗ₀_C = x    
    range = [(0.0001,0.1),
    (0.1,0.5),
    (0.001,1.0),
    (0.0004,0.1),
    (1.0,15.0),
    (10.0,25.0),
    (0.0001,5.0),
    (0.0001,3.0),
    (0.0001,0.1),
    (0.0004,0.1)]      

    res = BlackBoxOptim.bboptimize(x->opt_par_obj(x,
    raw_input_F,
    GPP_data_F,
    Ec_data_F,
    train_set,    
    raw_input_C,
    GPP_data_C,
    Ec_data_C,
    train_set;
    weight_GPP=weight_GPP); SearchRange = range,NThreads=Threads.nthreads()-1,MaxFuncEvals=MaxFuncEvals)
    x_opt = BlackBoxOptim.best_candidate(res) 
    
    x_opt_F = x_opt[1:8]
    x_opt_C = x_opt[[9,2,3,10,5,6,7,8]]
    
    par_opt_F = ModelPar(x_opt_F;stand_type=stand_type_F)
    par_opt_C = ModelPar(x_opt_C;stand_type=stand_type_C)
    
    JLD.save("output/"*folder_name*"/"*filename*"_F.jld",
    "x_opt",x_opt_F,
    "par_opt",par_opt_F,
    "stand_type",stand_type_F,
    "weight_GPP",weight_GPP,
    "train_set",train_set,
    "val_set",val_set,
    "train_set_weekly",train_set_weekly,
    "val_set_weekly",val_set_weekly)

    JLD.save("output/"*folder_name*"/"*filename*"_C.jld",
    "x_opt",x_opt_C,
    "par_opt",par_opt_C,
    "stand_type",stand_type_C,
    "weight_GPP",weight_GPP,
    "train_set",train_set,
    "val_set",val_set,
    "train_set_weekly",train_set_weekly,
    "val_set_weekly",val_set_weekly)

    return nothing
end  

function run_crossval_work_list()
    
    #run_crossval_opt("crossval_20240306_F_W_1_5_run_2",:Fertilized,1.5)
    #run_crossval_opt("crossval_20240306_F_W_1_5_run_3",:Fertilized,1.5)    
    #run_crossval_opt("crossval_20240306_F_W_1_5_run_4",:Fertilized,1.5)
    #run_crossval_opt("crossval_20240306_F_W_1_5_run_5",:Fertilized,1.5)
    #run_crossval_opt("crossval_20240306_F_W_1_5_run_6",:Fertilized,1.5)
    #run_crossval_opt("crossval_20240306_F_W_1_5_run_7",:Fertilized,1.5)
    #run_crossval_opt("crossval_20240306_F_W_1_5_run_8",:Fertilized,1.5)
    #run_crossval_opt("crossval_20240306_F_W_1_5_run_9",:Fertilized,1.5)
    run_crossval_opt("crossval_20240306_F_W_1_5_run_10",:Fertilized,1.5)

    run_crossval_opt("crossval_20240306_C_W_1_5_run_1",:Control,1.5)
    run_crossval_opt("crossval_20240306_C_W_1_5_run_2",:Control,1.5)
    run_crossval_opt("crossval_20240306_C_W_1_5_run_3",:Control,1.5)
    run_crossval_opt("crossval_20240306_C_W_1_5_run_4",:Control,1.5)
    run_crossval_opt("crossval_20240306_C_W_1_5_run_5",:Control,1.5)
    run_crossval_opt("crossval_20240306_C_W_1_5_run_6",:Control,1.5)
    run_crossval_opt("crossval_20240306_C_W_1_5_run_7",:Control,1.5)
    run_crossval_opt("crossval_20240306_C_W_1_5_run_8",:Control,1.5)
    run_crossval_opt("crossval_20240306_C_W_1_5_run_9",:Control,1.5)
    run_crossval_opt("crossval_20240306_C_W_1_5_run_10",:Control,1.5)
    
    run_crossval_opt("crossval_20240306_shared_W_1_5_run_1",1.5)
    run_crossval_opt("crossval_20240306_shared_W_1_5_run_2",1.5)
    run_crossval_opt("crossval_20240306_shared_W_1_5_run_3",1.5)
    run_crossval_opt("crossval_20240306_shared_W_1_5_run_4",1.5)
    run_crossval_opt("crossval_20240306_shared_W_1_5_run_5",1.5)
    run_crossval_opt("crossval_20240306_shared_W_1_5_run_6",1.5)
    run_crossval_opt("crossval_20240306_shared_W_1_5_run_7",1.5)
    run_crossval_opt("crossval_20240306_shared_W_1_5_run_8",1.5)
    run_crossval_opt("crossval_20240306_shared_W_1_5_run_9",1.5)
    run_crossval_opt("crossval_20240306_shared_W_1_5_run_10",1.5)
end