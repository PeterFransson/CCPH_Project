function run_week(data_week::Vector{CCPH.WeatherDataStruct},
    Xₜ_week::Vector{Real},
    treesize::CCPH.TreeSize,
    par::ModelPar,
    mean_Nₘ_f::Real)
    
    model,envfun,daylength,kinetic,cons,treepar,treesize,hydPar = intitiate_model(data_week,Xₜ_week,treesize,par)
    
    optval_weekly,output_weekly = get_model_output(daylength,model,kinetic,envfun)
    
    output_day_vec = Vector{CCPH.CCPHOutput}(undef,7)
    optval_day_vec =Vector{OptVal}(undef,7)

    for i = 1:7

        model_day,envfun_day,daylength_day = intitiate_model(data_week[i],Xₜ_week[i],kinetic,cons,treepar,treesize,hydPar)
                
        optval_day,output_day = get_model_output(mean_Nₘ_f,
        optval_weekly.gₛ₁,
        optval_weekly.gₛ₂,
        daylength_day,
        model_day,
        kinetic,
        envfun_day)  
        
        output_day_vec[i] = output_day
        optval_day_vec[i] = optval_day
    end
    
    return ModelOutput(output_weekly,optval_weekly,output_day_vec,optval_day_vec)
end

function run_week(input_data::RawInputData,
    Xₜ::Vector{Real},
    par::ModelPar,
    mean_Nₘ_f::Real)

    treesize = input_data.treesize    

    modeloutput = [run_week(input_data.weather_growth[idx[1]:idx[2]],Xₜ[idx[1]:idx[2]],treesize,par,mean_Nₘ_f) for idx in input_data.growth_indices_weekly]

    return modeloutput
end

function run_model(par::ModelPar,raw_input::Vector{RawInputData},mean_Nₘ_f::Real;stand_type::Symbol=:Fertilized)
    Xₜ = Xₜ_fun.(raw_input,Ref(par))
    modeloutput = run_week.(raw_input,Xₜ,Ref(par),mean_Nₘ_f)    
    GPP_model = get_GPP_model.(modeloutput)
    Ec_model = get_Ec_model.(modeloutput)
    Nₘ_f_model = get_Nₘ_f_model.(modeloutput)
    return (GPP_model,Ec_model,Nₘ_f_model)
end

function test_static_N()
    fld = "crossval_20240306_shared_W_1_5_run_4"
    
    #--Fertilized--
    println("--Fertilized--")
    stand_type = JLD.load("output/"*fld*"/result_F.jld","stand_type")
    x_opt = JLD.load("output/"*fld*"/result_F.jld","x_opt")  

    raw_input = RawInputData(;stand_type=stand_type)
    Ec_data = calc_Ec_data.(raw_input)
    GPP_data = get_GPP_data.(raw_input;stand_type=stand_type) 
   
    par = ModelPar(x_opt;stand_type=stand_type)

    GPP_model,Ec_model,Nₘ_f_model = run_model(par,raw_input;stand_type=stand_type)

    GPP_R2,GPP_RMSE,GPP_MAPE,GPP_cor = get_sum_stat(GPP_data,[GPP_model[i]*raw_input[i].ζ for i in 1:4])
    Ec_R2,Ec_RMSE,Ec_MAPE,Ec_cor = get_sum_stat(Ec_data,Ec_model)

    println("--Dynamic Nₘ_f")
    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)")

    @show mean_Nₘ_f = sum(sum.(Nₘ_f_model))/sum(length.(Nₘ_f_model))

    GPP_model,Ec_model,Nₘ_f_model = run_model(par,raw_input,mean_Nₘ_f;stand_type=stand_type)

    GPP_R2,GPP_RMSE,GPP_MAPE,GPP_cor = get_sum_stat(GPP_data,[GPP_model[i]*raw_input[i].ζ for i in 1:4])
    Ec_R2,Ec_RMSE,Ec_MAPE,Ec_cor = get_sum_stat(Ec_data,Ec_model)

    println("--Static Nₘ_f")
    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)")  
    
    #--Control--
    println("--Control--")
    stand_type = JLD.load("output/"*fld*"/result_C.jld","stand_type")
    x_opt = JLD.load("output/"*fld*"/result_C.jld","x_opt")  

    raw_input = RawInputData(;stand_type=stand_type)
    Ec_data = calc_Ec_data.(raw_input)
    GPP_data = get_GPP_data.(raw_input;stand_type=stand_type) 
    
    par = ModelPar(x_opt;stand_type=stand_type)

    GPP_model,Ec_model,Nₘ_f_model = run_model(par,raw_input;stand_type=stand_type)

    GPP_R2,GPP_RMSE,GPP_MAPE,GPP_cor = get_sum_stat(GPP_data,[GPP_model[i]*raw_input[i].ζ for i in 1:4])
    Ec_R2,Ec_RMSE,Ec_MAPE,Ec_cor = get_sum_stat(Ec_data,Ec_model)

    println("--Dynamic Nₘ_f")
    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)")

    @show mean_Nₘ_f = sum(sum.(Nₘ_f_model))/sum(length.(Nₘ_f_model))

    GPP_model,Ec_model,Nₘ_f_model = run_model(par,raw_input,mean_Nₘ_f;stand_type=stand_type)

    GPP_R2,GPP_RMSE,GPP_MAPE,GPP_cor = get_sum_stat(GPP_data,[GPP_model[i]*raw_input[i].ζ for i in 1:4])
    Ec_R2,Ec_RMSE,Ec_MAPE,Ec_cor = get_sum_stat(Ec_data,Ec_model)

    println("--Static Nₘ_f")
    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)")    
end