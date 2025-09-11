function test_static_swc()
    fld = "crossval_20241105_shared_W_1_5_run_9"
    
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

    println("--Dynamic θₛ")
    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)")

    θₛ_growth = []
    θₛ_growth_2015 = [weather.θₛ for weather in raw_input[1].weather_growth]
    θₛ_growth_2016 = [weather.θₛ for weather in raw_input[2].weather_growth]
    θₛ_growth_2017 = [weather.θₛ for weather in raw_input[3].weather_growth]
    θₛ_growth_2018 = [weather.θₛ for weather in raw_input[4].weather_growth]
    append!(θₛ_growth,θₛ_growth_2015,θₛ_growth_2016,θₛ_growth_2017,θₛ_growth_2018)
    θₛ_mean = mean(θₛ_growth)
    println("Mean θₛ=$(θₛ_mean)")
    
    raw_input_mean = deepcopy(raw_input)
    
    #Change SWC (θₛ)
    for i in 1:4
        for j in eachindex(raw_input_mean[i].weather_raw)
            raw_input_mean[i].weather_raw[j].θₛ = θₛ_mean
        end
        for j in eachindex(raw_input_mean[i].weather_growth)
            raw_input_mean[i].weather_raw[j].θₛ = θₛ_mean
        end
    end
    
    GPP_model,Ec_model,Nₘ_f_model = run_model(par,raw_input_mean;stand_type=stand_type)

    GPP_R2,GPP_RMSE,GPP_MAPE,GPP_cor = get_sum_stat(GPP_data,[GPP_model[i]*raw_input_mean[i].ζ for i in 1:4])
    Ec_R2,Ec_RMSE,Ec_MAPE,Ec_cor = get_sum_stat(Ec_data,Ec_model)

    println("--Static θₛ")
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

    println("--Dynamic θₛ")
    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)")  
    
    θₛ_growth = []
    θₛ_growth_2015 = [weather.θₛ for weather in raw_input[1].weather_growth]
    θₛ_growth_2016 = [weather.θₛ for weather in raw_input[2].weather_growth]
    θₛ_growth_2017 = [weather.θₛ for weather in raw_input[3].weather_growth]
    θₛ_growth_2018 = [weather.θₛ for weather in raw_input[4].weather_growth]
    append!(θₛ_growth,θₛ_growth_2015,θₛ_growth_2016,θₛ_growth_2017,θₛ_growth_2018)
    θₛ_mean = mean(θₛ_growth)
    println("Mean θₛ=$(θₛ_mean)")
    
    raw_input_mean = deepcopy(raw_input)

    #Change SWC (θₛ)
    for i in 1:4
        for j in eachindex(raw_input_mean[i].weather_raw)
            raw_input_mean[i].weather_raw[j].θₛ = θₛ_mean
        end
        for j in eachindex(raw_input_mean[i].weather_growth)
            raw_input_mean[i].weather_raw[j].θₛ = θₛ_mean
        end
    end

    GPP_model,Ec_model,Nₘ_f_model = run_model(par,raw_input_mean;stand_type=stand_type)

    GPP_R2,GPP_RMSE,GPP_MAPE,GPP_cor = get_sum_stat(GPP_data,[GPP_model[i]*raw_input_mean[i].ζ for i in 1:4])
    Ec_R2,Ec_RMSE,Ec_MAPE,Ec_cor = get_sum_stat(Ec_data,Ec_model)

    println("--Static θₛ")
    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)")    
end