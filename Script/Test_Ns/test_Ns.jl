function test_Ns()
    fld = "crossval_20241105_shared_W_1_5_run_9"

    #First test
    println("First test")

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

    println("--Default C Nᵤ=$(par.Nₛ-0.0056)")
    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)")

    par.Nₛ = 0.0116#0.0056+0.012 #Nᵤ=0.017

    GPP_model,Ec_model,Nₘ_f_model = run_model(par,raw_input;stand_type=stand_type)

    GPP_R2,GPP_RMSE,GPP_MAPE,GPP_cor = get_sum_stat(GPP_data,[GPP_model[i]*raw_input[i].ζ for i in 1:4])
    Ec_R2,Ec_RMSE,Ec_MAPE,Ec_cor = get_sum_stat(Ec_data,Ec_model)

    println("--Set F Nᵤ=$(par.Nₛ-0.0056)")
    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)") 
    
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

    println("--Default C Nᵤ=$(par.Nₛ-0.0056)")
    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)")  

    par.Nₛ = 0.0116#0.0056 #Nᵤ=0.0

    GPP_model,Ec_model,Nₘ_f_model = run_model(par,raw_input;stand_type=stand_type)

    GPP_R2,GPP_RMSE,GPP_MAPE,GPP_cor = get_sum_stat(GPP_data,[GPP_model[i]*raw_input[i].ζ for i in 1:4])
    Ec_R2,Ec_RMSE,Ec_MAPE,Ec_cor = get_sum_stat(Ec_data,Ec_model)

    println("--Set C Nᵤ to F Nᵤ=$(par.Nₛ-0.0056)")
    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)")  

    println("Second test")

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

    println("--Default C Nᵤ=$(par.Nₛ-0.0056)")
    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)")

    par.Nₛ = 0.0056+0.012 #Nᵤ=0.012

    GPP_model,Ec_model,Nₘ_f_model = run_model(par,raw_input;stand_type=stand_type)

    GPP_R2,GPP_RMSE,GPP_MAPE,GPP_cor = get_sum_stat(GPP_data,[GPP_model[i]*raw_input[i].ζ for i in 1:4])
    Ec_R2,Ec_RMSE,Ec_MAPE,Ec_cor = get_sum_stat(Ec_data,Ec_model)

    println("--Set F Nᵤ=$(par.Nₛ-0.0056)")
    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)") 
    
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

    println("--Default C Nᵤ=$(par.Nₛ-0.0056)")
    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)") 

    par.Nₛ = 0.0056+0.0 #Nᵤ=0.0

    GPP_model,Ec_model,Nₘ_f_model = run_model(par,raw_input;stand_type=stand_type)

    GPP_R2,GPP_RMSE,GPP_MAPE,GPP_cor = get_sum_stat(GPP_data,[GPP_model[i]*raw_input[i].ζ for i in 1:4])
    Ec_R2,Ec_RMSE,Ec_MAPE,Ec_cor = get_sum_stat(Ec_data,Ec_model)

    println("--Set C Nᵤ=$(par.Nₛ-0.0056)")
    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)")     
end