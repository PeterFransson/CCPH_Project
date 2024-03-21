using Pkg
Pkg.activate(".")
using Revise, CCPH, Statistics, Distributions, Plots, JLD, Random, StatsPlots
import Dates, CSV, DataFrames, BlackBoxOptim

#=
include("./Weather_Rosinedal_Struct/create_weather_struct_RO.jl")
include("./Script/Fit_Model_2_RO_Data/model_weather_ts_tuning_C_F_Auxiliary.jl")
include("Script/Fit_Model_2_RO_Data/result_plots_script.jl")
include("Script/Fit_Model_2_RO_Data/result_output_script.jl")
include("Script/Fit_Model_2_RO_Data/model_tuning_C_F_Opt.jl")
include("Script/Fit_Model_2_RO_Data/model_tuning_C_F_Separate_Opt.jl")
include("Script/Fit_Model_2_RO_Data/run_validation.jl")
include("CrossValidation/RO_CrossValidation.jl")
include("Script/Fit_Model_2_RO_Data/work_list.jl")
include("Script/Fit_Model_2_RO_Data/create_result_plots.jl")
include("Script/Fit_Model_2_RO_Data/plot_work_list.jl")
=#

#include("Script/Test_gs_I/test_gs_I.jl")
#include("Script/Fit_Model_2_RO_Data/test_A_N.jl")

include("./Weather_Rosinedal_Struct/create_weather_struct_RO.jl")
include("./Script/Fit_Model_2_RO_Data/run_par_fit.jl")
include("./Script/Fit_Model_2_RO_Data/run_crossval.jl")
include("./Script/paper_results/weather_plot.jl")
include("./Script/paper_results/model_data.jl")
include("./Script/paper_results/plant_var_vs_weather.jl")
#include("Script/Test_static_N/test_static_N.jl")

function run_opt_test(folder_name::String,stand_type::Symbol,weight_GPP::Real)

    #folder_name = "test_new_opt_20240207_C_W_1_5"
    filename = "result"

    isdir("./output/"*folder_name)|| mkdir("./output/"*folder_name)

    #stand_type=:Control
    raw_input = RawInputData(;stand_type=stand_type)
    Ec_data = calc_Ec_data.(raw_input)
    GPP_data = get_GPP_data.(raw_input;stand_type=stand_type)    
    
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
    Ec_data;
    stand_type=stand_type,
    weight_GPP=weight_GPP); SearchRange = range,NThreads=Threads.nthreads()-1)
    x_opt = BlackBoxOptim.best_candidate(res) 
    par_opt = ModelPar(x_opt;stand_type=stand_type)
    
    JLD.save("output/"*folder_name*"/"*filename*".jld","x_opt",x_opt,"par_opt",par_opt,"stand_type",stand_type,"weight_GPP",weight_GPP)
    return nothing
end  

function draw_opt_test(folder_name::String)

    #folder_name = "test_new_opt_20240207_C_W_1_5"
    filename = "result"
    
    stand_type = JLD.load("output/"*folder_name*"/"*filename*".jld","stand_type")

    raw_input = RawInputData(;stand_type=stand_type)
    Ec_data = calc_Ec_data.(raw_input)
    GPP_data = get_GPP_data.(raw_input;stand_type=stand_type)  

    #Load opt val 
    @show x_opt = JLD.load("output/"*folder_name*"/"*filename*".jld","x_opt")   
        
    par = ModelPar(x_opt;stand_type=stand_type)

    GPP_model,Ec_model,Nₘ_f_model = run_model(par,raw_input;stand_type=stand_type)

    GPP_R2,GPP_RMSE,GPP_MAPE,GPP_cor = get_sum_stat(GPP_data,GPP_model)
    Ec_R2,Ec_RMSE,Ec_MAPE,Ec_cor = get_sum_stat(Ec_data,Ec_model)

    println("GPP: R²:$(GPP_R2), RMSE:$(GPP_RMSE), MAPE:$(GPP_MAPE), corr:$(GPP_cor)")
    println("Ec: R²:$(Ec_R2), RMSE:$(Ec_RMSE), MAPE:$(Ec_MAPE), corr:$(Ec_cor)")
    
    plot([data.date for data in raw_input[1].weather_growth],GPP_data[1],linecolor=:blue,xlabel="Date",ylabel="GPP")
    plot!([data.date for data in raw_input[1].weather_growth],GPP_model[1]*raw_input[1].ζ,linecolor=:red) 
       
    for i in 2:4 
        plot!([data.date for data in raw_input[i].weather_growth],GPP_data[i],linecolor=:blue)
        plot!([data.date for data in raw_input[i].weather_growth],GPP_model[i]*raw_input[i].ζ,linecolor=:red)
    end
    
    pl1 = plot!(legends=false)    
    
    plot([data.date for data in raw_input[1].weather_growth],Ec_data[1],linecolor=:blue,xlabel="Date",ylabel="Ec",label="Data")
    plot!([data.date for data in raw_input[1].weather_growth],Ec_model[1],linecolor=:red,label="Model") 
      
    for i in 2:4 
        plot!([data.date for data in raw_input[i].weather_growth],Ec_data[i],linecolor=:blue,label="")
        plot!([data.date for data in raw_input[i].weather_growth],Ec_model[i],linecolor=:red,label="")
    end
    
    pl2 = plot!() 
    

    plot([data.date for data in raw_input[1].weather_growth],Nₘ_f_model[1]*100,linecolor=:red,xlabel="Date",ylabel="Nₘ (%)")
         
    for i in 2:4 
        plot!([data.date for data in raw_input[i].weather_growth],Nₘ_f_model[i]*100,linecolor=:red)        
    end

    pl3 = plot!(legends=false) 
    
    plot(pl1,pl2,pl3,layout=(3,1))    
end

function run_work_list()
    #=
    run_opt_test("test_new_opt_20240207_F_W_1_0",:Fertilized,1.0)
    run_opt_test("test_new_opt_20240207_C_W_1_0",:Control,1.0)

    run_opt_test("test_new_opt_20240207_F_W_1_3",:Fertilized,1.3)
    run_opt_test("test_new_opt_20240207_C_W_1_3",:Control,1.3)

    run_opt_test("test_new_opt_20240207_F_W_1_7",:Fertilized,1.7)
    run_opt_test("test_new_opt_20240207_C_W_1_7",:Control,1.7)
    =#

    run_opt_test("test_new_opt_20240207_F_W_1_5",:Fertilized,1.5)
    run_opt_test("test_new_opt_20240207_C_W_1_5",:Control,1.5)
end

function run_draw_work_list()

    println("--Weight: 1.0--")
    println("--Fetrilized")
    draw_opt_test("test_new_opt_20240207_F_W_1_0")
    println("--Control")
    draw_opt_test("test_new_opt_20240207_C_W_1_0")

    println("--Weight: 1.3--")
    println("--Fetrilized")
    draw_opt_test("test_new_opt_20240207_F_W_1_3")
    println("--Control")
    draw_opt_test("test_new_opt_20240207_C_W_1_3")

    println("--Weight: 1.5--")
    println("--Fetrilized")
    draw_opt_test("test_new_opt_20240207_F_W_1_5")
    println("--Control")
    draw_opt_test("test_new_opt_20240207_C_W_1_5")
   
    println("--Weight: 1.7--")
    println("--Fetrilized")
    draw_opt_test("test_new_opt_20240207_F_W_1_7")
    println("--Control")
    draw_opt_test("test_new_opt_20240207_C_W_1_7")
end

function test_train_val()
    stand_type=:Control
    raw_input = RawInputData(;stand_type=stand_type)
    Ec_data = calc_Ec_data.(raw_input)
    GPP_data = get_GPP_data.(raw_input;stand_type=stand_type)  

    @show raw_input[1].growth_indices_weekly
    @show raw_input[2].growth_indices_weekly
    @show raw_input[3].growth_indices_weekly
    @show raw_input[4].growth_indices_weekly
     
    (train_set,val_set,train_set_weekly,val_set_weekly) = CreateTrainValSet(raw_input) 
        
    GPP_data_train = [GPP_data[i][train_set[i]] for i in 1:4]   
    GPP_data_val = [GPP_data[i][val_set[i]] for i in 1:4]     
end

function run_get_data_work_list()  
    fld_F = "crossval_20240306_F_W_1_5_run"
    fld_C = "crossval_20240306_C_W_1_5_run"

    for i = 1:10
        println("----Run_$(i)---")
        println("--Fertilized")
        get_sum_stat(fld_F*"_$(i)","result")  
        println("--Control")     
        get_sum_stat(fld_C*"_$(i)","result")        
    end
end   

#test_train_val()

#run_work_list()
#run_draw_work_list()
#draw_opt_test("test_new_opt_20240207_C_W_1_7")
#run_opt_test()
#draw_opt_test()
#test_static_N()
#run_crossval_work_list()
#run_get_data_work_list() 
#create_weather_plot()
#draw_shared_model()
#draw_nonshared_model()
plant_var_vs_weather()