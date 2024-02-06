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

function run_opt_test()

    folder_name = "test_new_opt_20240206"
    filename = "result"

    isdir("./output/"*folder_name)|| mkdir("./output/"*folder_name)

    stand_type=:Fertilized
    raw_input_F = RawInputData(;stand_type=stand_type)
    Ec_data_F = calc_Ec_data.(raw_input_F)
    GPP_data_F = get_GPP_data.(raw_input_F;stand_type=stand_type)    
    
    #Nₛ,α_max,a_Jmax,Kₓₗ₀,τ,ΔS,a_GPP,b_GPP = x

    
    range = [(0.0001,0.1),
    (0.1,0.5),
    (0.01,1.0),
    (0.0004,0.1),
    (1.0,15.0),
    (10.0,25.0),
    (0.0001,5.0),
    (0.0001,3.0)]

    x0 = [0.015,
    0.13,
    0.010,
    0.00059,
    14.6,
    17.82,
    0.53,
    0.042]
    

    #=
    range = [(0.005,0.1),
    (0.1,0.3),
    (0.008,0.2),
    (0.0004,0.001),
    (7.0,15.0),
    (10.0,25.0),
    (0.0001,5.0),
    (0.0001,3.0)]
    =#

    res = BlackBoxOptim.bboptimize(x->opt_par_obj(x,
    raw_input_F,
    GPP_data_F,
    Ec_data_F;
    stand_type=stand_type); SearchRange = range,NThreads=Threads.nthreads()-1)
    x_opt = BlackBoxOptim.best_candidate(res) 
    par_opt = ModelPar(x_opt;stand_type=stand_type)
    
    JLD.save("output/"*folder_name*"/"*filename*".jld","x_opt",x_opt,"par_opt",par_opt,"stand_type",stand_type)
end  

function draw_opt_test()

    folder_name = "test_new_opt_20240206"
    filename = "result"
    
    stand_type = JLD.load("output/"*folder_name*"/"*filename*".jld","stand_type")

    raw_input_F = RawInputData(;stand_type=stand_type)
    Ec_data_F = calc_Ec_data.(raw_input_F)
    GPP_data_F = get_GPP_data.(raw_input_F;stand_type=stand_type)  

    #Load opt val 
    x_opt = JLD.load("output/"*folder_name*"/"*filename*".jld","x_opt")

    par = ModelPar(x_opt;stand_type=stand_type)

    GPP_model,Ec_model,Nₘ_f_model = run_model(par,raw_input_F;stand_type=stand_type)
    
    plot([data.date for data in raw_input_F[1].weather_growth],GPP_data_F[1],linecolor=:blue)
    plot!([data.date for data in raw_input_F[1].weather_growth],GPP_model[1]*raw_input_F[1].ζ,linecolor=:red)    
    for i in 2:4 
        plot!([data.date for data in raw_input_F[i].weather_growth],GPP_data_F[i],linecolor=:blue)
        plot!([data.date for data in raw_input_F[i].weather_growth],GPP_model[i]*raw_input_F[i].ζ,linecolor=:red)
    end
    pl1 = plot!()    
    
    plot([data.date for data in raw_input_F[1].weather_growth],Ec_data_F[1],linecolor=:blue)
    plot!([data.date for data in raw_input_F[1].weather_growth],Ec_model[1],linecolor=:red)    
    for i in 2:4 
        plot!([data.date for data in raw_input_F[i].weather_growth],Ec_data_F[i],linecolor=:blue)
        plot!([data.date for data in raw_input_F[i].weather_growth],Ec_model[i],linecolor=:red)
    end
    pl2 = plot!() 
    
    plot(pl1,pl2,layout=(2,1),legends=false)
end
run_opt_test()
draw_opt_test()