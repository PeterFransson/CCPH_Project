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

function run_test()
    raw_input_F = RawInputData(;stand_type=:Fertilized)
    Ec_data_F = calc_Ec_data.(raw_input_F)
    GPP_data_F = get_GPP_data(stand_type=:Fertilized)
    #plot([first(data).date for data in raw_input_F[1].weather],GPP_data_F[1])
    par = ModelPar(0.013,0.14,0.011,0.00054,14.8,17.3,0.57,0.035)
    Xₜ_F = Xₜ_fun.(raw_input_F,Ref(par))

    modeloutput = run_week.(raw_input_F,Xₜ_F,Ref(par))

    GPP_model = get_GPP_model(modeloutput[1])
   
    #plot([first(data).date for data in raw_input_F[1].weather],Xₜ_F[1])
end

run_test()