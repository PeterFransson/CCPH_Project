using Pkg
Pkg.activate(".")
using Revise, CCPH, Statistics, Distributions, Plots, JLD, Random, StatsPlots
import Dates, CSV, Interpolations, AdaptiveMCMC, BlackBoxOptim, Optim, DiffEqSensitivity

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

#include("Script/Test_gs_I/test_gs_I.jl")
#include("Script/Fit_Model_2_RO_Data/test_A_N.jl")