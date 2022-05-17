using Pkg
Pkg.activate(".")
using Revise, CCPH, Statistics, Distributions, Plots, JLD, Random, StatsPlots
import Dates, CSV, Interpolations, AdaptiveMCMC, BlackBoxOptim, Optim, DiffEqSensitivity

include("./Weather_Rosinedal_Struct/create_weather_struct_RO.jl")
include("./Script/Fit_Model_2_RO_Data/model_weather_ts_tuning_C_F_Auxiliary.jl")
include("Script/Fit_Model_2_RO_Data/model_tuning_C_F_Opt.jl")
include("Script/Fit_Model_2_RO_Data/model_tuning_C_F_Separate_Opt.jl")
include("Script/Fit_Model_2_RO_Data/model_tuning_C_F_Separate_MCMC.jl")
include("Script/Fit_Model_2_RO_Data/run_validation.jl")
include("CrossValidation/RO_CrossValidation.jl")
#include("Script/error_plot.jl")
include("Script/Fit_Model_2_RO_Data/work_list.jl")
#include("Script/Fit_Model_2_RO_Data/model_C_F_Sensitivity_Analysis.jl")

#=
parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
        ranges = [(0.0001,0.1),(10.0,40.0),(0.01,1.0),(0.0005,0.1),
        (0.1,10.0),(0.1,0.5),(1.0,15.0),(10.0,25.0),
        (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
        file_name = "RO_Opt_Separate_GPP_EC_Nm_f_20220428"
run_opt_par(file_name,ranges,parasym;ParaDictInit_C=ParaDictInit_C) 
run_validation_RO_2019(file_name,ranges,parasym;ParaDictInit_C=ParaDictInit_C)
=#