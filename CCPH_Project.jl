using Pkg
Pkg.activate(".")
using Revise, CCPH, Statistics, Distributions, Plots, JLD
import Dates, CSV, Interpolations, AdaptiveMCMC, BlackBoxOptim, Optim

include("./mcmc/metropolis_log_P.jl")
include("./mcmc/analysis.jl")
include("./Weather_Rosinedal_Struct/create_weather_struct_RO.jl")
include("./Script/Fit_Model_2_RO_Data/model_weather_ts_tuning_C_F_Auxiliary.jl")
include("./Transpiration/transpiration.jl")
include("Script/Fit_Model_2_RO_Data/model_tuning_C_F_Separate_Opt.jl")