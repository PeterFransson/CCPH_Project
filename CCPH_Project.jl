using Pkg
Pkg.activate(".")
using Revise, CCPH, Statistics, Distributions, Plots, JLD
import Dates, CSV, Interpolations

include("./mcmc/metropolis.jl")
include("./mcmc/analysis.jl")
include("./Weather_Rosinedal_Struct/create_weather_struct_RO.jl")
include("./Script/model_weather_ts_tuning_C_F_Auxiliary.jl")
include("Script/model_weather_ts_tuning_C_F_GPP_only_Few_para_NoRoot.jl")
#include("./Transpiration/transpiration.jl")
#include("./Transpiration/test_transpiration.jl")
#include("Script/model_weather_ts_tuning_C_F_GPP_only_w_ET.jl")
#include("./Script/model_weather_ts_tuning_C_F.jl")
#include("Script/model_weather_ts_tuning_C_F_min_par.jl")
#include("Script/model_weather_ts_tuning_C_F_GPP_only.jl")
#include("./Script/model_weather_ts_tuning_C_F_with_Wf.jl")