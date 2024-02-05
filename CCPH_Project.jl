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
    GPP_data_F = get_GPP_data.(raw_input_F;stand_type=:Fertilized)    
    par = ModelPar(0.013,0.14,0.011,0.00054,14.8,17.3,0.57,0.035)
    Xₜ_F = Xₜ_fun.(raw_input_F,Ref(par))

    modeloutput = run_week.(raw_input_F,Xₜ_F,Ref(par))
    
    GPP_model = get_GPP_model.(modeloutput)
    Ec_model = get_Ec_model.(modeloutput)   
    
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

run_test()