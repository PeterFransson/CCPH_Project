function plot_work_list()
    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.0018,:a_Ec=>0.0433,:b_Ec=>0.179)   
    ParaDictInit_F = Dict(:μ_Nₘ_f=>0.018,:b_Nₘ_f=>0.0018,:a_Ec=>0.0688,:b_Ec=>0.146)
    file_name = "RO_Opt_GPP_Ec_Nm_f_20220619_2"

    calibparavec = CalibParaVec(
    (:Nₛ,0.0001,0.1,true),
    (:a_Jmax,0.01,1.0),
    (:Kₓₗ₀,0.0005,0.1,true),
    (:α_max,0.1,0.5),
    (:τ,1.0,15.0),
    (:Smax,10.0,25.0),                        
    (:a_GPP,0.0001,5.0),
    (:b_GPP,0.0001,3.0)) 

    ranges,para2ind = CreateOptVar(calibparavec)     

    run_create_result_plots(file_name,
    para2ind,
    10;
    ParaDictInit_F=ParaDictInit_F,
    ParaDictInit_C=ParaDictInit_C)
end

plot_work_list()