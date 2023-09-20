function plot_work_list()
    #=
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
    =#

    
    #---Run non-sharing parameter case, Start---
    parasym = [:Nₛ,
    :a_Jmax,
    :Kₓₗ₀,    
    :α_max,
    :τ,
    :Smax,
    :a_GPP,
    :b_GPP]
    ranges = [(0.0001,0.1),
    (0.01,1.0),
    (0.0005,0.1),    
    (0.1,0.5),
    (1.0,15.0),
    (10.0,25.0),
    (0.0001,5.0),
    (0.0001,3.0)]

    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.0018,:a_Ec=>0.0433,:b_Ec=>0.179)   
    ParaDictInit_F = Dict(:μ_Nₘ_f=>0.018,:b_Nₘ_f=>0.0018,:a_Ec=>0.0688,:b_Ec=>0.146)
    file_name = "RO_Opt_Separate_GPP_Ec_Nm_f_eco_scaling_New_Wf_20230829"

    run_create_result_plots(file_name,
    ranges,
    parasym,
    10;
    ParaDictInit_F=ParaDictInit_F,
    ParaDictInit_C=ParaDictInit_C)
    #---Run non-sharing parameter case, Done---
    

    #---Run sharing parameter case, Start---
    calibparavec = CalibParaVec(
        (:Nₛ,0.0001,0.1,true),(:a_Jmax,0.01,1.0),
        (:Kₓₗ₀,0.0005,0.1,true),
        (:α_max,0.1,0.5),
        (:τ,1.0,15.0),(:Smax,10.0,25.0),                        
        (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0)) 

    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.0018,:a_Ec=>0.0433,:b_Ec=>0.179)   
    ParaDictInit_F = Dict(:μ_Nₘ_f=>0.018,:b_Nₘ_f=>0.0018,:a_Ec=>0.0688,:b_Ec=>0.146) 
    ranges,para2ind = CreateOptVar(calibparavec) 
    file_name = "RO_Opt_GPP_Ec_Nm_f_eco_scaling_New_Wf_20230829"

    run_create_result_plots(file_name,
    ranges,
    para2ind,
    10;
    ParaDictInit_F=ParaDictInit_F,
    ParaDictInit_C=ParaDictInit_C)
    #---Run sharing parameter case, Done---
end

plot_work_list()