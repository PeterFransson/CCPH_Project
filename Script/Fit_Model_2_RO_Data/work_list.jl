function work_list()
    #=
    Threads.@threads for ix = 1:3
    if ix == 1    
    parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:X0,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
    ranges = [(0.001,0.1),(10.0,40.0),(0.01,1.0),(0.001,0.1),
        (0.1,6.0),(0.1,0.5),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),
        (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
    x0_F = [0.0010018365156108655, 13.055483612590724, 0.3002418081993492,
        0.07147806011034544, 2.3366352738114964, 0.48788007386624604, -1.1099087213234984,
        14.921141416504296, 23.16827082107623, 0.6507597997899615, 0.0451402105520481, 0.12222273101753175, 0.0687942434461461]
    x0_C = [0.035852741991519714, 11.62745193905442, 0.9627222534724391,
        0.0012322987818022921, 4.147640506056484, 0.1681112431600152, -3.775281121258116,
        3.9015734019298396, 22.54869224534012, 0.19973955332162818, 0.3334662663349886, 0.19660296597291332, 0.07248679732356901] 
    file_name = "RO_MCMC_Separate_GPP_EC_20220326"
    run_sampler(file_name,ranges,parasym,x0_F,x0_C;n_samples=10000,n_burn_in=3000)
    plot_samples(file_name,parasym)
    end

    if ix == 2
    parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:X0,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
    ranges = [(0.001,0.1),(10.0,40.0),(0.01,1.0),(0.001,0.1),
    (0.1,6.0),(0.1,0.5),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),
    (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
    file_name = "RO_Opt_Separate_GPP_EC_Nm_f_20220326" 
    run_calibration(file_name,ranges,parasym;
    PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C)
    end
    end 
    =#

    #=
    parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:X0,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
    ranges = [(0.0001,0.1),(10.0,40.0),(0.01,1.0),(0.0005,0.1),
    (0.1,10.0),(0.1,0.5),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),
    (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
    file_name = "RO_Opt_Separate_GPP_EC_Nm_f_20220421" 
    run_calibration(file_name,ranges,parasym;
    PopulationSize=100,MaxSteps=10000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C)
    =#  
      
    #=
    Threads.@threads for ix = 1:6
    if ix == 1
    calibparavec = CalibParaVec(
                (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                (:Kₓₗ₀,0.0005,0.1),(:i,0.1,10.0),(:α_max,0.1,0.5),
                (:τ,1.0,15.0),(:Smax,10.0,25.0),
                (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
            ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
            ranges,para2ind = CreateOptVar(calibparavec) 
            file_name = "RO_Opt_GPP_EC_Nm_f_202204028_1" 
            run_calibration(file_name,ranges,para2ind;
            PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C)

    end
    if ix == 2
        calibparavec = CalibParaVec(
                    (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                    (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5),
                    (:τ,1.0,15.0),(:Smax,10.0,25.0),
                    (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
                ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
                ranges,para2ind = CreateOptVar(calibparavec) 
                file_name = "RO_Opt_GPP_EC_Nm_f_202204028_2" 
                run_calibration(file_name,ranges,para2ind;
                PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C)
    
    end
    if ix == 3
        calibparavec = CalibParaVec(
                    (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                    (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5),
                    (:τ,1.0,15.0,true),(:Smax,10.0,25.0),
                    (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
                ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
                ranges,para2ind = CreateOptVar(calibparavec) 
                file_name = "RO_Opt_GPP_EC_Nm_f_202204028_3"
                run_calibration(file_name,ranges,para2ind;
                PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C)
    
    end
    if ix == 4
        calibparavec = CalibParaVec(
                    (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                    (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5),
                    (:τ,1.0,15.0,true),(:Smax,10.0,25.0,true),
                    (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
                ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
                ranges,para2ind = CreateOptVar(calibparavec) 
                file_name = "RO_Opt_GPP_EC_Nm_f_202204028_4"
                run_calibration(file_name,ranges,para2ind;
                PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C)
    
    end
    if ix == 5
        calibparavec = CalibParaVec(
        (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
        (:Kₓₗ₀,0.0005,0.1),(:i,0.1,10.0),(:α_max,0.1,0.5),
        (:τ,1.0,15.0),(:Smax,10.0,25.0,true),
        (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
        ranges,para2ind = CreateOptVar(calibparavec) 
        file_name = "RO_Opt_GPP_EC_Nm_f_202204028_5" 
        run_calibration(file_name,ranges,para2ind;
        PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C)
    end
    if ix == 6
        parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
        ranges = [(0.0001,0.1),(10.0,40.0),(0.01,1.0),(0.0005,0.1),
        (0.1,10.0),(0.1,0.5),(1.0,15.0),(10.0,25.0),
        (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
        file_name = "RO_Opt_Separate_GPP_EC_Nm_f_20220428" 
        run_calibration(file_name,ranges,parasym;
        PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C)
    end
    end  
    =# 
    #= 
    parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:τ,:Smax,:a_GPP,:b_GPP]
        ranges = [(0.0001,0.1),(10.0,40.0),(0.01,1.0),(0.0005,0.1),
        (0.1,10.0),(0.1,0.5),(1.0,15.0),(10.0,25.0),
        (0.0001,5.0),(0.0001,3.0)]
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
        file_name = "RO_Opt_Separate_GPP_Nm_f_20220429" 
        run_calibration(file_name,ranges,parasym;
        PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Nm_f,ParaDictInit_C=ParaDictInit_C)

        parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:τ,:Smax,:a_GPP,:b_GPP]
        ranges = [(0.0001,0.1),(10.0,40.0),(0.01,1.0),(0.0005,0.1),
        (0.1,10.0),(0.1,0.5),(1.0,15.0),(10.0,25.0),
        (0.0001,5.0),(0.0001,3.0)]
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
        file_name = "RO_Opt_Separate_GPP_20220429" 
        run_calibration(file_name,ranges,parasym;
        PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP,ParaDictInit_C=ParaDictInit_C)
        =#
        #=
        parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
        ranges = [(0.0001,0.1),(10.0,40.0),(0.01,1.0),(0.0005,0.1),
        (0.1,10.0),(0.1,0.5),(1.0,15.0),(10.0,25.0),
        (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
        file_name = "RO_Opt_Separate_GPP_EC_Nm_f_20220501" 
        run_calibration(file_name,ranges,parasym;
        PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=append!(collect(1:21),collect(45:84)))
        run_validation_RO_2019(file_name,ranges,parasym;ParaDictInit_C=ParaDictInit_C,ind=collect(22:44))
        =#
    #=
    Threads.@threads for ix = 1:3
    if ix == 1
    calibparavec = CalibParaVec(
                (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                (:Kₓₗ₀,0.0005,0.1),(:i,0.1,10.0),(:α_max,0.1,0.5),
                (:τ,1.0,15.0),(:Smax,10.0,25.0),
                (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
            ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
            ranges,para2ind = CreateOptVar(calibparavec) 
            file_name = "RO_Opt_GPP_EC_Nm_f_20220502_1"
            run_calibration(file_name,ranges,para2ind;
            PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=append!(collect(1:21),collect(45:84)))
            run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(22:44))
    end
    if ix == 2
        calibparavec = CalibParaVec(
                    (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                    (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5),
                    (:τ,1.0,15.0,true),(:Smax,10.0,25.0),
                    (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
                ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
                ranges,para2ind = CreateOptVar(calibparavec) 
                file_name = "RO_Opt_GPP_EC_Nm_f_20220502_2"
                run_calibration(file_name,ranges,para2ind;
                PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=append!(collect(1:21),collect(45:84)))
                run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(22:44))
    end
    if ix == 3
        calibparavec = CalibParaVec(
                    (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                    (:Kₓₗ₀,0.0005,0.1),(:i,0.1,10.0),(:α_max,0.1,0.5),
                    (:τ,1.0,15.0),(:Smax,10.0,25.0,true),
                    (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
                ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
                ranges,para2ind = CreateOptVar(calibparavec) 
                file_name = "RO_Opt_GPP_EC_Nm_f_20220502_3"
                run_calibration(file_name,ranges,para2ind;
                PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=append!(collect(1:21),collect(45:84)))
                run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(22:44))
    end
    end
    =#
    
    #=
    Threads.@threads for ix = 1:4
    if ix == 1
        calibparavec = CalibParaVec(
            (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
            (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5),
            (:τ,1.0,15.0,true),
            (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
        ranges,para2ind = CreateOptVar(calibparavec) 
        file_name = "RO_Opt_GPP_EC_Nm_f_20220506_1"
        run_calibration(file_name,ranges,para2ind;
        PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=collect(22:84))
        run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(1:21))
    end
    if ix == 2
        calibparavec = CalibParaVec(
            (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
            (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5),
            (:τ,1.0,15.0,true),
            (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
        ranges,para2ind = CreateOptVar(calibparavec) 
        file_name = "RO_Opt_GPP_EC_Nm_f_20220506_2"
        run_calibration(file_name,ranges,para2ind;
        PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=append!(collect(1:21),collect(45:84)))
        run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(22:44))
    end
    if ix == 3
    calibparavec = CalibParaVec(
        (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
        (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5),
        (:τ,1.0,15.0,true),
        (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
    ranges,para2ind = CreateOptVar(calibparavec) 
    file_name = "RO_Opt_GPP_EC_Nm_f_20220506_3"
    run_calibration(file_name,ranges,para2ind;
    PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=append!(collect(1:44),collect(65:84)))
    run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(45:64))
    end
    if ix == 4    
        calibparavec = CalibParaVec(
            (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
            (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5),
            (:τ,1.0,15.0,true),
            (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
        ranges,para2ind = CreateOptVar(calibparavec) 
        file_name = "RO_Opt_GPP_EC_Nm_f_20220506_4"
        run_calibration(file_name,ranges,para2ind;
        PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=collect(1:64))
        run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(65:84))
    end
    end
    =#
    
    #=
    Threads.@threads for ix = 1:4
    if ix == 1
    parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:τ,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
        ranges = [(0.0001,0.1),(10.0,40.0),(0.01,1.0),(0.0005,0.1),
        (0.1,10.0),(0.1,0.5),(1.0,15.0),
        (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
        file_name = "RO_Opt_Separate_GPP_EC_Nm_f_20220505_1" 
        run_calibration(file_name,ranges,parasym;
        PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=collect(22:84))
        run_validation_RO(file_name,ranges,parasym;ParaDictInit_C=ParaDictInit_C,ind=collect(1:21))
    end
    if ix == 2
        parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:τ,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
        ranges = [(0.0001,0.1),(10.0,40.0),(0.01,1.0),(0.0005,0.1),
        (0.1,10.0),(0.1,0.5),(1.0,15.0),
        (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
        file_name = "RO_Opt_Separate_GPP_EC_Nm_f_20220505_2" 
        run_calibration(file_name,ranges,parasym;
        PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=append!(collect(1:21),collect(45:84)))
        run_validation_RO(file_name,ranges,parasym;ParaDictInit_C=ParaDictInit_C,ind=collect(22:44))
    end
    if ix == 3
        parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:τ,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
        ranges = [(0.0001,0.1),(10.0,40.0),(0.01,1.0),(0.0005,0.1),
        (0.1,10.0),(0.1,0.5),(1.0,15.0),
        (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
        file_name = "RO_Opt_Separate_GPP_EC_Nm_f_20220505_3" 
        run_calibration(file_name,ranges,parasym;
        PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=append!(collect(1:44),collect(65:84)))
        run_validation_RO(file_name,ranges,parasym;ParaDictInit_C=ParaDictInit_C,ind=collect(45:64))
    end
    if ix == 4
        parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:τ,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
        ranges = [(0.0001,0.1),(10.0,40.0),(0.01,1.0),(0.0005,0.1),
        (0.1,10.0),(0.1,0.5),(1.0,15.0),
        (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
        file_name = "RO_Opt_Separate_GPP_EC_Nm_f_20220505_4" 
        run_calibration(file_name,ranges,parasym;
        PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=collect(1:64))
        run_validation_RO(file_name,ranges,parasym;ParaDictInit_C=ParaDictInit_C,ind=collect(65:84))
    end
    end
    =#    
    #=
    Threads.@threads for ix = 1:4   
        if ix == 1
            calibparavec = CalibParaVec(
            (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
            (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5),
            (:τ,1.0,15.0,true),
            (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
            ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
            ranges,para2ind = CreateOptVar(calibparavec) 
            file_name = "RO_Opt_GPP_EC_Nm_f_20220506_4_weight_1"
            run_calibration(file_name,ranges,para2ind;
            PopulationSize=200,MaxSteps=30000,
            Calc_logP=(model::ModelResult,data::RoData,ParaDict::Dict{Symbol,Float64};ind::Array{Int64,1})->Calc_logP_GPP_Ec_Nm_f_weight(model,data,ParaDict;λ_gpp=1.4,ind=ind),
            ParaDictInit_C=ParaDictInit_C,ind=collect(1:64))
            run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(65:84))
        end
        if ix == 2
            calibparavec = CalibParaVec(
            (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
            (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5),
            (:τ,1.0,15.0,true),
            (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
            ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
            ranges,para2ind = CreateOptVar(calibparavec) 
            file_name = "RO_Opt_GPP_EC_Nm_f_20220506_4_weight_2"
            run_calibration(file_name,ranges,para2ind;
            PopulationSize=200,MaxSteps=30000,
            Calc_logP=(model::ModelResult,data::RoData,ParaDict::Dict{Symbol,Float64};ind::Array{Int64,1})->Calc_logP_GPP_Ec_Nm_f_weight(model,data,ParaDict;λ_gpp=1.8,ind=ind),
            ParaDictInit_C=ParaDictInit_C,ind=collect(1:64))
            run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(65:84))
        end
        if ix == 3
            calibparavec = CalibParaVec(
            (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
            (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5),
            (:τ,1.0,15.0,true),
            (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
            ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
            ranges,para2ind = CreateOptVar(calibparavec) 
            file_name = "RO_Opt_GPP_EC_Nm_f_20220506_4_weight_3"
            run_calibration(file_name,ranges,para2ind;
            PopulationSize=200,MaxSteps=30000,
            Calc_logP=(model::ModelResult,data::RoData,ParaDict::Dict{Symbol,Float64};ind::Array{Int64,1})->Calc_logP_GPP_Ec_Nm_f_weight(model,data,ParaDict;λ_gpp=2.0,ind=ind),
            ParaDictInit_C=ParaDictInit_C,ind=collect(1:64))
            run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(65:84))
        end
        if ix == 4
            calibparavec = CalibParaVec(
            (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
            (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5),
            (:τ,1.0,15.0,true),
            (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
            ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
            ranges,para2ind = CreateOptVar(calibparavec) 
            file_name = "RO_Opt_GPP_EC_Nm_f_20220506_4_weight_4"
            run_calibration(file_name,ranges,para2ind;
            PopulationSize=200,MaxSteps=30000,
            Calc_logP=(model::ModelResult,data::RoData,ParaDict::Dict{Symbol,Float64};ind::Array{Int64,1})->Calc_logP_GPP_Ec_Nm_f_weight(model,data,ParaDict;λ_gpp=2.4,ind=ind),
            ParaDictInit_C=ParaDictInit_C,ind=collect(1:64))
            run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(65:84))
        end
        end
        =#
        #=
        Threads.@threads for ix = 1:4
            if ix == 1
            parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:τ,:a_GPP,:b_GPP]
                ranges = [(0.0001,0.1),(10.0,40.0),(0.01,1.0),(0.0005,0.1),
                (0.1,10.0),(0.1,0.5),(1.0,15.0),
                (0.0001,5.0),(0.0001,3.0)]
                ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
                file_name = "RO_Opt_Separate_GPP_Nm_f_20220508_1" 
                run_calibration(file_name,ranges,parasym;
                PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=collect(22:84))
                run_validation_RO(file_name,ranges,parasym;ParaDictInit_C=ParaDictInit_C,ind=collect(1:21))
            end
            if ix == 2
                parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:τ,:a_GPP,:b_GPP]
                ranges = [(0.0001,0.1),(10.0,40.0),(0.01,1.0),(0.0005,0.1),
                (0.1,10.0),(0.1,0.5),(1.0,15.0),
                (0.0001,5.0),(0.0001,3.0)]
                ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
                file_name = "RO_Opt_Separate_GPP_Nm_f_20220506_2" 
                run_calibration(file_name,ranges,parasym;
                PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=append!(collect(1:21),collect(45:84)))
                run_validation_RO(file_name,ranges,parasym;ParaDictInit_C=ParaDictInit_C,ind=collect(22:44))
            end
            if ix == 3
                parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:τ,:a_GPP,:b_GPP]
                ranges = [(0.0001,0.1),(10.0,40.0),(0.01,1.0),(0.0005,0.1),
                (0.1,10.0),(0.1,0.5),(1.0,15.0),
                (0.0001,5.0),(0.0001,3.0)]
                ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
                file_name = "RO_Opt_Separate_GPP_Nm_f_20220508_3" 
                run_calibration(file_name,ranges,parasym;
                PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=append!(collect(1:44),collect(65:84)))
                run_validation_RO(file_name,ranges,parasym;ParaDictInit_C=ParaDictInit_C,ind=collect(45:64))
            end
            if ix == 4
                parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:τ,:a_GPP,:b_GPP]
                ranges = [(0.0001,0.1),(10.0,40.0),(0.01,1.0),(0.0005,0.1),
                (0.1,10.0),(0.1,0.5),(1.0,15.0),
                (0.0001,5.0),(0.0001,3.0)]
                ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
                file_name = "RO_Opt_Separate_GPP_Nm_f_20220508_4" 
                run_calibration(file_name,ranges,parasym;
                PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=collect(1:64))
                run_validation_RO(file_name,ranges,parasym;ParaDictInit_C=ParaDictInit_C,ind=collect(65:84))
            end
            end
            =#

            #=
            Threads.@threads for ix = 1:8
                if ix == 1
                    calibparavec = CalibParaVec(
                        (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                        (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5,true),
                        (:τ,1.0,15.0,true),                        
                        (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,5.0),(:b_Ec,0.0001,3.0)) 
                    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
                    ranges,para2ind = CreateOptVar(calibparavec) 
                    file_name = "RO_Opt_GPP_Ec_Nm_f_20220511_1"
                    run_calibration(file_name,ranges,para2ind;
                    PopulationSize=300,MaxSteps=40000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=collect(22:84))
                    run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(1:21))
                end
                if ix == 2
                    calibparavec = CalibParaVec(
                        (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                        (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5,true),
                        (:τ,1.0,15.0,true),
                        (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,5.0),(:b_Ec,0.0001,3.0)) 
                    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
                    ranges,para2ind = CreateOptVar(calibparavec) 
                    file_name = "RO_Opt_GPP_Ec_Nm_f_20220511_2"
                    run_calibration(file_name,ranges,para2ind;
                    PopulationSize=300,MaxSteps=40000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=append!(collect(1:21),collect(45:84)))
                    run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(22:44))
                end
                if ix == 3
                calibparavec = CalibParaVec(
                    (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                    (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5,true),  
                    (:τ,1.0,15.0,true),                  
                    (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,5.0),(:b_Ec,0.0001,3.0)) 
                ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
                ranges,para2ind = CreateOptVar(calibparavec) 
                file_name = "RO_Opt_GPP_Ec_Nm_f_20220511_3"
                run_calibration(file_name,ranges,para2ind;
                PopulationSize=300,MaxSteps=40000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=append!(collect(1:44),collect(65:84)))
                run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(45:64))
                end
                if ix == 4    
                    calibparavec = CalibParaVec(
                        (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                        (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5,true), 
                        (:τ,1.0,15.0,true),                       
                        (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,5.0),(:b_Ec,0.0001,3.0)) 
                    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
                    ranges,para2ind = CreateOptVar(calibparavec) 
                    file_name = "RO_Opt_GPP_Ec_Nm_f_20220511_4"
                    run_calibration(file_name,ranges,para2ind;
                    PopulationSize=300,MaxSteps=40000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=collect(1:64))
                    run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(65:84))
                end

                if ix == 5
                    calibparavec = CalibParaVec(
                        (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                        (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5,true),
                        (:τ,1.0,15.0,true),                        
                        (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0)) 
                    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
                    ranges,para2ind = CreateOptVar(calibparavec) 
                    file_name = "RO_Opt_GPP_Nm_f_20220511_1"
                    run_calibration(file_name,ranges,para2ind;
                    PopulationSize=300,MaxSteps=40000,Calc_logP=Calc_logP_GPP_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=collect(22:84))
                    run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(1:21))
                end
                if ix == 6
                    calibparavec = CalibParaVec(
                        (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                        (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5,true),
                        (:τ,1.0,15.0,true),
                        (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0)) 
                    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
                    ranges,para2ind = CreateOptVar(calibparavec) 
                    file_name = "RO_Opt_GPP_Nm_f_20220511_2"
                    run_calibration(file_name,ranges,para2ind;
                    PopulationSize=300,MaxSteps=40000,Calc_logP=Calc_logP_GPP_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=append!(collect(1:21),collect(45:84)))
                    run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(22:44))
                end
                if ix == 7
                calibparavec = CalibParaVec(
                    (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                    (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5,true),  
                    (:τ,1.0,15.0,true),                  
                    (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0)) 
                ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
                ranges,para2ind = CreateOptVar(calibparavec) 
                file_name = "RO_Opt_GPP_Nm_f_20220511_3"
                run_calibration(file_name,ranges,para2ind;
                PopulationSize=300,MaxSteps=40000,Calc_logP=Calc_logP_GPP_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=append!(collect(1:44),collect(65:84)))
                run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(45:64))
                end
                if ix == 8    
                    calibparavec = CalibParaVec(
                        (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                        (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5,true), 
                        (:τ,1.0,15.0,true),                       
                        (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0)) 
                    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
                    ranges,para2ind = CreateOptVar(calibparavec) 
                    file_name = "RO_Opt_GPP_Nm_f_20220511_4"
                    run_calibration(file_name,ranges,para2ind;
                    PopulationSize=300,MaxSteps=40000,Calc_logP=Calc_logP_GPP_Nm_f,ParaDictInit_C=ParaDictInit_C,ind=collect(1:64))
                    run_validation_RO(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind=collect(65:84))
                end
                end
                =#
    
    #=            
    calibparavec = CalibParaVec(
        (:Nₛ,0.0001,0.1,true),(:a_Jmax,0.01,1.0),
        (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),
        (:α_max,0.1,0.5),
        (:τ,1.0,15.0),(:Smax,10.0,25.0),                        
        (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,5.0),(:b_Ec,0.0001,3.0)) 
    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.0018)   
    ranges,para2ind = CreateOptVar(calibparavec) 
    file_name = "RO_Opt_GPP_Ec_Nm_f_20220604"
    
    run_CrossValidation(file_name,
    10,
    ranges,
    para2ind;
    PopulationSize=200,
    MaxSteps=30000,
    Calc_logP=Calc_logP_GPP_Ec_Nm_f,
    ParaDictInit_C=ParaDictInit_C)   
    =#
    
    #=
    calibparavec = CalibParaVec(
        (:Nₛ,0.0001,0.1,true),(:a_Jmax,0.01,1.0),
        (:Kₓₗ₀,0.0005,0.1,true),
        (:α_max,0.1,0.5),
        (:τ,1.0,15.0),(:Smax,10.0,25.0),                        
        (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0))#,(:a_Ec,0.0001,5.0),(:b_Ec,0.0001,3.0)) 
    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.0018,:a_Ec=>0.0433,:b_Ec=>0.179)   
    ParaDictInit_F = Dict(:μ_Nₘ_f=>0.018,:b_Nₘ_f=>0.0018,:a_Ec=>0.0688,:b_Ec=>0.146) 
    ranges,para2ind = CreateOptVar(calibparavec) 
    file_name = "RO_Opt_GPP_Ec_Nm_f_20230113"
    
    run_CrossValidation(file_name,
    10,
    ranges,
    para2ind;
    PopulationSize=350,
    MaxSteps=35000,
    Calc_logP=Calc_logP_GPP_Ec_Nm_f,
    ParaDictInit_C=ParaDictInit_C,
    ParaDictInit_F=ParaDictInit_F)   
    =# 

    #=
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
    file_name = "RO_Opt_Separate_GPP_Ec_Nm_f_20230114"
    run_CrossValidation(file_name,
    10,
    ranges,
    parasym;
    PopulationSize=350,
    MaxSteps=35000,
    Calc_logP=Calc_logP_GPP_Ec_Nm_f,
    ParaDictInit_C=ParaDictInit_C,
    ParaDictInit_F=ParaDictInit_F)  
    =#

    #=
    parasym = [:Nₛ,
    :a_Jmax,
    :Kₓₗ₀,
    :i,
    :α_max,
    :τ,
    :a_GPP,
    :b_GPP,
    :a_Ec,
    :b_Ec]
    ranges = [(0.0001,0.1),
    (0.01,1.0),
    (0.0005,0.1),
    (0.1,10.0),
    (0.1,0.5),
    (1.0,15.0),
    (0.0001,5.0),
    (0.0001,3.0),
    (0.0001,5.0),
    (0.0001,3.0)]
    
    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
    file_name = "RO_Opt_Separate_GPP_Ec_Nm_f_20220602"
    run_CrossValidation(file_name,
    10,
    ranges,
    parasym;
    PopulationSize=200,
    MaxSteps=30000,
    Calc_logP=Calc_logP_GPP_Ec_Nm_f,
    ParaDictInit_C=ParaDictInit_C)  
    =#

    #=
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
    file_name = "RO_Opt_Separate_GPP_Ec_Nm_f_eco_scaling_20230731"

    run_CrossValidation(file_name,
    10,
    ranges,
    parasym;
    PopulationSize=350,
    MaxSteps=35000,
    Calc_logP=Calc_logP_GPP_Ec_Nm_f_eco_scaling,
    ParaDictInit_C=ParaDictInit_C,
    ParaDictInit_F=ParaDictInit_F)
    #---Run non-sharing parameter case, Done--- 
    =#

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
    file_name = "RO_Opt_GPP_Ec_Nm_f_eco_scaling_20230801"
    
    run_CrossValidation(file_name,
    10,
    ranges,
    para2ind;
    PopulationSize=400,
    MaxSteps=45000,
    Calc_logP=Calc_logP_GPP_Ec_Nm_f_eco_scaling,
    ParaDictInit_C=ParaDictInit_C,
    ParaDictInit_F=ParaDictInit_F)
    #---Run sharing parameter case, Done---
    
    return nothing
end

work_list()