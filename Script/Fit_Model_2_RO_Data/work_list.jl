function work_list()
    #=
    Threads.@threads for ix = 1:12 

    if ix == 1
        calibparavec = CalibParaVec(
            (:Nₛ,0.001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
            (:Kₓₗ₀,0.001,0.1),(:i,0.1,6.0),(:α_max,0.1,0.5),
            (:X0,-18.0,-0.9),(:τ,1.0,15.0),(:Smax,10.0,25.0),
            (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0))        
        ranges,para2ind = CreateOptVar(calibparavec) 
        file_name = "RO_Opt_GPP_EC_20220315" 
        run_calibration(file_name,ranges,para2ind;
        PopulationSize=100,MaxSteps=20000,Calc_logP=Calc_logP_GPP_Ec) 
    end   
    
    if ix == 2
        calibparavec = CalibParaVec(
            (:Nₛ,0.001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
            (:Kₓₗ₀,0.001,0.1),(:i,0.1,6.0),(:α_max,0.1,0.5),
            (:X0,-18.0,-0.9),(:τ,1.0,15.0,true),(:Smax,10.0,25.0),
            (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0))        
        ranges,para2ind = CreateOptVar(calibparavec) 
        file_name = "RO_Opt_GPP_EC_20220315_2" 
        run_calibration(file_name,ranges,para2ind;
        PopulationSize=100,MaxSteps=20000,Calc_logP=Calc_logP_GPP_Ec) 
    end

    if ix == 3
        calibparavec = CalibParaVec(
            (:Nₛ,0.001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
            (:Kₓₗ₀,0.001,0.1,true),(:i,0.1,6.0),(:α_max,0.1,0.5),
            (:X0,-18.0,-0.9),(:τ,1.0,15.0,true),(:Smax,10.0,25.0),
            (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0))        
        ranges,para2ind = CreateOptVar(calibparavec) 
        file_name = "RO_Opt_GPP_EC_20220315_3" 
        run_calibration(file_name,ranges,para2ind;
        PopulationSize=100,MaxSteps=20000,Calc_logP=Calc_logP_GPP_Ec) 
    end

    if ix == 4
        calibparavec = CalibParaVec(
            (:Nₛ,0.001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
            (:Kₓₗ₀,0.001,0.1),(:i,0.1,6.0),(:α_max,0.1,0.5),
            (:X0,-18.0,-0.9),(:τ,1.0,15.0),(:Smax,10.0,25.0),
            (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
        ranges,para2ind = CreateOptVar(calibparavec) 
        file_name = "RO_Opt_GPP_EC_Nm_f_20220315" 
        run_calibration(file_name,ranges,para2ind;
        PopulationSize=100,MaxSteps=20000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C) 
    end  

    if ix == 5
        calibparavec = CalibParaVec(
            (:Nₛ,0.001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
            (:Kₓₗ₀,0.001,0.1),(:i,0.1,6.0),(:α_max,0.1,0.5),
            (:X0,-18.0,-0.9),(:τ,1.0,15.0,true),(:Smax,10.0,25.0),
            (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
        ranges,para2ind = CreateOptVar(calibparavec) 
        file_name = "RO_Opt_GPP_EC_Nm_f_20220315_2" 
        run_calibration(file_name,ranges,para2ind;
        PopulationSize=100,MaxSteps=20000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C) 
    end  

    if ix == 6
        calibparavec = CalibParaVec(
            (:Nₛ,0.001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
            (:Kₓₗ₀,0.001,0.1,true),(:i,0.1,6.0),(:α_max,0.1,0.5),
            (:X0,-18.0,-0.9),(:τ,1.0,15.0,true),(:Smax,10.0,25.0),
            (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
        ranges,para2ind = CreateOptVar(calibparavec) 
        file_name = "RO_Opt_GPP_EC_Nm_f_20220315_3" 
        run_calibration(file_name,ranges,para2ind;
        PopulationSize=100,MaxSteps=20000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C) 
    end  
        
    if ix == 7
        parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:X0,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
        ranges = [(0.001,0.1),(10.0,40.0),(0.01,1.0),(0.001,0.1),
        (0.1,6.0),(0.1,0.5),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),
        (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]        
        file_name = "RO_Opt_Separate_GPP_EC_20220315" 
        run_calibration(file_name,ranges,parasym;
        PopulationSize=100,MaxSteps=20000,Calc_logP=Calc_logP_GPP_Ec)
    end  

    if ix == 8
        parasym = [:Nₛ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:X0,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
        ranges = [(0.001,0.1),(0.01,1.0),(0.001,0.1),
        (0.1,6.0),(0.1,0.5),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),
        (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]        
        file_name = "RO_Opt_Separate_GPP_EC_20220315_2" 
        run_calibration(file_name,ranges,parasym;
        PopulationSize=100,MaxSteps=20000,Calc_logP=Calc_logP_GPP_Ec)
    end 

    if ix == 9
        parasym = [:Nₛ,:a_Jmax,:Kₓₗ₀,:α_max,:X0,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
        ranges = [(0.001,0.1),(0.01,1.0),(0.001,0.1),
        (0.1,0.5),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),
        (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]        
        file_name = "RO_Opt_Separate_GPP_EC_20220315_3" 
        run_calibration(file_name,ranges,parasym;
        PopulationSize=100,MaxSteps=20000,Calc_logP=Calc_logP_GPP_Ec)
    end 

    if ix == 10
        parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:X0,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
        ranges = [(0.001,0.1),(10.0,40.0),(0.01,1.0),(0.001,0.1),
        (0.1,6.0),(0.1,0.5),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),
        (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
        file_name = "RO_Opt_Separate_GPP_EC_Nm_f_20220315" 
        run_calibration(file_name,ranges,parasym;
        PopulationSize=100,MaxSteps=20000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C)
    end  
    
    if ix == 11
        parasym = [:Nₛ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:X0,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
        ranges = [(0.001,0.1),(0.01,1.0),(0.001,0.1),
        (0.1,6.0),(0.1,0.5),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),
        (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
        file_name = "RO_Opt_Separate_GPP_EC_Nm_f_20220315_2" 
        run_calibration(file_name,ranges,parasym;
        PopulationSize=100,MaxSteps=20000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C)
    end 
     
    if ix == 12
        parasym = [:Nₛ,:a_Jmax,:Kₓₗ₀,:α_max,:X0,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
        ranges = [(0.001,0.1),(0.01,1.0),(0.001,0.1),
        (0.1,0.5),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),
        (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
        ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
        file_name = "RO_Opt_Separate_GPP_EC_Nm_f_20220315_3" 
        run_calibration(file_name,ranges,parasym;
        PopulationSize=100,MaxSteps=20000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C)
    end  

    end   
    =#
    
    #=
    Threads.@threads for ix = 1:3
        if ix == 1
            parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:X0,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
            ranges = [(0.001,0.1),(10.0,40.0),(0.01,1.0),(0.001,0.1),
            (0.1,6.0),(0.1,0.5),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),
            (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
            ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
            file_name = "RO_Opt_Separate_GPP_EC_Nm_f_20220316" 
            run_calibration(file_name,ranges,parasym;
            PopulationSize=200,MaxSteps=30000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C)
        end
        if ix == 2
            parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:X0,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
            ranges = [(0.001,0.1),(10.0,40.0),(0.01,1.0),(0.001,0.1),
            (0.1,6.0),(0.1,0.5),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),
            (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
            ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
            file_name = "RO_Opt_Separate_GPP_EC_20220316" 
            run_calibration(file_name,ranges,parasym;
            PopulationSize=200,MaxSteps=20000,Calc_logP=Calc_logP_GPP_Ec,ParaDictInit_C=ParaDictInit_C)
        end
        if ix == 3
            calibparavec = CalibParaVec(
                (:Nₛ,0.001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                (:Kₓₗ₀,0.001,0.1),(:i,0.1,6.0),(:α_max,0.1,0.5),
                (:X0,-18.0,-0.9),(:τ,1.0,15.0),(:Smax,10.0,25.0,true),
                (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
            ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
            ranges,para2ind = CreateOptVar(calibparavec) 
            file_name = "RO_Opt_GPP_EC_Nm_f_20220316" 
            run_calibration(file_name,ranges,para2ind;
            PopulationSize=200,MaxSteps=20000,Calc_logP=Calc_logP_GPP_Ec_Nm_f,ParaDictInit_C=ParaDictInit_C)
        end
        if ix == 4
            calibparavec = CalibParaVec(
                (:Nₛ,0.001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                (:Kₓₗ₀,0.001,0.1),(:i,0.1,6.0),(:α_max,0.1,0.5),
                (:X0,-18.0,-0.9),(:τ,1.0,15.0),(:Smax,10.0,25.0,true),
                (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
            ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
            ranges,para2ind = CreateOptVar(calibparavec) 
            file_name = "RO_Opt_GPP_EC_20220316" 
            run_calibration(file_name,ranges,para2ind;
            PopulationSize=200,MaxSteps=20000,Calc_logP=Calc_logP_GPP_Ec,ParaDictInit_C=ParaDictInit_C)
        end
    end 
    =#   
    Threads.@threads for ix = 1:3
        if ix == 1
            calibparavec = CalibParaVec(
                (:Nₛ,0.001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                (:Kₓₗ₀,0.001,0.1),(:i,0.1,6.0),(:α_max,0.1,0.5),
                (:X0,-18.0,-0.9),(:τ,1.0,15.0),(:Smax,10.0,25.0),
                (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
            ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
            ranges,para2ind = CreateOptVar(calibparavec) 
            file_name = "RO_Opt_GPP_EC_20220316" 
            run_calibration(file_name,ranges,para2ind;
            PopulationSize=100,MaxSteps=7000,Calc_logP=Calc_logP_GPP_Ec,ParaDictInit_C=ParaDictInit_C)                
        end
        if ix == 2
            calibparavec = CalibParaVec(
                (:Nₛ,0.001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                (:Kₓₗ₀,0.001,0.1),(:i,0.1,6.0),(:α_max,0.1,0.5),
                (:X0,-18.0,-0.9),(:τ,1.0,15.0,true),(:Smax,10.0,25.0),
                (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
            ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
            ranges,para2ind = CreateOptVar(calibparavec) 
            file_name = "RO_Opt_GPP_EC_20220316_2" 
            run_calibration(file_name,ranges,para2ind;
            PopulationSize=100,MaxSteps=7000,Calc_logP=Calc_logP_GPP_Ec,ParaDictInit_C=ParaDictInit_C)  
        end
        if ix == 3
            calibparavec = CalibParaVec(
                (:Nₛ,0.001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
                (:Kₓₗ₀,0.001,0.1,true),(:i,0.1,6.0),(:α_max,0.1,0.5),
                (:X0,-18.0,-0.9),(:τ,1.0,15.0,true),(:Smax,10.0,25.0),
                (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,3.0),(:b_Ec,0.0001,3.0)) 
            ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
            ranges,para2ind = CreateOptVar(calibparavec) 
            file_name = "RO_Opt_GPP_EC_20220316_3" 
            run_calibration(file_name,ranges,para2ind;
            PopulationSize=100,MaxSteps=7000,Calc_logP=Calc_logP_GPP_Ec,ParaDictInit_C=ParaDictInit_C)  
        end
    end
end

work_list()