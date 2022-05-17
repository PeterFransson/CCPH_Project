function log_post_distri_RO_CCPH_sep(para::Array{Float64,1},parasym::Array{Symbol,1},ranges::Array{Tuple{Float64,Float64},1},
    RO_data::RO_raw_data;stand_type::String="Fertilized",Calc_logP::Function=Calc_logP_GPP_Ec,
    ParaDictInit=nothing)
    if any(x->x==false,first.(ranges).<para.<last.(ranges)) 
        println("Boundary Error")      
        post = [-Inf,-Inf,-Inf]
        return post
    else        
        try             
            ParaDict = CreateParaDict(parasym,para;ParaDictInit=ParaDictInit)

            model,weatherts,data = Get_Result_RO_CCPH(ParaDict,RO_data;stand_type=stand_type)
                
            logP = Calc_logP(model,data,ParaDict) 
            
            return logP      
        catch err
            println("Parameters Error: ", err)
            println("Parameter=$(para)")
            
            post = [-10^4,-10^4,-10^4]
            return post
        end
    end           
end

function log_post_distri_RO_CCPH_combined(para::Array{Float64,1},parasym::Array{Symbol,1},ranges::Array{Tuple{Float64,Float64},1},
    RO_data::RO_raw_data;stand_type::String="Fertilized",Calc_logP::Function=Calc_logP_GPP_Ec,
    ParaDictInit=nothing)
    if any(x->x==false,first.(ranges).<para.<last.(ranges)) 
        println("Boundary Error")      
        post = [-Inf,-Inf,-Inf]
        return post
    else        
        try             
            ParaDict = CreateParaDict(parasym,para;ParaDictInit=ParaDictInit)

            model,weatherts,data = Get_Result_RO_CCPH(ParaDict,RO_data;stand_type=stand_type)
                
            logP = Calc_logP(model,data,ParaDict) 
            
            return logP      
        catch err
            println("Parameters Error: ", err)
            println("Parameter=$(para)")
            
            post = -10^4
            return post
        end
    end           
end

function run_sensitivity_analysis()
    file_name = "RO_Opt_Separate_GPP_EC_Nm_f_sep_20220429" 

    parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:X0,:τ,:Smax,:a_GPP,:b_GPP,:a_Ec,:b_Ec]
    ranges = [(0.001,0.1),(10.0,40.0),(0.01,1.0),(0.0006,0.1),
    (0.1,8.0),(0.1,0.5),(-18.0,-0.9),(1.0,15.0),(10.0,25.0),
    (0.0001,5.0),(0.0001,3.0),(0.0001,5.0),(0.0001,3.0)]
    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
    RO_data = Load_RO_data() 
    
    Threads.@threads for ix = 1:2
    if ix == 1
    m = DiffEqSensitivity.gsa(x->log_post_distri_RO_CCPH_sep(x,parasym,ranges,RO_data;
    stand_type="Fertilized",Calc_logP=Calc_logP_GPP_Ec_Nm_f_sep),DiffEqSensitivity.Sobol(),
    ranges,N=2000)  
    
    @show m.ST

    bar(m.ST[1,:],title="Total Order Indices GPP",legend=false)
    p1 = bar(string.(parasym),m.ST[1,:],title="Total Order Indices GPP",legend=false)
    p2 = bar(string.(parasym),m.ST[2,:],title="Total Order Indices Ec",legend=false)
    p3 = bar(string.(parasym),m.ST[3,:],title="Total Order Indices Nm_f",legend=false)
    p1_ = bar(string.(parasym),m.S1[1,:],title="First Order Indices GPP",legend=false) 
    p2_ = bar(string.(parasym),m.S1[2,:],title="First Order Indices Ec",legend=false)
    p3_ = bar(string.(parasym),m.S1[3,:],title="First Order Indices Nm_f",legend=false)
    plot(p1,p2,p3,p1_,p2_,p3_)
    savefig("./plots/"*file_name*"_SA_F.svg")
    end  
    if ix == 2  
    m = DiffEqSensitivity.gsa(x->log_post_distri_RO_CCPH_sep(x,parasym,ranges,RO_data;
    stand_type="Control",Calc_logP=Calc_logP_GPP_Ec_Nm_f_sep,ParaDictInit=ParaDictInit_C),DiffEqSensitivity.Sobol(),
    ranges,N=2000)   

    @show m.ST
    
    p1 = bar(string.(parasym),m.ST[1,:],title="Total Order Indices GPP",legend=false)
    p2 = bar(string.(parasym),m.ST[2,:],title="Total Order Indices Ec",legend=false)
    p3 = bar(string.(parasym),m.ST[3,:],title="Total Order Indices Nm_f",legend=false)
    p1_ = bar(string.(parasym),m.S1[1,:],title="First Order Indices GPP",legend=false) 
    p2_ = bar(string.(parasym),m.S1[2,:],title="First Order Indices Ec",legend=false)
    p3_ = bar(string.(parasym),m.S1[3,:],title="First Order Indices Nm_f",legend=false)
    plot(p1,p2,p3,p1_,p2_,p3_)
    savefig("./plots/"*file_name*"_SA_C.svg")
    end
    end
end

run_sensitivity_analysis()