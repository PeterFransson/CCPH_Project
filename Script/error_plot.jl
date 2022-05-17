mutable struct ModelStat{T<:Float64}
    R²::T
    PearsonCor::T
    MAPE::T
    MSE::T
    RMSE::T
end

mutable struct ModelStatSum{T<:ModelStat}
    GPP::T
    Ec::T
end

function Get_Model_Stat(ParaDict::Dict{Symbol,Float64};
    stand_type::String="Fertilized",ind_val::Array{Int64,1}=collect(1:21))
    
    ind_train = setdiff(1:84, ind_val)
    RO_data = Load_RO_data()
    model,weatherts,data = Get_Result_RO_CCPH(ParaDict,RO_data;stand_type=stand_type)

    modeltrain_GPP = ModelStat(calcR²(data.GPP[ind_train],model.GPP[ind_train]),
    cor(data.GPP[ind_train],model.GPP[ind_train]),
    calcMAPE(data.GPP[ind_train],model.GPP[ind_train]),
    calcMSE(data.GPP[ind_train],model.GPP[ind_train]),
    calcRMSE(data.GPP[ind_train],model.GPP[ind_train]))
    
    modelval_GPP = ModelStat(calcR²(data.GPP[ind_val],model.GPP[ind_val]),
    cor(data.GPP[ind_val],model.GPP[ind_val]),
    calcMAPE(data.GPP[ind_val],model.GPP[ind_val]),
    calcMSE(data.GPP[ind_val],model.GPP[ind_val]),
    calcRMSE(data.GPP[ind_val],model.GPP[ind_val]))

    modeltrain_Ec = ModelStat(calcR²(data.Ec[ind_train],model.Ec[ind_train]),
    cor(data.Ec[ind_train],model.Ec[ind_train]),
    calcMAPE(data.Ec[ind_train],model.Ec[ind_train]),
    calcMSE(data.Ec[ind_train],model.Ec[ind_train]),
    calcRMSE(data.Ec[ind_train],model.Ec[ind_train]))
    
    modelval_Ec = ModelStat(calcR²(data.Ec[ind_val],model.Ec[ind_val]),
    cor(data.Ec[ind_val],model.Ec[ind_val]),
    calcMAPE(data.Ec[ind_val],model.Ec[ind_val]),
    calcMSE(data.Ec[ind_val],model.Ec[ind_val]),
    calcRMSE(data.Ec[ind_val],model.Ec[ind_val]))

    modeltrain = ModelStatSum(modeltrain_GPP,modeltrain_Ec)
    modelval = ModelStatSum(modelval_GPP,modelval_Ec)    

    return (modeltrain,modelval)
end

function Get_Model_Stat(file_name::String,ranges::Array{Tuple{Float64, Float64},1},
    parasym::Array{Symbol,1};ParaDictInit_F=nothing,ParaDictInit_C=nothing,ind_val::Array{Int64,1}=collect(1:21))
    file_name_F = file_name*"_F"
    file_name_C = file_name*"_C"

    par_F = load("./output/"*file_name_F*".jld","xopt")     
    par_C = load("./output/"*file_name_C*".jld","xopt") 

    ParaDict_F = CreateParaDict(parasym,par_F;ParaDictInit=ParaDictInit_F)   
    modeltrain_F,modelval_F = Get_Model_Stat(ParaDict_F;stand_type="Fertilized",ind_val=ind_val)

    ParaDict_C = CreateParaDict(parasym,par_C;ParaDictInit=ParaDictInit_C)     
    modeltrain_C,modelval_C = Get_Model_Stat(ParaDict_C;stand_type="Control",ind_val=ind_val)

    return (modeltrain_F,modelval_F,modeltrain_C,modelval_C)
end

function Get_Model_Stat(file_name::String,ranges::Array{Tuple{Float64, Float64},1},
    para2ind::Dict{Symbol,Any};ParaDictInit_F=nothing,ParaDictInit_C=nothing,ind_val::Array{Int64,1}=collect(1:21))
    
    par = load("./output/"*file_name*".jld","xopt")   

    ParaDict_F,ParaDict_C = CreateParaDict(para2ind,par;
    ParaDictInit_F=ParaDictInit_F,ParaDictInit_C=ParaDictInit_C)
      
    modeltrain_F,modelval_F = Get_Model_Stat(ParaDict_F;stand_type="Fertilized",ind_val=ind_val)       
    modeltrain_C,modelval_C = Get_Model_Stat(ParaDict_C;stand_type="Control",ind_val=ind_val)

    return (modeltrain_F,modelval_F,modeltrain_C,modelval_C)
end

function no_shared_parameters()
    modeltrain_F = Array{ModelStatSum,1}(undef,4)
    modelval_F = Array{ModelStatSum,1}(undef,4)
    modeltrain_C = Array{ModelStatSum,1}(undef,4)
    modelval_C = Array{ModelStatSum,1}(undef,4)

    parasym = [:Nₛ,:rₘ,:a_Jmax,:Kₓₗ₀,:i,:α_max,:τ,:a_GPP,:b_GPP]
    ranges = [(0.0001,0.1),(10.0,40.0),(0.01,1.0),(0.0005,0.1),
    (0.1,10.0),(0.1,0.5),(1.0,15.0),
    (0.0001,5.0),(0.0001,3.0)]
    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)
    
    file_name = "RO_Opt_Separate_GPP_Nm_f_20220508_1" 
    modeltrain_F[1],modelval_F[1],modeltrain_C[1],modelval_C[1] = Get_Model_Stat(file_name,ranges,parasym;ParaDictInit_C=ParaDictInit_C,ind_val=collect(1:21))

    file_name = "RO_Opt_Separate_GPP_Nm_f_20220506_2" 
    modeltrain_F[2],modelval_F[2],modeltrain_C[2],modelval_C[2] = Get_Model_Stat(file_name,ranges,parasym;ParaDictInit_C=ParaDictInit_C,ind_val=collect(22:44))

    file_name = "RO_Opt_Separate_GPP_Nm_f_20220508_3" 
    modeltrain_F[3],modelval_F[3],modeltrain_C[3],modelval_C[3] = Get_Model_Stat(file_name,ranges,parasym;ParaDictInit_C=ParaDictInit_C,ind_val=collect(45:64))

    file_name = "RO_Opt_Separate_GPP_Nm_f_20220508_4" 
    modeltrain_F[4],modelval_F[4],modeltrain_C[4],modelval_C[4] = Get_Model_Stat(file_name,ranges,parasym;ParaDictInit_C=ParaDictInit_C,ind_val=collect(65:84))

    return (modeltrain_F,modelval_F,modeltrain_C,modelval_C)
end

function shared_parameters()
    modeltrain_F = Array{ModelStatSum,1}(undef,4)
    modelval_F = Array{ModelStatSum,1}(undef,4)
    modeltrain_C = Array{ModelStatSum,1}(undef,4)
    modelval_C = Array{ModelStatSum,1}(undef,4)

    calibparavec = CalibParaVec(
        (:Nₛ,0.0001,0.1,true),(:rₘ,10.0,40.0),(:a_Jmax,0.01,1.0),
        (:Kₓₗ₀,0.0005,0.1,true),(:i,0.1,10.0),(:α_max,0.1,0.5,true),        
        (:a_GPP,0.0001,5.0),(:b_GPP,0.0001,3.0),(:a_Ec,0.0001,5.0),(:b_Ec,0.0001,3.0)) 
    ParaDictInit_C = Dict(:μ_Nₘ_f=>0.0113,:b_Nₘ_f=>0.000634)       
    ranges,para2ind = CreateOptVar(calibparavec) 

    file_name = "RO_Opt_GPP_Nm_f_20220511_1"
    modeltrain_F[1],modelval_F[1],modeltrain_C[1],modelval_C[1] = Get_Model_Stat(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind_val=collect(1:21))

    file_name = "RO_Opt_GPP_Nm_f_20220511_2"
    modeltrain_F[2],modelval_F[2],modeltrain_C[2],modelval_C[2] = Get_Model_Stat(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind_val=collect(22:44))

    file_name = "RO_Opt_GPP_Nm_f_20220511_3"
    modeltrain_F[3],modelval_F[3],modeltrain_C[3],modelval_C[3] = Get_Model_Stat(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind_val=collect(45:64))

    file_name = "RO_Opt_GPP_Nm_f_20220511_4"
    modeltrain_F[4],modelval_F[4],modeltrain_C[4],modelval_C[4] = Get_Model_Stat(file_name,ranges,para2ind;ParaDictInit_C=ParaDictInit_C,ind_val=collect(65:84))

    return (modeltrain_F,modelval_F,modeltrain_C,modelval_C)
end

function create_statplot(modeltrain::T,modelval::T) where {T<:Array{ModelStatSum,1}}
    mn_R²_GPP = [modeltrain[1].GPP.R², modeltrain[2].GPP.R², modeltrain[3].GPP.R²,modeltrain[4].GPP.R²,
    modelval[1].GPP.R², modelval[2].GPP.R², modelval[3].GPP.R², modelval[4].GPP.R²]
    
    mn_MAPE_GPP = [modeltrain[1].GPP.MAPE, modeltrain[2].GPP.MAPE, modeltrain[3].GPP.MAPE,modeltrain[4].GPP.MAPE,
    modelval[1].GPP.MAPE, modelval[2].GPP.MAPE, modelval[3].GPP.MAPE, modelval[4].GPP.MAPE]

    mn_RMSE_GPP = [modeltrain[1].GPP.RMSE, modeltrain[2].GPP.RMSE, modeltrain[3].GPP.RMSE,modeltrain[4].GPP.RMSE,
    modelval[1].GPP.RMSE, modelval[2].GPP.RMSE, modelval[3].GPP.RMSE, modelval[4].GPP.RMSE]

    mn_cor_GPP = [modeltrain[1].GPP.PearsonCor, modeltrain[2].GPP.PearsonCor,
    modeltrain[3].GPP.PearsonCor,modeltrain[4].GPP.PearsonCor,
    modelval[1].GPP.PearsonCor, modelval[2].GPP.PearsonCor,
    modelval[3].GPP.PearsonCor, modelval[4].GPP.PearsonCor]

    mn_R²_Ec = [modeltrain[1].Ec.R², modeltrain[2].Ec.R², modeltrain[3].Ec.R²,modeltrain[4].Ec.R²,
    modelval[1].Ec.R², modelval[2].Ec.R², modelval[3].Ec.R², modelval[4].Ec.R²]
    
    mn_MAPE_Ec = [modeltrain[1].Ec.MAPE, modeltrain[2].Ec.MAPE, modeltrain[3].Ec.MAPE,modeltrain[4].Ec.MAPE,
    modelval[1].Ec.MAPE, modelval[2].Ec.MAPE, modelval[3].Ec.MAPE, modelval[4].Ec.MAPE]

    mn_RMSE_Ec = [modeltrain[1].Ec.RMSE, modeltrain[2].Ec.RMSE, modeltrain[3].Ec.RMSE,modeltrain[4].Ec.RMSE,
    modelval[1].Ec.RMSE, modelval[2].Ec.RMSE, modelval[3].Ec.RMSE, modelval[4].Ec.RMSE]

    mn_cor_Ec = [modeltrain[1].Ec.PearsonCor, modeltrain[2].Ec.PearsonCor,
    modeltrain[3].Ec.PearsonCor,modeltrain[4].Ec.PearsonCor,
    modelval[1].Ec.PearsonCor, modelval[2].Ec.PearsonCor,
    modelval[3].Ec.PearsonCor, modelval[4].Ec.PearsonCor]
    
    sx = repeat(["Training", "Validation"], inner = 4)    
    nam = repeat(string.(2015:2018), outer = 2)

    p_R²_GPP = groupedbar(nam, mn_R²_GPP, group = sx, ylabel = "R²",title="GPP")
    p_MAPE_GPP = groupedbar(nam, mn_MAPE_GPP, group = sx, ylabel = "MAPE")
    p_RMSE_GPP = groupedbar(nam, mn_RMSE_GPP, group = sx, ylabel = "RMSE")  
    p_cor_GPP = groupedbar(nam, mn_cor_GPP, group = sx, ylabel = "cor")  
    
    p_R²_Ec = groupedbar(nam, mn_R²_Ec, group = sx, ylabel = "R²",title="Ec")
    p_MAPE_Ec = groupedbar(nam, mn_MAPE_Ec, group = sx, ylabel = "MAPE")
    p_RMSE_Ec = groupedbar(nam, mn_RMSE_Ec, group = sx, ylabel = "RMSE")  
    p_cor_Ec = groupedbar(nam, mn_cor_Ec, group = sx, ylabel = "cor") 

    statplot = plot(p_R²_GPP,p_R²_Ec,p_MAPE_GPP,p_MAPE_Ec,
    p_RMSE_GPP,p_RMSE_Ec,p_cor_GPP,p_cor_Ec,layout=(4,2),legends=false)

    return statplot
end

function create_statplot(file_name::String,modeltrain::T,modelval::T) where {T<:Array{ModelStatSum,1}}
    GPPR² = [modeltrain[1].GPP.R², modeltrain[2].GPP.R², modeltrain[3].GPP.R²,modeltrain[4].GPP.R²,
    modelval[1].GPP.R², modelval[2].GPP.R², modelval[3].GPP.R², modelval[4].GPP.R²]
    GPPR² = round.(GPPR²,digits=2)

    GPPMAPE = [modeltrain[1].GPP.MAPE, modeltrain[2].GPP.MAPE, modeltrain[3].GPP.MAPE,modeltrain[4].GPP.MAPE,
    modelval[1].GPP.MAPE, modelval[2].GPP.MAPE, modelval[3].GPP.MAPE, modelval[4].GPP.MAPE]
    GPPMAPE = round.(GPPMAPE,digits=2)

    GPPRMSE = [modeltrain[1].GPP.RMSE, modeltrain[2].GPP.RMSE, modeltrain[3].GPP.RMSE,modeltrain[4].GPP.RMSE,
    modelval[1].GPP.RMSE, modelval[2].GPP.RMSE, modelval[3].GPP.RMSE, modelval[4].GPP.RMSE]
    GPPRMSE = round.(GPPRMSE,digits=2)

    GPPcor = [modeltrain[1].GPP.PearsonCor, modeltrain[2].GPP.PearsonCor,
    modeltrain[3].GPP.PearsonCor,modeltrain[4].GPP.PearsonCor,
    modelval[1].GPP.PearsonCor, modelval[2].GPP.PearsonCor,
    modelval[3].GPP.PearsonCor, modelval[4].GPP.PearsonCor]
    GPPcor = round.(GPPcor,digits=2)

    EcR² = [modeltrain[1].Ec.R², modeltrain[2].Ec.R², modeltrain[3].Ec.R²,modeltrain[4].Ec.R²,
    modelval[1].Ec.R², modelval[2].Ec.R², modelval[3].Ec.R², modelval[4].Ec.R²]
    EcR² = round.(EcR²,digits=2)
    
    EcMAPE = [modeltrain[1].Ec.MAPE, modeltrain[2].Ec.MAPE, modeltrain[3].Ec.MAPE,modeltrain[4].Ec.MAPE,
    modelval[1].Ec.MAPE, modelval[2].Ec.MAPE, modelval[3].Ec.MAPE, modelval[4].Ec.MAPE]
    EcMAPE = round.(EcMAPE,digits=2)

    EcRMSE = [modeltrain[1].Ec.RMSE, modeltrain[2].Ec.RMSE, modeltrain[3].Ec.RMSE,modeltrain[4].Ec.RMSE,
    modelval[1].Ec.RMSE, modelval[2].Ec.RMSE, modelval[3].Ec.RMSE, modelval[4].Ec.RMSE]
    EcRMSE = round.(EcRMSE,digits=2)

    Eccor = [modeltrain[1].Ec.PearsonCor, modeltrain[2].Ec.PearsonCor,
    modeltrain[3].Ec.PearsonCor,modeltrain[4].Ec.PearsonCor,
    modelval[1].Ec.PearsonCor, modelval[2].Ec.PearsonCor,
    modelval[3].Ec.PearsonCor, modelval[4].Ec.PearsonCor]
    Eccor = round.(Eccor,digits=2)

    open("./output/"*file_name*".txt","w") do io
        println(io,"--Statistics--")
        println(io,"Info: validation values in parentheses")
        println(io,"");println(io,"")
        println(io,"GPP:")
        println(io,"")
        println(io,"Year: 2015, 2016, 2017, 2018")
        println(io,"R²: $(GPPR²[1])($(GPPR²[5])), $(GPPR²[2])($(GPPR²[6])), $(GPPR²[3])($(GPPR²[7])), $(GPPR²[4])($(GPPR²[8]))")
        println(io,"MAPE: $(GPPMAPE[1])($(GPPMAPE[5])), $(GPPMAPE[2])($(GPPMAPE[6])), $(GPPMAPE[3])($(GPPMAPE[7])), $(GPPMAPE[4])($(GPPMAPE[8]))")
        println(io,"RMSE: $(GPPRMSE[1])($(GPPRMSE[5])), $(GPPRMSE[2])($(GPPRMSE[6])), $(GPPRMSE[3])($(GPPRMSE[7])), $(GPPRMSE[4])($(GPPRMSE[8]))")
        println(io,"Cor: $(GPPcor[1])($(GPPcor[5])), $(GPPcor[2])($(GPPcor[6])), $(GPPcor[3])($(GPPcor[7])), $(GPPcor[4])($(GPPcor[8]))")
        println(io,"")
        println(io,"Mean MAPE: $(round(mean(GPPMAPE[1:4]),digits=2))($(round(mean(GPPMAPE[5:8]),digits=2)))")
        println(io,"Mean RMSE: $(round(mean(GPPRMSE[1:4]),digits=2))($(round(mean(GPPRMSE[5:8]),digits=2)))")
        println(io,"");println(io,"")
        println(io,"Ec:")
        println(io,"")
        println(io,"Year: 2015, 2016, 2017, 2018")
        println(io,"R²: $(EcR²[1])($(EcR²[5])), $(EcR²[2])($(EcR²[6])), $(EcR²[3])($(EcR²[7])), $(EcR²[4])($(EcR²[8]))")
        println(io,"MAPE: $(EcMAPE[1])($(EcMAPE[5])), $(EcMAPE[2])($(EcMAPE[6])), $(EcMAPE[3])($(EcMAPE[7])), $(EcMAPE[4])($(EcMAPE[8]))")
        println(io,"RMSE: $(EcRMSE[1])($(EcRMSE[5])), $(EcRMSE[2])($(EcRMSE[6])), $(EcRMSE[3])($(EcRMSE[7])), $(EcRMSE[4])($(EcRMSE[8]))")
        println(io,"Cor: $(Eccor[1])($(Eccor[5])), $(Eccor[2])($(Eccor[6])), $(Eccor[3])($(Eccor[7])), $(Eccor[4])($(Eccor[8]))")
        println(io,"")
        println(io,"Mean MAPE: $(round(mean(EcMAPE[1:4]),digits=2))($(round(mean(EcMAPE[5:8]),digits=2)))")
        println(io,"Mean RMSE: $(round(mean(EcRMSE[1:4]),digits=2))($(round(mean(EcRMSE[5:8]),digits=2)))")
    end 
end

function create_plot()
    
    modeltrain_F,modelval_F,modeltrain_C,modelval_C = shared_parameters()
    
    F_stat = create_statplot(modeltrain_F,modelval_F)
    C_stat = create_statplot(modeltrain_C,modelval_C)
    
    file_name = "RO_Opt_Shared_Parameters_GPP_Nm_f_20220511"

    savefig(F_stat,"./plots/"*file_name*"_Stat_Fertilized.svg")
    create_statplot(file_name*"_Stat_Fertilized",modeltrain_F,modelval_F)
    savefig(C_stat,"./plots/"*file_name*"_Stat_Control.svg")
    create_statplot(file_name*"_Stat_Control",modeltrain_C,modelval_C) 
        
    #=    
    modeltrain_F,modelval_F,modeltrain_C,modelval_C = no_shared_parameters()

    F_stat = create_statplot(modeltrain_F,modelval_F)
    C_stat = create_statplot(modeltrain_C,modelval_C)
    
    
    file_name = "RO_Opt_non_Shared_Parameters_20220508"

    savefig(F_stat,"./plots/"*file_name*"_Stat_Fertilized.svg")
    savefig(C_stat,"./plots/"*file_name*"_Stat_Control.svg")
    =#
end

create_plot()