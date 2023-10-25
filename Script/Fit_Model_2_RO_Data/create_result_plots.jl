function create_result_plots(file_name::String,
    parainfo::Union{Dict{Symbol,Any},Array{Symbol,1}};
    ParaDictInit_F=nothing,
    ParaDictInit_C=nothing)

    RO_data = Load_RO_data()          
    
    if typeof(parainfo)==Dict{Symbol,Any}
        ind = load("./output/"*file_name*".jld","ind_train") 
        par = load("./output/"*file_name*".jld","xopt") 
        para2ind = parainfo
        ParaDict_F,ParaDict_C = CreateParaDict(para2ind,par;
        ParaDictInit_F=ParaDictInit_F,ParaDictInit_C=ParaDictInit_C)
    elseif typeof(parainfo)==Array{Symbol,1}  
        file_name_F = file_name*"_F"
        file_name_C = file_name*"_C"
        ind = load("./output/"*file_name_F*".jld","ind_train")
        par_F = load("./output/"*file_name_F*".jld","xopt")
        par_C = load("./output/"*file_name_C*".jld","xopt")
        parasym = parainfo 
        ParaDict_F = CreateParaDict(parasym,par_F;ParaDictInit=ParaDictInit_F)
        ParaDict_C = CreateParaDict(parasym,par_C;ParaDictInit=ParaDictInit_C)
    else
        error("parainfo Wrong type")
    end

    model_F,weatherts_F,data_F = Get_Result_RO_CCPH(ParaDict_F,RO_data;stand_type="Fertilized")
    model_C,weatherts_C,data_C = Get_Result_RO_CCPH(ParaDict_C,RO_data;stand_type="Control")
    
    a_GPP_F,b_GPP_F,a_Ec_F,b_Ec_F = ParaDict_F[:a_GPP],ParaDict_F[:b_GPP],ParaDict_F[:a_Ec],ParaDict_F[:b_Ec]
    a_GPP_C,b_GPP_C,a_Ec_C,b_Ec_C = ParaDict_C[:a_GPP],ParaDict_C[:b_GPP],ParaDict_C[:a_Ec],ParaDict_C[:b_Ec]

    σ_GPP_F = a_GPP_F.+model_F.GPP*b_GPP_F
    σ_GPP_C = a_GPP_C.+model_C.GPP*b_GPP_C

    σ_EC_F = a_Ec_F.+b_Ec_F*model_F.Ec
    σ_EC_C = a_Ec_C.+b_Ec_C*model_C.Ec

    c = -log(0.025*2) #Use for calculating 95% credible intervals
        
    plot(weatherts_F.date[ind],data_F.GPP[ind],label="Data",ylabel="GPP",ylims=(0.0,12.0))
    plot!(weatherts_F.date[ind],model_F.GPP[ind]+c*σ_GPP_F[ind])
    plot!(weatherts_F.date[ind],model_F.GPP[ind]-c*σ_GPP_F[ind])
    pl1 = plot!(weatherts_F.date[ind],model_F.GPP[ind],label="Model")
    plot(weatherts_C.date[ind],data_C.GPP[ind],label="Data",ylabel="GPP",ylims=(0.0,12.0))
    plot!(weatherts_C.date[ind],model_C.GPP[ind]+c*σ_GPP_C[ind])
    plot!(weatherts_C.date[ind],model_C.GPP[ind]-c*σ_GPP_C[ind])
    pl2 = plot!(weatherts_C.date[ind],model_C.GPP[ind],label="Model")

    plot(weatherts_F.date[ind],data_F.Ec[ind],label="Data",ylabel="E_C",ylims=(0.0,3.0))
    plot!(weatherts_F.date[ind],model_F.Ec[ind]+c*σ_EC_F[ind])
    plot!(weatherts_F.date[ind],model_F.Ec[ind]-c*σ_EC_F[ind])
    pl3 = plot!(weatherts_F.date[ind],model_F.Ec[ind],label="Model")
    plot(weatherts_C.date[ind],data_C.Ec[ind],label="Data",ylabel="E_C",ylims=(0.0,3.0))
    plot!(weatherts_C.date[ind],model_C.Ec[ind]+c*σ_EC_C[ind])
    plot!(weatherts_C.date[ind],model_C.Ec[ind]-c*σ_EC_C[ind])
    pl4 = plot!(weatherts_C.date[ind],model_C.Ec[ind],label="Model")

    plot(pl1,pl2,pl3,pl4,layout=(2,2),legends=false)
    savefig("./plots/"*file_name*"_result.svg")  
    
    savefig(pl1,"./plots/"*file_name*"_result_GPP_F.svg")
    savefig(pl2,"./plots/"*file_name*"_result_GPP_C.svg")
    savefig(pl3,"./plots/"*file_name*"_result_E_C_F.svg")
    savefig(pl4,"./plots/"*file_name*"_result_E_C_C.svg")

    plot(weatherts_F.date[ind],data_F.gₛ[ind],label="Data",ylabel="gₛ")
    pl1 = plot!(weatherts_F.date[ind],model_F.gₛ[ind],label="Model")
    plot(weatherts_C.date[ind],data_C.gₛ[ind],label="Data",ylabel="gₛ")
    pl2 = plot!(weatherts_C.date[ind],model_C.gₛ[ind],label="Model")

    plot(pl1,pl2,layout=(1,2),legends=false)
    savefig("./plots/"*file_name*"_result_gs.svg")

    pl1 = plot(weatherts_F.date[ind],model_F.Nₘ_f[ind]*100,label="Model",ylabel="Nₘ_f (%)")    
    pl2 = plot(weatherts_C.date[ind],model_C.Nₘ_f[ind]*100,label="Model",ylabel="Nₘ_f (%)")
    
    plot(pl1,pl2,layout=(1,2),legends=false)
    savefig("./plots/"*file_name*"_result_Nm_f.svg")
    
    ind_val = setdiff(1:84, ind)
    plot([0.0,15.0],[0.0,15.0],xlabel="Data GPP",ylabel="Model GPP")
    plot!(data_F.GPP[ind],model_F.GPP[ind],seriestype=:scatter)
    pl1 = plot!(data_F.GPP[ind_val],model_F.GPP[ind_val],seriestype=:scatter)
    plot([0.0,15.0],[0.0,15.0],xlabel="Data GPP",ylabel="Model GPP")
    plot!(data_C.GPP[ind],model_C.GPP[ind],seriestype=:scatter)
    pl2 = plot!(data_C.GPP[ind_val],model_C.GPP[ind_val],seriestype=:scatter)

    plot([0.0,2.0],[0.0,2.0],xlabel="Data Ec",ylabel="Model Ec")
    plot!(data_F.Ec[ind],model_F.Ec[ind],seriestype=:scatter)
    pl3 = plot!(data_F.Ec[ind_val],model_F.Ec[ind_val],seriestype=:scatter)
    plot([0.0,2.0],[0.0,2.0],xlabel="Data Ec",ylabel="Model Ec")
    plot!(data_C.Ec[ind],model_C.Ec[ind],seriestype=:scatter)
    pl4 = plot!(data_C.Ec[ind_val],model_C.Ec[ind_val],seriestype=:scatter)

    plot(pl1,pl2,pl3,pl4,layout=(2,2),legends=false)    
    savefig("./plots/"*file_name*"_annual_error.svg")

    plot([0.0,15.0],[0.0,0.0],xlabel="Data GPP",ylabel="Res GPP")
    plot!(data_F.GPP[ind],model_F.GPP[ind]-data_F.GPP[ind],seriestype=:scatter)
    pl1 = plot!(data_F.GPP[ind_val],model_F.GPP[ind_val]-data_F.GPP[ind_val],seriestype=:scatter)
    plot([0.0,15.0],[0.0,0.0],xlabel="Data GPP",ylabel="Res GPP")
    plot!(data_C.GPP[ind],model_C.GPP[ind]-data_C.GPP[ind],seriestype=:scatter)
    pl2 = plot!(data_C.GPP[ind_val],model_C.GPP[ind_val]-data_C.GPP[ind_val],seriestype=:scatter)

    plot([0.0,2.0],[0.0,0.0],xlabel="Data Ec",ylabel="Res Ec")
    plot!(data_F.Ec[ind],model_F.Ec[ind]-data_F.Ec[ind],seriestype=:scatter)
    pl3 = plot!(data_F.Ec[ind_val],model_F.Ec[ind_val]-data_F.Ec[ind_val],seriestype=:scatter)
    plot([0.0,2.0],[0.0,0.0],xlabel="Data Ec",ylabel="Res Ec")
    plot!(data_C.Ec[ind],model_C.Ec[ind]-data_C.Ec[ind],seriestype=:scatter)
    pl4 = plot!(data_C.Ec[ind_val],model_C.Ec[ind_val]-data_C.Ec[ind_val],seriestype=:scatter)

    plot(pl1,pl2,pl3,pl4,layout=(2,2),legends=false)    
    savefig("./plots/"*file_name*"_residuals.svg")

    create_wue_plot("./plots/"*file_name*"_wue_F",model_F,weatherts_F)  
    create_wue_plot("./plots/"*file_name*"_wue_C",model_C,weatherts_C) 
    create_wue_stat("./output/"*file_name*"_wue_F",model_F,weatherts_F)
    create_wue_stat("./output/"*file_name*"_wue_C",model_C,weatherts_C)

    create_trait_plots("./plots/"*file_name*"_trait_F",model_F,weatherts_F)
    create_trait_plots("./plots/"*file_name*"_trait_C",model_C,weatherts_C)
    create_trait_stat("./output/"*file_name*"_trait_F",model_F,weatherts_F)
    create_trait_stat("./output/"*file_name*"_trait_C",model_C,weatherts_C)

    create_A_plot("./plots/"*file_name*"_A_F",model_F,weatherts_F)
    create_A_plot("./plots/"*file_name*"_A_C",model_C,weatherts_C)

    create_cᵢ_plot("./plots/"*file_name*"_cᵢ_F",model_F,weatherts_F)
    create_cᵢ_plot("./plots/"*file_name*"_cᵢ_C",model_C,weatherts_C)

    create_trait_plots("./plots/"*file_name*"_trait_F_C",
    model_F,
    weatherts_F,    
    model_C,
    weatherts_C)

    create_wue_plot("./plots/"*file_name*"_wue_F_C",
    model_F,
    weatherts_F,    
    model_C,
    weatherts_C)
end
function run_create_result_plots(file_name::String,
    ranges::Array{Tuple{Float64, Float64}, 1},
    parainfo::Union{Dict{Symbol,Any},Array{Symbol,1}},
    n_runs::Integer;
    ParaDictInit_F=nothing,
    ParaDictInit_C=nothing)

    isdir("./output/"*file_name)|| mkdir("./output/"*file_name)
    isdir("./plots/"*file_name) || mkdir("./plots/"*file_name)

    modeltrain_F = Array{ModelStatSum,1}(undef,n_runs)
    modelval_F = Array{ModelStatSum,1}(undef,n_runs)
    modeltrain_C = Array{ModelStatSum,1}(undef,n_runs)
    modelval_C = Array{ModelStatSum,1}(undef,n_runs)
    
    for ix = 1:n_runs
        file_name_save = file_name*"/"*file_name*"_$(ix)"

        run_opt_par(file_name_save,
        ranges,
        parainfo;
        ParaDictInit_F=ParaDictInit_F,
        ParaDictInit_C=ParaDictInit_C)        

        run_validation_RO(file_name_save,
        ranges,
        parainfo;
        ParaDictInit_F=ParaDictInit_F,
        ParaDictInit_C=ParaDictInit_C)

        modeltrain_F[ix],modelval_F[ix],modeltrain_C[ix],modelval_C[ix] = Get_Model_Stat(file_name_save,
        ranges,
        parainfo;
        ParaDictInit_F=ParaDictInit_F,
        ParaDictInit_C=ParaDictInit_C)

        create_result_plots(file_name_save,
        parainfo;
        ParaDictInit_F=ParaDictInit_F,
        ParaDictInit_C=ParaDictInit_C)
    end

    writestat2file(file_name*"/"*file_name,
    modeltrain_F,
    modelval_F,
    modeltrain_C,
    modelval_C)    
    
    CreateStatPlot(file_name*"/"*file_name,
    n_runs,
    modeltrain_F,
    modelval_F,
    modeltrain_C,
    modelval_C)
end

