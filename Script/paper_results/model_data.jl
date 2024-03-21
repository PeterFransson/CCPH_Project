mutable struct TypeStat
    R2::Real
    RMSE::Real
    MAPE::Real
    cor::Real
end
mutable struct ModelStat
    GPP::TypeStat
    Ec::TypeStat
    log_L::Real #Log likelihood
end

function get_sum_stat_crossval(folder_name::String,filename::String)   
    
    stand_type = JLD.load("output/"*folder_name*"/"*filename*".jld","stand_type")
    weight_GPP = JLD.load("output/"*folder_name*"/"*filename*".jld","weight_GPP")

    x_opt = JLD.load("output/"*folder_name*"/"*filename*".jld","x_opt")   
    train_set = JLD.load("output/"*folder_name*"/"*filename*".jld","train_set")
    val_set = JLD.load("output/"*folder_name*"/"*filename*".jld","val_set") 

    raw_input = RawInputData(;stand_type=stand_type)
    Ec_data = calc_Ec_data.(raw_input)
    GPP_data = get_GPP_data.(raw_input;stand_type=stand_type)      
        
    par = ModelPar(x_opt;stand_type=stand_type)

    GPP_model,Ec_model,Nₘ_f_model = run_model(par,raw_input;stand_type=stand_type)

    log_L = Calc_logP_GPP_Ec_Nm_f(GPP_model,Ec_model,Nₘ_f_model,GPP_data,Ec_data,par,raw_input;weight_GPP=weight_GPP)

    GPP_R2_train,GPP_RMSE_train,GPP_MAPE_train,GPP_cor_train = get_sum_stat([GPP_data[i][train_set[i]] for i in 1:4],[GPP_model[i][train_set[i]]*raw_input[i].ζ for i in 1:4])
    Ec_R2_train,Ec_RMSE_train,Ec_MAPE_train,Ec_cor_train = get_sum_stat([Ec_data[i][train_set[i]] for i in 1:4],[Ec_model[i][train_set[i]] for i in 1:4])

    GPP_R2_val,GPP_RMSE_val,GPP_MAPE_val,GPP_cor_val = get_sum_stat([GPP_data[i][val_set[i]] for i in 1:4],[GPP_model[i][val_set[i]]*raw_input[i].ζ for i in 1:4])
    Ec_R2_val,Ec_RMSE_val,Ec_MAPE_val,Ec_cor_val = get_sum_stat([Ec_data[i][val_set[i]] for i in 1:4],[Ec_model[i][val_set[i]] for i in 1:4])

    result_stat_train = ModelStat(TypeStat(GPP_R2_train,GPP_RMSE_train,GPP_MAPE_train,GPP_cor_train),
    TypeStat(Ec_R2_train,Ec_RMSE_train,Ec_MAPE_train,Ec_cor_train),
    log_L)

    result_stat_val = ModelStat(TypeStat(GPP_R2_val,GPP_RMSE_val,GPP_MAPE_val,GPP_cor_val),
    TypeStat(Ec_R2_val,Ec_RMSE_val,Ec_MAPE_val,Ec_cor_val),
    log_L)

    return result_stat_train,result_stat_val
end

function write_par2file(filename::String,par::ModelPar,stand_type::Symbol)
    open(filename,"w") do io
        println(io,"Stand type: $(stand_type)")    
        println(io,"Nₛ: $(round(par.Nₛ,sigdigits=2))")
        println(io,"α_max: $(round(par.α_max,digits=2))")
        println(io,"a_Jmax: $(round(par.a_Jmax,digits=2))")
        println(io,"ΔS: $(round(par.ΔS,digits=2))")
        println(io,"τ: $(round(par.τ,digits=2))")
        println(io,"Kₓₗ₀: $(round(par.Kₓₗ₀,sigdigits=2))")
        println(io,"a_GPP: $(round(par.a_GPP,digits=2))")
        println(io,"b_GPP: $(round(par.b_GPP,digits=2))")
    end
end
function write_stat2file(filename::String,
    GPP_data::Vector{Vector{R}},
    Ec_data::Vector{Vector{T}},
    GPP_model::Vector{Vector{W}},
    Ec_model::Vector{Vector{S}},
    train_set::Vector{Vector{Bool}},
    val_set::Vector{Vector{Bool}},
    stand_type::Symbol) where {R<:Real,T<:Real,W<:Real,S<:Real}

    GPP_data_train = [GPP_data[i][train_set[i]] for i in 1:4] 
    Ec_data_train = [Ec_data[i][train_set[i]] for i in 1:4] 

    GPP_model_train = [GPP_model[i][train_set[i]] for i in 1:4]
    Ec_model_train = [Ec_model[i][train_set[i]] for i in 1:4]

    GPP_data_val = [GPP_data[i][val_set[i]] for i in 1:4] 
    Ec_data_val = [Ec_data[i][val_set[i]] for i in 1:4] 

    GPP_model_val = [GPP_model[i][val_set[i]] for i in 1:4]
    Ec_model_val = [Ec_model[i][val_set[i]] for i in 1:4]

    GPP_R2,GPP_RMSE,GPP_MAPE,GPP_cor = get_sum_stat(GPP_data,GPP_model)
    Ec_R2,Ec_RMSE,Ec_MAPE,Ec_cor = get_sum_stat(Ec_data,Ec_model)

    GPP_R2_train,GPP_RMSE_train,GPP_MAPE_train,GPP_cor_train = get_sum_stat(GPP_data_train,GPP_model_train)
    Ec_R2_train,Ec_RMSE_train,Ec_MAPE_train,Ec_cor_train = get_sum_stat(Ec_data_train,Ec_model_train)

    GPP_R2_val,GPP_RMSE_val,GPP_MAPE_val,GPP_cor_val = get_sum_stat(GPP_data_val,GPP_model_val)
    Ec_R2_val,Ec_RMSE_val,Ec_MAPE_val,Ec_cor_val = get_sum_stat(Ec_data_val,Ec_model_val)

    open(filename,"w") do io
        println(io,"Stand type: $(stand_type)") 
        println(io,"--Total--")
        println(io,"GPP: R²:$(round(GPP_R2, digits=2)), RMSE:$(round(GPP_RMSE, digits=2)), MAPE:$(round(GPP_MAPE, digits=2)), corr:$(round(GPP_cor, digits=2))")
        println(io,"Ec: R²:$(round(Ec_R2, digits=2)), RMSE:$(round(Ec_RMSE, digits=2)), MAPE:$(round(Ec_MAPE, digits=2)), corr:$(round(Ec_cor, digits=2))")
        println(io,"--Train--")
        println(io,"GPP: R²:$(round(GPP_R2_train, digits=2)), RMSE:$(round(GPP_RMSE_train, digits=2)), MAPE:$(round(GPP_MAPE_train, digits=2)), corr:$(round(GPP_cor_train, digits=2))")
        println(io,"Ec: R²:$(round(Ec_R2_train, digits=2)), RMSE:$(round(Ec_RMSE_train, digits=2)), MAPE:$(round(Ec_MAPE_train, digits=2)), corr:$(round(Ec_cor_train, digits=2))")
        println(io,"--Validation--")
        println(io,"GPP: R²:$(round(GPP_R2_val, digits=2)), RMSE:$(round(GPP_RMSE_val, digits=2)), MAPE:$(round(GPP_MAPE_val, digits=2)), corr:$(round(GPP_cor_val, digits=2))")
        println(io,"Ec: R²:$(round(Ec_R2_val, digits=2)), RMSE:$(round(Ec_RMSE_val, digits=2)), MAPE:$(round(Ec_MAPE_val, digits=2)), corr:$(round(Ec_cor_val, digits=2))")
    end
end

function draw_shared_model()
    fld = "crossval_20240306_shared_W_1_5_run_4"

    stand_type_F = JLD.load("output/"*fld*"/result_F.jld","stand_type")
    raw_input_F = RawInputData(;stand_type=stand_type_F)
    Ec_data_F = calc_Ec_data.(raw_input_F)
    GPP_data_F = get_GPP_data.(raw_input_F;stand_type=stand_type_F)   
    x_opt_F = JLD.load("output/"*fld*"/result_F.jld","x_opt")   
    train_set_F = JLD.load("output/"*fld*"/result_F.jld","train_set")
    val_set_F = JLD.load("output/"*fld*"/result_F.jld","val_set")       
    par_F = ModelPar(x_opt_F;stand_type=stand_type_F)
    GPP_model_F,Ec_model_F,Nₘ_f_model_F = run_model(par_F,raw_input_F;stand_type=stand_type_F)
    write_par2file("./output/paper_results/shared_par/par_F.txt",par_F,stand_type_F)
    write_stat2file("./output/paper_results/shared_par/stat_F.txt",
    GPP_data_F,
    Ec_data_F,
    [GPP_model_F[i]*raw_input_F[i].ζ for i in 1:4],
    Ec_model_F,
    train_set_F,
    val_set_F,
    stand_type_F)
    
    stand_type_C = JLD.load("output/"*fld*"/result_C.jld","stand_type")
    raw_input_C = RawInputData(;stand_type=stand_type_C)
    Ec_data_C = calc_Ec_data.(raw_input_C)
    GPP_data_C = get_GPP_data.(raw_input_C;stand_type=stand_type_C)   
    x_opt_C = JLD.load("output/"*fld*"/result_C.jld","x_opt") 
    train_set_C = JLD.load("output/"*fld*"/result_C.jld","train_set")
    val_set_C = JLD.load("output/"*fld*"/result_C.jld","val_set")          
    par_C = ModelPar(x_opt_C;stand_type=stand_type_C)
    GPP_model_C,Ec_model_C,Nₘ_f_model_C = run_model(par_C,raw_input_C;stand_type=stand_type_C)
    write_par2file("./output/paper_results/shared_par/par_C.txt",par_C,stand_type_C)
    write_stat2file("./output/paper_results/shared_par/stat_C.txt",
    GPP_data_C,
    Ec_data_C,
    [GPP_model_C[i]*raw_input_C[i].ζ for i in 1:4],
    Ec_model_C,
    train_set_C,
    val_set_C,
    stand_type_C)
       
    c = -log(0.025*2) #Use for calculating 95% credible intervals

    pl1 = [plot(xlabel="2015",ylabel="GPP (g C day⁻¹ m⁻²)",legends=false, ylims = (0,13),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (0,13),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (0,13),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (0,13),yaxis=false,guidefontsize=12)]

    pl2 = [plot(xlabel="2015",ylabel="Ec (mm day⁻¹)",legends=false, ylims = (0,3),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (0,3),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (0,3),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (0,3),yaxis=false,guidefontsize=12)]

    pl3 = [plot(xlabel="2015",ylabel="GPP (g C day⁻¹ m⁻²)",legends=false, ylims = (0,13),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (0,13),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (0,13),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (0,13),yaxis=false,guidefontsize=12)]

    pl4 = [plot(xlabel="2015",ylabel="Ec (mm day⁻¹)",legends=false, ylims = (0,3),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (0,3),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (0,3),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (0,3),yaxis=false,guidefontsize=12)]

    pl5 = plot(xlabel="Data GPP",ylabel="Res GPP",legends=false,guidefontsize=12,ytickfontsize=12)

    pl6 = plot(xlabel="Data Ec",ylabel="Res Ec",legends=false,guidefontsize=12,ytickfontsize=12)

    pl7 = plot(xlabel="Data GPP",ylabel="Res GPP",legends=false,guidefontsize=12,ytickfontsize=12)

    pl8 = plot(xlabel="Data Ec",ylabel="Res Ec",legends=false,guidefontsize=12,ytickfontsize=12)

    for i = 1:4
        date = [weather.date for weather in raw_input_F[i].weather_growth]

        start_tick = ""
        end_tick = ""

        GPP_model_F_temp = GPP_model_F[i]*raw_input_F[i].ζ
        Ec_model_F_temp = Ec_model_F[i]

        σ_GPP_F = par_F.a_GPP.+GPP_model_F_temp*par_F.b_GPP
        σ_EC_F = par_F.a_Ec.+Ec_model_F_temp*par_F.b_Ec        

        xs, ys = [date; reverse(date)], [GPP_model_F_temp.+c*σ_GPP_F; reverse(GPP_model_F_temp.-c*σ_GPP_F)]
        plot!(pl1[i],xs, ys, st=:shape, fc=:blues, lw=0)

        plot!(pl1[i],date,GPP_model_F_temp,linecolor=:blue)
        plot!(pl1[i],date,GPP_data_F[i],linecolor=:green)
        plot!(pl1[i],xticks=([date[1],date[end]],[start_tick,end_tick]))

        xs, ys = [date; reverse(date)], [Ec_model_F_temp.+c*σ_EC_F; reverse(Ec_model_F_temp.-c*σ_EC_F)]
        plot!(pl2[i],xs, ys, st=:shape, fc=:blues, lw=0)

        plot!(pl2[i],date,Ec_model_F_temp,linecolor=:blue)
        plot!(pl2[i],date,Ec_data_F[i],linecolor=:green)
        plot!(pl2[i],xticks=([date[1],date[end]],[start_tick,end_tick]))

        GPP_model_C_temp = GPP_model_C[i]*raw_input_C[i].ζ
        Ec_model_C_temp = Ec_model_C[i]

        σ_GPP_C = par_C.a_GPP.+GPP_model_C_temp*par_C.b_GPP
        σ_EC_C = par_C.a_Ec.+Ec_model_C_temp*par_C.b_Ec        

        xs, ys = [date; reverse(date)], [GPP_model_C_temp.+c*σ_GPP_C; reverse(GPP_model_C_temp.-c*σ_GPP_C)]
        plot!(pl3[i],xs, ys, st=:shape, fc=:blues, lw=0)

        plot!(pl3[i],date,GPP_model_C_temp,linecolor=:blue)
        plot!(pl3[i],date,GPP_data_C[i],linecolor=:green)
        plot!(pl3[i],xticks=([date[1],date[end]],[start_tick,end_tick]))

        xs, ys = [date; reverse(date)], [Ec_model_C_temp.+c*σ_EC_C; reverse(Ec_model_C_temp.-c*σ_EC_C)]
        plot!(pl4[i],xs, ys, st=:shape, fc=:blues, lw=0)

        plot!(pl4[i],date,Ec_model_C_temp,linecolor=:blue)
        plot!(pl4[i],date,Ec_data_C[i],linecolor=:green)
        plot!(pl4[i],xticks=([date[1],date[end]],[start_tick,end_tick]))

        GPP_data_F_train = GPP_data_F[i][train_set_F[i]]
        Ec_data_F_train = Ec_data_F[i][train_set_F[i]]

        GPP_model_F_train = GPP_model_F_temp[train_set_F[i]]
        Ec_model_F_train = Ec_model_F_temp[train_set_F[i]]

        GPP_data_F_val = GPP_data_F[i][val_set_F[i]]
        Ec_data_F_val = Ec_data_F[i][val_set_F[i]]

        GPP_model_F_val = GPP_model_F_temp[val_set_F[i]]
        Ec_model_F_val = Ec_model_F_temp[val_set_F[i]]

        plot!(pl5,GPP_data_F_train,GPP_data_F_train-GPP_model_F_train,seriestype=:scatter,markercolor=:red)
        plot!(pl5,GPP_data_F_val,GPP_data_F_val-GPP_model_F_val,seriestype=:scatter,markercolor=:green)
        plot!(pl5,[minimum(GPP_data_F_train),maximum(GPP_data_F_train)],[0,0])

        plot!(pl6,Ec_data_F_train,Ec_data_F_train-Ec_model_F_train,seriestype=:scatter,markercolor=:red)        
        plot!(pl6,Ec_data_F_val,Ec_data_F_val-Ec_model_F_val,seriestype=:scatter,markercolor=:green)
        plot!(pl6,[minimum(Ec_data_F_train),maximum(Ec_data_F_train)],[0,0])

        GPP_data_C_train = GPP_data_C[i][train_set_C[i]]
        Ec_data_C_train = Ec_data_C[i][train_set_C[i]]

        GPP_model_C_train = GPP_model_C_temp[train_set_C[i]]
        Ec_model_C_train = Ec_model_C_temp[train_set_C[i]]

        GPP_data_C_val = GPP_data_C[i][val_set_C[i]]
        Ec_data_C_val = Ec_data_C[i][val_set_C[i]]

        GPP_model_C_val = GPP_model_C_temp[val_set_C[i]]
        Ec_model_C_val = Ec_model_C_temp[val_set_C[i]]

        plot!(pl7,GPP_data_C_train,GPP_data_C_train-GPP_model_C_train,seriestype=:scatter,markercolor=:red)        
        plot!(pl7,GPP_data_C_val,GPP_data_C_val-GPP_model_C_val,seriestype=:scatter,markercolor=:green)
        plot!(pl7,[minimum(GPP_data_C_train),maximum(GPP_data_C_train)],[0,0])

        plot!(pl8,Ec_data_C_train,Ec_data_C_train-Ec_model_C_train,seriestype=:scatter,markercolor=:red)        
        plot!(pl8,Ec_data_C_val,Ec_data_C_val-Ec_model_C_val,seriestype=:scatter,markercolor=:green)
        plot!(pl8,[minimum(Ec_data_C_train),maximum(Ec_data_C_train)],[0,0])
    end

    save_fld = "./plots/paper_results/shared_par"
    isdir(save_fld)|| mkdir(save_fld)

    #Data vs Model output
    pl1_fin = plot(pl1[1],pl1[2],pl1[3],pl1[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)    
    savefig(pl1_fin,save_fld*"/GPP_F.svg")

    pl2_fin = plot(pl2[1],pl2[2],pl2[3],pl2[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl2_fin,save_fld*"/Ec_F.svg")

    pl3_fin = plot(pl3[1],pl3[2],pl3[3],pl3[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl3_fin,save_fld*"/GPP_C.svg")

    pl4_fin = plot(pl4[1],pl4[2],pl4[3],pl4[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl4_fin,save_fld*"/Ec_C.svg")
    
    savefig(pl5,save_fld*"/Error_GPP_F.svg")
    savefig(pl6,save_fld*"/Error_Ec_F.svg")

    savefig(pl7,save_fld*"/Error_GPP_C.svg")
    savefig(pl8,save_fld*"/Error_Ec_C.svg")
    
    sx = repeat(["Training", "Validation"], inner = 10)    
    nam = repeat(string.(1:10), outer = 2)

    plot(pl5,pl7,pl6,pl8,layout=(2,2),size=(900,900))
    savefig(save_fld*"/Error.svg")

    fld = "crossval_20240306_shared_W_1_5_run"
    
    RMSE_GPP_F = zeros(20)
    MAPE_GPP_F = zeros(20)
    cor_GPP_F = zeros(20)

    RMSE_Ec_F = zeros(20)
    MAPE_Ec_F = zeros(20)
    cor_Ec_F = zeros(20)

    RMSE_GPP_C = zeros(20)
    MAPE_GPP_C = zeros(20)
    cor_GPP_C = zeros(20)

    RMSE_Ec_C = zeros(20)
    MAPE_Ec_C = zeros(20)
    cor_Ec_C = zeros(20)

    for i = 1:10
        result_stat_F_train,result_stat_F_val = get_sum_stat_crossval(fld*"_$(i)","result_F")

        RMSE_GPP_F[i] = result_stat_F_train.GPP.RMSE
        RMSE_GPP_F[i+10] = result_stat_F_val.GPP.RMSE
        MAPE_GPP_F[i] = result_stat_F_train.GPP.MAPE
        MAPE_GPP_F[i+10] = result_stat_F_val.GPP.MAPE
        cor_GPP_F[i] = result_stat_F_train.GPP.cor
        cor_GPP_F[i+10] = result_stat_F_val.GPP.cor

        RMSE_Ec_F[i] = result_stat_F_train.Ec.RMSE
        RMSE_Ec_F[i+10] = result_stat_F_val.Ec.RMSE
        MAPE_Ec_F[i] = result_stat_F_train.Ec.MAPE
        MAPE_Ec_F[i+10] = result_stat_F_val.Ec.MAPE
        cor_Ec_F[i] = result_stat_F_train.Ec.cor
        cor_Ec_F[i+10] = result_stat_F_val.Ec.cor

        result_stat_C_train,result_stat_C_val = get_sum_stat_crossval(fld*"_$(i)","result_C")

        RMSE_GPP_C[i] = result_stat_C_train.GPP.RMSE
        RMSE_GPP_C[i+10] = result_stat_C_val.GPP.RMSE
        MAPE_GPP_C[i] = result_stat_C_train.GPP.MAPE
        MAPE_GPP_C[i+10] = result_stat_C_val.GPP.MAPE
        cor_GPP_C[i] = result_stat_C_train.GPP.cor
        cor_GPP_C[i+10] = result_stat_C_val.GPP.cor

        RMSE_Ec_C[i] = result_stat_C_train.Ec.RMSE
        RMSE_Ec_C[i+10] = result_stat_C_val.Ec.RMSE
        MAPE_Ec_C[i] = result_stat_C_train.Ec.MAPE
        MAPE_Ec_C[i+10] = result_stat_C_val.Ec.MAPE
        cor_Ec_C[i] = result_stat_C_train.Ec.cor
        cor_Ec_C[i+10] = result_stat_C_val.Ec.cor
    end
    
    pl1_1 = groupedbar(nam, RMSE_GPP_F, group = sx, ylabel = "RMSE",legends=false,title="GPP")
    pl1_2 = groupedbar(nam, MAPE_GPP_F, group = sx, ylabel = "MAPE",legends=false)    
    pl1_3 = groupedbar(nam, cor_GPP_F, group = sx, ylabel = "COR",legends=false)
    pl1_fin = plot(pl1_1,pl1_2,pl1_3,layout=(3,1),left_margin = 3Plots.mm,guidefontsize=12,ytickfontsize=12,xtickfontsize=12)
    savefig(save_fld*"/Stat_GPP_F.svg")

    pl2_1 = groupedbar(nam, RMSE_Ec_F, group = sx, ylabel = "RMSE",legends=false,title="Ec")
    pl2_2 = groupedbar(nam, MAPE_Ec_F, group = sx, ylabel = "MAPE",legends=false)    
    pl2_3 = groupedbar(nam, cor_Ec_F, group = sx, ylabel = "COR",legends=false)
    pl2_fin = plot(pl2_1,pl2_2,pl2_3,layout=(3,1),left_margin = 3Plots.mm,guidefontsize=12,ytickfontsize=12,xtickfontsize=12)
    savefig(save_fld*"/Stat_Ec_F.svg")

    pl3_1 = groupedbar(nam, RMSE_GPP_C, group = sx, ylabel = "RMSE",legends=false,title="GPP")
    pl3_2 = groupedbar(nam, MAPE_GPP_C, group = sx, ylabel = "MAPE",legends=false)    
    pl3_3 = groupedbar(nam, cor_GPP_C, group = sx, ylabel = "COR",legends=false)
    pl3_fin = plot(pl3_1,pl3_2,pl3_3,layout=(3,1),left_margin = 3Plots.mm,guidefontsize=12,ytickfontsize=12,xtickfontsize=12)
    savefig(save_fld*"/Stat_GPP_C.svg")

    pl4_1 = groupedbar(nam, RMSE_Ec_C, group = sx, ylabel = "RMSE",legends=false,title="Ec")
    pl4_2 = groupedbar(nam, MAPE_Ec_C, group = sx, ylabel = "MAPE",legends=false)    
    pl4_3 = groupedbar(nam, cor_Ec_C, group = sx, ylabel = "COR",legends=false)
    pl4_fin = plot(pl4_1,pl4_2,pl4_3,layout=(3,1),left_margin = 3Plots.mm,guidefontsize=12,ytickfontsize=12,xtickfontsize=12)
    savefig(save_fld*"/Stat_Ec_C.svg")

    plot(pl1_fin,pl2_fin,pl3_fin,pl4_fin,layout=(1,4),size=(1800,900))
    savefig(save_fld*"/Stat.svg")
end

function draw_nonshared_model()
    fld_F = "crossval_20240306_F_W_1_5_run_4"
    fld_C = "crossval_20240306_C_W_1_5_run_4"

    stand_type_F = JLD.load("output/"*fld_F*"/result.jld","stand_type")
    raw_input_F = RawInputData(;stand_type=stand_type_F)
    Ec_data_F = calc_Ec_data.(raw_input_F)
    GPP_data_F = get_GPP_data.(raw_input_F;stand_type=stand_type_F)   
    x_opt_F = JLD.load("output/"*fld_F*"/result.jld","x_opt")   
    train_set_F = JLD.load("output/"*fld_F*"/result.jld","train_set")
    val_set_F = JLD.load("output/"*fld_F*"/result.jld","val_set")       
    par_F = ModelPar(x_opt_F;stand_type=stand_type_F)
    GPP_model_F,Ec_model_F,Nₘ_f_model_F = run_model(par_F,raw_input_F;stand_type=stand_type_F)
    write_par2file("./output/paper_results/nonshared_par/par_F.txt",par_F,stand_type_F)
    write_stat2file("./output/paper_results/nonshared_par/stat_F.txt",
    GPP_data_F,
    Ec_data_F,
    [GPP_model_F[i]*raw_input_F[i].ζ for i in 1:4],
    Ec_model_F,
    train_set_F,
    val_set_F,
    stand_type_F)
    
    stand_type_C = JLD.load("output/"*fld_C*"/result.jld","stand_type")
    raw_input_C = RawInputData(;stand_type=stand_type_C)
    Ec_data_C = calc_Ec_data.(raw_input_C)
    GPP_data_C = get_GPP_data.(raw_input_C;stand_type=stand_type_C)   
    x_opt_C = JLD.load("output/"*fld_C*"/result.jld","x_opt") 
    train_set_C = JLD.load("output/"*fld_C*"/result.jld","train_set")
    val_set_C = JLD.load("output/"*fld_C*"/result.jld","val_set")          
    par_C = ModelPar(x_opt_C;stand_type=stand_type_C)
    GPP_model_C,Ec_model_C,Nₘ_f_model_C = run_model(par_C,raw_input_C;stand_type=stand_type_C)
    write_par2file("./output/paper_results/nonshared_par/par_C.txt",par_C,stand_type_C)
    write_stat2file("./output/paper_results/nonshared_par/stat_C.txt",
    GPP_data_C,
    Ec_data_C,
    [GPP_model_C[i]*raw_input_C[i].ζ for i in 1:4],
    Ec_model_C,
    train_set_C,
    val_set_C,
    stand_type_C)
       
    c = -log(0.025*2) #Use for calculating 95% credible intervals

    pl1 = [plot(xlabel="2015",ylabel="GPP (g C day⁻¹ m⁻²)",legends=false, ylims = (0,13),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (0,13),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (0,13),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (0,13),yaxis=false,guidefontsize=12)]

    pl2 = [plot(xlabel="2015",ylabel="Ec (mm day⁻¹)",legends=false, ylims = (0,3),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (0,3),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (0,3),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (0,3),yaxis=false,guidefontsize=12)]

    pl3 = [plot(xlabel="2015",ylabel="GPP (g C day⁻¹ m⁻²)",legends=false, ylims = (0,13),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (0,13),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (0,13),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (0,13),yaxis=false,guidefontsize=12)]

    pl4 = [plot(xlabel="2015",ylabel="Ec (mm day⁻¹)",legends=false, ylims = (0,3),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (0,3),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (0,3),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (0,3),yaxis=false,guidefontsize=12)]

    pl5 = plot(xlabel="Data GPP",ylabel="Res GPP",legends=false,guidefontsize=12,ytickfontsize=12)

    pl6 = plot(xlabel="Data Ec",ylabel="Res Ec",legends=false,guidefontsize=12,ytickfontsize=12)

    pl7 = plot(xlabel="Data GPP",ylabel="Res GPP",legends=false,guidefontsize=12,ytickfontsize=12)

    pl8 = plot(xlabel="Data Ec",ylabel="Res Ec",legends=false,guidefontsize=12,ytickfontsize=12)

    for i = 1:4
        date = [weather.date for weather in raw_input_F[i].weather_growth]

        start_tick = ""
        end_tick = ""

        GPP_model_F_temp = GPP_model_F[i]*raw_input_F[i].ζ
        Ec_model_F_temp = Ec_model_F[i]

        σ_GPP_F = par_F.a_GPP.+GPP_model_F_temp*par_F.b_GPP
        σ_EC_F = par_F.a_Ec.+Ec_model_F_temp*par_F.b_Ec        

        xs, ys = [date; reverse(date)], [GPP_model_F_temp.+c*σ_GPP_F; reverse(GPP_model_F_temp.-c*σ_GPP_F)]
        plot!(pl1[i],xs, ys, st=:shape, fc=:blues, lw=0)

        plot!(pl1[i],date,GPP_model_F_temp,linecolor=:blue)
        plot!(pl1[i],date,GPP_data_F[i],linecolor=:green)
        plot!(pl1[i],xticks=([date[1],date[end]],[start_tick,end_tick]))

        xs, ys = [date; reverse(date)], [Ec_model_F_temp.+c*σ_EC_F; reverse(Ec_model_F_temp.-c*σ_EC_F)]
        plot!(pl2[i],xs, ys, st=:shape, fc=:blues, lw=0)

        plot!(pl2[i],date,Ec_model_F_temp,linecolor=:blue)
        plot!(pl2[i],date,Ec_data_F[i],linecolor=:green)
        plot!(pl2[i],xticks=([date[1],date[end]],[start_tick,end_tick]))

        GPP_model_C_temp = GPP_model_C[i]*raw_input_C[i].ζ
        Ec_model_C_temp = Ec_model_C[i]

        σ_GPP_C = par_C.a_GPP.+GPP_model_C_temp*par_C.b_GPP
        σ_EC_C = par_C.a_Ec.+Ec_model_C_temp*par_C.b_Ec        

        xs, ys = [date; reverse(date)], [GPP_model_C_temp.+c*σ_GPP_C; reverse(GPP_model_C_temp.-c*σ_GPP_C)]
        plot!(pl3[i],xs, ys, st=:shape, fc=:blues, lw=0)

        plot!(pl3[i],date,GPP_model_C_temp,linecolor=:blue)
        plot!(pl3[i],date,GPP_data_C[i],linecolor=:green)
        plot!(pl3[i],xticks=([date[1],date[end]],[start_tick,end_tick]))

        xs, ys = [date; reverse(date)], [Ec_model_C_temp.+c*σ_EC_C; reverse(Ec_model_C_temp.-c*σ_EC_C)]
        plot!(pl4[i],xs, ys, st=:shape, fc=:blues, lw=0)

        plot!(pl4[i],date,Ec_model_C_temp,linecolor=:blue)
        plot!(pl4[i],date,Ec_data_C[i],linecolor=:green)
        plot!(pl4[i],xticks=([date[1],date[end]],[start_tick,end_tick]))

        GPP_data_F_train = GPP_data_F[i][train_set_F[i]]
        Ec_data_F_train = Ec_data_F[i][train_set_F[i]]

        GPP_model_F_train = GPP_model_F_temp[train_set_F[i]]
        Ec_model_F_train = Ec_model_F_temp[train_set_F[i]]

        GPP_data_F_val = GPP_data_F[i][val_set_F[i]]
        Ec_data_F_val = Ec_data_F[i][val_set_F[i]]

        GPP_model_F_val = GPP_model_F_temp[val_set_F[i]]
        Ec_model_F_val = Ec_model_F_temp[val_set_F[i]]

        plot!(pl5,GPP_data_F_train,GPP_data_F_train-GPP_model_F_train,seriestype=:scatter,markercolor=:red)
        plot!(pl5,GPP_data_F_val,GPP_data_F_val-GPP_model_F_val,seriestype=:scatter,markercolor=:green)
        plot!(pl5,[minimum(GPP_data_F_train),maximum(GPP_data_F_train)],[0,0])

        plot!(pl6,Ec_data_F_train,Ec_data_F_train-Ec_model_F_train,seriestype=:scatter,markercolor=:red)        
        plot!(pl6,Ec_data_F_val,Ec_data_F_val-Ec_model_F_val,seriestype=:scatter,markercolor=:green)
        plot!(pl6,[minimum(Ec_data_F_train),maximum(Ec_data_F_train)],[0,0])

        GPP_data_C_train = GPP_data_C[i][train_set_C[i]]
        Ec_data_C_train = Ec_data_C[i][train_set_C[i]]

        GPP_model_C_train = GPP_model_C_temp[train_set_C[i]]
        Ec_model_C_train = Ec_model_C_temp[train_set_C[i]]

        GPP_data_C_val = GPP_data_C[i][val_set_C[i]]
        Ec_data_C_val = Ec_data_C[i][val_set_C[i]]

        GPP_model_C_val = GPP_model_C_temp[val_set_C[i]]
        Ec_model_C_val = Ec_model_C_temp[val_set_C[i]]

        plot!(pl7,GPP_data_C_train,GPP_data_C_train-GPP_model_C_train,seriestype=:scatter,markercolor=:red)        
        plot!(pl7,GPP_data_C_val,GPP_data_C_val-GPP_model_C_val,seriestype=:scatter,markercolor=:green)
        plot!(pl7,[minimum(GPP_data_C_train),maximum(GPP_data_C_train)],[0,0])

        plot!(pl8,Ec_data_C_train,Ec_data_C_train-Ec_model_C_train,seriestype=:scatter,markercolor=:red)        
        plot!(pl8,Ec_data_C_val,Ec_data_C_val-Ec_model_C_val,seriestype=:scatter,markercolor=:green)
        plot!(pl8,[minimum(Ec_data_C_train),maximum(Ec_data_C_train)],[0,0])
    end

    save_fld = "./plots/paper_results/nonshared_par"
    isdir(save_fld)|| mkdir(save_fld)

    #Data vs Model output
    pl1_fin = plot(pl1[1],pl1[2],pl1[3],pl1[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)    
    savefig(pl1_fin,save_fld*"/GPP_F.svg")

    pl2_fin = plot(pl2[1],pl2[2],pl2[3],pl2[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl2_fin,save_fld*"/Ec_F.svg")

    pl3_fin = plot(pl3[1],pl3[2],pl3[3],pl3[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl3_fin,save_fld*"/GPP_C.svg")

    pl4_fin = plot(pl4[1],pl4[2],pl4[3],pl4[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl4_fin,save_fld*"/Ec_C.svg")
    
    savefig(pl5,save_fld*"/Error_GPP_F.svg")
    savefig(pl6,save_fld*"/Error_Ec_F.svg")

    savefig(pl7,save_fld*"/Error_GPP_C.svg")
    savefig(pl8,save_fld*"/Error_Ec_C.svg")
    
    sx = repeat(["Training", "Validation"], inner = 10)    
    nam = repeat(string.(1:10), outer = 2)

    plot(pl5,pl7,pl6,pl8,layout=(2,2),size=(900,900))
    savefig(save_fld*"/Error.svg")

    fld_F = "crossval_20240306_F_W_1_5_run"
    fld_C = "crossval_20240306_C_W_1_5_run"
    
    RMSE_GPP_F = zeros(20)
    MAPE_GPP_F = zeros(20)
    cor_GPP_F = zeros(20)

    RMSE_Ec_F = zeros(20)
    MAPE_Ec_F = zeros(20)
    cor_Ec_F = zeros(20)

    RMSE_GPP_C = zeros(20)
    MAPE_GPP_C = zeros(20)
    cor_GPP_C = zeros(20)

    RMSE_Ec_C = zeros(20)
    MAPE_Ec_C = zeros(20)
    cor_Ec_C = zeros(20)

    for i = 1:10
        result_stat_F_train,result_stat_F_val = get_sum_stat_crossval(fld_F*"_$(i)","result")

        RMSE_GPP_F[i] = result_stat_F_train.GPP.RMSE
        RMSE_GPP_F[i+10] = result_stat_F_val.GPP.RMSE
        MAPE_GPP_F[i] = result_stat_F_train.GPP.MAPE
        MAPE_GPP_F[i+10] = result_stat_F_val.GPP.MAPE
        cor_GPP_F[i] = result_stat_F_train.GPP.cor
        cor_GPP_F[i+10] = result_stat_F_val.GPP.cor

        RMSE_Ec_F[i] = result_stat_F_train.Ec.RMSE
        RMSE_Ec_F[i+10] = result_stat_F_val.Ec.RMSE
        MAPE_Ec_F[i] = result_stat_F_train.Ec.MAPE
        MAPE_Ec_F[i+10] = result_stat_F_val.Ec.MAPE
        cor_Ec_F[i] = result_stat_F_train.Ec.cor
        cor_Ec_F[i+10] = result_stat_F_val.Ec.cor

        result_stat_C_train,result_stat_C_val = get_sum_stat_crossval(fld_C*"_$(i)","result")

        RMSE_GPP_C[i] = result_stat_C_train.GPP.RMSE
        RMSE_GPP_C[i+10] = result_stat_C_val.GPP.RMSE
        MAPE_GPP_C[i] = result_stat_C_train.GPP.MAPE
        MAPE_GPP_C[i+10] = result_stat_C_val.GPP.MAPE
        cor_GPP_C[i] = result_stat_C_train.GPP.cor
        cor_GPP_C[i+10] = result_stat_C_val.GPP.cor

        RMSE_Ec_C[i] = result_stat_C_train.Ec.RMSE
        RMSE_Ec_C[i+10] = result_stat_C_val.Ec.RMSE
        MAPE_Ec_C[i] = result_stat_C_train.Ec.MAPE
        MAPE_Ec_C[i+10] = result_stat_C_val.Ec.MAPE
        cor_Ec_C[i] = result_stat_C_train.Ec.cor
        cor_Ec_C[i+10] = result_stat_C_val.Ec.cor
    end
    
    pl1_1 = groupedbar(nam, RMSE_GPP_F, group = sx, ylabel = "RMSE",legends=false,title="GPP")
    pl1_2 = groupedbar(nam, MAPE_GPP_F, group = sx, ylabel = "MAPE",legends=false)    
    pl1_3 = groupedbar(nam, cor_GPP_F, group = sx, ylabel = "COR",legends=false)
    pl1_fin = plot(pl1_1,pl1_2,pl1_3,layout=(3,1),left_margin = 3Plots.mm,guidefontsize=12,ytickfontsize=12,xtickfontsize=12)
    savefig(save_fld*"/Stat_GPP_F.svg")

    pl2_1 = groupedbar(nam, RMSE_Ec_F, group = sx, ylabel = "RMSE",legends=false,title="Ec")
    pl2_2 = groupedbar(nam, MAPE_Ec_F, group = sx, ylabel = "MAPE",legends=false)    
    pl2_3 = groupedbar(nam, cor_Ec_F, group = sx, ylabel = "COR",legends=false)
    pl2_fin = plot(pl2_1,pl2_2,pl2_3,layout=(3,1),left_margin = 3Plots.mm,guidefontsize=12,ytickfontsize=12,xtickfontsize=12)
    savefig(save_fld*"/Stat_Ec_F.svg")

    pl3_1 = groupedbar(nam, RMSE_GPP_C, group = sx, ylabel = "RMSE",legends=false,title="GPP")
    pl3_2 = groupedbar(nam, MAPE_GPP_C, group = sx, ylabel = "MAPE",legends=false)    
    pl3_3 = groupedbar(nam, cor_GPP_C, group = sx, ylabel = "COR",legends=false)
    pl3_fin = plot(pl3_1,pl3_2,pl3_3,layout=(3,1),left_margin = 3Plots.mm,guidefontsize=12,ytickfontsize=12,xtickfontsize=12)
    savefig(save_fld*"/Stat_GPP_C.svg")

    pl4_1 = groupedbar(nam, RMSE_Ec_C, group = sx, ylabel = "RMSE",legends=false,title="Ec")
    pl4_2 = groupedbar(nam, MAPE_Ec_C, group = sx, ylabel = "MAPE",legends=false)    
    pl4_3 = groupedbar(nam, cor_Ec_C, group = sx, ylabel = "COR",legends=false)
    pl4_fin = plot(pl4_1,pl4_2,pl4_3,layout=(3,1),left_margin = 3Plots.mm,guidefontsize=12,ytickfontsize=12,xtickfontsize=12)
    savefig(save_fld*"/Stat_Ec_C.svg")

    plot(pl1_fin,pl2_fin,pl3_fin,pl4_fin,layout=(1,4),size=(1800,900))
    savefig(save_fld*"/Stat.svg")
end