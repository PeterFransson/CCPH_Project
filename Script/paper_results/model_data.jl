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
    end

    save_fld = "./plots/paper_results/shared_par"
    isdir(save_fld)|| mkdir(save_fld)

    pl1_fin = plot(pl1[1],pl1[2],pl1[3],pl1[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl1_fin,save_fld*"/GPP_F.svg")

    pl2_fin = plot(pl2[1],pl2[2],pl2[3],pl2[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl2_fin,save_fld*"/Ec_F.svg")

    pl3_fin = plot(pl3[1],pl3[2],pl3[3],pl3[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl3_fin,save_fld*"/GPP_C.svg")

    pl4_fin = plot(pl4[1],pl4[2],pl4[3],pl4[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm)
    savefig(pl4_fin,save_fld*"/Ec_C.svg")
end