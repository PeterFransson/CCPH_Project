function draw_shared_model()
    fld = "crossval_20240306_shared_W_1_5_run_4"

    stand_type_F = JLD.load("output/"*fld*"/result_F.jld","stand_type")
    raw_input_F = RawInputData(;stand_type=stand_type_F)
    Ec_data_F = calc_Ec_data.(raw_input_F)
    GPP_data_F = get_GPP_data.(raw_input_F;stand_type=stand_type_F)   
    x_opt_F = JLD.load("output/"*fld*"/result_F.jld","x_opt")          
    par_F = ModelPar(x_opt_F;stand_type=stand_type_F)
    GPP_model_F,Ec_model_F,Nₘ_f_model_F = run_model(par_F,raw_input_F;stand_type=stand_type_F)

    stand_type_C = JLD.load("output/"*fld*"/result_C.jld","stand_type")
    raw_input_C = RawInputData(;stand_type=stand_type_C)
    Ec_data_C = calc_Ec_data.(raw_input_C)
    GPP_data_C = get_GPP_data.(raw_input_C;stand_type=stand_type_C)   
    x_opt_C = JLD.load("output/"*fld*"/result_C.jld","x_opt")          
    par_C = ModelPar(x_opt_C;stand_type=stand_type_C)
    GPP_model_C,Ec_model_C,Nₘ_f_model_C = run_model(par_C,raw_input_C;stand_type=stand_type_C)
       
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