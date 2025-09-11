function test_gs_vs_I()
    fld = "crossval_20241105_shared_W_1_5_run_9"

    #Default
    stand_type_F = JLD.load("output/"*fld*"/result_F.jld","stand_type")
    raw_input_F = RawInputData(;stand_type=stand_type_F)     
    x_opt_F = JLD.load("output/"*fld*"/result_F.jld","x_opt")         
    par_F = ModelPar(x_opt_F;stand_type=stand_type_F)
    gₛ_F,Nₘ_f_F,wue_F = get_plant_var(par_F,raw_input_F;stand_type=stand_type_F)

    stand_type_C = JLD.load("output/"*fld*"/result_C.jld","stand_type")
    raw_input_C = RawInputData(;stand_type=stand_type_C)      
    x_opt_C = JLD.load("output/"*fld*"/result_C.jld","x_opt")         
    par_C = ModelPar(x_opt_C;stand_type=stand_type_C)
    gₛ_C,Nₘ_f_C,wue_C = get_plant_var(par_C,raw_input_C;stand_type=stand_type_C)

    pl1 = [plot(xlabel="2015",ylabel="Δgₛ (mol s⁻¹ m⁻²)",legends=false, ylims = (-0.001,0.005),guidefontsize=12,ytickfontsize=12), 
    plot(xlabel="2016",ylabel="",legends=false, ylims = (-0.001,0.005),yaxis=false,guidefontsize=12),
    plot(xlabel="2017",ylabel="",legends=false, ylims = (-0.001,0.005),yaxis=false,guidefontsize=12),
    plot(xlabel="2018",ylabel="",legends=false, ylims = (-0.001,0.005),yaxis=false,guidefontsize=12)]
        
    #Increased I₀
    facI = 1.05

    raw_input_F_inc_I = deepcopy(raw_input_F)
    raw_input_C_inc_I = deepcopy(raw_input_C)
    
    #Change I₀
    for i in 1:4
        for j in eachindex(raw_input_F_inc_I[i].weather_raw)
            raw_input_F_inc_I[i].weather_raw[j].Radₜₒ = facI*raw_input_F_inc_I[i].weather_raw[j].Radₜₒ
            raw_input_C_inc_I[i].weather_raw[j].Radₜₒ = facI*raw_input_C_inc_I[i].weather_raw[j].Radₜₒ
        end
        for j in eachindex(raw_input_F_inc_I[i].weather_growth)
            raw_input_F_inc_I[i].weather_raw[j].Radₜₒ = facI*raw_input_F_inc_I[i].weather_raw[j].Radₜₒ
            raw_input_C_inc_I[i].weather_raw[j].Radₜₒ = facI*raw_input_C_inc_I[i].weather_raw[j].Radₜₒ
        end
    end

    gₛ_F_inc_I,Nₘ_f_F_inc_I,wue_F_inc_I = get_plant_var(par_F,raw_input_F_inc_I;stand_type=stand_type_F)
    gₛ_C_inc_I,Nₘ_f_C_inc_I,wue_C_inc_I = get_plant_var(par_C,raw_input_C_inc_I ;stand_type=stand_type_C)
    

    for i = 1:4
        date = [weather.date for weather in raw_input_F[i].weather_growth]
        
        start_tick = ""
        end_tick = ""

        plot!(pl1[i],date,gₛ_F_inc_I[i].-gₛ_F[i],linecolor=:blue)
        plot!(pl1[i],date,gₛ_C_inc_I[i].-gₛ_C[i],linecolor=:red)
        plot!(pl1[i],xticks=([date[1],date[end]],[start_tick,end_tick]))        
    end 

    save_fld = "./plots/paper_results/plant_vs_weather"
    isdir(save_fld)|| mkdir(save_fld)

    pl1_fin = plot(pl1[1],pl1[2],pl1[3],pl1[4],layout=(1,4),size=(900,300),left_margin = 4Plots.mm) 
    savefig(pl1_fin,save_fld*"/gs_I_date.svg")
end