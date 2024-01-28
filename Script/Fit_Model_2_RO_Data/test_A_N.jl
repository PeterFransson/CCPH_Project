function calc_A_N(Tₐ,Nₘ_f_vec)    
    y_vec = zeros(size(Nₘ_f_vec))
    A_vec = zeros(size(Nₘ_f_vec))
    cost_vec = zeros(size(Nₘ_f_vec))
    
    for i in eachindex(Nₘ_f_vec)
        #Parameters
        α_max = 0.14
        a_Jmax = 0.011
        k = 0.52
        m = 0.05
        r_gₛ = 0.42
        Nₛ = 0.041 #Control
        #Traits
        gₛ = 0.1
        Nₘ_f = Nₘ_f_vec[i]
        #Weather variables
        Tₐ = Tₐ
        I₀ = 600*10^-6
        θₛ = 0.175
        VPD = 250.0

        env = CCPH.EnvironmentStruct(θₛ=θₛ,I₀=I₀,Tₐ=Tₐ,VPD=VPD)
        photo_kinetic = CCPH.PhotoKineticRates()
        photopar = CCPH.PhotoPar(photo_kinetic,env.Tₐ;α=α_max)

        #Calculate Jₘₐₓ
        Jₘₐₓ = CCPH.Calc_Jₘₐₓ(Nₘ_f,a_Jmax,0.0,photopar.b_Jmax,1.0) 
        #Quantum yield
        α = CCPH.Calc_α(1.0,α_max)
        #Irradiance incident on a leaf at canopy top
        Iᵢ = env.I₀*k/(1-m)    
        #calcualte total conductance
        gₜ = gₛ*r_gₛ   
        #calculate electron transport
        J = CCPH.Calc_J(Iᵢ,Jₘₐₓ,α,photopar.θ)
        #Calcualte intercellular carbon dioxide concentration
        cᵢ = CCPH.calc_opt_cᵢ(gₜ,J,photopar,env)  
        #Calculate leaf C assimilation
        A = CCPH.calc_Assimilation(gₜ,cᵢ,env.P,env.Cₐ)
        #Calculate N cost
        N_cost = Jₘₐₓ*Nₛ

        A_vec[i] = A
        cost_vec[i] = N_cost
        y_vec[i] = A-N_cost        
    end
    return (y_vec,A_vec,cost_vec)
end


function test_A_N()
    Nₘ_f_vec = range(0.005,stop=0.035,length=200)

    y_vec_5,A_vec_5,cost_vec_5 = calc_A_N(5.0,Nₘ_f_vec)
    y_vec_10,A_vec_10,cost_vec_10 = calc_A_N(10.0,Nₘ_f_vec)

    pl1 = plot(Nₘ_f_vec*100,y_vec_5*10^6,ylabel="A-Cost (μmol s⁻¹ m⁻²)",xlabel="Nₘ_f (%)",label="5 °C")
    plot!(pl1,Nₘ_f_vec*100,y_vec_10*10^6,label="10 °C")

    pl2 = plot(Nₘ_f_vec,A_vec_5*10^6,ylabel="(μmol s⁻¹ m⁻²)",xlabel="Nₘ_f (%)",label="A 5 °C")
    plot!(pl2,Nₘ_f_vec,cost_vec_5*10^6,label="Cost 5 °C")
    plot!(pl2,Nₘ_f_vec,A_vec_10*10^6,label="A 10 °C")
    plot!(pl2,Nₘ_f_vec,cost_vec_10*10^6,label="Cost 10 °C")

    plot(pl1,pl2,layout=(2,1))

    savefig("./plots/A_Cost_plot.svg")    

    
    pl3 = plot(Nₘ_f_vec*100,A_vec_5*10^6,ylabel="(μmol s⁻¹ m⁻²)",xlabel="Nₘ_f (%)",label="A 5 °C")
    plot!(pl3,Nₘ_f_vec*100,cost_vec_5*10^6,label="Cost 5 °C")
    plot!(pl3,Nₘ_f_vec*100,A_vec_10*10^6,label="A 10 °C")
    plot!(pl3,Nₘ_f_vec*100,cost_vec_10*10^6,label="Cost 10 °C")

    savefig(pl3,"./plots/A_Cost_plot_2.svg")    
end

test_A_N()