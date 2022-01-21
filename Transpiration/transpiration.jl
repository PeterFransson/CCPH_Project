#Parameters used to estimate canopy transpiration
mutable struct TranspirationParameter{T<:Float64}
    ρₐ::T #density of dry air (Kg m⁻³)
    cₚₐ::T #Specific heat capacity of dry air (J kg⁻¹ K⁻¹)
    γ::T #Psychrometric constant (Pa K⁻¹)
    λ::T #Latent heat of water vaporization (J kg⁻¹)       
end
TranspirationParameter(;ρₐ::T=1.204,cₚₐ::T=1004.0,γ::T=66.1,λ::T=2.454e6) where {T<:Float64} = TranspirationParameter(ρₐ,cₚₐ,γ,λ)

#derivative of saturation Vapor Pressure (Pa) with respect to temperature at T degree celsius
#Derived from Equation 21 in Alduchov and Eskridge (1996)
dSVPₜ(T::AbstractFloat) = 610.94*17.625*243.04*exp(17.625*T/(243.04+T))/(243.04+T)^2

function CalcE_eq(ϕₙₐ::T,Tₐ::T;par::TranspirationParameter=TranspirationParameter()) where {T<:Float64}
    #Equation for calcualting the equilibrium evaporation rate E_eq (Kg s⁻¹ m⁻²)
    #ϕₙₐ net radiation absorbed by the canopy (J s⁻¹ m⁻²)
    #Tₐ air temperature (°C)    
    s = dSVPₜ(Tₐ)    
    E_eq = s*ϕₙₐ/((s+par.γ)*par.λ)
    return E_eq
end

function CalcE_imp(g_C::T,VPD::T;par::TranspirationParameter=TranspirationParameter()) where {T<:Float64}
    #Equation for calcualting the imposed evaporation rate E_imp (Kg s⁻¹ m⁻²)
    #g_C canopy conductance (water based conductance) (m³ s⁻¹ m⁻²)
    #VPD vapour-pressure deficit (Pa)
    E_imp = par.ρₐ*par.cₚₐ*g_C*VPD/(par.γ*par.λ)
    return E_imp
end

function CalcE(ϕₙₐ::T,Tₐ::T,g_C::T,VPD::T,Ω::T;par::TranspirationParameter=TranspirationParameter()) where {T<:Float64}
    #Equaiton for estimating the canopy tranpiraiton E (Kg s⁻¹ m⁻²)
    #ϕₙₐ net radiation absorbed by the canopy (J s⁻¹ m⁻²)
    #Tₐ air temperature (°C)
    #g_C canopy conductance (water based conductance) (m³ s⁻¹ m⁻²)
    #VPD vapour-pressure deficit (Pa)
    #Ω a dimentionless number in the range 0.0-1.0
    E_eq = CalcE_eq(ϕₙₐ,Tₐ;par=par)
    E_imp = CalcE_imp(g_C,VPD;par=par)
    E = Ω*E_eq+(1- Ω)*E_imp
    return E
end

function Estg_C(E::T,ϕₙₐ::T,Tₐ::T,VPD::T,Ω::T;par::TranspirationParameter=TranspirationParameter()) where {T<:Float64}
    #Equaiton for estimating the canopy conductance g_C (m³ s⁻¹ m⁻²)
    #E canopy tranpiraiton (Kg s⁻¹ m⁻²)
    #ϕₙₐ net radiation absorbed by the canopy (J s⁻¹ m⁻²)
    #Tₐ air temperature (°C)    
    #VPD vapour-pressure deficit (Pa)
    #Ω a dimentionless number in the range 0.0-1.0
    E_eq = CalcE_eq(ϕₙₐ,Tₐ;par)
    g_C = (E-Ω*E_eq)/(1-Ω)*(par.γ*par.λ)/(par.ρₐ*par.cₚₐ*VPD)
    return g_C
end

function PAR2Rad(PAR::T;conv_coeff::T=10^6/2.3) where {T<:Float64}
    #Funciton for converting From PAR (mol s⁻¹ m⁻²) to global
    #radiation (Rad, J s⁻¹ m⁻²) 
    Rad = PAR*conv_coeff
    return Rad
end 

function gₛ2g_C(gₛ::T,PAR::T,LAI::T;cons::CCPH.Constants=CCPH.Constants(),treepar::CCPH.TreePar=CCPH.TreePar(),
    I_gₛ::T=50.0e-6) where {T<:Float64}
    #Estimates canopy conductance g_C (m³ s⁻¹ m⁻²)
    #gₛ stomatal conductanceg (mol CO₂ s⁻¹ m⁻²)
    #PAR at canopy top (mol s⁻¹ m⁻²) 
    #LAI Canopy leaf area (m² leaf area m⁻² ground area)
    #I_gₛ Irradiace level at which gₛ is half of canopy top gₛ (mol s⁻¹ m⁻²) 
       
    g_C = cons.r*cons.M_H2O/cons.ρ_H2O*gₛ*log((I_gₛ+PAR)/(I_gₛ+PAR*exp(-treepar.k*LAI)))/treepar.k    
    
    return g_C
end

function Calctranpiration(VPD_z::T,Sₑ::T;E_cm_ref::T=1.812,s_VPD_z::T=3.121,
    s_Sₑ::T=18.342) where {T<:Float64}
    #Tor-Ngern et al. 2017 model for estimating canopy transpiration E_c (mm/d)
    #Standard paramater values (E_cm_ref,s_VPD_z,s_Sₑ) are taken from the 
    #original paper for the Rosinedal scots pine stand
    #--Input--
    #VPD_z vapor pressure deficit (Pa)
    #Sₑ effective saturation (-)    

    E_cm = E_cm_ref*(1-exp(-s_VPD_z*VPD_z*10^-3))
    E_c = E_cm*(1-exp(-s_Sₑ*Sₑ))

    return E_c
end

function Est_gₛ(LAI::T,daylight::T;cons::CCPH.Constants=CCPH.Constants(),
    env::CCPH.EnvironmentStruct=CCPH.EnvironmentStruct()) where {T<:Float64}
    # Estimate stomatal conductance gₛ (mol C m⁻² leaf area s⁻¹) based on 
    # Weather data and leaf area index, using the canopy transpiration from
    # Tor-Ngern et al. 2017
    #--Input-- 
    # LAI leaf area index (m² leaf area m⁻² ground area)   
    # daylight lenght of a single day (s)
    VPD = env.VPD
    θₛ = env.θₛ
    P = env.P    
    VPD_z = VPD.*daylight/(24*3600)    
    Sₑ = CCPH.CalcSₑ.(θₛ)
    E_c = Calctranpiration.(VPD_z,Sₑ) #mm/day
    E_c *= 10^-3*cons.ρ_H2O/cons.M_H2O/daylight #mol m⁻² s⁻¹
    gₛ = E_c*P/(VPD*cons.r*LAI)
    return gₛ
end