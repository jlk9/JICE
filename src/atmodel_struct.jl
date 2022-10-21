# Written by Joseph Kump (josek97@utexas.edu)
# Object that provides important adjacent variables for the sea ice model,
# particularly for getting atmospheric energy flux
# Initially, we'll assume fluxes are constant on our time scale

const emissivity = 0.985            # Emissivity of snow or ice
const sbc        = 5.670374419e-8   # Stefan-Boltzman Constant (FIX)
const q_1        = 1.16378e7        # kg / m^3
const q_2        = 5897.8           # K
const κ          = 0.40             # von Karman constant
const g          = 9.8              # gravitational acceleration (positive downward)
const z_ref      = 10.0             # ice reference height
const U_dmin     = 1.0              # minimum allowable value for |U_a| (since we don't use high frequency)
const z_deg      = 1.0              # Assumed level height

#=
Properties:
    nt          Number of time steps
    F_s         Sensible heat flux
    F_l         Latent heat flux
    F_Ld        Downward longwave flux
    F_Lu        Upward longwave flux
    F_sw        Shortwave flux    (Common is 120 W/m^2)
    T_a         Atmosphere temperature
    θ_a         Air potential temperature
    ρ_a         Surface air density
    Q_a         Air specific humidity
    c_p         Specific heat of air
    L_vap       Latent heat of vaporization
    L_ice       Latent heat of fusion
    c_u         exchange coefficient
    c_θ         exchange coefficient
    c_q         exchange coefficient
    Υ           "stability"
    C_l         nonlinear turbulent heat transfer coefficient
    C_s         nonlinear turbulent heat transfer coefficient
    U_a         Wind velocity
=#
mutable struct ATModel

    nt::Int64

    F_s::Float64
    F_l::Float64
    F_Ld::Float64
    F_Lu::Float64
    F_sw::Float64

    T_a::Float64
    Θ_a::Float64
    ρ_a::Float64
    Q_a::Float64
    c_p::Float64

    L_vap::Float64
    L_ice::Float64

    c_u::Float64
    c_Θ::Float64
    c_q::Float64

    Y::Float64

    C_l::Float64
    C_s::Float64

    U_a::Vector{Float64}

end

# Creates an atmodel object and initial parameters
function initialize_ATModel(nt, F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, L_vap, L_ice, U_a)

    model = ATModel(nt, 0.0, 0.0, F_Ld, 0.0, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, L_vap, L_ice, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, U_a)
    return model
end

#= Sets initial flux values for model run. This is only done once
    at the start of the model run, to provide (constant) flux terms
=#
function set_atm_helper_values(atmodel, T_sfc, z_ice, natmiter)

    λ = log(z_deg / z_ref)

    # Compute initial exchange coefficients:
    atmodel.c_u = κ / log(z_ref / z_ice)
    atmodel.c_Θ = atmodel.c_u
    atmodel.c_q = atmodel.c_u

    Q_sfc = (q_1/atmodel.ρ_a)*exp(-q_2 / (T_sfc + 273.15))

    # Iterate and update exchange coefficients:
    for k = 1:natmiter

        # Update turbulent scales
        u_star = atmodel.c_u * max(U_dmin, (atmodel.U_a[1]^2+atmodel.U_a[2]^2+atmodel.U_a[3]^2)^.5)
        Θ_star = atmodel.c_Θ * (atmodel.Θ_a - T_sfc)
        Q_star = atmodel.c_q * (atmodel.Q_a - Q_sfc)

        # Update Y and compute χ
        Y = (κ*g*z_deg)*(Θ_star / (atmodel.Θ_a*(1+0.606atmodel.Q_a)) + Q_star / (1/0.606 + atmodel.Q_a))/(u_star^2)

        # Compute ψ_m and ψ_s TODO: add alternative if Y is unstable? (< 0)
        ψ_m = -(0.7Y + 0.75*(Y-14.3)*exp(-0.35Y) + 10.7)
        ψ_s = ψ_m

        if Y < 0
            χ   = (1 - 16Y)^0.25
            ψ_m = 2*log(0.5*(1+χ)) + log(0.5*(1+χ^2)) - 2*atan(χ) + π/2
            ψ_s = 2*log(0.5*(1+χ^2))
        end

        # Update exchange coefficients
        atmodel.c_u = atmodel.c_u / (1 + atmodel.c_u*(λ-ψ_m)/κ)
        atmodel.c_Θ = atmodel.c_Θ / (1 + atmodel.c_Θ*(λ-ψ_s)/κ)
        atmodel.c_q = atmodel.c_Θ

    end

    # Compute heatflux coefficients. TODO: wind stress?
    atmodel.C_l = atmodel.ρ_a*(atmodel.L_vap + atmodel.L_ice)*u_star*atmodel.c_q
    atmodel.C_s = atmodel.ρ_a*atmodel.c_p*u_star*atmodel.c_Θ + 1

    return nothing
end

# Computes values of fluxes affected by T_sf
function set_atm_flux_values(atmodel, T_sf)

    # Outgoing longwave flux:
    atmodel.F_Lu = emissivity*sbc*(T_sf + 273.15)^4

    # Sensible heat flux:
    atmodel.F_s = atmodel.C_s*(atmodel.Θ_a - (T_sf + 273.15))

    # Latent heat flux:
    Q_sf        = (q_1/atmodel.ρ_a)*exp(-q_2 / (T_sf + 273.15))
    atmodel.F_l = atmodel.C_l*(atmodel.Q_a - Q_sf)

end

# Computes derivatives of fluxes affected by T_sf, based on analytical formula for above
function set_atm_dflux_values(atmodel, T_sf)

    # Outgoing longwave flux:
    dF_Lu = 4.0*emissivity*sbc*(T_sf + 273.15)^3

    # Sensible heat flux:
    dF_s  = -atmodel.C_s

    # Latent heat flux:
    dQ_sf = (q_1/atmodel.ρ_a)*exp(-q_2 / (T_sf + 273.15))*(-q_2 / (T_sf + 273.15)^2)
    dF_l = -atmodel.C_l*dQ_sf

    return dF_Lu + dF_s + dF_l

end