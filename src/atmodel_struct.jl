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
    N_t          Number of time steps
    F_s         Sensible heat flux
    F_l         Latent heat flux
    F_Ld        Downward longwave flux
    F_Lu        Upward longwave flux
    F_Lu        Derivative of upward longwave flux relative to T_sf
    F_s         Derivative of sensible heat flux relative to T_sf
    F_l         Derivative of latent heat flux relative to T_sf
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

    N_t::Int64

    F_s::Float64
    F_l::Float64
    F_Ld::Float64
    F_Lu::Float64
    F_sw::Float64

    dF_Lu::Float64
    dF_s::Float64
    dF_l::Float64

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

    u_star::Float64
    Θ_star::Float64
    Q_star::Float64
    Y::Float64
    χ::Float64
    ψ_m::Float64
    ψ_s::Float64

    C_l::Float64
    C_s::Float64

    U_a::Vector{Float64}

end

# Creates an atmodel object and initial parameters
function initialize_ATModel(N_t, F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, L_vap, L_ice, U_a)

    model = ATModel(N_t, 0.0, 0.0, F_Ld, 0.0, F_sw, 0.0, 0.0, 0.0, T_a, Θ_a, ρ_a, Q_a, c_p, L_vap, L_ice,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, U_a)
    return model
end
