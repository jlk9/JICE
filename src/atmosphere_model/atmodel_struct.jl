# Written by Joseph Kump (josek97@utexas.edu)
# Object that provides important adjacent variables for the sea ice model,
# particularly for getting atmospheric energy flux
# Initially, we'll assume fluxes are constant on our time scale

const emissivity = 0.985                # Emissivity of snow or ice
const sbc        = 5.670374419e-8       # Stefan-Boltzman Constant (FIX)
const q_1        = 1.16378e7            # kg / m^3
const q_2        = 5897.8               # K
const κ          = 0.40                 # von Karman constant
const g          = 9.8                  # gravitational acceleration (positive downward)
const z_ref      = 10.0                 # ice reference height
const U_dmin     = 1.0                  # minimum allowable value for |U_a| (since we don't use high frequency)
const z_deg      = 1.0                  # Assumed level height
const C_to_K     = 273.15               # converts Celcius to Kelvin
const λ          = -2.3025850929940455  # log(z_deg / z_ref)
const L_vap      = 2260.0               # Latent heat of vaporization
const L_ice      = 334.0                # Latent heat of fusion
#=
Properties:
    N_t         Number of time steps
    F_sw        Shortwave flux    (Common is 120 W/m^2)
    T_a         Atmosphere temperature
    θ_a         Air potential temperature
    ρ_a         Surface air density
    Q_a         Air specific humidity
    c_p         Specific heat of air
    c_u         exchange coefficient
    c_θ         exchange coefficient
    c_q         exchange coefficient
    U_a         Wind velocity
=#
mutable struct ATModel

    N_t::Int64

    F_SWvdr::Float64
    F_SWidr::Float64
    F_SWvdf::Float64
    F_SWidf::Float64

    F_Ld::Float64

    T_a::Float64
    Θ_a::Float64
    ρ_a::Float64
    Q_a::Float64
    c_p::Float64

    c_u::Vector{Float64}
    c_Θ::Vector{Float64}
    c_q::Vector{Float64}

    atm_u_star::Vector{Float64}

    U_a::Vector{Float64}

end

# Creates an atmodel object and initial parameters
function initialize_ATModel(N_t, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a)
    
    model = ATModel(N_t, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p,
                    zeros(Float64, 1), zeros(Float64, 1), zeros(Float64, 1), zeros(Float64, 1), U_a)
    
    return model
end
