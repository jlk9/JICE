# Written by Joseph Kump (josek97@utexas.edu)
# Object that provides important adjacent variables for the sea ice model,
# particularly for getting atmospheric energy flux

const emissivity = 0.985     # Emissivity of snow or ice
const sbc        = 0.0       # Stefan-Boltzman Constant (FIX)
const q_1        = 1.16378e7 # kg / m^3
const q_2        = 5897.8    # K
const κ          = 0.40      # von Karman constant

#=
Properties:
    nt          Number of time steps
    F_s         Sensible heat flux
    F_l         Latent heat flux
    F_Ld        Downward longwave flux
    F_Lu        Upward longwave flux
    F_sw        Shortwave flux
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
    C_l         nonlinear turbulent heat transfer coefficient
    C_s         nonlinear turbulent heat transfer coefficient
    U_a         Wind velocity
=#
mutable struct Atmodel

    nt::Int64

    F_s::Float64
    F_l::Float64
    F_Ld::Float64
    F_Lu::Float64
    F_sw::Float64

    T_a::Float64
    θ_a::Float64
    ρ_a::Float64
    Q_a::Float64
    c_p::Float64

    L_vap::Float64
    L_ice::Float64

    c_u::Float64
    c_θ::Float64
    c_q::Float64

    C_l::Float64
    C_s::Float64

    U_a::Vector{Float64}

end

# Creates an atmodel object and initial parameters
function initialize_atmodel(F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a)

    model = Atmodel(0.0, 0.0, F_Ld, 0.0, F_sw, T_a, Θ_a, ρ_a, Q_a)
    return model
end