# Written by Joseph Kump (josek97@utexas.edu)
# Object that provides important adjacent variables for the sea ice model,
# particularly for getting atmospheric energy flux
# Initially, we'll assume fluxes are constant on our time scale

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
