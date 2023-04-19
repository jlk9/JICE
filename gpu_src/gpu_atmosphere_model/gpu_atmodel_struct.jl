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
mutable struct ATModelArrays

    N_t::Int64
    N_c::Int64

    F_SWvdr::AbstractArray{Float64, 1}
    F_SWidr::AbstractArray{Float64, 1}
    F_SWvdf::AbstractArray{Float64, 1}
    F_SWidf::AbstractArray{Float64, 1}

    F_Ld::AbstractArray{Float64, 1}

    T_a::AbstractArray{Float64, 1}
    Θ_a::AbstractArray{Float64, 1}
    ρ_a::AbstractArray{Float64, 1}
    Q_a::AbstractArray{Float64, 1}
    c_p::AbstractArray{Float64, 1}

    c_u::AbstractArray{Float64, 1}
    c_Θ::AbstractArray{Float64, 1}
    c_q::AbstractArray{Float64, 1}

    atm_u_star::AbstractArray{Float64, 1}

    U_a::AbstractArray{Float64, 1}

    Q_sfc::AbstractArray{Float64, 1}
    F_SWsfc::AbstractArray{Float64, 1}
    F_SWpen::AbstractArray{Float64, 1}

end

# Creates an atmodel object and initial parameters
function initialize_ATModelArrays(N_t, N_c, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onGPU)

    if onGPU
        F_SWvdr_d = CuArray(F_SWvdr)
        F_SWidr_d = CuArray(F_SWidr)
        F_SWvdf_d = CuArray(F_SWvdf)
        F_SWidf_d = CuArray(F_SWidf)
        F_Ld_d    = CuArray(F_Ld)
        T_a_d     = CuArray(T_a)
        Θ_a_d     = CuArray(Θ_a)
        ρ_a_d     = CuArray(ρ_a)
        Q_a_d     = CuArray(Q_a)
        c_p_d     = CuArray(c_p)
        U_a_d     = CuArray(U_a)

        c_u_d = CUDA.zeros(Float64, N_c)
        c_Θ_d = CUDA.zeros(Float64, N_c)
        c_q_d = CUDA.zeros(Float64, N_c)

        atm_u_star_d = CUDA.zeros(Float64, N_c)
        Q_sfc_d      = CUDA.zeros(Float64, N_c)
        F_SWsfc_d    = CUDA.zeros(Float64, N_c)
        F_SWpen_d    = CUDA.zeros(Float64, N_c)
    else
        F_SWvdr_d = deepcopy(F_SWvdr)
        F_SWidr_d = deepcopy(F_SWidr)
        F_SWvdf_d = deepcopy(F_SWvdf)
        F_SWidf_d = deepcopy(F_SWidf)
        F_Ld_d    = deepcopy(F_Ld)
        T_a_d     = deepcopy(T_a)
        Θ_a_d     = deepcopy(Θ_a)
        ρ_a_d     = deepcopy(ρ_a)
        Q_a_d     = deepcopy(Q_a)
        c_p_d     = deepcopy(c_p)
        U_a_d     = deepcopy(U_a)

        c_u_d = zeros(Float64, N_c)
        c_Θ_d = zeros(Float64, N_c)
        c_q_d = zeros(Float64, N_c)

        atm_u_star_d = zeros(Float64, N_c)
        Q_sfc_d      = zeros(Float64, N_c)
        F_SWsfc_d    = zeros(Float64, N_c)
        F_SWpen_d    = zeros(Float64, N_c)
    end
    
    model = ATModelArrays(N_t, N_c, F_SWvdr_d, F_SWidr_d, F_SWvdf_d, F_SWidf_d, F_Ld_d, T_a_d, Θ_a_d, ρ_a_d, Q_a_d, c_p_d,
                          c_u_d, c_Θ_d, c_q_d, atm_u_star_d, U_a_d, Q_sfc_d, F_SWsfc_d, F_SWpen_d)
    
    return model
end
