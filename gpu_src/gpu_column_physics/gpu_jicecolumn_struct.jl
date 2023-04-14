# Written by Joseph Kump (josek97@utexas.edu)
# Outlines JICEColumn, a Julia structure for single-column ice thermodynamics

using CUDA

include("../../src/constants/jice_constants.jl")

#= JICEColumnArrays struct
Properties:
    N_t    (dim'less)   number of time steps, int
    N_i    (dim'less)   number of ice layers, int
    N_s    (dim'less)   number of snow layers, int
    H_i    (m)          total ice thickness, float
    H_s    (m)          total snow thickness, float
    T_frz  (K)          freezing point of water at salinity S, float
    Δt     (sec)        size of time step, float
    u_star ()           friction velocity, float
    T_w    (K)          temperature of sea surface
    α_vdr  ()           visible direct albedo, ice and snow
    α_idr  ()           infrared direct albedo, ice and snow
    α_vdf  ()           visible diffuse albedo, ice and snow
    α_idf  ()           infrared diffuse albedo, ice and snow
    T_0    (C)          initial temps, array of floats (length N_i+1, including skin layer)
    F_0    (W/m^2)      total heat flux at surface, array of floats
    dF_0   (W/m^2 C)    derivative total heat flux at surface, array of floats

    Δh          (m)     thickness of each layer, AbstractArray{Float64, 1}
    Δh̄          (m)     averaged thickness of each layer, AbstractArray{Float64, 1}
    S           (ppt)   salinity, AbstractArray{Float64, 1}
    c_i         (J/kg)  specific heat of ice, AbstractArray{Float64, 1}
    K           (W/m K) thermal conductivity of sea ice, AbstractArray{Float64, 1}
    K̄           (W/m K) pairwise-averaged thermal conductivity of sea ice, AbstractArray{Float64, 1}
    I_pen       (W m^2) penetrating solar flux, AbstractArray{Float64, 1}
    maindiag    ()      main diagonal of tridiagonal system to solve, AbstractArray{Float64, 1}
    subdiag     ()      sub diagonal of tridiagonal system to solve, AbstractArray{Float64, 1}
    supdiag     ()      super diagonal of tridiagonal system to solve, AbstractArray{Float64, 1}

    F_Lu        (W/m^2)     Upward longwave flux
    F_s         (W/m^2)     Sensible heat flux
    F_l         (W/m^2)     Latent heat flux
    F_Ld        (W/m^2)     Downward longwave flux
    dF_Lu       (W/m^2 C)   Derivative of upward longwave flux relative to T_sf
    dF_s        (W/m^2 C)   Derivative of sensible heat flux relative to T_sf
    dF_l        (W/m^2 C)   Derivative of latent heat flux relative to T_sf

    T_array     (C)         array of temperatures stored at each timestep, Matrix{T}
    Δh_array    (m)         array of layer thicknesses stored at each timestep, Matrix{T}
=#
struct JICEColumnArrays

    # Variables that must be provided to initialize the model
    N_t::Int64
    N_c::Int64
    N_i::Int64
    N_s::Int64
    N_layers::Int64

    T_frz::Float64
    Δt::Float64
    u_star::Float64
    T_w::AbstractArray{Float64, 1}

    T_n::AbstractArray{Float64, 1}

    # Variables that are created based on the above:
    H_i::AbstractArray{Float64, 1}
    H_iold::AbstractArray{Float64, 1}
    H_s::AbstractArray{Float64, 1}

    rside::AbstractArray{Float64, 1}
    fside::AbstractArray{Float64, 1}
    f_bot::AbstractArray{Float64, 1}

    α_vdr::AbstractArray{Float64, 1}
    α_idr::AbstractArray{Float64, 1}
    α_vdf::AbstractArray{Float64, 1}
    α_idf::AbstractArray{Float64, 1}

    T_nplus::AbstractArray{Float64, 1}
    F_0::AbstractArray{Float64, 1}
    dF_0::AbstractArray{Float64, 1}

    Δh::AbstractArray{Float64, 1}
    S::AbstractArray{Float64, 1}
    c_i::AbstractArray{Float64, 1}
    K::AbstractArray{Float64, 1}
    K̄::AbstractArray{Float64, 1}
    I_pen::AbstractArray{Float64, 1}
    q::AbstractArray{Float64, 1}
    q_new::AbstractArray{Float64, 1}
    z_old::AbstractArray{Float64, 1}
    z_new::AbstractArray{Float64, 1}
    maindiag::AbstractArray{Float64, 1}
    subdiag::AbstractArray{Float64, 1}
    supdiag::AbstractArray{Float64, 1}

    F_Lu::AbstractArray{Float64, 1}
    F_s::AbstractArray{Float64, 1}
    F_l::AbstractArray{Float64, 1}

    dF_Lu::AbstractArray{Float64, 1}
    dF_s::AbstractArray{Float64, 1}
    dF_l::AbstractArray{Float64, 1}

end

struct JICEColumnArrays_Device

    T_w::CuArray{Float64}
    T_n::CuArray{Float64}
    H_i::CuArray{Float64}
    H_iold::CuArray{Float64}
    H_s::CuArray{Float64}
    T_nplus::CuArray{Float64}
    F_0::CuArray{Float64}
    dF_0::CuArray{Float64}
    Δh::CuArray{Float64}
    S::CuArray{Float64}
    c_i::CuArray{Float64}
    K::CuArray{Float64}
    K̄::CuArray{Float64}
    I_pen::CuArray{Float64}
    maindiag::CuArray{Float64}
    subdiag::CuArray{Float64}
    supdiag::CuArray{Float64}

end


# Constructs a JICEColumn object given the initial parameters
function initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onDevice)

    H_iold, rside, fside, f_bot, α_vdr, α_idr, α_vdf, α_idf, T_nplus, F_0, dF_0, Δh, S, c_i, K, K̄, I_pen, q, q_new, z_old, z_new, maindiag, subdiag, supdiag, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l = allocate_memory(N_c, N_i, N_i+N_s+1, N_t, onDevice)

    if !onDevice
        T_w_d = deepcopy(T_w)
        T_0_d = deepcopy(T_0)
        H_i_d = deepcopy(H_i)
        H_s_d = deepcopy(H_s)
    else
        T_w_d = CuArray(T_w)
        T_0_d = CuArray(T_0)
        H_i_d = CuArray(H_i)
        H_s_d = CuArray(H_s)
    end

    copyto!(T_nplus, T_0_d)

    jcolumn = JICEColumnArrays(N_t, N_c, N_i, N_s, N_i+N_s+1, T_frz, Δt, u_star, T_w_d, T_0_d,
                               H_i_d, H_iold, H_s_d, rside, fside, f_bot, α_vdr, α_idr, α_vdf, α_idf,
                               T_nplus, F_0, dF_0, Δh, S, c_i, K, K̄, I_pen, q, q_new, z_old, z_new, maindiag, subdiag, supdiag,
                               F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l)

    # Some preliminary work before running the model:
    for k in 1:(jcolumn.N_layers*jcolumn.N_c)
        jcolumn.T_nplus[k] = jcolumn.T_n[k]
    end

    generate_S(jcolumn.S, jcolumn.N_i, jcolumn.N_s, jcolumn.N_c, jcolumn.N_layers)
    
    # Get thicknesses of each snow and ice layer:
    for n in 1:jcolumn.N_c
        Δh_i = jcolumn.H_i[n] / jcolumn.N_i
        Δh_s = jcolumn.H_s[n] / jcolumn.N_s
        index = (n-1) * jcolumn.N_layers
        for k in 2:(jcolumn.N_s+1)
            jcolumn.Δh[index + k] = Δh_s
        end
        for k in (jcolumn.N_s+2):jcolumn.N_layers
            jcolumn.Δh[index + k] = Δh_i
        end
    end

    return jcolumn
end

# Allocates all necessary memory for intermediate variables in the model
function allocate_memory(N_c, N_i, N_layers, N_t, onDevice)

    rside = zeros(Float64, N_c)
    fside = zeros(Float64, N_c)
    f_bot = zeros(Float64, N_c)

    α_vdr = zeros(Float64, 2N_c)
    α_idr = zeros(Float64, 2N_c)
    α_vdf = zeros(Float64, 2N_c)
    α_idf = zeros(Float64, 2N_c)

    # Other intermediate data to keep:
    q      = zeros(Float64, N_layers*N_c)
    q_new  = zeros(Float64, N_layers*N_c)
    z_old  = zeros(Float64, (N_layers+1)*N_c)
    z_new  = zeros(Float64, (N_layers+1)*N_c)

    F_Lu = zeros(Float64, N_t*N_c)
    F_s  = zeros(Float64, N_t*N_c)
    F_l  = zeros(Float64, N_t*N_c)

    dF_Lu = zeros(Float64, N_t*N_c)
    dF_s  = zeros(Float64, N_t*N_c)
    dF_l  = zeros(Float64, N_t*N_c)

    # Either on host or device:
    if !onDevice
        H_iold  = zeros(Float64, N_c)
        T_nplus = zeros(Float64, N_layers*N_c)

        F_0  = zeros(Float64, N_t*N_c)
        dF_0 = zeros(Float64, N_t*N_c)
        Δh   = zeros(Float64, N_layers*N_c)

        S      = zeros(Float64, N_layers*N_c)
        c_i    = zeros(Float64, N_layers*N_c)
        K      = zeros(Float64, N_layers*N_c)
        K̄      = zeros(Float64, N_layers*N_c)
        I_pen  = zeros(Float64, N_i*N_c)

        maindiag = zeros(Float64, N_layers*N_c)
        subdiag  = zeros(Float64, N_layers*N_c)
        supdiag  = zeros(Float64, N_layers*N_c)
    else
        H_iold  = CUDA.zeros(Float64, N_c)
        T_nplus = CUDA.zeros(Float64, N_layers*N_c)

        F_0  = CUDA.zeros(Float64, N_t*N_c)
        dF_0 = CUDA.zeros(Float64, N_t*N_c)
        Δh   = CUDA.zeros(Float64, N_layers*N_c)

        S      = CUDA.zeros(Float64, N_layers*N_c)
        c_i    = CUDA.zeros(Float64, N_layers*N_c)
        K      = CUDA.zeros(Float64, N_layers*N_c)
        K̄      = CUDA.zeros(Float64, N_layers*N_c)
        I_pen  = CUDA.zeros(Float64, N_i*N_c)

        maindiag = CUDA.zeros(Float64, N_layers*N_c)
        subdiag  = CUDA.zeros(Float64, N_layers*N_c)
        supdiag  = CUDA.zeros(Float64, N_layers*N_c)
    end

    return H_iold, rside, fside, f_bot, α_vdr, α_idr, α_vdf, α_idf, T_nplus, F_0, dF_0, Δh, S, c_i, K, K̄, I_pen, q, q_new, z_old, z_new, maindiag, subdiag, supdiag, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l

end

# Gets the salinity profile for this column of sea ice. Since this obeys the BL99 model the salinity
# is constant
@inline function generate_S(S, N_i, N_s, N_c, N_layers)

    for k in 1:N_i
        z = k / N_i
        for l in 0:(N_c-1)
            S[l*N_layers + N_s + 1 + k] = 0.5S_max * (1 - cos(pi*z^(0.407/(z+0.573))))
        end
    end

    return nothing
end