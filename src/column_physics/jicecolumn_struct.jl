# Written by Joseph Kump (josek97@utexas.edu)
# Outlines JICEColumn, a Julia structure for single-column ice thermodynamics

const V_a    = 0.025
const K_a    = 0.03
const T_mlt0 = 273.15   # K         melting point of freshwater
const ρ_0    = 917.0    # kg/m^3    density of fresh (pure) ice
const ρ_i    = 917.0    # kg/m^3    density of sea ice
const ρ_w    = 1025.0   # kg/m^3    density of seawater (based on common estimate)
const ρ_s    = 330.0    # kg/m^3    density of snow
const c_0    = 2106.0   # J/kg/K    specific heat of fresh ice
const c_w    = 3900.0   # J kg K    specific heat of seawater (estimate)
const c_h    = 0.006    #           heat transfer coefficient
const K_s    = 0.3      # W/m/deg   thermal conductivity of snow
const L_0    = 334000.0 # J / kg    latent heat of fusion of fresh ice
const κ_i    = 1.4      #           extinction coefficient
const μ      = 0.054    # deg/ppt   liquidus ratio between the freezing temperature and salinity of brine
const S_max  = 3.2      # ppt       maximum salinity of sea ice
const ahmax  = 0.3      # m         thickness above which albedo is constant
const α_o    = 0.06     #           ocean albedo

const α_icev    = 0.78  # visible ice albedo for h > ahmax
const α_icei    = 0.36  # near-ir ice albedo for h > ahmax
const α_snowv   = 0.98  # cold snow albedo, visible
const α_snowi   = 0.70  # cold snow albedo, near IR

const dα_mlt    = -0.075 # albedo change for temp change of 1 degree for ice
const dα_mltv   = -0.1   # albedo change for temp change of 1 degree for snow, visible
const dα_mlti   = -0.15  # albedo change for temp change of 1 degree for snow, infrared

const i0vis     = 0.7   # fraction of penetrating solar radiation at top surface

const puny       = 1.0e-11  # For numerical tests
const Tsf_errmax = 0.01     # For numerical test of convergence

#= JICEColumn struct
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

    Δh          (m)     thickness of each layer, Vector{Float64}
    Δh̄          (m)     averaged thickness of each layer, Vector{Float64}
    S           (ppt)   salinity, Vector{Float64}
    c_i         (J/kg)  specific heat of ice, Vector{Float64}
    K           (W/m K) thermal conductivity of sea ice, Vector{Float64}
    K̄           (W/m K) pairwise-averaged thermal conductivity of sea ice, Vector{Float64}
    I_pen       (W m^2) penetrating solar flux, Vector{Float64}
    maindiag    ()      main diagonal of tridiagonal system to solve, Vector{Float64}
    subdiag     ()      sub diagonal of tridiagonal system to solve, Vector{Float64}
    supdiag     ()      super diagonal of tridiagonal system to solve, Vector{Float64}

    F_Lu        (W/m^2)     Upward longwave flux
    F_s         (W/m^2)     Sensible heat flux
    F_l         (W/m^2)     Latent heat flux
    F_Ld        (W/m^2)     Downward longwave flux
    dF_Lu       (W/m^2 C)   Derivative of upward longwave flux relative to T_sf
    dF_s        (W/m^2 C)   Derivative of sensible heat flux relative to T_sf
    dF_l        (W/m^2 C)   Derivative of latent heat flux relative to T_sf

    T_array     (C)         array of temperatures stored at each timestep, Matrix{Float64}
    Δh_array    (M)         array of layer thicknesses stored at each timestep, Matrix{Float64}
=#
mutable struct JICEColumn

    # Variables that must be provided to initialize the model
    N_t::Int64
    N_i::Int64
    N_s::Int64

    H_i::Float64
    H_s::Float64
    T_frz::Float64
    Δt::Float64
    u_star::Float64
    T_w::Float64

    T_n::Vector{Float64}

    # Variables that are created based on the above:
    α_vdr::Vector{Float64}
    α_idr::Vector{Float64}
    α_vdf::Vector{Float64}
    α_idf::Vector{Float64}

    T_nplus::Vector{Float64}
    F_0::Vector{Float64}
    dF_0::Vector{Float64}

    Δh::Vector{Float64}
    Δh̄::Vector{Float64}
    S::Vector{Float64}
    c_i::Vector{Float64}
    K::Vector{Float64}
    K̄::Vector{Float64}
    I_pen::Vector{Float64}
    q_i::Vector{Float64}
    q_inew::Vector{Float64}
    z_old::Vector{Float64}
    z_new::Vector{Float64}
    maindiag::Vector{Float64}
    subdiag::Vector{Float64}
    supdiag::Vector{Float64}

    F_Lu::Vector{Float64}
    F_s::Vector{Float64}
    F_l::Vector{Float64}

    dF_Lu::Vector{Float64}
    dF_s::Vector{Float64}
    dF_l::Vector{Float64}

    T_array::Matrix{Float64}
    Δh_array::Matrix{Float64}

end

# Constructs a JICEColumn object given the initial parameters
function initialize_JICEColumn(N_t, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0)

    T_nplus, F_0, dF_0, Δh, Δh̄, S, c_i, K, K̄, I_pen, q_i, q_inew, z_old, z_new, maindiag, subdiag, supdiag, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l, T_array, Δh_array = allocate_memory(N_i, N_s, N_t)

    # Get initial thicknesses of each snow and ice layer:
    for k in 1:N_s
        Δh[k+1] = H_s / N_s
    end
    for k in 1:N_i
        Δh[k+N_s+1] = H_i / N_i
    end

    jcolumn = JICEColumn(N_t, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0,
                        zeros(Float64,2), zeros(Float64,2), zeros(Float64,2), zeros(Float64,2), T_nplus, F_0, dF_0,
                        Δh, Δh̄, S, c_i, K, K̄, I_pen, q_i, q_inew, z_old, z_new, maindiag,
                        subdiag, supdiag, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l, T_array, Δh_array)

    # Some preliminary work before running the model:
    jcolumn.T_nplus[:] = jcolumn.T_n

    generate_S(jcolumn.S, jcolumn.N_i)
    
    # TODO: more detailed implementation of α
    #generate_α(jcolumn.H_i, jcolumn.α_vdr, jcolumn.α_idr, jcolumn.α_vdf, jcolumn.α_idf, jcolumn.T_n[1])
                
    jcolumn.T_array[:, 1] = jcolumn.T_n
    jcolumn.Δh_array[:,1] = jcolumn.Δh

    return jcolumn
end

# Allocates all necessary memory for intermediate variables in the model
function allocate_memory(N_i, N_s, N_t)

    T_nplus = zeros(Float64, N_i+N_s+1)

    F_0  = zeros(Float64, N_t)
    dF_0 = zeros(Float64, N_t)

    # Getting Δh̄ (averages of adjacent thicknesses, length K)
    Δh = zeros(Float64, N_i+N_s+1)
    Δh̄ = zeros(Float64, N_i+N_s)

    # Other intermediate data to keep:
    S      = zeros(Float64, N_i)
    c_i    = zeros(Float64, N_i+N_s+1)
    K      = zeros(Float64, N_i+N_s+1)
    K̄      = zeros(Float64, N_i+N_s)
    I_pen  = zeros(Float64, N_i)
    q_i    = zeros(Float64, N_i)
    q_inew = zeros(Float64, N_i)
    z_old  = zeros(Float64, N_i+1)
    z_new  = zeros(Float64, N_i+1)

    maindiag = zeros(Float64, N_i+N_s+1)
    subdiag  = zeros(Float64, N_i+N_s)
    supdiag  = zeros(Float64, N_i+N_s)

    F_Lu = zeros(Float64, N_t)
    F_s  = zeros(Float64, N_t)
    F_l  = zeros(Float64, N_t)

    dF_Lu = zeros(Float64, N_t)
    dF_s  = zeros(Float64, N_t)
    dF_l  = zeros(Float64, N_t)

    T_array  = zeros(Float64, N_i+N_s+1, N_t+1)
    Δh_array = zeros(Float64, N_i+N_s+1, N_t+1)

    return T_nplus, F_0, dF_0, Δh, Δh̄, S, c_i, K, K̄, I_pen, q_i, q_inew, z_old, z_new, maindiag, subdiag, supdiag, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l, T_array, Δh_array

end

# Gets the salinity profile for this column of sea ice. Since this obeys the BL99 model the salinity
# is constant
@inline function generate_S(S, N_i)

    for k in 1:N_i
        z      = k / N_i
        S[k] = 0.5 * S_max * (1 - cos(pi*z^(0.407/(z+0.573))))
    end

    return nothing
end
