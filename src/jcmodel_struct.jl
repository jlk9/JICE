# Written by Joseph Kump (josek97@utexas.edu)
# Outlines JCModel, a Julia structure for single-column ice thermodynamics

const V_a    = 0.025
const K_a    = 0.03
const T_mlt0 = 273.15   # K      melting point of freshwater
const ρ_0    = 917.0    # kg/m^3 density of fresh (pure) ice
const ρ_i    = 917.0    # kg/m^3 density of sea ice
const ρ_w    = 1025.0   # kg/m^3 density of seawater (based on common estimate)
const c_0    = 2106.0   # J/kg/K specific heat of fresh ice
const c_w    = 3900.0   # J kg K specific heat of seawater (estimate)
const c_h    = 0.006    #        heat transfer coefficient
const L_0    = 334000.0 # J / kg latent heat of fusion of fresh ice
const S_max  = 3.2      # ppt    maximum salinity of sea ice
const ahmax  = 0.3      # m      thickness above which albedo is constant

const puny       = 1.0e-11  # For numerical tests
const Tsf_errmax = 0.01     # For numerical test of convergence

#= JCModel struct
Properties:
    N_i    (dim'less)   number of ice layers, int (we consider the surface "skin layer")
    nt     (dim'less)   number of time steps, int
    H      (m)          total ice thickness, float
    L      (J kg^-1)    latent heat of sublimation, float
    T_frz  (K)          freezing point of water at salinity S, float
    i_0    (W m^-2)     the penetrating solar flux reduction at the top surface, float
    κ_i    ()           extinction coefficient, float
    Δt     (sec)        size of time step, float
    u_star ()           friction velocity, float
    T_w    (K)          temperature of sea surface
    α      ()           albedo, float
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

    T_array     ()      array of temperatures stored at each timestep, Matrix{Float64}
    Δh_array    ()      array of layer thicknesses stored at each timestep, Matrix{Float64}
=#
mutable struct JCModel

    # Variables that must be provided to initialize the model
    N_i::Int64
    nt::Int64

    H::Float64
    L::Float64
    T_frz::Float64
    i_0::Float64
    κ_i::Float64
    Δt::Float64
    u_star::Float64
    T_w::Float64
    α::Float64

    T_0::Vector{Float64}
    # Top layer fluxes
    # TODO: derive from lines ~98-217 in icepack_therm_shared.F90
    F_0::Vector{Float64}
    dF_0::Vector{Float64}

    # Variables that are created based on the above:
    Δh::Vector{Float64}
    Δh̄::Vector{Float64}
    S::Vector{Float64}
    c_i::Vector{Float64}
    K::Vector{Float64}
    K̄::Vector{Float64}
    I_pen::Vector{Float64}
    maindiag::Vector{Float64}
    subdiag::Vector{Float64}
    supdiag::Vector{Float64}

    T_array::Matrix{Float64}
    Δh_array::Matrix{Float64}

end

# Constructs a JCModel object given the initial parameters
function initialize_JCModel(N_i, nt, H, L, T_frz, i_0, κ_i, Δt, u_star, T_w, T_0, F_0, dF_0)

    Δh, Δh̄, S, c_i, K, K̄, I_pen, maindiag, subdiag, supdiag, T_array, Δh_array = allocate_memory(H, N_i)

    model = JCModel(N_i, nt, H, L, T_frz, i_0, κ_i, Δt, u_star, T_w, 0.0, T_0, F_0, dF_0,
                   Δh, Δh̄, S, c_i, K, K̄, I_pen, maindiag, subdiag, supdiag, T_array, Δh_array)

    return model
end

# Allocates all necessary memory for intermediate variables in the model
function allocate_memory(H, N_i)

    # Get the initial thicknesses of each layer
    h_star  = min(0.5, 0.25*H)
    Δh_1    = min(h_star, (H - h_star)/(N_i-1))
    Δh_k    = (H - Δh_1) / (N_i-1)

    Δh    = Δh_k .+ zeros(Float64, N_i+1)
    Δh[1] = 0.0
    Δh[2] = Δh_1

    # Getting Δz̄ (averages of adjacent thicknesses, length K)
    Δh̄ = zeros(Float64, N_i)

    # Other intermediate data to keep:
    S     = zeros(Float64, N_i+1)
    c_i   = zeros(Float64, N_i+1)
    K     = zeros(Float64, N_i+1)
    K̄     = zeros(Float64, N_i)
    I_pen = zeros(Float64, N_i+1)

    maindiag = zeros(Float64, N_i+1)
    subdiag  = zeros(Float64, N_i)
    supdiag  = zeros(Float64, N_i)
    T_array  = zeros(Float64, N_i+1, nt+1)
    Δh_array = zeros(Float64, N_i+1, nt+1)

    return Δh, Δh̄, S, c_i, K, K̄, I_pen, maindiag, subdiag, supdiag, T_array, Δh_array

end