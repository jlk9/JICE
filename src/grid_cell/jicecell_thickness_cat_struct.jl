# Written by Joseph Kump (josek97@utexas.edu)
# Outlines a single thickness category of the JICECell object, JICECell_thickness_cat.
# This resembles the JICEColumn struct in column_physics, except that it contains subarrays

include("../constants/jice_constants.jl")

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
    Δh_array    (m)         array of layer thicknesses stored at each timestep, Matrix{Float64}
=#
mutable struct JICECell_Cat

    # Variables that must be provided to initialize the model
    T_n::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}

    # Variables that are created based on the above:
    H_i::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    H_iold::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    H_s::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    α_vdr::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    α_idr::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    α_vdf::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    α_idf::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}

    T_nplus::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    F_0::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    dF_0::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}

    Δh::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    S::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    c_i::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    K::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    K̄::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    I_pen::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    q::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    q_new::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    z_old::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    z_new::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    maindiag::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    subdiag::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    supdiag::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}

    F_Lu::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    F_s::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    F_l::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}

    dF_Lu::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    dF_s::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}
    dF_l::SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}

end

# Constructs a JICEColumn object given the initial parameters
function initialize_JICECell_Cat(N_t, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0)

    H_i_array, H_s_array, T_nplus, F_0, dF_0, Δh, S, c_i, K, K̄, I_pen, q_i, q_inew, z_old, z_new, maindiag, subdiag, supdiag, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l, T_array, Δh_array = allocate_memory(N_i, N_s, N_t)

    # Get initial thicknesses of each snow and ice layer:
    for k in 1:N_s
        Δh[k+1] = H_s / N_s
    end
    for k in 1:N_i
        Δh[k+N_s+1] = H_i / N_i
    end

    jcolumn = JICECell_Cat(deepcopy(T_0), H_i .+ zeros(Float64,1), zeros(Float64,1), H_s .+ zeros(Float64,1),
                        zeros(Float64,2), zeros(Float64,2), zeros(Float64,2), zeros(Float64,2), T_nplus, F_0, dF_0,
                        Δh, S, c_i, K, K̄, I_pen, q_i, q_inew, z_old, z_new, maindiag,
                        subdiag, supdiag, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l)

    # Some preliminary work before running the model:
    jcolumn.T_nplus[:] = jcolumn.T_n

    generate_S(jcolumn.S, jcolumn.N_i)
    
    # Set initial temperatures and thicknesses:
    jcolumn.H_i_array[1]  = H_i
    jcolumn.H_s_array[1]  = H_s
    jcolumn.T_array[:, 1] = jcolumn.T_n
    jcolumn.Δh_array[:,1] = jcolumn.Δh

    return jcolumn
end

# Allocates all necessary memory for intermediate variables in the model
function allocate_memory(N_i, N_s, N_t)

    H_i_array = zeros(Float64, N_t+1)
    H_s_array = zeros(Float64, N_t+1)
    T_nplus   = zeros(Float64, N_i+N_s+1)

    F_0  = zeros(Float64, N_t)
    dF_0 = zeros(Float64, N_t)
    Δh   = zeros(Float64, N_i+N_s+1)

    # Other intermediate data to keep:
    S      = zeros(Float64, N_i)
    c_i    = zeros(Float64, N_i+N_s+1)
    K      = zeros(Float64, N_i+N_s+1)
    K̄      = zeros(Float64, N_i+N_s)
    I_pen  = zeros(Float64, N_i)
    q      = zeros(Float64, N_i+N_s+1)
    q_new = zeros(Float64, N_i+N_s+1)
    z_old  = zeros(Float64, N_i+N_s+2)
    z_new  = zeros(Float64, N_i+N_s+2)

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

    return H_i_array, H_s_array, T_nplus, F_0, dF_0, Δh, S, c_i, K, K̄, I_pen, q, q_new, z_old, z_new, maindiag, subdiag, supdiag, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l, T_array, Δh_array

end