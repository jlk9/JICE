# Written by Joseph Kump (josek97@utexas.edu)
# Outlines JICEColumn, a Julia structure for single-column ice thermodynamics

using CUDA, KernelAbstractions

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

    Δh          (m)     thickness of each layer, VFT
    Δh̄          (m)     averaged thickness of each layer, VFT
    S           (ppt)   salinity, VFT
    c_i         (J/kg)  specific heat of ice, VFT
    K           (W/m K) thermal conductivity of sea ice, VFT
    K̄           (W/m K) pairwise-averaged thermal conductivity of sea ice, VFT
    I_pen       (W m^2) penetrating solar flux, VFT
    maindiag    ()      main diagonal of tridiagonal system to solve, VFT
    subdiag     ()      sub diagonal of tridiagonal system to solve, VFT
    supdiag     ()      super diagonal of tridiagonal system to solve, VFT

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
mutable struct JICEColumnArrays{A, I, FT, VFT}

    onGPU::Bool
    architecture::A

    # Variables that must be provided to initialize the model
    N_t::I
    N_c::I
    N_i::I
    N_s::I
    N_layers::I

    T_frz::FT
    Δt::FT
    u_star::FT
    T_w::VFT

    T_n::VFT
    T_sfc::VFT

    # Variables that are created based on the above:
    H_i::VFT
    H_iold::VFT
    H_s::VFT

    rside::VFT
    fside::VFT
    f_bot::VFT

    α_vdr_i::VFT
    α_idr_i::VFT
    α_vdf_i::VFT
    α_idf_i::VFT

    α_vdr_s::VFT
    α_idr_s::VFT
    α_vdf_s::VFT
    α_idf_s::VFT

    T_nplus::VFT
    F_0::VFT
    dF_0::VFT

    Δh::VFT
    S::VFT
    c_i::VFT
    K::VFT
    K̄::VFT
    I_pen::VFT
    q::VFT
    q_new::VFT
    z_old::VFT
    z_new::VFT
    maindiag::VFT
    subdiag::VFT
    supdiag::VFT

    F_Lu::VFT
    F_s::VFT
    F_l::VFT

    dF_Lu::VFT
    dF_s::VFT
    dF_l::VFT

end


# Constructs a JICEColumn object given the initial parameters
function initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onGPU)

    # TODO: make Architecture able to check if CUDA is available here:
    if onGPU
        arch = GPU()
    else
        arch = CPU()
    end

    H_iold, rside, fside, f_bot, α_vdr_i, α_idr_i, α_vdf_i, α_idf_i, α_vdr_s, α_idr_s, α_vdf_s, α_idf_s, T_nplus, T_sfc, F_0, dF_0, Δh, S, c_i, K, K̄, I_pen, q, q_new, z_old, z_new, maindiag, subdiag, supdiag, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l = allocate_memory(N_c, N_i, N_i+N_s+1, arch)

    if arch == CPU()
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

    jcolumn = JICEColumnArrays(onGPU, arch, N_t, N_c, N_i, N_s, N_i+N_s+1, T_frz, Δt, u_star, T_w_d, T_0_d, T_sfc,
                               H_i_d, H_iold, H_s_d, rside, fside, f_bot, α_vdr_i, α_idr_i, α_vdf_i, α_idf_i, α_vdr_s, α_idr_s, α_vdf_s, α_idf_s,
                               T_nplus, F_0, dF_0, Δh, S, c_i, K, K̄, I_pen, q, q_new, z_old, z_new, maindiag, subdiag, supdiag,
                               F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l)
   
    if arch == CPU()
        kernel_generate_S! = generate_S!(arch,4)
        kernel_generate_initial_thicknesses! = generate_initial_thicknesses!(arch,4)
    else
        #numblocks = ceil(Int, jcolumn.N_layers*jcolumn.N_c/256)
        #@cuda threads=256 blocks=numblocks gpu_generate_S(jcolumn.S, jcolumn.N_i, jcolumn.N_s, jcolumn.N_c, jcolumn.N_layers)
        #@cuda threads=256 blocks=numblocks gpu_generate_initial_thicknesses(jcolumn.Δh, jcolumn.H_i, jcolumn.H_s, jcolumn.N_i, jcolumn.N_s, jcolumn.N_c, jcolumn.N_layers)
        kernel_generate_S! = generate_S!(arch,256)
        kernel_generate_initial_thicknesses! = generate_initial_thicknesses!(arch,256)
    end

    kernel_generate_S!(jcolumn.S, jcolumn.N_i, jcolumn.N_s, jcolumn.N_c, jcolumn.N_layers, ndrange=jcolumn.N_c*jcolumn.N_layers)
    kernel_generate_initial_thicknesses!(jcolumn.Δh, jcolumn.H_i, jcolumn.H_s, jcolumn.N_i, jcolumn.N_s, jcolumn.N_layers, ndrange=jcolumn.N_c*jcolumn.N_layers)

    return jcolumn
end

# Allocates all necessary memory for intermediate variables in the model
function allocate_memory(N_c, N_i, N_layers, arch)
    
    rside = zeros(Float64, 1) #N_c)
    fside = zeros(Float64, 1) #N_c)
    f_bot = zeros(Float64, 1) #N_c)

    # Other intermediate data to keep:
    q      = zeros(Float64, 1) #N_layers*N_c)
    q_new  = zeros(Float64, 1) #N_layers*N_c)
    z_old  = zeros(Float64, 1) #(N_layers+1)*N_c)
    z_new  = zeros(Float64, 1) #(N_layers+1)*N_c)
    

    # Either on host or device:
    if arch == CPU()
        F_Lu = zeros(Float64, N_c)
        F_s  = zeros(Float64, N_c)
        F_l  = zeros(Float64, N_c)

        dF_Lu = zeros(Float64, N_c)
        dF_s  = zeros(Float64, N_c)
        dF_l  = zeros(Float64, N_c)

        α_vdr_i = zeros(Float64, N_c)
        α_idr_i = zeros(Float64, N_c)
        α_vdf_i = zeros(Float64, N_c)
        α_idf_i = zeros(Float64, N_c)

        α_vdr_s = zeros(Float64, N_c)
        α_idr_s = zeros(Float64, N_c)
        α_vdf_s = zeros(Float64, N_c)
        α_idf_s = zeros(Float64, N_c)

        H_iold  = zeros(Float64, N_c)
        T_nplus = zeros(Float64, N_layers*N_c)
        T_sfc   = zeros(Float64, N_c)

        F_0  = zeros(Float64, N_c)
        dF_0 = zeros(Float64, N_c)
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
        F_Lu = CuArray(zeros(Float64, N_c))
        F_s  = CuArray(zeros(Float64, N_c))
        F_l  = CuArray(zeros(Float64, N_c))

        dF_Lu = CuArray(zeros(Float64, N_c))
        dF_s  = CuArray(zeros(Float64, N_c))
        dF_l  = CuArray(zeros(Float64, N_c))

        α_vdr_i = CuArray(zeros(Float64, N_c))
        α_idr_i = CuArray(zeros(Float64, N_c))
        α_vdf_i = CuArray(zeros(Float64, N_c))
        α_idf_i = CuArray(zeros(Float64, N_c))

        α_vdr_s = CuArray(zeros(Float64, N_c))
        α_idr_s = CuArray(zeros(Float64, N_c))
        α_vdf_s = CuArray(zeros(Float64, N_c))
        α_idf_s = CuArray(zeros(Float64, N_c))

        H_iold  = CuArray(zeros(Float64, N_c))
        T_nplus = CuArray(zeros(Float64, N_layers*N_c))
        T_sfc   = CuArray(zeros(Float64, N_c))

        F_0  = CuArray(zeros(Float64, N_c))
        dF_0 = CuArray(zeros(Float64, N_c))
        Δh   = CuArray(zeros(Float64, N_layers*N_c))

        S      = CuArray(zeros(Float64, N_layers*N_c))
        c_i    = CuArray(zeros(Float64, N_layers*N_c))
        K      = CuArray(zeros(Float64, N_layers*N_c))
        K̄      = CuArray(zeros(Float64, N_layers*N_c))
        I_pen  = CuArray(zeros(Float64, N_i*N_c))

        maindiag = CuArray(zeros(Float64, N_layers*N_c))
        subdiag  = CuArray(zeros(Float64, N_layers*N_c))
        supdiag  = CuArray(zeros(Float64, N_layers*N_c))
    end

    return H_iold, rside, fside, f_bot, α_vdr_i, α_idr_i, α_vdf_i, α_idf_i, α_vdr_s, α_idr_s, α_vdf_s, α_idf_s, T_nplus, T_sfc, F_0, dF_0, Δh, S, c_i, K, K̄, I_pen, q, q_new, z_old, z_new, maindiag, subdiag, supdiag, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l

end

# Gets the salinity profile for this column of sea ice. Since this obeys the BL99 model the salinity
# is constant
@kernel function generate_S!(S, N_i, N_s, N_c, N_layers)

    index = @index(Global)
    k = (index - 1) % N_layers + 1
    z = (k - N_s - 1) / N_i
    if k > N_s + 1
        S[index] = 0.5S_max * (1 - cos(pi*z^(0.407/(z+0.573))))
    end
end
#=
# Gets the salinity profile for this column of sea ice. Since this obeys the BL99 model the salinity
# is constant
@inline function generate_S(S, N_i, N_s, N_c, N_layers)
    
    for index in 1:(N_layers*N_c)
        k = (index - 1) % N_layers + 1
        z = (k - N_s - 1) / N_i
        if k > N_s + 1
            S[index] = 0.5S_max * (1 - cos(pi*z^(0.407/(z+0.573))))
        end
    end

    return nothing
end

# Gets the salinity profile for this column of sea ice. Since this obeys the BL99 model the salinity
# is constant
@inline function gpu_generate_S(S, N_i, N_s, N_c, N_layers)
    
    index_start  = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride       = gridDim().x * blockDim().x

    for index in index_start:stride:(N_c*N_layers)
        k = (index - 1) % N_layers + 1
        z = (k - N_s - 1) / N_i
        if k > N_s + 1
            S[index] = 0.5S_max * (1 - cos(pi*z^(0.407/(z+0.573))))
        end
    end

    return nothing
end
=#
# Generates initial thicknesses of ice columns.
@kernel function generate_initial_thicknesses!(Δh, H_i, H_s, N_i, N_s, N_layers)
    
    index = @index(Global)
    col = ((index - 1) ÷ N_layers) + 1
    k   = (index - 1) % N_layers + 1
    if k > 1 && k < N_s + 2
        Δh[index] = H_s[col] / N_s
    elseif k > N_s + 1
        Δh[index] = H_i[col] / N_i
    end
end
#=
# Generates initial thicknesses of ice columns.
@inline function generate_initial_thicknesses(Δh, H_i, H_s, N_i, N_s, N_c, N_layers)
    
    for index in 1:(N_layers*N_c)
        col = ((index - 1) ÷ N_layers) + 1
        k   = (index - 1) % N_layers + 1
        if k > 1 && k < N_s + 2
            Δh[index] = H_s[col] / N_s
        elseif k > N_s + 1
            Δh[index] = H_i[col] / N_i
        end
    end
end

# Generates initial thicknesses of ice columns.
@inline function gpu_generate_initial_thicknesses(Δh, H_i, H_s, N_i, N_s, N_c, N_layers)
    
    index_start  = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride       = gridDim().x * blockDim().x

    for index in index_start:stride:(N_c*N_layers)
        col = ((index - 1) ÷ N_layers) + 1
        k   = (index - 1) % N_layers + 1
        if k > 1 && k < N_s + 2
            Δh[index] = H_s[col] / N_s
        elseif k > N_s + 1
            Δh[index] = H_i[col] / N_i
        end
    end
end
=#