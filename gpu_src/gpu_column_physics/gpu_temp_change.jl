# Written by Joseph Kump (josek97@utexas.edu)
# Computes the thermodynamic temperature changes in once column
# of sea ice

using CUDA

# Runs a single time step of the temperature changes using FDM
function step_temp_change(N_c, N_i, N_s, N_layers, H_s, S, T_frz, Δh, T_old, T_new, c_i, K, K̄, I_pen, F_0, dF_0, maindiag, subdiag, supdiag, Δt, onGPU)

    if onGPU

        numblocks = ceil(Int, N_layers*N_c/256)
        # Ice thermal conductivity (length N_i+1)
        @cuda threads=256 blocks=numblocks gpu_generate_K(K, N_c, N_s, N_layers, S, T_old)
        #gpu_generate_K(K, N_c, N_s, N_layers, H_s, S, T_old)
    
        # Specific heat (length N_i+1), generated by formula
        @cuda threads=256 blocks=numblocks gpu_generate_c_i(c_i, N_c, N_s, N_layers, S, T_old, T_new)
        
        # Get the Matrix and RHS:
        @cuda threads=256 blocks=numblocks gpu_generate_matrix_rhs(N_c, N_i, N_s, N_layers, Δh, c_i, K, K̄, dF_0, F_0, T_frz, I_pen, maindiag, subdiag, supdiag, T_new, T_old, Δt)
	    #=
        handle     = CUSPARSE.cusparseCreate()
        #bufferTemp = zeros(UInt64, 1)
	    #bufferSize = 0
	    bufferSize = zeros(UInt64, 1)
        CUSPARSE.cusparseDgtsv2StridedBatch_bufferSizeExt(handle, N_layers, subdiag, maindiag, supdiag, T_new, N_c, N_layers, pointer(bufferSize))
	    #println(bufferSize)
	    pbuffer = CUDA.CuPtr{Nothing}
	    CUDA.cudaMalloc(pointer(pbuffer), bufferSize)
	    #CUDA.Mem.alloc(pbuffer, bufferSize)
	    CUSPARSE.cusparseDgtsv2StridedBatch(handle, N_layers, subdiag, maindiag, supdiag, T_new, N_c, N_layers, pbuffer)
	    CUSPARSE.gtsv2!(subdiag, maindiag, supdiag, T_new, 'O', false)
        =#
	    #CUSPARSE.gtsv2!(subdiag, maindiag, supdiag, T_new)
        numblocks = ceil(Int, N_c/256)
        @cuda threads=256 blocks=numblocks gpu_batched_tridiagonal_solve(T_new, N_layers-1, N_c, maindiag, subdiag, supdiag)
    else
        # Ice thermal conductivity (length N_i+1)
        generate_K(K, N_c, N_s, N_layers, S, T_old)
        #gpu_generate_K(K, N_c, N_s, N_layers, H_s, S, T_old)
    
        # Specific heat (length N_i+1), generated by formula
        generate_c_i(c_i, N_c, N_s, N_layers, S, T_old, T_new)
        
        # Get the Matrix and RHS:
        generate_matrix_rhs(N_c, N_i, N_s, N_layers, Δh, c_i, K, K̄, dF_0, F_0, T_frz, I_pen, maindiag, subdiag, supdiag, T_new, T_old, Δt)

        batched_tridiagonal_solve(T_new, N_layers-1, N_c, maindiag, subdiag, supdiag)
    end


    return nothing
end

# Produces the main diagonal of the matrices for the tridiagonal solve and temperature changes
@inline function generate_matrix_rhs(N_c, N_i, N_s, N_layers, Δh, c_i, K, K̄, dF_0, F_0, T_frz, I_pen, maindiag, subdiag, supdiag, rhs, T_old, Δt)

    # First get K̄, then fill in entries of maindiag::
    for index in 1:(N_c*N_layers)
        col = ((index-1) ÷ N_layers) + 1
        k   = (index-1) % N_layers

        # First get K̄:
        if index % N_layers == 0
            K̄[index] = 0.5K[index]
        else
            K̄[index] = (2*K[index]*K[index+1]) / (Δh[index+1]*K[index] + Δh[index]*K[index+1])
        end
        
        # Initialize maindiag as 1's, off diags and rhs as 0's:
        maindiag[index] = 1.0
        rhs[index]      = 0.0
        subdiag[index]  = 0.0
        supdiag[index]  = 0.0
        
        # Compute η_k (Δh index special for lower boundary), then multiply by snow or ice density
        η_k  = Δt / (c_i[index]*Δh[min(index+1, col*N_layers)])
        η_k /= ρ_s + (k > N_s) * (ρ_i - ρ_s) # conditional is for ice density
        if k == 0
            maindiag[index] = dF_0[col] - K̄[index]
            rhs[index]      = dF_0[col]*T_old[index] - F_0[col]
        else
            maindiag[index] = 1 + η_k*(K̄[index-1] + K̄[index])
            rhs[index]      = T_old[index] #- (k > N_s) * 
        end

        if k > N_s
            rhs[index] += η_k * I_pen[(col-1)*N_i + k - N_s]
            if k == N_layers - 1
                rhs[index] += η_k * 0.5K[index]*T_frz
            end
        end

        # Setting values of subdiag and supdiag
        if k != 0
            subdiag[index] = -η_k * K̄[index-1]
        end
        if k != N_layers-1
            supdiag[index] = K̄[index]
            if k != 0
                supdiag[index] *= -η_k
            end
        end
    end
    return nothing
end

# Produces the main diagonal of the matrices for the tridiagonal solve and temperature changes
@inline function gpu_generate_matrix_rhs(N_c, N_i, N_s, N_layers, Δh, c_i, K, K̄, dF_0, F_0, T_frz, I_pen, maindiag, subdiag, supdiag, rhs, T_old, Δt)

    index_start  = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride       = gridDim().x * blockDim().x

    for index in index_start:stride:(N_c*N_layers)
        col = ((index-1) ÷ N_layers) + 1
        k   = (index-1) % N_layers

        # First get K̄:
        if index % N_layers == 0
            K̄[index] = 0.5K[index]
        else
            K̄[index] = (2*K[index]*K[index+1]) / (Δh[index+1]*K[index] + Δh[index]*K[index+1])
        end
        
        # Initialize maindiag as 1's, off diags and rhs as 0's:
        maindiag[index] = 1.0
        rhs[index]      = 0.0
        subdiag[index]  = 0.0
        supdiag[index]  = 0.0
        
        # Compute η_k (Δh index special for lower boundary), then multiply by snow or ice density
        η_k  = Δt / (c_i[index]*Δh[min(index+1, col*N_layers)])
        η_k /= ρ_s + (k > N_s) * (ρ_i - ρ_s) # conditional is for ice density
        if k == 0
            maindiag[index] = dF_0[col] - K̄[index]
            rhs[index]      = dF_0[col]*T_old[index] - F_0[col]
        else
            maindiag[index] = 1 + η_k*(K̄[index-1] + K̄[index])
            rhs[index]      = T_old[index] #- (k > N_s) * 
        end

        if k > N_s
            rhs[index] += η_k * I_pen[(col-1)*N_i + k - N_s]
            if k == N_layers - 1
                rhs[index] += η_k * 0.5K[index]*T_frz
            end
        end

        # Setting values of subdiag and supdiag
        if k != 0
            subdiag[index] = -η_k * K̄[index-1]
        end
        if k != N_layers-1
            supdiag[index] = K̄[index]
            if k != 0
                supdiag[index] *= -η_k
            end
        end
    end
    return nothing
end

# Solves a large number of tridiagonal systems in batch on a GPU
function gpu_batched_tridiagonal_solve(x, N, count, maindiag, subdiag, supdiag)

    index_start  = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride       = gridDim().x * blockDim().x

    for c in index_start:stride:count

        # First entry for the tridiagonal system this thread solves:
        system_start = (c-1) * (N+1) + 1
        
        for i in (system_start+1):(system_start+N)
            w           = (subdiag[i-1]/maindiag[i-1])
            maindiag[i] = maindiag[i] - w*supdiag[i-1]
            #x[i]        = x[i] - w*x[i-1]
        end
        #=
        x[N+system_start] = x[N+system_start]/maindiag[N+system_start]

        for i in N+system_start-1:-1:system_start
            x[i] = (x[i] - supdiag[i]*x[i+1])/maindiag[i]
        end
        =#
    end
end

# Solves a large number of tridiagonal systems in batch
function batched_tridiagonal_solve(x, N, count, maindiag, subdiag, supdiag)

    for c in 1:count

        x_n        = view(x, ((c-1)*(N+1)+1):c*(N+1))
        maindiag_n = view(maindiag, ((c-1)*(N+1)+1):c*(N+1))
        subdiag_n  = view(subdiag, ((c-1)*(N+1)+2):c*(N+1))
        supdiag_n  = view(supdiag, ((c-1)*(N+1)+1):c*(N+1)-1)

        tridiagonal_solve(x_n, N, maindiag_n, subdiag_n, supdiag_n)
    end
end

# Solves a tridiagonal system using the tridiagonal matrix algorithm (aka Thomas algorithm)
# x contains the original rhs values which are permuted
# the subdiag and supdiag entries are also permuted
@inline function tridiagonal_solve(x, N, maindiag, subdiag, supdiag)

    @inbounds for i in 2:(N+1)

        w           = (subdiag[i-1]/maindiag[i-1])
        maindiag[i] = maindiag[i] - w*supdiag[i-1]
        #x[i]        = x[i] - w*x[i-1]
    end
    #=
    x[N+1] = x[N+1]/maindiag[N+1]

    @inbounds for i in N:-1:1

        x[i] = (x[i] - supdiag[i]*x[i+1])/maindiag[i]
    end
    =#
    return nothing

end

# Gets the thermal conductivity at this time step. (Maykut,1971)
# 2.03 is K_0, conductivity of fresh ice
# 0.13 is β, an empirical constant
@inline function generate_K(K, N_c, N_s, N_layers, S, T)
    
    for index in 1:(N_c*N_layers)
        col      = ((index-1) ÷ N_layers) + 1
        k        = (index-1) % N_layers
        K[index] = 0.3
        if k > N_s
            K[index] = 2.03 + ((0.13S[index])/T[index])
        end
    end
    
    return nothing
end

# Generating K on the GPU:
function gpu_generate_K(K, N_c, N_s, N_layers, S, T)
    
    index_start  = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride       = gridDim().x * blockDim().x

    for index in index_start:stride:(N_c*N_layers)
        col      = ((index-1) ÷ N_layers) + 1
        k        = (index-1) % N_layers
        K[index] = 0.3
        if k > N_s
            K[index] = 2.03 + ((0.13S[index])/T[index])
        end
    end
    
    return nothing
end

# gets the specific heat of sea ice
@inline function generate_c_i(c_i, N_c, N_s, N_layers, S, T_old, T_new)
    
    for index in 1:(N_c*N_layers)
        k = (index-1) % N_layers
        # snow/fresh ice     conditional for interior ice
        c_i[index] = c_0
        if k > N_s
            c_i[index] += (L_0*μ*S[index])/(T_old[index]*T_new[index])
        end
    end

    return nothing
end

# gets the specific heat of sea ice
@inline function gpu_generate_c_i(c_i, N_c, N_s, N_layers, S, T_old, T_new)
    
    index_start  = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride       = gridDim().x * blockDim().x

    for index in index_start:stride:(N_c*N_layers)
        k = (index-1) % N_layers
        # snow/fresh ice     conditional for interior ice
        c_i[index] = c_0
        if k > N_s
            c_i[index] += (L_0*μ*S[index])/(T_old[index]*T_new[index])
        end
    end

    return nothing
end