# Written by Joseph Kump (josek97@utexas.edu)
# Computes the thermodynamic temperature changes in once column
# of sea ice

#using LinearAlgebra

# Runs a single time step of the temperature changes using FDM
@inline function step_temp_change(N_i, N_s, S, T_frz, Δh, Δh̄, T_old, T_new, c_i, K, K̄, I_pen, F_0, dF_0, maindiag, subdiag, supdiag, Δt)

    # Getting the average thicknesses of each layer:
    for i in 1:(N_i+N_s)
        Δh̄[i] = (Δh[i] + Δh[i+1]) ./ 2
    end

    # Ice thermal conductivity (length N_i+1)
    generate_K(K, N_i, N_s, S, T_old)

    # Specific heat (length N_i+1), generated by formula
    generate_c_i(c_i, N_i, N_s, c_0, S, T_old, T_new)
    # Getting the right hand side. First coefficients for the current temps T_old:
    generate_coeffs(T_new, T_old, N_i, N_s, dF_0)
    # Get the Matrix and complete RHS:
    generate_matrix_rhs(N_i, N_s, Δh, c_i, K, K̄, dF_0, F_0, T_frz, I_pen, maindiag, subdiag, supdiag, T_new, Δt)

    #mat = Tridiagonal(subdiag, maindiag, supdiag)

    #println(c_i)
    #println(K)
    #println(T_old)
    #println(mat)
    #println(T_new)
    # Calculate new time step
    tridiagonal_solve(T_new, N_i+N_s, maindiag, subdiag, supdiag)

    return nothing
end

# Generates the explicit coefficients for time step T_n, which are used for the RHS
# of the tridiagonal system
@inline function generate_coeffs(coeffs, old, N_i, N_s, dF_0)

    for i in 1:(N_i+N_s+1)
        coeffs[i] = old[i]
    end
    coeffs[1] *= dF_0

    return nothing
end

# Generates the implicit (tridiagonal) matrix for time step T_nplus
# I think this matrix is stable in the Thomas algorithm because it is
# diagonally dominant
@inline function generate_matrix_rhs(N_i, N_s, Δh, c_i, K, K̄, dF_0, F_0, T_frz, I_pen, maindiag, subdiag, supdiag, rhs, Δt)

    # Getting K̄ (length N_i):
    for k in 1:(N_i+N_s)
        K̄[k] = (2*K[k]*K[k+1]) / (Δh[k+1]*K[k] + Δh[k]*K[k+1])
    end

    # Top surface, based on (50) in science guide
    maindiag[1] = (dF_0 - K̄[1])
    supdiag[1]  = K̄[1]
    rhs[1]     -= F_0

    # Interior snow layers, based on (52) in science guide
    for k in 2:(N_s+1)

        η_k          = Δt / (ρ_s*c_i[k]*Δh[k+1])
        maindiag[k]  = 1 + η_k*(K̄[k-1] + K̄[k])
        subdiag[k-1] = -η_k*K̄[k-1]
        supdiag[k]   = -η_k*K̄[k]
    end
    # Interior ice layers, based on (52) in science guide
    for k in (N_s+2):(N_i+N_s)

        η_k          = Δt / (ρ_i*c_i[k]*Δh[k+1])
        maindiag[k]  = 1 + η_k*(K̄[k-1] + K̄[k])
        subdiag[k-1] = -η_k*K̄[k-1]
        supdiag[k]   = -η_k*K̄[k]
        rhs[k]      += η_k*I_pen[k-N_s-1]
    end

    # Bottom surface
    η_end               = Δt / (ρ_i*c_i[N_i+N_s+1]*Δh[N_i+N_s+1])
    maindiag[N_i+N_s+1] = 1 + η_end*(K̄[N_i+N_s] + 0.5*K[N_i+N_s+1]) #1 + η_Ni*(K[N_i] + K[N_i+1])
    subdiag[N_i+N_s]    = -η_end*K̄[N_i+N_s]
    rhs[N_i+N_s+1]     += η_end*(I_pen[N_i] + 0.5*K[N_i+N_s+1]*T_frz) #η_Ni*(I_pen[N_i] + K[N_i+1]*T_frz)

    return nothing
end

# Solves a tridiagonal system using the tridiagonal matrix algorithm (aka Thomas algorithm)
# x contains the original rhs values which are permuted
# the subdiag and supdiag entries are also permuted
@inline function tridiagonal_solve(x, N, maindiag, subdiag, supdiag)

    @inbounds for i in 2:(N+1)

        w           = (subdiag[i-1]/maindiag[i-1])
        maindiag[i] = maindiag[i] - w*supdiag[i-1]
        x[i]        = x[i] - w*x[i-1]
    end

    x[N+1] = x[N+1]/maindiag[N+1]

    @inbounds for i in N:-1:1

        x[i] = (x[i] - supdiag[i]*x[i+1])/maindiag[i]
    end

    return nothing

end

# Gets the melting point as a function of salinity
@inline function t_mlt(S)

    return 273.15 - 0.0592.*S - 9.37*(10^(-6)).*(S.^2) - 5.33*(10^(-7)).*(S.^3)

end

# Gets the thermal conductivity at this time step. (Maykut,1971)
# 2.03 is K_0, conductivity of fresh ice
# 0.13 is β, an empirical constant
@inline function generate_K(K, N_i, N_s, S, T)

    # Surface, either snow or ice
    K[1] = 2.03
    if N_s > 0
        K[1] = 0.3
    end
    # Snow layers
    for k in 2:(N_s+1)
        K[k] = 0.3
    end
    # Ice layers
    for i in (N_s+2):(N_i+N_s+1)
    
        K[i] = 2.03 + (0.13*S[i-N_s-1])/T[i]
    end

    return nothing
end

# Older, more complicated function for thermal conductivity
@inline function generate_K_old(K, K_b, N_i, N_s, θ, T_mltS)

    for i in 1:(N_i+1)
    
        K_b[i] = 1.162 * (0.45  + 0.0108*θ[i]  + 0.0000504*(θ[i]^2))
        K[i]   = 1.162 * (1.905 + 0.00866*θ[i] + 0.0000297*(θ[i]^2))
        K[i]  *= (2*K[i] + K_a - (2*V_a)*(K[i]-K_a)) / (2*K[i] + K_a - (2*V_a)*(K[i]-K_a))
        K[i]  -= (K[i]-K_b[i])*((T_mltS-T_mlt0)/θ[i])
    end

    return nothing
end

# gets the specific heat of sea ice
@inline function generate_c_i(c_i, N_i, N_s, c_0, S, T_old, T_new)

    # Snow and ice surface both match fresh ice
    for k in 1:(N_s+1)
        c_i[k] = c_0
    end
    # Deeper layers of ice
    for k in (N_s+2):(N_i+N_s+1)
        c_i[k] = c_0 + (L_0*μ*S[k-N_s-1])/(T_old[k]*T_new[k])
    end

    return nothing
end