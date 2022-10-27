# Written by Joseph Kump (josek97@utexas.edu)
# Basic code that is just ice temperature in column, with no formulation or melting.
# This version eliminates vectorized operations, enabling use of Enzyme.jl
# It also reduces the number of memory allocations substantially, increasing runtime
# and lowering memory footprint.

include("./jcmodel_struct.jl")
include("./surface_flux.jl")

#= Model Function
Input:
    jcmodel         an initialized JCModel object
    atmodel         an initialized ATModel object for surface flux computations

Output:
    Technically nothing, but the function updates jcmodel's T_array and Δh_arrays to have a log
    of updated ice temperatures and thicknesses
=#
function run_ice(jcmodel, atmodel)

    T_n     = deepcopy(jcmodel.T_0)
    T_nplus = deepcopy(jcmodel.T_0)

    generate_S(jcmodel.S, jcmodel.N_i)
    compute_surface_flux(jcmodel, atmodel)
    generate_I_pen(jcmodel.I_pen, jcmodel.i_0*(1-jcmodel.α)*atmodel.F_sw, jcmodel.κ_i, jcmodel.H, jcmodel.N_i)

    jcmodel.T_array[:, 1] = T_n
    jcmodel.Δh_array[:,1] = jcmodel.Δh

    # Main loop of temperature modifications:
    for step in 1:jcmodel.N_t

        run_ice_step(jcmodel.N_i, jcmodel.S, jcmodel.L, jcmodel.T_frz, jcmodel.κ_i, jcmodel.Δh, jcmodel.Δh̄, T_n, T_nplus,
                    jcmodel.c_i, jcmodel.K, jcmodel.K̄, jcmodel.I_pen, jcmodel.q_i, jcmodel.F_0[step], jcmodel.dF_0[step], 
                    jcmodel.maindiag, jcmodel.subdiag, jcmodel.supdiag, jcmodel.Δt, jcmodel.u_star, jcmodel.T_w)

        # Update T_n
        T_n[:] = T_nplus

        jcmodel.T_array[:, step+1] = T_n
        jcmodel.Δh_array[:,step+1] = jcmodel.Δh

    end

    return nothing
end

# Runs a single time step of the thermodynamic model for easier AD implementation
@inline function run_ice_step(N_i, S, L, T_frz, κ_i, Δh, Δh̄, T_old, T_new, c_i, K, K̄,
                              I_pen, q_i, F_0, dF_0, maindiag, subdiag, supdiag, Δt, u_star, T_w)

    # Getting the average thicknesses of each layer:
    for i in 1:N_i
        Δh̄[i] = (Δh[i] + Δh[i+1]) ./ 2
    end
    # Specific heat (length N_i+1), generated by formula
    generate_c_i(c_i, N_i, c_0, S, T_old, T_new)
    # Ice thermal conductivity (length N_i+1)
    generate_K(K, N_i, S, T_old)

    # Getting the right hand side. First coefficients for the current temps T_old:
    generate_coeffs(T_new, N_i, dF_0)
    # Multiply these with T_old to get rhs (stored inside T_new)
    for i in 1:(N_i+1)
        T_new[i] *= T_old[i]
    end

    # Get the Matrix and complete RHS:
    generate_matrix_rhs(N_i, Δh, c_i, K, K̄, dF_0, F_0, T_frz, I_pen, maindiag, subdiag, supdiag, T_new, Δt)

    # Calculate new time step
    tridiagonal_solve(T_new, N_i, maindiag, subdiag, supdiag)

    # Updating ice enthalpy:
    generate_q_i_from_T(q_i, N_i, T_new, S)

    # Updating bottom ice thickness:
    bottom_ice_change(N_i, Δh, T_new, q_i, K, T_frz, Δt, u_star, T_w)

    # Restore the sea ice thicknesses so that they are uniform again, and calculate new
    # enthalpies to preserve conservation of energy:
    #(Since K is no longer needed at this time step, we use it for intermediate steps)
    #rebalance_ice_layers(Δh, q_i, N_i, K, z_old, z_new)

    # Finally, recompute layer tempreatures to match new enthalpies:
    #generate_T_from_q_i(T_nplus, N_i, q_i, S)

    # TODO: add z thickness fields to jcmodel and compute them for thickness rebalancing

    return nothing

end

# Generates the explicit coefficients for time step T_n, which are used for the RHS
# of the tridiagonal system
@inline function generate_coeffs(coeffs, N_i, dF_0)

    coeffs[1] = dF_0

    for i in 2:(N_i+1)

        coeffs[i] = 1

    end

    return nothing

end

# Generates the implicit (tridiagonal) matrix for time step T_nplus
# I think this matrix is stable in the Thomas algorithm because it is
# diagonally dominant
@inline function generate_matrix_rhs(N_i, Δh, c_i, K, K̄, dF_0, F_0, T_frz, I_pen, maindiag, subdiag, supdiag, rhs, Δt)

    # Getting K̄ (length N_i):
    for k in 1:N_i
        K̄[k] = (2*K[k]*K[k+1]) / (Δh[k]*K[k] + Δh[k+1]*K[k+1])
    end

    # Interior ice layers, based on (52) in science guide
    for k in 2:N_i

        η_k = Δt / (ρ_i*c_i[k+1]*Δh[k+1])
        #= Pretty sure it's supposed to use K̄ but I'll follow documentation here
        # this is for 1:(N_i-1)
        maindiag[k]  = 1 + η_k*(K[k] + K[k+1])
        subdiag[k-1] = -η_k*K[k]
        supdiag[k]   = -η_k*K[k+1]
        rhs[k]      += η_k*I_pen[k]
        =#
        maindiag[k]  = 1 + η_k*(K̄[k-1] + K̄[k])
        subdiag[k-1] = -η_k*K̄[k-1]
        supdiag[k]   = -η_k*K̄[k]
        rhs[k]      += η_k*I_pen[k]
    end

    # Top surface, based on (50) in science guide
    maindiag[1] = (dF_0 - K̄[1])
    supdiag[1]  = K̄[1]
    rhs[1]     -= F_0

    # Bottom surface
    η_Ni            = Δt / (ρ_i*c_i[N_i+1]*Δh[N_i+1])
    maindiag[N_i+1] = 1 + η_Ni*(K̄[N_i] + 0.5*K[N_i+1]) #1 + η_Ni*(K[N_i] + K[N_i+1])
    subdiag[N_i]    = -η_Ni*K̄[N_i]
    rhs[N_i]       += η_Ni*(I_pen[N_i] + 0.5*K[N_i+1]*T_frz) #η_Ni*(I_pen[N_i] + K[N_i+1]*T_frz)

    #= OLD APPROACH
    # Interior ice layers, based on (52) in science guide
    for k in 2:N_i

        η_k = Δt / (ρ_i*c_i[k+1]*Δh[k+1])
        #= Pretty sure it's supposed to use K̄ but I'll follow documentation here
        # this is for 1:(N_i-1)
        maindiag[k+1] = 1 + η_k*(K̄[k] + K̄[k+1])
        subdiag[k]    = -η_k*K̄[k]
        supdiag[k+1]  = -η_k*K̄[k+1]
        =#
        maindiag[k]  = 1 + η_k*(K[k] + K[k+1])
        subdiag[k-1] = -η_k*K[k]
        supdiag[k]   = -η_k*K[k+1]
        rhs[k]      += η_k*I_pen[k]
    end

    # Top surface, based on (50) in science guide
    maindiag[1] = (dF_0 - K̄[1])
    supdiag[1]  = K̄[1]
    rhs[1]     -= F_0

    # Bottom surface
    η_Ni            = Δt / (ρ_i*c_i[N_i+1]*Δh[N_i+1])
    maindiag[N_i+1] = 1 + η_Ni*(K[N_i] + K[N_i+1])
    subdiag[N_i]    = -η_Ni*K[N_i]
    rhs[N_i]       += η_Ni*(I_pen[N_i] + K[N_i+1]*T_frz)
    =#
    return nothing

end

# Solves a tridiagonal system using the tridiagonal matrix algorithm (aka Thomas algorithm)
# x contains the original rhs values which are permuted
# the subdiag and supdiag entries are also permuted
@inline function tridiagonal_solve(x, N_i, maindiag, subdiag, supdiag)

    @inbounds for i in 2:(N_i+1)

        w           = (subdiag[i-1]/maindiag[i-1])
        maindiag[i] = maindiag[i] - w*supdiag[i-1]
        x[i]        = x[i] - w*x[i-1]
    end

    x[N_i+1] = x[N_i+1]/maindiag[N_i+1]

    @inbounds for i in N_i:-1:1

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
@inline function generate_K(K, N_i, S, T)

    for i in 1:(N_i+1)
    
        K[i] = 2.03 + (0.13*S[i])/T[i]

    end

    return nothing
end

# Older, more complicated function for thermal conductivity
@inline function generate_K_old(K, K_b, N_i, θ, T_mltS)

    for i in 1:(N_i+1)
    
        K_b[i] = 1.162 * (0.45  + 0.0108*θ[i]  + 0.0000504*(θ[i]^2))
        K[i]   = 1.162 * (1.905 + 0.00866*θ[i] + 0.0000297*(θ[i]^2))
        K[i]  *= (2*K[i] + K_a - (2*V_a)*(K[i]-K_a)) / (2*K[i] + K_a - (2*V_a)*(K[i]-K_a))
        K[i]  -= (K[i]-K_b[i])*((T_mltS-T_mlt0)/θ[i])

    end

    return nothing
end

# Gets the penetrating solar radiation for this column of sea ice
# NOTE: since we don't update this at each time step, we're doing just
# initial depth of each layer
@inline function generate_I_pen(I_pen, I_0, κ_i, H, N_i)

    for k in 0:N_i
        z      = k / N_i
        I_pen[k+1] = I_0*exp(-κ_i*z*H)
    end

    return nothing
end

# Gets the salinity profile for this column of sea ice
@inline function generate_S(S, N_i)

    for k in 0:N_i
        z      = k / N_i
        S[k+1] = 0.5 * S_max * (1 - cos(pi*z^(0.407/(z+0.573))))
    end

    return nothing
end

# gets the specific heat of sea ice
@inline function generate_c_i(c_i, N_i, c_0, S, T_old, T_new)

    for k in 1:N_i+1
        c_i[k] = c_0 + (L_0*μ*S[k])/(T_old[k]*T_new[k])
    end

    return nothing

end

# Computes the (constant) atmospheric flux affecting the model
function compute_surface_flux(jcmodel, atmodel)

    # Compute atmospheric fluxes dependent on ice:
    set_atm_helper_values(atmodel, jcmodel.T_0[1], jcmodel.H, 5)
    set_atm_flux_values(atmodel, jcmodel.T_0[1])
    set_atm_dflux_values(atmodel, jcmodel.T_0[1])

    # Reduce shortwave flux with albedo
    jcmodel.α = 0.7 # CHANGE TO NOT PRESET later

    # Now compute total surface flux:
    jcmodel.F_0[1] = (1-jcmodel.α)*jcmodel.i_0*atmodel.F_sw + atmodel.F_Ld + atmodel.F_Lu + atmodel.F_l + atmodel.F_s
    # And now compute derivative of flux:
    jcmodel.dF_0[1] = atmodel.dF_Lu + atmodel.dF_s + atmodel.dF_l

    # Apply to all time values since we (for now) assume constant flux:
    for k in 2:jcmodel.N_t

        jcmodel.F_0[k]  = jcmodel.F_0[1]
        jcmodel.dF_0[k] = jcmodel.dF_0[1]
    end

    return nothing
end

# Computes the enthalpy q from the temerature T
@inline function generate_q_i_from_T(q_i, N_i, T, S)

    for k in 1:N_i+1
        T_m  = μ*S[k]
        q_i[k] = -ρ_i*(c_0*(T_m - T[k]) + L_0*(1 - T_m/T[k]) - c_w*T_m)
    end

    return nothing
end

# Computes the temperature T from the enthalpy q
@inline function generate_T_from_q_i(T, N_i, q_i, S)

    for k in 1:N_i+1
        T_m  = μ*S[k]
        b    = (c_w-c_0)*T_m - q_i[k]/ρ_i - L_0
        T[k] = (-b - sqrt(b^2 - 4c_0*L_0*T_m))/(2c_0)
    end

end

# Computes the change in thickness of the bottom ice layer due to growth / melting
# This function is based on CICE documentation
@inline function bottom_ice_change(N_i, Δh, T_n, q_i, K, T_frz, Δt, u_star, T_w)

    # Get conductive heat flux at bottom surface:
    F_cb = (T_n[N_i+1] - T_frz) * K[N_i+1]/Δh[N_i+1]

    # Get net downward heat flux:
    F_bot = -ρ_w*c_w*c_h*u_star*(T_w - T_frz)

    # Now calculate the change in ice thickness in the bottom layer:
    δz = (F_cb - F_bot)*Δt / q_i[N_i+1]

    # And finally update the bottom later thickness:
    Δh[N_i+1] -= δz

    return nothing
end

# Redistributes the ice thicknesses so they are uniform, then updates the
# enthalpy q_i to match.
@inline function rebalance_ice_layers(Δh, q_i, N_i, q_new, z_old, z_new)

    # Get new total thickness
    H_new = 0.0
    for k in 1:(N_i+1)
        H_new += Δh[k]
    end

    # Split evenly
    Δh_new = H_new / N_i

    # Update enthalpies for each layer:
    for k in 2:(N_i+1)
        q_new[k] = 0.0
        for m in 2:(N_i+1)
        
            η_km      = min(z_old[m], z_new[k]) - max(z_old[m-1], z_new[k-1])
            q_new[k] += η_km*q_i[m]
        end
        q_new[k] /= Δh_new
    end

    # Now that we have q_new and Δh_new, we can replace current q_i and Δh with these:
    for k in 2:(N_i+1)
        q_i[k] = q_new[k]
        Δh[k]  = Δh_new
    end

end