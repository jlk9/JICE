# Written by Joseph Kump (josek97@utexas.edu)
# Basic code that is just ice temperature in column, with no formulation or melting.
# This version eliminates vectorized operations, enabling use of Enzyme.jl
# It also reduces the number of memory allocations substantially, increasing runtime
# and lowering memory footprint.

include("./jicecolumn_struct.jl")
include("./surface_flux.jl")
include("./temp_change.jl")
include("./growth_melt.jl")

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

    jcmodel.α = 0.7

    jcmodel.T_array[:, 1] = T_n
    jcmodel.Δh_array[:,1] = jcmodel.Δh

    # Main loop of temperature modifications:
    for step in 1:jcmodel.N_t
        
        run_ice_step(jcmodel.N_i, jcmodel.N_t, jcmodel.H, jcmodel.T_frz, jcmodel.i_0, jcmodel.κ_i, jcmodel.Δt,
                    jcmodel.u_star, jcmodel.T_w, jcmodel.α, jcmodel.F_0, jcmodel.dF_0,
                    jcmodel.Δh, jcmodel.Δh̄, jcmodel.S, jcmodel.c_i, jcmodel.K, jcmodel.K̄, jcmodel.I_pen, jcmodel.q_i,
                    jcmodel.q_inew, jcmodel.z_old, jcmodel.z_new, jcmodel.maindiag,
                    jcmodel.subdiag, jcmodel.supdiag, jcmodel.F_Lu, jcmodel.F_s, jcmodel.F_l, jcmodel.dF_Lu, jcmodel.dF_s, jcmodel.dF_l,
                    atmodel.F_sw, atmodel.F_Ld, atmodel.T_a, atmodel.Θ_a, atmodel.ρ_a, atmodel.Q_a, atmodel.c_p,
                    atmodel.c_u, atmodel.c_Θ, atmodel.c_q, atmodel.atm_u_star, atmodel.U_a,
                    T_n, T_nplus, step)
        
        # Update T_n
        T_n[:] = T_nplus

        jcmodel.T_array[:, step+1] = T_n
        jcmodel.Δh_array[:,step+1] = jcmodel.Δh

    end

    return nothing
end

# Runs one step of ice process:
@inline function run_ice_step(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w, α, F_0, dF_0,
                                Δh, Δh̄, S, c_i, K, K̄, I_pen, q_i, q_inew, z_old, z_new, maindiag,
                                subdiag, supdiag, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l,
                                F_sw, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p,
                                c_u, c_Θ, c_q, atm_u_star, U_a,
                                T_n, T_nplus, step)

    # Computes the surface fluxes at this time step
    step_surface_flux(α, i_0, T_n[1], H, F_0, dF_0, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l,
                        F_sw, F_Ld, c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, ρ_a, c_p, step)
        
    # Gets the penetrating shortwave radiation
    generate_I_pen(I_pen, i_0*(1-α)*F_sw, κ_i, H, N_i)
        
    # Computes the temperature changes at this step
    step_temp_change(N_i, S, T_frz, Δh, Δh̄, T_n, T_nplus, c_i,
                    K, K̄, I_pen, F_0[step], dF_0[step], maindiag,
                    subdiag, supdiag, Δt)

    # Gets the growth/melt and rebalances
    step_growth_melt(N_i, S, T_frz, Δh, T_n, T_nplus, K, q_i, q_inew,
                    z_old, z_new, Δt, u_star, T_w)

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
