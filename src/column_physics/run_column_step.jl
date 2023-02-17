# Written by Joseph Kump (josek97@utexas.edu)
# Basic code that is just ice temperature in column, with no formulation or melting.
# This version eliminates vectorized operations, enabling use of Enzyme.jl
# It also reduces the number of memory allocations substantially, increasing runtime
# and lowering memory footprint.

include("./surface_flux.jl")
include("./temp_change.jl")
include("./growth_melt.jl")

# Runs one step of ice process. 
@inline function run_column_step(N_i, N_s,
                                N_layers, N_t,
                                T_frz, Δt, u_star,
                                T_w, T_n,
                                H_i, H_iold, H_s,
                                α_vdr, α_idr,
                                α_vdf, α_idf,
                                T_nplus,
                                F_0, dF_0,
                                Δh, S, c_i,
                                K, K̄,
                                I_pen,
                                q, q_new,
                                z_old, z_new,
                                maindiag, subdiag, supdiag,
                                F_Lu, F_s, F_l,
                                dF_Lu, dF_s, dF_l,
                                F_SWvdr, F_SWidr,
                                F_SWvdf, F_SWidf,
                                F_Ld,
                                Θ_a, ρ_a,
                                Q_a, c_p,
                                c_u, c_Θ, c_q,
                                atm_u_star, U_a,
                                step, index)
    
    # Computes the current albedo
    generate_α(H_i[1], α_vdr, α_idr, α_vdf, α_idf, T_n[N_layers*index+1], index)
    
    # Computes the surface fluxes at this time step
    step_surface_flux(N_i, α_vdr, α_idr, α_vdf, α_idf, T_n[N_layers*index+1], H_i[1], H_s[1],
                         F_0, dF_0, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l,
                         F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, I_pen,
                         c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star,
                         ρ_a, c_p, step, index)
    
    # Computes the temperature changes at this step
    step_temp_change(N_i, N_s, H_s[1], S, T_frz, Δh, T_n,
                     T_nplus, c_i, K, K̄, I_pen, F_0[N_t*index+step], dF_0[N_t*index+step],
                     maindiag, subdiag, supdiag, Δt, index)

    # Gets the growth/melt and rebalances
    step_growth_melt(N_i, N_s, H_s[1], S, T_frz, Δh, T_nplus, K, K̄,
                      q, q_new, z_old, z_new, Δt, u_star, T_w, F_0[N_t*index+step], index)

    # Add up the layer thicknesses to get the new total thickness
    readd_total_thickness(N_i, N_s, Δh, H_i, H_iold, H_s, step, index)

end

# Julia's sum() operation is memory inefficient for slices of arrays, this is oddly much faster
function readd_total_thickness(N_i, N_s, Δh, H_i, H_iold, H_s, step, index)

    H_iold[1] = H_i[1]
    H_i[1]    = 0.0
    H_s[1]    = 0.0

    for k in 1:(N_s+1)
        H_s[1] += Δh[k]
    end
    for k in (N_s+2):(N_i+N_s+1)
        H_i[1] += Δh[k]
    end
end