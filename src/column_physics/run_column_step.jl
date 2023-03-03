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
                                step)


    # Redo of column steps to more closely follow step_therm1
    # First there are a bunch of optional tags

    # Then we need to adjust the frzmlt for ice-ocean heat fluxes (frzmlt_bottom_lateral, line 2535)
    #frzmlt_bottom_lateral(Δt, N_i, N_s, H_i, H_s, T_w, T_frz, frzmlt, area, u_star_min, rside, fside, f_bot)

    # Then we might need to set a number of helper arrays to 0 (line 2550)


    # Then calculate atmosphere boundary since we normally calculate T_sfc (line 2582):


    # Calculate thermo vertical (line 2655)


    # Increment area-weighted fluxes (line 2827)



    
    # Computes the current albedo
    generate_α(H_i[1], α_vdr, α_idr, α_vdf, α_idf, T_n[1])
    
    # Computes the surface fluxes at this time step
    step_surface_flux(N_i, α_vdr, α_idr, α_vdf, α_idf, T_n[1], H_i[1], H_s[1],
                         F_0, dF_0, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l,
                         F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, I_pen,
                         c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star,
                         ρ_a, c_p, step)
    
    # Computes the temperature changes at this step
    step_temp_change(N_i, N_s, H_s[1], S, T_frz, Δh, T_n,
                     T_nplus, c_i, K, K̄, I_pen, F_0[step], dF_0[step],
                     maindiag, subdiag, supdiag, Δt)

    # Gets the growth/melt and rebalances
    step_growth_melt(N_i, N_s, H_s[1], S, T_frz, Δh, T_nplus, K, K̄,
                      q, q_new, z_old, z_new, Δt, u_star, T_w, F_0[step])

    # Add up the layer thicknesses to get the new total thickness
    readd_total_thickness(N_i, N_s, Δh, H_i, H_iold, H_s)

end

# Julia's sum() operation is memory inefficient for slices of arrays, this is oddly much faster
function readd_total_thickness(N_i, N_s, Δh, H_i, H_iold, H_s)

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