# Written by Joseph Kump (josek97@utexas.edu)
# Basic code that is just ice temperature in column, with no formulation or melting.
# This version eliminates vectorized operations, enabling use of Enzyme.jl
# It also reduces the number of memory allocations substantially, increasing runtime
# and lowering memory footprint.

include("./jicecolumn_struct.jl")
include("../atmosphere_model/atmodel_struct.jl")
include("./surface_flux.jl")
include("./temp_change.jl")
include("./growth_melt.jl")

#= Model Function
Input:
    jcolumn         an initialized JICEColumn object
    atmodel         an initialized ATModel object for surface flux computations

Output:
    Technically nothing, but the function updates jcolumn's T_array and Δh_arrays to have a log
    of updated ice temperatures and thicknesses
=#
@inline function run_ice_column(jcolumn, atmodel)

    # Main loop of temperature modifications:
    for step in 1:jcolumn.N_t

        #= # Old, unused version of time step
        run_column_step_old(jcolumn.N_i, jcolumn.N_s, jcolumn.N_t, jcolumn.H_i, jcolumn.H_s, jcolumn.T_frz,
                        jcolumn.Δt, jcolumn.u_star, jcolumn.T_w,
                        jcolumn.α_vdr, jcolumn.α_idr, jcolumn.α_vdf, jcolumn.α_idf, jcolumn.F_0, jcolumn.dF_0,
                        jcolumn.Δh, jcolumn.Δh̄, jcolumn.S, jcolumn.c_i, jcolumn.K, jcolumn.K̄, jcolumn.I_pen, jcolumn.q_i,
                        jcolumn.q_inew, jcolumn.z_old, jcolumn.z_new, jcolumn.maindiag,
                        jcolumn.subdiag, jcolumn.supdiag, jcolumn.F_Lu, jcolumn.F_s, jcolumn.F_l, jcolumn.dF_Lu, jcolumn.dF_s, jcolumn.dF_l,
                        atmodel.F_SWvdr, atmodel.F_SWidr, atmodel.F_SWvdf, atmodel.F_SWidf,
                        atmodel.F_Ld, atmodel.T_a, atmodel.Θ_a, atmodel.ρ_a, atmodel.Q_a, atmodel.c_p,
                        atmodel.c_u, atmodel.c_Θ, atmodel.c_q, atmodel.atm_u_star, atmodel.U_a,
                        jcolumn.T_n, jcolumn.T_nplus, step)
        =#

        run_column_step(jcolumn, atmodel, step)
        
        # Update T_n and store current temps and thicknesses:
        jcolumn.T_n[:] = jcolumn.T_nplus
        

        jcolumn.T_array[:, step+1] = jcolumn.T_n
        jcolumn.Δh_array[:,step+1] = jcolumn.Δh
    end

    return nothing
end

# Runs one step of ice process. This approach uses the struct within each step:
@inline function run_column_step(jcolumn, atmodel, step)
    
    # Computes the current albedo
    generate_α(jcolumn)
    
    # Computes the surface fluxes at this time step
    step_surface_flux(jcolumn, atmodel, step)
    
    # Computes the temperature changes at this step
    step_temp_change(jcolumn, step)

    # Gets the growth/melt and rebalances
    step_growth_melt(jcolumn, step)

    # Add up the layer thicknesses to get the new total thickness
    readd_total_thickness(jcolumn)
end

# Julia's sum() operation is memory inefficient for slices of arrays, this is much faster
function readd_total_thickness(jcolumn)

    jcolumn.H_iold = jcolumn.H_i

    jcolumn.H_s = 0.0
    jcolumn.H_i = 0.0

    for k in 1:(jcolumn.N_s+1)
        jcolumn.H_s += jcolumn.Δh[k]
    end
    for k in (jcolumn.N_s+2):(jcolumn.N_i+jcolumn.N_s+1)
        jcolumn.H_i += jcolumn.Δh[k]
    end
end

#=
# Runs one step of ice process. This field-by-field approach was meant for old AD functions
@inline function run_column_step_old(N_i, N_s, N_t, H_i, H_s, T_frz, Δt, u_star, T_w,
                                α_vdr, α_idr, α_vdf, α_idf, F_0, dF_0,
                                Δh, Δh̄, S, c_i, K, K̄, I_pen, q_i, q_inew, z_old, z_new, maindiag,
                                subdiag, supdiag, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l,
                                F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p,
                                c_u, c_Θ, c_q, atm_u_star, U_a,
                                T_n, T_nplus, step)

    # Computes the current albedo
    generate_α(H_i, α_vdr, α_idr, α_vdf, α_idf, T_n[1])

    # Computes the surface fluxes at this time step
    step_surface_flux(N_i, α_vdr, α_idr, α_vdf, α_idf, T_n[1], H_i, H_s, F_0, dF_0, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l,
                        F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, I_pen, c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, ρ_a, c_p, step)
    
    # Computes the temperature changes at this step
    step_temp_change(N_i, N_s, S, T_frz, Δh, Δh̄, T_n, T_nplus, c_i,
                    K, K̄, I_pen, F_0[step], dF_0[step], maindiag,
                    subdiag, supdiag, Δt)

    # Gets the growth/melt and rebalances
    step_growth_melt(N_i, N_s, S, T_frz, Δh, T_nplus, K, q_i, q_inew,
                    z_old, z_new, Δt, u_star, T_w)

end
=#