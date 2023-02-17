# Written by Joseph Kump (josek97@utexas.edu)
# Basic code that is just ice temperature in column, with no formulation or melting.
# This version eliminates vectorized operations, enabling use of Enzyme.jl
# It also reduces the number of memory allocations substantially, increasing runtime
# and lowering memory footprint.

include("./jicecolumn_struct.jl")
include("../atmosphere_model/atmodel_struct.jl")
include("./run_column_step.jl")

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
        
        run_column_step(jcolumn.N_i, jcolumn.N_s,
                        jcolumn.N_i+jcolumn.N_s+1, jcolumn.N_t,
                        jcolumn.T_frz, jcolumn.Δt, jcolumn.u_star,
                        jcolumn.T_w, jcolumn.T_n,
                        jcolumn.H_i, jcolumn.H_iold, jcolumn.H_s,
                        jcolumn.α_vdr, jcolumn.α_idr,
                        jcolumn.α_vdf, jcolumn.α_idf,
                        jcolumn.T_nplus,
                        jcolumn.F_0, jcolumn.dF_0,
                        jcolumn.Δh, jcolumn.S, jcolumn.c_i,
                        jcolumn.K, jcolumn.K̄,
                        jcolumn.I_pen,
                        jcolumn.q, jcolumn.q_new,
                        jcolumn.z_old, jcolumn.z_new,
                        jcolumn.maindiag, jcolumn.subdiag, jcolumn.supdiag,
                        jcolumn.F_Lu, jcolumn.F_s, jcolumn.F_l,
                        jcolumn.dF_Lu, jcolumn.dF_s, jcolumn.dF_l,
                        atmodel.F_SWvdr, atmodel.F_SWidr,
                        atmodel.F_SWvdf, atmodel.F_SWidf,
                        atmodel.F_Ld,
                        atmodel.Θ_a, atmodel.ρ_a,
                        atmodel.Q_a, atmodel.c_p,
                        atmodel.c_u, atmodel.c_Θ, atmodel.c_q,
                        atmodel.atm_u_star, atmodel.U_a,
                        step, 0)

        # Update T_n and store current temps and thicknesses:
        jcolumn.T_n[:] = jcolumn.T_nplus
        
        jcolumn.H_i_array[step+1] = jcolumn.H_i[1]
        jcolumn.H_s_array[step+1] = jcolumn.H_s[1]

        jcolumn.T_array[:, step+1] = jcolumn.T_n
        jcolumn.Δh_array[:,step+1] = jcolumn.Δh
    end

    return nothing
end

#=
# Runs one step of ice process. This approach uses the struct within each step:
@inline function run_column_step(jcolumn, atmodel, step)
    
    # Computes the current albedo
    generate_α(jcolumn.H_i_array[step], jcolumn.α_vdr, jcolumn.α_idr, jcolumn.α_vdf, jcolumn.α_idf, jcolumn.T_n[1])
    
    # Computes the surface fluxes at this time step
    step_surface_flux(jcolumn.N_i, jcolumn.α_vdr, jcolumn.α_idr, jcolumn.α_vdf, jcolumn.α_idf, jcolumn.T_n[1], jcolumn.H_i_array[step], jcolumn.H_s_array[step],
                         jcolumn.F_0, jcolumn.dF_0, jcolumn.F_Lu, jcolumn.F_s, jcolumn.F_l, jcolumn.dF_Lu, jcolumn.dF_s, jcolumn.dF_l,
                         atmodel.F_SWvdr, atmodel.F_SWidr, atmodel.F_SWvdf, atmodel.F_SWidf, atmodel.F_Ld, jcolumn.I_pen,
                         atmodel.c_u, atmodel.c_Θ, atmodel.c_q, atmodel.U_a, atmodel.Θ_a, atmodel.Q_a, atmodel.atm_u_star,
                         atmodel.ρ_a, atmodel.c_p, step)
    
    # Computes the temperature changes at this step
    step_temp_change(jcolumn.N_i, jcolumn.N_s, jcolumn.H_s_array[step], jcolumn.S, jcolumn.T_frz, jcolumn.Δh, jcolumn.T_n,
                     jcolumn.T_nplus, jcolumn.c_i, jcolumn.K, jcolumn.K̄, jcolumn.I_pen, jcolumn.F_0[step], jcolumn.dF_0[step],
                     jcolumn.maindiag, jcolumn.subdiag, jcolumn.supdiag, jcolumn.Δt)

    # Gets the growth/melt and rebalances
    step_growth_melt(jcolumn.N_i, jcolumn.N_s, jcolumn.H_s_array[step], jcolumn.S, jcolumn.T_frz, jcolumn.Δh, jcolumn.T_nplus, jcolumn.K, jcolumn.K̄,
                      jcolumn.q, jcolumn.q_new, jcolumn.z_old, jcolumn.z_new, jcolumn.Δt, jcolumn.u_star, jcolumn.T_w, jcolumn.F_0[step])

    # Add up the layer thicknesses to get the new total thickness
    readd_total_thickness(jcolumn.N_i, jcolumn.N_s, jcolumn.Δh, jcolumn.H_i_array, jcolumn.H_s_array, step)

end
=#