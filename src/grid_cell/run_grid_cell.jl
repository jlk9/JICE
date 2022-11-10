# Written by Joseph Kump (josek97@utexas.edu)
# Code that runs the time step for a single grid cell object in
# our sea ice model, first running the model for each column, then
# applying horizontal transport amongst columns

include("../column_physics/run_ice_column.jl")
include("./jicecell_struct.jl")

#= Runs the model of a JICE Cell object

=#
function run_ice_cell(jcell)

    for step in 1:jcell.N_t
        run_cell_step(jcell, step)
    end

end

#= Runs one step of the model for a JICE Cell object

=#
function run_cell_step(jcell, step)

    # Run one step for each column
    for k in 1:length(jcell.areas)

        jcolumn = jcell.columns[k]

        run_column_step(jcolumn.N_i, jcolumn.N_t, jcolumn.H_i, jcolumn.T_frz, jcolumn.i_0, jcolumn.κ_i, jcolumn.Δt,
                        jcolumn.u_star, jcolumn.T_w, jcolumn.α_vdr, jcolumn.F_0, jcolumn.dF_0,
                        jcolumn.Δh, jcolumn.Δh̄, jcolumn.S, jcolumn.c_i, jcolumn.K, jcolumn.K̄, jcolumn.I_pen, jcolumn.q_i,
                        jcolumn.q_inew, jcolumn.z_old, jcolumn.z_new, jcolumn.maindiag,
                        jcolumn.subdiag, jcolumn.supdiag, jcolumn.F_Lu, jcolumn.F_s, jcolumn.F_l, jcolumn.dF_Lu, jcolumn.dF_s, jcolumn.dF_l,
                        jcell.atm.F_SWvdr, jcell.atm.F_Ld, jcell.atm.T_a, jcell.atm.Θ_a, jcell.atm.ρ_a, jcell.atm.Q_a, jcell.atm.c_p,
                        jcell.atm.c_u, jcell.atm.c_Θ, jcell.atm.c_q, jcell.atm.atm_u_star, jcell.atm.U_a,
                        jcolumn.T_n, jcolumn.T_nplus, step)

        # Update T_n and store current temps and thicknesses:
        jcolumn.T_n[:] = jcolumn.T_nplus
        jcolumn.H_i    = sum(jcolumn.Δh)

        jcolumn.T_array[:, step+1] = jcolumn.T_n
        jcolumn.Δh_array[:,step+1] = jcolumn.Δh

    end

    # Apply horizontal transport between columns
    

end
