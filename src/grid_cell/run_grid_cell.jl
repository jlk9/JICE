# Written by Joseph Kump (josek97@utexas.edu)
# Code that runs the time step for a single grid cell object in
# our sea ice model, first running the model for each column, then
# applying horizontal transport amongst columns

include("../column_physics/run_ice_inline.jl")
include("./jicecell_struct.jl")

#= Runs the model of a JICE Cell object

=#
function run_ice_cell(jice_cell)


end

#= Runs one step of the model for a JICE Cell object

=#
function run_cell_step(jice_cell, step)

    # Run one step for each column
    for k in length(jice_cell.areas)

        jcolumn = jice_cell.columns[k]

        run_column_step(jcolumn.N_i, jcolumn.N_t, jcolumn.H, jcolumn.T_frz, jcolumn.i_0, jcolumn.κ_i, jcolumn.Δt,
                        jcolumn.u_star, jcolumn.T_w, jcolumn.α, jcolumn.F_0, jcolumn.dF_0,
                        jcolumn.Δh, jcolumn.Δh̄, jcolumn.S, jcolumn.c_i, jcolumn.K, jcolumn.K̄, jcolumn.I_pen, jcolumn.q_i,
                        jcolumn.q_inew, jcolumn.z_old, jcolumn.z_new, jcolumn.maindiag,
                        jcolumn.subdiag, jcolumn.supdiag, jcolumn.F_Lu, jcolumn.F_s, jcolumn.F_l, jcolumn.dF_Lu, jcolumn.dF_s, jcolumn.dF_l,
                        jice_cell.atm.F_sw, jice_cell.atm.F_Ld, jice_cell.atm.T_a, jice_cell.atm.Θ_a, jice_cell.atm.ρ_a, jice_cell.atm.Q_a, jice_cell.atm.c_p,
                        jice_cell.atm.c_u, jice_cell.atm.c_Θ, jice_cell.atm.c_q, jice_cell.atm.atm_u_star, jice_cell.atm.U_a,
                        jcolumn.T_n, jcolumn.T_nplus, step)

        # Update T_n and store current temps and thicknesses:
        jcolumn.T_n[:] = jcolumn.T_nplus
        jcolumn.H      = sum(jcolumn.Δh)

        jcolumn.T_array[:, step+1] = jcolumn.T_n
        jcolumn.Δh_array[:,step+1] = jcolumn.Δh

    end

    # Apply horizontal transport between columns
    

end
