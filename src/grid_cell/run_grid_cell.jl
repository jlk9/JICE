# Written by Joseph Kump (josek97@utexas.edu)
# Code that runs the time step for a single grid cell object in
# our sea ice model, first running the model for each column, then
# applying horizontal transport amongst columns

include("../column_physics/run_column_step.jl")
include("./jicecell_struct.jl")
include("./linear_itd_change.jl")
include("./add_new_ice.jl")
include("./lateral_melt.jl")

#= Runs the model of a JICE Cell object

=#
@inline function run_ice_cell(jcell)

    for step in 1:jcell.N_t
        run_cell_step(jcell, step)
    end

    return nothing
end

#= Runs one step of the model for a JICE Cell object

=#
@inline function run_cell_step(jcell, step)

    # Run one step for each column
    for n in 1:jcell.N_cat

        jcolumn = jcell.columns[n]

        run_column_step(jcell.N_i, jcell.N_s,
                        jcell.N_i+jcell.N_s+1, jcell.N_t,
                        jcell.T_frz, jcell.Δt, jcell.u_star,
                        jcell.T_w, jcolumn.T_n,
                        jcolumn.H_i_array, jcolumn.H_s_array,
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
                        jcell.atm.F_SWvdr, jcell.atm.F_SWidr,
                        jcell.atm.F_SWvdf, jcell.atm.F_SWidf,
                        jcell.atm.F_Ld,
                        jcell.atm.Θ_a, jcell.atm.ρ_a,
                        jcell.atm.Q_a, jcell.atm.c_p,
                        jcell.atm.c_u, jcell.atm.c_Θ, jcell.atm.c_q,
                        jcell.atm.atm_u_star, jcell.atm.U_a,
                        step, 0)

        # Update T_n:
        jcolumn.T_n[:] = jcolumn.T_nplus
    end
    
    # Apply horizontal transport between columns vis a linear map
    # This is done with function linear_itd in icepack_therm_itd
    # Also check that conservation of energy is conserved
    linear_itd_change(jcell, step)
    conservation_check_itd(jcell)

    # Add new ice growing in ocean run_cell_step
    # This is done with add_new_ice in icepack_therm_itd
    add_new_ice(jcell, step)
    conservation_check_new_ice(jcell)

    # Melt ice laterally
    # This is done with lateral_melt in icepack_therm_itd
    lateral_melt(jcell)
    
    # Store current temps and thicknesses:
    for n in 1:jcell.N_cat

        jcolumn = jcell.columns[n]

        jcolumn.T_array[:,step+1]  = jcolumn.T_n
        jcolumn.Δh_array[:,step+1] = jcolumn.Δh
    end
    
    return nothing
end

#= Checks that ice and snow volume and energy are conserved through linear remapping and other processes

=#
@inline function conservation_check_itd(jcell)

    total_i_vol_old    = 0.0
    total_s_vol_old    = 0.0
    total_i_energy_old = 0.0
    total_s_energy_old = 0.0

    total_i_vol    = 0.0
    total_s_vol    = 0.0
    total_i_energy = 0.0
    total_s_energy = 0.0

    for n in 1:jcell.N_cat
        total_i_vol_old    += jcell.vol_i_old[n]
        total_s_vol_old    += jcell.vol_s_old[n]
        total_i_energy_old += jcell.i_energy_old[n]
        total_s_energy_old += jcell.s_energy_old[n]

        total_i_vol    += jcell.vol_i[n]
        total_s_vol    += jcell.vol_s[n]
        total_i_energy += jcell.i_energy[n]
        total_s_energy += jcell.s_energy[n]
    end

    i_vol_change    = abs(total_i_vol_old - total_i_vol) / abs(total_i_vol_old)
    s_vol_change    = abs(total_s_vol_old - total_s_vol) / abs(total_s_vol_old)

    i_energy_change = abs(total_i_energy_old - total_i_energy) / abs(total_i_energy_old)
    s_energy_change = abs(total_s_energy_old - total_s_energy) / abs(total_s_energy_old)

    if i_vol_change > puny
        jcell.i_vol_itd_change = true
    end
    if s_vol_change > puny
        jcell.s_vol_itd_change = true
    end
    if i_energy_change > puny
        jcell.i_energy_itd_change = true
    end
    if s_energy_change > puny
        jcell.s_energy_itd_change = true
    end

    return nothing
end


#= Checks that ice volume and energy are conserved through adding new ice

=#
@inline function conservation_check_new_ice(jcell)

    total_i_vol_old    = 0.0
    total_i_energy_old = 0.0

    total_i_vol    = 0.0
    total_i_energy = 0.0

    for n in 1:jcell.N_cat
        total_i_vol_old    += jcell.vol_i_old[n]
        total_i_energy_old += jcell.i_energy_old[n]

        total_i_vol    += jcell.vol_i[n]
        total_i_energy += jcell.i_energy[n]
    end

    i_vol_change    = abs(total_i_vol_old - total_i_vol) / abs(total_i_vol_old)
    i_energy_change = abs(total_i_energy_old - total_i_energy) / abs(total_i_energy_old)

    if i_vol_change > puny
        jcell.i_vol_new_change = true
    end
    if i_energy_change > puny
        jcell.i_energy_new_change = true
    end

    return nothing
end
