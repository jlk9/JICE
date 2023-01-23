# Written by Joseph Kump (josek97@utexas.edu)
# Code that runs the time step for a single grid cell object in
# our sea ice model, first running the model for each column, then
# applying horizontal transport amongst columns

include("../column_physics/run_ice_column.jl")
include("./jicecell_struct.jl")
include("./linear_itd_change.jl")
include("./add_new_ice.jl")

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
    for k in 1:jcell.N_cat

        jcolumn = jcell.columns[k]

        run_column_step(jcolumn, jcell.atm, step)

        # Update T_n and store current temps and thicknesses:
        jcolumn.T_n[:] = jcolumn.T_nplus

        jcolumn.T_array[:, step+1] = jcolumn.T_n
        jcolumn.Δh_array[:,step+1] = jcolumn.Δh

    end

    # Apply horizontal transport between columns vis a linear map
    # This is done with function linear_itd in icepack_therm_itd
    linear_itd_change(jcell)

    # Check that conservation of energy is conserved
    conservation_check_itd(jcell)

    # Add new ice growing in ocean run_cell_step
    # This is done with add_new_ice in icepack_therm_itd
    add_new_ice(jcell)
    conservation_check_new_ice(jcell)

    # Melt ice laterally
    # This is done with lateral_melt in icepack_therm_itd
    
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
        println("Oh no! The total change in ice volume from horizontal transport is above acceptable machine precision.")
    end
    if s_vol_change > puny
        println("Oh no! The total change in snow volume from horizontal transport is above acceptable machine precision.")
    end
    if i_energy_change > puny
        println("Oh no! The total change in ice energy from horizontal transport is above acceptable machine precision.")
    end
    if s_energy_change > puny
        println("Oh no! The total change in snow energy from horizontal transport is above acceptable machine precision.")
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
        println("Oh no! The total change in ice volume from adding new ice (when factoring in added component) is above acceptable machine precision.")
    end
    if i_energy_change > puny
        println("Oh no! The total change in ice energy from adding new ice (when factoring in added component) is above acceptable machine precision.")
    end

    return nothing
end
