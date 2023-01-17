# Written by Joseph Kump (josek97@utexas.edu)
# Code that runs the time step for a single grid cell object in
# our sea ice model, first running the model for each column, then
# applying horizontal transport amongst columns

include("../column_physics/run_ice_column.jl")
include("./jicecell_struct.jl")
include("linear_itd_change.jl")

#= Runs the model of a JICE Cell object

=#
@inline function run_ice_cell(jcell)

    for step in 1:jcell.N_t
        run_cell_step(jcell, step)
    end

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


    # Add new ice growing in ocean run_cell_step
    # This is done with add_new_ice in icepack_therm_itd


    # Melt ice laterally
    # This is done with lateral_melt in icepack_therm_itd
    

end
