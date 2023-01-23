# Written by Joseph Kump (josek97@utexas.edu)
# Code that computes the addition of new ice into a grid cell after remapping

include("./jicecell_struct.jl")

#= Based on function of the same name in icepack_therm_itd, line 1284
=#
@inline function add_new_ice(jcell)

    # Initialize new variables (line 1493)
    h_surp  = 0.0
    hi0_new = 0.0
    ai0_new = 0.0
    v_surp  = 0.0

    #=
    d_an_latg(:) = 0.0
    area2(:) = aicen(:)
    
    lead_area    = 0.0
    latsurf_area = 0.0
    G_radial     = 0.0
    tot_latg     = 0.0
    =#

    hi0_max = 1000000.0
    if jcell.N_cat > 1
        hi0_max = 0.9jcell.H_bds[2]
    end

    # conditional for fsd

    # We will reuse old area and ice volume fields here since we are done with itd shift (line 1526):
    for n in 1:jcell.N_cat

        jcell.areas_old[n] = jcell.areas[n]
        jcell.vol_i_old[n] = jcell.vol_i[n]
        jcell.vol_s_old[n] = jcell.vol_s[n]

        jcell.i_energy_old[n] = 0.0
        jcell.i_energy[n]     = 0.0

        jcell.d_area_i_new[n]  = 0.0
        jcell.d_area_total[n]  = 0.0
        jcell.vi0_new_array[n] = 0.0
    end

    # Conditional for conservation of energy (line 1537)
    for n in 1:jcell.N_cat

        # Get ice energy:
        jcolumn = jcell.columns[n]
        w_i     = (jcolumn.H_i * jcell.areas[n]) / jcolumn.N_i
        for k in (jcolumn.N_s+2):(jcolumn.N_s+jcolumn.N_i+1)
            jcell.i_energy_old[n] += jcolumn.q[k] * w_i
        end
    end

    # Compute average enthalpy of new ice (line 1553), using B-L thermodynamics of course
    # This assumes salinity profile is same for all columns, so we pick the first column arbitrarily
    #salinity_prof = jcell.columns[1].S
    qi0_new       = -ρ_i*L_0

    # Compute the volume, area, and thickness of new ice (line 1579)
    # FOR NOW we're assuming freezing/melting potential is 0
    f_new    = max(jcell.frzmlt, 0.0)
    vi0_new  = -f_new * jcell.columns[1].Δt / qi0_new
    #vi0_init = vi0_new

    # increment ice volume and energy (line 1586)
    jcell.vol_i_old[1]    += vi0_new
    jcell.i_energy_old[1] += vi0_new * qi0_new

    # Update freshwater and salt fluxes (line 1601)
    # For now it doesnt look like we do anything, unless we want to update ocen flux (nothing else for BL thermodynamics)

    # Decide how to distribute new ice (line 1633)
    if vi0_new > 0.0

        # Note there are extra conditionals for fsd
        area_water = 1.0
        for n in 1:jcell.N_cat
            area_water -= jcell.areas[n]
        end

        # new ice area and thickness
        # H_bds[1] < new ice thickness < H_bds[2]
        if area_water > puny
            hi0_new = max(vi0_new/area_water, hfrazilmin)
            # Distribute new volume over all categories if necessary
            if (hi0_new > hi0_max) && (area_water + puny < 1.0)
                hi0_new = hi0_max
                ai0_new = area_water
                v_surp  = vi0_new - ai0_new * hi0_new
                h_surp  = v_surp / (1.0-area_water)
                vi0_new = ai0_new * hi0_new
            else # Otherwise just in category 1
                ai0_new = vi0_new / hi0_new
            end
        else
            h_surp  = vi0_new / (1.0-area_water)
            vi0_new = 0.0
        end

        # Volume added to each category from lateral growth (line 1676), currently ignored since
        # it's only for fsd

        # combine area changes from new ice growth and lateral growth (lateral growth is fsd only):
        jcell.d_area_i_new[1]  = ai0_new
        jcell.d_area_total[1]  = jcell.d_area_i_new[1]
        jcell.vi0_new_array[1] = vi0new

    end

    # Distribute excess ice volume among categories by increasing thickness (line 1694)
    if h_surp > 0.0

        for n in 1:jcell.N_cat

            jcolumn = jcell.columns[n]

            v_surp = h_surp * jcell.areas[n]
            v_tmp  = jcolumn.H_i * jcell.areas[n] + v_surp

            # number of conditionals for things we don't trace

            # update volume:
            jcell.vol_i[n] = v_tmp
            v_tmp         -= v_surp # now temp is old volume

            # We get new ice thickness
            jcolumn     = jcell.columns[n]
            jcolumn.H_i = jcell.vol_i[n] / jcell.areas[n]
            new_Δh      = jcolumn.H_i / jcolumn.N_i
            for k in (jcolumn.N_s+2):(jcolumn.N_i+jcolumn.N_s+1)
                jcolumn.Δh[k] = new_Δh
                # Now we add corresponding enthalpy into new ice growth
                if jcell.vol_i[n] > 0.0
                    jcolumn.q[k] = (jcolumn.q[k]*v_tmp + qi0_new*v_surp) / jcell.vol_i[n]
                end
            end

            # And then get the new temperatures from these enthalpies:
            generate_T_from_q(jcolumn.T_n, jcolumn.N_i, jcolumn.N_s, jcolumn.H_s, jcolumn.q, jcolumn.S)
            # TODO: consider modifying salinity(?)
        end
    end

    # Combine new ice into cat 1 since we don't track FSD yet (line 1792)
    # NOTE: this will be more involved with FSD factored
    N_cats_to_add_ice = 1

    for n in 1:N_cats_to_add_ice
        if (jcell.d_area_total[n] > 0.0) && (jcell.vi0_new_array[n] > 0.0)
            
            jcell.areas[n] += jcell.d_area_total[n]
            jcell.vol_i[n] += jcell.vi0_new_array[n]

            # TODO: trace change in surface temperature(?)

            # a lot of conditionals not yet implemented, next is line 1893

            # We get new ice thickness
            jcolumn     = jcell.columns[n]
            jcolumn.H_i = jcell.vol_i[n] / jcell.areas[n]
            new_Δh      = jcolumn.H_i / jcolumn.N_i
            for k in (jcolumn.N_s+2):(jcolumn.N_i+jcolumn.N_s+1)
                jcolumn.Δh[k] = new_Δh
                # Now we add corresponding enthalpy into new ice growth
                if jcell.vol_i[n] > 0.0
                    jcolumn.q[k] = (jcolumn.q[k]*v_tmp + qi0_new*v_surp) / jcell.vol_i[n]
                end
            end

            # And then get the new temperatures from these enthalpies:
            generate_T_from_q(jcolumn.T_n, jcolumn.N_i, jcolumn.N_s, jcolumn.H_s, jcolumn.q, jcolumn.S)

        end
    end

    # Conservation check (line 1910)

    # Biogeochemistry (ignored for now)

    return nothing
end