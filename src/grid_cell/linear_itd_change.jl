# Written by Joseph Kump (josek97@utexas.edu)
# Code that computes the remapping of ice column areas
# between thickness categories within one grid cell

include("./jicecell_struct.jl")
#=
    MAJOR ASSUMPTION: here we assume each thickness category has the same number of ice/snow layers
=#
@inline function linear_itd_change(jcell)

    # Compute energy sums remapping should conserve:
    sum_total_energy(jcell.i_energy_old, jcell.s_energy_old, jcell)

    #println(jcell.i_energy_old)
    #println(jcell.s_energy_old)

    # Compute thickness changes in each category
    for n in 1:jcell.N_cat
        jcell.dH[n] = jcell.columns[n].H_i - jcell.columns[n].H_iold
    end

    # Compute new category boundaries
    jcell.H_bnew[1] = jcell.H_bds[1]
    for n in 1:(jcell.N_cat-1)

        # To clean up indexing
        nplus1 = n+1

        H_ioldn  = jcell.columns[n].H_iold
        H_ioldnp = jcell.columns[nplus1].H_iold

        if (H_ioldn > puny) && (H_ioldnp > puny)

            if H_ioldnp > H_ioldn
                # Interpolate between adjacent categories in this case:
                slope                = (jcell.dH[nplus1] - jcell.dH[n]) / (H_ioldnp - H_ioldn)
                jcell.H_bnew[nplus1] = jcell.H_bds[nplus1] + jcell.dH[n] + slope*(jcell.H_bds[nplus1] - H_ioldn)
            else
                println("Oh no! Category bounds crossed streams!")
                return nothing
            end
        elseif H_ioldn > puny # implies H_ioldnp ~= 0
            jcell.H_bnew[nplus1] = jcell.H_bds[nplus1] + jcell.dH[n]
        elseif H_ioldnp > puny # implies H_ioldn ~= 0
            jcell.H_bnew[nplus1] = jcell.H_bds[nplus1] + jcell.dH[nplus1]
        else
            jcell.H_bnew[nplus1] = jcell.H_bds[nplus1]
        end

        # OPTION: can error check here like CICE does (icepack_therm_itd.F90 lines ~340-400)
    end

    # For last thickness category bound:
    N_catplus1 = jcell.N_cat+1
    if jcell.areas[jcell.N_cat] > puny
        jcell.H_bnew[N_catplus1] = 3.0jcell.columns[jcell.N_cat].H_i - 2.0jcell.H_bnew[jcell.N_cat]
    else
        jcell.H_bnew[N_catplus1] = jcell.H_bds[N_catplus1]
    end
    jcell.H_bnew[N_catplus1] = max(jcell.H_bnew[N_catplus1], jcell.H_bds[jcell.N_cat])

    # Boundary check (each one lies between adjacent values of H_i)
    for n in 1:jcell.N_cat
        if jcell.H_bnew[n] > jcell.H_bds[n+1]
            println("Oh no! One of the new catagory bounds exceeds the next (old) bound up.")
            return nothing
        end
        if jcell.H_bnew[n+1] < jcell.H_bds[n]
            println("Oh no! One of the new catagory bounds is smaller than the previous old bound.")
            return nothing
        end
    end

    # Compute g0, g1, hL, hR for remapping:
    fit_line(jcell)

    # Find area lost due to melting of thin ice (line ~439) (NOTE: we keep track of all areas including open water)
    if jcell.areas[1] > puny

        dh0 = jcell.dH[1]
        if dh0 < 0.0

            dh0 = min(-dh0,jcell.H_bds[2])
            # Integrate g(1) from 0 to dh0 to estimate area melted
            # First the right integration limit:
            η_max = min(dh0, jcell.hR[1]) - jcell.hL[1]

            if η_max > 0.0
                x1  = η_max
                x2  = 0.5η_max^2
                da0 = jcell.g1[1]*x2 + jcell.g0[1]*x1 # total ice area removed

                # Now constrain the new thickness so it is <= H_iold
                damax = jcell.areas[1] * (1.0 - jcell.columns[1].H_i / jcell.columns[1].H_iold)
                da0   = min(da0, damax)

                # Remove this area, while conserving volume
                jcell.columns[1].H_i = jcell.columns[1].H_i * jcell.areas[1] / (jcell.areas[1] - da0)
                jcell.areas[1]      -= da0
            end
        else #dh0 >= 0.0
            jcell.H_bnew[1] = min(dh0, jcell.H_bds[2]) # shift lower boundary to the right, as thin ice grew
        end
    end

    # Compute the area and volume to be shifted across each boundary (line ~493)
    for n in 1:(jcell.N_cat-1)

        nplus1 = n + 1
        nplus2 = n + 2
        η_min  = 0.0
        η_max  = 0.0

        if jcell.H_bnew[nplus1] > jcell.H_bds[nplus1] # transfer from n+1 to n+2
            η_min          = max(jcell.H_bds[nplus1], jcell.hL[nplus1]) - jcell.hL[nplus1]
            η_max          = max(jcell.H_bnew[nplus1], jcell.hR[nplus1]) - jcell.hL[nplus1]
            jcell.donor[n] = n
        else # transfer down a category, η_min = 0
            η_max          = min(jcell.H_bds[nplus1], jcell.hR[nplus2]) - jcell.hL[nplus2]
            jcell.donor[n] = nplus1
        end

        if η_max > η_min
            x1   = η_max - η_min
            wk1  = η_min^2
            wk2  = η_max^2
            x2   = 0.5*(wk2-wk1)
            wk1 *= η_min
            wk2 *= η_max
            x3   = (wk2 - wk1) / 3.0
            nd   = jcell.donor[n] + 1 # increment by 1 in our implementation vs CICE
            
            jcell.dareas[n] = jcell.g1[nd]*x2 + jcell.g0[nd]*x1
            jcell.dvol_i[n] = jcell.g1[nd]*x3 + jcell.g0[nd]*x2 + jcell.dareas[n]*jcell.hL[nd]
        end

        # Shift 0 ice if dareas[n] too small:
        nd = jcell.donor[n]
        if (jcell.dareas[n] < jcell.areas[nd]*puny) || (jcell.dvol_i[n] < jcell.columns[nd].H_i*jcell.areas[n]*puny)
            jcell.dareas[n] = 0.0
            jcell.donor[n]  = 0
        end

        # Shift entire category if dareas[n] is close to areas[n+1]:
        if (jcell.dareas[n] > jcell.areas[nd]*(1.0-puny)) || (jcell.dvol_i[n] > jcell.columns[nd].H_i*jcell.areas[n]*(1.0-puny))
            jcell.dareas[n] = jcell.areas[nd]
        end
    end
    #=
    # Shift ice between categories as necessary (line ~564) NOT SURE WHY CICE DOES THIS
    # Maintain negative definiteness of snow enthalpy:
    for n in 1:jcell.N_cat
        for k in 2:(jcell.columns[n].N_s+1)
            jcell.columns[n].q[k] += ρ_s*L_0
        end
    end
    =#
    shift_ice(jcell)
    #=
    # Maintain negative definiteness of snow enthalpy:
    for n in 1:jcell.N_cat
        for k in 2:(jcell.columns[n].N_s+1)
            jcell.columns[n].q[k] -= ρ_s*L_0
        end
    end
    =#
    # Compute new energy sums remapping should still conserve:
    sum_total_energy(jcell.i_energy, jcell.s_energy, jcell)

    #println(jcell.i_energy)
    #println(jcell.s_energy)

    return nothing
end

#=  Fit g(h) with a line, satisfying area and volume constraints.
    To reduce roundoff errors caused by large values of g0 and g1,
    we actually compute g(eta), where eta = h - hL, and hL is the
    left boundary.

    Compared to version in CICE, here we want to do it in one array
=#
@inline function fit_line(jcell)

    # First for first column (different rules, line ~431 in icepack_therm_itd.F90)
    if (jcell.areas[1] > puny) && (jcell.H_bds[2] - jcell.H_bnew[1] > puny)

        # Initialize hL and hR
        hL   = jcell.H_bnew[1]
        hR   = jcell.H_bds[2]
        hice = jcell.columns[1].H_iold

        # Change hL or hR if hicen(n) falls outside central third of range
        h13 = (2.0hL + hR) / 3.0
        h23 = (hL + 2.0hR) / 3.0
        if (hice < h13)
            hR = 3.0hice - 2.0hL
        elseif (hice > h23)
            hL = 3.0hice - 2.0hR
        end

        # Compute coefficients of g(eta) = g0 + g1*eta
        dhr = 1.0 / (hR - hL)
        wk1 = 6.0jcell.areas[1] * dhr
        wk2 = (hice - hL) * dhr

        jcell.g0[1] = wk1 * (2.0/3.0 - wk2)
        jcell.g1[1] = 2.0dhr * wk1 * (wk2 - 0.5)
        jcell.hL[1] = hL
        jcell.hR[1] = hR
    else
        jcell.g0[1] = 0.0
        jcell.g1[1] = 0.0
        jcell.hL[1] = 0.0
        jcell.hR[1] = 0.0
    end

    # Then for subsequent columns (line ~485 in icepack_therm_itd.F90)
    for n = 1:jcell.N_cat
        nplus1 = n + 1

        if (jcell.areas[n] > puny) && (jcell.H_bnew[nplus1] - jcell.H_bnew[n] > puny)
            # Initialize hL and hR
            hL   = jcell.H_bnew[n]
            hR   = jcell.H_bnew[nplus1]
            hice = jcell.columns[n].H_i

            # Change hL or hR if hicen(n) falls outside central third of range
            h13 = (2.0hL + hR) / 3.0
            h23 = (hL + 2.0hR) / 3.0
            if (hice < h13)
                hR = 3.0hice - 2.0hL
            elseif (hice > h23)
                hL = 3.0hice - 2.0hR
            end

            # Compute coefficients of g(eta) = g0 + g1*eta
            dhr = 1.0 / (hR - hL)
            wk1 = 6.0jcell.areas[n] * dhr
            wk2 = (hice - hL) * dhr

            jcell.g0[nplus1] = wk1 * (2.0/3.0 - wk2)
            jcell.g1[nplus1] = 2.0dhr * wk1 * (wk2 - 0.5)
            jcell.hL[nplus1] = hL
            jcell.hR[nplus1] = hR
        else
            jcell.g0[nplus1] = 0.0
            jcell.g1[nplus1] = 0.0
            jcell.hL[nplus1] = 0.0
            jcell.hR[nplus1] = 0.0
        end

        # Also reset change variables:
        jcell.donor[n]  = 0
        jcell.dareas[n] = 0.0
        jcell.dvol_i[n] = 0.0
        jcell.dvol_s[n] = 0.0
    end

    return nothing
end


# Shifts ice between layers, as well as ice enthalpy for conservation of energy
# Based on shift_ice in Icepack/icepack_itd.F90, line ~356
@inline function shift_ice(jcell)

    # First set old areas and volumes:
    for n in 1:jcell.N_cat
        jcell.areas_old[n] = jcell.areas[n]
        jcell.vol_i[n]     = jcell.columns[n].H_i * jcell.areas[n]
        jcell.vol_s[n]     = jcell.columns[n].H_s * jcell.areas[n]
        jcell.vol_i_old[n] = jcell.vol_i[n]
        jcell.vol_s_old[n] = jcell.vol_s[n]
    end
    
    # Define variables equal to (area or volume) * tracers (line ~440 in icepack_itd)
    for n in 1:jcell.N_cat

        # idea: atrcrn is the tracers (ice/snow enthalpy) * the column's area
        # the conditional statement handles multiple layers
        # it looks like atrcrn just takes the product of all layers in a category if there are multiple
        jcolumn = jcell.columns[n]

        # NOTE: we are temporalily modifying jcolumn's enthalpy with its area/volume
        # first store q values in q_new as well:
        #for k in 1:(jcolumn.N_s+jcolumn.N_i+1)
        #    jcolumn.q_new[k] = jcolumn.q[k]
        #end

        # get area-adjusted surface enthalpy:
        jcolumn.q[1] *= jcell.areas_old[n]

        # we'll get volume-adjusted snow enthalpy:
        w_s = jcell.vol_s_old[n] / jcolumn.N_s
        for k in 2:(jcolumn.N_s+1)
            jcolumn.q[k] *= w_s
        end

        # then we'll get volume-adjusted ice enthalpy:
        w_i = jcell.vol_i_old[n] / jcolumn.N_i
        for k in (jcolumn.N_s+2):(jcolumn.N_s+jcolumn.N_i+1)
            jcolumn.q[k] *= w_i
        end
    end

    # Check for daice or dvice out of range (line ~458)
    for n in 1:(jcell.N_cat-1)

        darea_neg = false
        dvoli_neg = false
        darea_gre = false
        dvoli_gre = false

        if jcell.donor[n] > 0

            nd = jcell.donor[n]

            # If change in area is negative, either round back to 0 or report an error
            if jcell.dareas[n] < 0.0
                if jcell.dareas[n] > -puny*jcell.areas[nd]
                    jcell.dareas[n] = 0.0
                    jcell.dvol_i[n] = 0.0
                else
                    darea_neg = true
                end
            end

            # If change in ice volume is negative, either round back to 0 or report an error
            if jcell.dvol_i[n] < 0.0
                if jcell.dvol_i[n] > -puny*jcell.vol_i[nd]
                    jcell.dareas[n] = 0.0
                    jcell.dvol_i[n] = 0.0
                else
                    dvoli_neg = true
                end
            end

            # If change in area is greater than source area, either round back to source or report error
            if jcell.dareas[n] > jcell.areas[nd]*(1.0-puny)
                if jcell.dareas[n] < jcell.areas[nd]*(1.0+puny)
                    jcell.dareas[n] = jcell.areas[nd]
                    jcell.dvol_i[n] = jcell.vol_i[nd]
                else
                    darea_gre = true
                end
            end

            # If change in ice volume is greater than source volume, either round back to source or report error
            if jcell.dvol_i[n] > jcell.vol_i[nd]*(1.0-puny)
                if jcell.dvol_i[n] < jcell.vol_i[nd]*(1.0+puny)
                    jcell.dareas[n] = jcell.areas[nd]
                    jcell.dvol_i[n] = jcell.vol_i[nd]
                else
                    dvoli_gre = true
                end
            end

        end

        # Error messages
        if darea_neg
            println("Oh no! The change in ice area for a column is negative, that is unphysical")
        end
        if dvoli_neg
            println("Oh no! The change in ice volumn for a column is negative, that is unphysical")
        end
        if darea_gre
            println("Oh no! The change in ice area for a column is greater than its source, that is numerically unstable")
        end
        if dvoli_gre
            println("Oh no! The change in ice volumn for a column is greater than its source, that is numerically unstable")
        end
    end

    # Transfer volume and energy between categories (line ~595)
    for n in 1:(jcell.N_cat-1)

        if jcell.dareas[n] > 0.0

            nd = jcell.donor[n]
            nr = n
            if nd == n
                nr += 1
            end
        

            # Modify the area
            jcell.areas[nd] -= jcell.dareas[n]
            jcell.areas[nr] += jcell.dareas[n]

            # Modify the ice volume
            jcell.vol_i[nd] -= jcell.dvol_i[n]
            jcell.vol_i[nr] += jcell.dvol_i[n]

            # Now modify the snow volume
            worka  = jcell.dareas[n] / jcell.areas_old[nd]
            dvsnow = jcell.vol_s_old[nd] * worka
            #workb  = dvsnow

            jcell.vol_s[nd] -= dvsnow
            jcell.vol_s[nr] += dvsnow

            # Where changes in enthalpy are traced (adapted from line ~621):
            # iterate over variables to trace (def enthalpy, maybe others too)
                # get nd and nr (necessary? we already get those)
                # compute change in tracer (datrcr)
                # modify datrcr further if n_trcr_strata(index of this tracer) > 0

                # finally modify area of tracer at indices nd and nr using datrcr

                # idea: atrcrn is the tracers (ice/snow enthalpy) * the column's area
                # the conditional statement handles multiple layers
                # it looks like atrcrn just takes the product of all layers in a category if there are multiple

            # we draw from the nd column to get the change in tracers
            jcolumn_d = jcell.columns[nd]
            jcolumn_r = jcell.columns[nr]

            # NOTE: we store change in enthalpy in jcolumn.q_new
            # NOTE: since we weighted enthalpy by area/volume, we have to divide q here by old areas/volumes

            # first we'll get surface change to weighted enthalpy:
            jcolumn_d.q_new[1] = jcolumn_d.q[1] * jcell.dareas[n] / jcell.areas_old[n]

            # next we'll get snow change to weighted enthalpy:
            w_s = dvsnow / jcell.vol_s_old[n]
            for k in 2:(jcolumn_d.N_s+1)
                jcolumn_d.q_new[k] = jcolumn_d.q[k] * w_s
            end

            # then we'll get ice change to weighted enthalpy:
            w_i = jcell.dvol_i[n] / jcell.vol_i_old[n]
            for k in (jcolumn_d.N_s+2):(jcolumn_d.N_s+jcolumn_d.N_i+1)
                jcolumn_d.q_new[k] = jcolumn_d.q[k] * w_i
            end

            # Now we modify the area/volume-adjusted enthalpies.
            # NOTE: since surface is sometimes snow and sometimes ice, we might need to
            # handle special cases there
            for k in 1:(jcolumn_d.N_s+jcolumn_d.N_i+1)

                jcolumn_d.q[k] -= jcolumn_d.q_new[k]
                jcolumn_r.q[k] += jcolumn_d.q_new[k]
            end
        end
    end

    # Update snow/ice thickness and tracers (line ~647)
    for n in 1:jcell.N_cat

        jcolumn = jcell.columns[n]

        if jcell.areas[n] > puny
            if jcell.vol_s[n] > puny
                jcolumn.H_s = jcell.vol_s[n] / jcell.areas[n]
                for k in 2:(jcolumn.N_s+1)
                    jcolumn.Δh[k] = jcolumn.H_s / jcolumn.N_s
                end
            end
            if jcell.vol_i[n] > puny
                jcolumn.H_i = jcell.vol_i[n] / jcell.areas[n]
                for k in (jcolumn.N_s+2):(jcolumn.N_s+jcolumn.N_i+1)
                    jcolumn.Δh[k] = jcolumn.H_i / jcolumn.N_i
                end
            end
        end
    end

    compute_new_enthalpy(jcell)

    return nothing
end

#=
    Computes the new enthalpy values based on the modified tracers
=#
@inline function compute_new_enthalpy(jcell)

    # We need to divide our enthalpies by the new weighted area/volume
    for n in 1:jcell.N_cat

        jcolumn = jcell.columns[n]
        w_s     = jcell.vol_s[n] / jcolumn.N_s
        w_i     = jcell.vol_i[n] / jcolumn.N_i

        jcolumn.q[1] /= jcell.areas[n]
        for k in 2:(jcolumn.N_s+1)
            jcolumn.q[k] /= w_s
        end
        for k in (jcolumn.N_s+2):(jcolumn.N_s+jcolumn.N_i+1)
            jcolumn.q[k] /= w_i
        end

        # And then get the new temperatures from these enthalpies:
        generate_T_from_q(jcolumn.T_n, jcolumn.N_i, jcolumn.N_s, jcolumn.H_s, jcolumn.q, jcolumn.S)
    end

    return nothing
end


#= Adds up sum for total ice/snow energy

    ice_energy is either i_energy or i_energy_old
    sno_energy is either s_energy or s_energy_old, before or after model run
=#
@inline function sum_total_energy(ice_energy, sno_energy, jcell)

    for n in 1:jcell.N_cat

        jcolumn = jcell.columns[n]

        ice_energy[n] = 0.0
        sno_energy[n] = 0.0
        
        # Get snow energy if it exists
        if jcolumn.H_s > 0.0
            w_s = jcell.areas[n]*jcolumn.H_s / jcolumn.N_s
            for k in 2:(jcolumn.N_s+1)
                sno_energy[n] += jcolumn.q[k] * w_s
            end
        end
        # Get ice energy
        w_i = jcell.areas[n]*jcolumn.H_i / jcolumn.N_i
        for k in (jcolumn.N_s+2):(jcolumn.N_s+jcolumn.N_i+1)
            ice_energy[n] += jcolumn.q[k] * w_i
        end

    end

    return nothing
end