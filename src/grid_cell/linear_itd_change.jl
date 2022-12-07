# Written by Joseph Kump (josek97@utexas.edu)
# Code that computes the remapping of ice column areas
# between thickness categories within one grid cell

include("./jicecell_struct.jl")

@inline function linear_itd_change(jcell)

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
                break
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
    if jcell.areas[N_catplus1] > puny
        jcell.H_bnew[N_catplus1] = 3.0jcell.columns[jcell.N_cat].H_i - 2.0jcell.H_bnew[jcell.N_cat]
    else
        jcell.H_bnew[N_catplus1] = jcell.H_bds[N_catplus1]
    end
    jcell.H_bnew[N_catplus1] = max(jcell.H_bnew[N_catplus1], jcell.H_bds[jcell.N_cat])

    # TODO: Boundary check (each one lies between adjacent values of H_i)

    # Compute g0, g1, hL, hR for remapping:
    fit_line(jcell)

    # Find area lost due to melting of thin ice (line ~439) (NOTE: we keep track of all areas including open water)
    if jcell.areas[2] > puny

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
                damax = jcell.areas[2] * (1.0 - jcell.columns[1].H_i / jcell.columns[1].H_iold)
                da0   = min(da0, damax)

                # Remove this area, while conserving volume
                jcell.columns[1].H_i = jcell.columns[1].H_i * jcell.areas[2] / (jcell.areas[2] - da0)
                jcell.areas[2]      -= da0
                jcell.areas[1]      += da0 # MAKE SURE THIS INDEXING IS RIGHT
            end
        else #dh0 >= 0.0
            jcell.H_bnew[1] = min(dh0, jcell.H_bds[2]) # shift lower boundary to the right, as thin ice grew
        end
    end

    # Compute the area and volume to be shifted across each boundary (line ~493)

    return nothing
end

#=  Fit g(h) with a line, satisfying area and volume constraints.
    To reduce roundoff errors caused by large values of g0 and g1,
    we actually compute g(eta), where eta = h - hL, and hL is the
    left boundary.

    Comapred to version in CICE, here we want to do it in one array
=#
@inline function fit_line(jcell)

    # First for first column (different rules, line ~431 in icepack_therm_itd.F90)
    if (jcell.areas[2] > puny) && (jcell.H_bds[2] - jcell.H_bnew[1] > puny)

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
        wk1 = 6.0jcell.areas[2] * dhr
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

        if (jcell.areas[nplus1] > puny) && (jcell.H_bnew[nplus1] - jcell.H_bnew[n] > puny)
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
            wk1 = 6.0jcell.areas[nplus1] * dhr
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
    end

    return nothing
end