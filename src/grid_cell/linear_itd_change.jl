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
    if areas[N_catplus1] > puny
        jcell.H_bnew[N_catplus1] = 3.0jcell.columns[jcell.N_cat].H_i - 2.0jcell.H_bnew[jcell.N_cat]
    else
        jcell.H_bnew[N_catplus1] = jcell.H_bds[N_catplus1]
    end
    jcell.H_bnew[N_catplus1] = max(jcell.H_bnew[N_catplus1], jcell.H_bds[jcell.N_cat])

    # Boundary check (each one lies between adjacent values of H_i)

    return nothing
end

#=  Fit g(h) with a line, satisfying area and volume constraints.
    To reduce roundoff errors caused by large values of g0 and g1,
    we actually compute g(eta), where eta = h - hL, and hL is the
    left boundary.

    Comapred to version in CICE, here we want to do it in one array
=#
@inline function fit_line()


end