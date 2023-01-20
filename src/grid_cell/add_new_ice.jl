# Written by Joseph Kump (josek97@utexas.edu)
# Code that computes the addition of new ice into a grid cell after remapping

include("./jicecell_struct.jl")

#= Based on function of the same name in icepack_therm_itd, line 1284
=#
@inline function add_new_ice(jcell)

    # Initialize new variables (line 1493)
    hsurp  = 0.0
    hi0new = 0.0
    ai0new = 0.0
    #=
    d_an_latg(:) = 0.0
    d_an_tot(:) = 0.0
    d_an_newi(:) = 0.0
    vin0new(:) = 0.0
    area2(:) = aicen(:)
    =#
    lead_area    = 0.0
    latsurf_area = 0.0
    G_radial     = 0.0
    tot_latg     = 0.0

    hi0max = 1000000.0
    if jcell.N_cat > 1
        hi0max = 0.9jcell.H_bds[2]
    end

    # conditional for fsd

    # We will reuse old area and ice volume fields here since we are done with itd shift (line 1526):
    for n in 1:jcell.N_cat

        jcell.areas_old[n] = jcell.areas[n]
        jcell.vol_i_old[n] = jcell.vol_i[n]
        jcell.vol_s_old[n] = jcell.vol_s[n]

        jcell.i_energy_old[n] = 0.0
        jcell.s_energy_old[n] = 0.0
        jcell.i_energy[n]     = 0.0
        jcell.s_energy[n]     = 0.0
    end

    # Conditional for conservation of energy (line 1537)


    # Compute average enthalpy of new ice (line 1553)


    # Compute the volume, area, and thickness of new ice (line 1579)


    # Update freshwater and salt fluxes (line 1601)


    # Decide how to distribute new ice (line 1633)


    # Distribute excess ice volume among categories by increasing thickness (line 1694)


    # Combine new ice into cat 1 since we don't track FSD yet (line 1792)
    # NOTE: this will be more involved with FSD factored


    # Biogeochemistry (ignored for now)

    return nothing
end