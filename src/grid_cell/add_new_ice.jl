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

    # We will reuse old area and ice volume fields here since we are done with itd shift:
    for n in 1:jcell.N_cat

        jcell.areas_old[n] = jcell.areas[n]
        jcell.vol_i_old[n] = jcell.vol_i[n]
        # TODO?: set energy = 0 for conservation of energy (could be good for itd as well) 
    end

end