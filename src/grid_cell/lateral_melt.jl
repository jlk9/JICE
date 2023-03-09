# Written by Joseph Kump (josek97@utexas.edu)
# Code that computes the lateral melt of ice in each grid cell column after remapping and adding of new ice

include("./jicecell_struct.jl")

#= TODO: add details here
=#
@inline function lateral_melt(jcell)

    # Set initial variables (line 994)
    dheat_flux    = 0.0
    dfresh_flux   = 0.0
    dsalt_flux    = 0.0
    #= Variables yet unused (maybe fsd?)
    dvssl         = 0.0
    dvint         = 0.0
    cat1_arealoss = 0.0
    tmp           = 0.0
    vicen_init    = 0.0
    G_radialn     = 0.0
    delta_an      = 0.0
    qin           = 0.0
    afsdn         = 1.0
    afsdn_init    = 0.0
    df_flx        = 0.0
    f_flx         = 0.0
    =#
    flag = false

    # Optional FSD calculations


    # Computing average enthalpy of ice (line 1028)
    # This is pretty simple since we're not tracing FSD
    if jcell.rside > 0.0
        flag = true
        for n in 1:jcell.N_cat
            jcell.rside_array[n] = rside
        end
    end

    # Melt ice and increment fluxes (line 1090)
    if flag

        for n in 1:jcell.N_cat

            jcolumn = jcell.columns[n]

            # Fluxes to coupler, assuming dfresh_flux, dsalt_flux both > 0

            dfresh_flux = ρ_i*jcell.vol_i[n] + ρ_s*jcell.vol_s[n]

            # Assuming saltflux option is not prognostic
            dsalt_flux = 0.001ρ_i*jcell.vol_i[n]*ice_ref_salinity*jcell.rside_array[n] / jcell.Δt

            jcell.fresh_flux += dfresh_flux
            jcell.salt_flux  += dsalt_flux

            # TODO(?) history diagnostics (line 1115)

            # State variables. First set old variables then modify new ones based on change in ice:
            jcell.areas_old[n] = jcell.areas[n]
            jcell.vol_i_old[n] = jcell.vol_i[n]
            jcell.vol_s_old[n] = jcell.vol_s[n]

            jcell.areas[n] *= (1.0 - jcell.rside_array[n])
            jcell.vol_i[n] *= (1.0 - jcell.rside_array[n])
            jcell.vol_s[n] *= (1.0 - jcell.rside_array[n])

            # fsd (ignore for now)

            # fluxes

            # enthalpy tracers don't change since energy/vol is constant, but we do need to track it for heat flux
            # first snow:
            # heat flux = enthalpy * volume of snow to be melted (rside*vol_s_old / N_s) per time = enthalpy * w_s here
            w_s = (jcell.rside_array[n] * jcell.vol_s_old[n]) / (jcell.Δt * jcell.N_s)
            for k in 2:(jcell.N_s+1)
                dheat_flux       = jcolumn.q[k] * w_s
                jcell.heat_flux += dheat_flux
            end
            # then ice:
            # heat flux = enthalpy * volume of ice to be melted (rside*vol_i_old / N_i) per time = enthalpy * w_i here
            w_i = (jcell.rside_array[n] * jcell.vol_i_old[n]) / (jcell.Δt * jcell.N_i)
            for k in (jcell.N_s+2):(jcell.N_s+jcell.N_i+1)
                # heat flux = enthalpy * volume of snow to be melted (rside*vol_s_old / N_s) per time
                dheat_flux       = jcolumn.q[k] * w_i
                jcell.heat_flux += dheat_flux
            end

        end
    end


    # Biogeochemistry (ignore for now)

    return nothing
end
