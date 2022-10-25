# Written by Joseph Kump (josek97@utexas.edu)
# Contains surface heat flux functions for running the sea ice model, featuring
# the ATModel struct

include("./atmodel_struct.jl")

#= Sets initial flux values for model run. This is only done once
    at the start of the model run, to provide (constant) flux terms
=#
function set_atm_helper_values(atmodel, T_sfc, z_ice, natmiter)

    λ = log(z_deg / z_ref)

    # Compute initial exchange coefficients:
    atmodel.c_u = κ / log(z_ref / z_ice)
    atmodel.c_Θ = atmodel.c_u
    atmodel.c_q = atmodel.c_u

    Q_sfc = (q_1/atmodel.ρ_a)*exp(-q_2 / (T_sfc + 273.15))

    # Iterate and update exchange coefficients:
    for k = 1:natmiter

        # Update turbulent scales
        atmodel.u_star = atmodel.c_u * max(U_dmin, (atmodel.U_a[1]^2+atmodel.U_a[2]^2+atmodel.U_a[3]^2)^.5)
        atmodel.Θ_star = atmodel.c_Θ * (atmodel.Θ_a - T_sfc)
        atmodel.Q_star = atmodel.c_q * (atmodel.Q_a - Q_sfc)

        # Update Y and compute χ
        atmodel.Y = (κ*g*z_deg)*(atmodel.Θ_star / (atmodel.Θ_a*(1+0.606atmodel.Q_a)) + atmodel.Q_star / (1.0/0.606 + atmodel.Q_a))/(atmodel.u_star^2)

        # Compute ψ_m and ψ_s TODO: add alternative if Y is unstable? (< 0)
        atmodel.ψ_m = -(0.7atmodel.Y + 0.75*(atmodel.Y-14.3)*exp(-0.35atmodel.Y) + 10.7)
        atmodel.ψ_s = atmodel.ψ_m

        if atmodel.Y < 0
            atmodel.χ   = (1 - 16atmodel.Y)^0.25
            atmodel.ψ_m = 2*log(0.5*(1+atmodel.χ)) + log(0.5*(1+atmodel.χ^2)) - 2*atan(atmodel.χ) + π/2
            atmodel.ψ_s = 2*log(0.5*(1+atmodel.χ^2))
        end

        # Update exchange coefficients
        atmodel.c_u = atmodel.c_u / (1 + atmodel.c_u*(λ-atmodel.ψ_m)/κ)
        atmodel.c_Θ = atmodel.c_Θ / (1 + atmodel.c_Θ*(λ-atmodel.ψ_s)/κ)
        atmodel.c_q = atmodel.c_Θ

    end

    # Compute heatflux coefficients. TODO: wind stress?
    atmodel.C_l = atmodel.ρ_a*(atmodel.L_vap + atmodel.L_ice)*atmodel.u_star*atmodel.c_q
    atmodel.C_s = atmodel.ρ_a*atmodel.c_p*atmodel.u_star*atmodel.c_Θ + 1

    return nothing
end

# Computes values of fluxes affected by T_sf
function set_atm_flux_values(atmodel, T_sf)

    # Outgoing longwave flux:
    atmodel.F_Lu = emissivity*sbc*(T_sf + 273.15)^4

    # Sensible heat flux:
    atmodel.F_s = atmodel.C_s*(atmodel.Θ_a - (T_sf + 273.15))

    # Latent heat flux:
    Q_sf        = (q_1/atmodel.ρ_a)*exp(-q_2 / (T_sf + 273.15))
    atmodel.F_l = atmodel.C_l*(atmodel.Q_a - Q_sf)

end

# Computes derivatives of fluxes affected by T_sf, based on analytical formula for above
function set_atm_dflux_values(atmodel, T_sf)

    # Outgoing longwave flux:
    atmodel.dF_Lu = 4.0*emissivity*sbc*(T_sf + 273.15)^3

    # Sensible heat flux:
    atmodel.dF_s  = -atmodel.C_s

    # Latent heat flux:
    dQ_sf = (q_1/atmodel.ρ_a)*exp(-q_2 / (T_sf + 273.15))*(-q_2 / (T_sf + 273.15)^2)
    atmodel.dF_l = -atmodel.C_l*dQ_sf

    return nothing

end