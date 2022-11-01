# Written by Joseph Kump (josek97@utexas.edu)
# Contains surface heat flux functions for running the sea ice model, featuring
# the ATModel struct

include("./atmodel_struct.jl")

# Computes the (constant) atmospheric flux affecting the model
@inline function step_surface_flux(α, i_0, T_sfc, H, F_0, dF_0, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l, F_sw, F_Ld, c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, ρ_a, c_p, step)

    # Compute atmospheric fluxes dependent on ice:
    set_atm_flux_values(F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l, c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, ρ_a, c_p, T_sfc, H, step)
    

    # Reduce shortwave flux with albedo
    #α = 0.7 # CHANGE TO NOT PRESET later

    # Now compute total surface flux:
    F_0[step] = (1-α)*i_0*F_sw + F_Ld + F_Lu[step] + F_l[step] + F_s[step]
    # And now compute derivative of flux:
    dF_0[step] = dF_Lu[step] + dF_s[step] + dF_l[step]

    return nothing
end

# Sets helper values needed to compute flux
@inline function set_atm_flux_values(F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l, c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, ρ_a, c_p, T_sfc, z_ice, step)

    # Compute initial exchange coefficients:
    c_u[1] = κ / log(z_ref / z_ice)
    c_Θ[1] = c_u[1]
    c_q[1] = c_u[1]

    Q_sfc = (q_1/ρ_a)*exp(-q_2 / (T_sfc + C_to_K))

    # Iterate and update exchange coefficients:
    #=
    for k in 1:natmiter
        set_atm_helper_values_step(c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, T_sfc, Q_sfc)
    end
    =#
    # Approach for work with AD
    set_atm_helper_values_step(c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, T_sfc, Q_sfc)
    set_atm_helper_values_step(c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, T_sfc, Q_sfc)
    set_atm_helper_values_step(c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, T_sfc, Q_sfc)
    set_atm_helper_values_step(c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, T_sfc, Q_sfc)
    set_atm_helper_values_step(c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, T_sfc, Q_sfc)

    # Compute heatflux coefficients. TODO: wind stress?
    C_l = ρ_a * (L_vap + L_ice) * atm_u_star[1] * c_q[1]
    C_s = ρ_a * c_p * atm_u_star[1] * c_Θ[1] + 1

    # Outgoing longwave flux:
    F_Lu[step]  = emissivity * sbc * (T_sfc + C_to_K)^4
    dF_Lu[step] = 4.0 * emissivity * sbc * (T_sfc + C_to_K)^3

    # Sensible heat flux:
    F_s[step]  = C_s * (Θ_a[1] - (T_sfc + C_to_K))
    dF_s[step] = -C_s

    # Latent heat flux:
    Q_sf    = (q_1/ρ_a) * exp(-q_2 / (T_sfc + C_to_K))
    dQ_sf   = (q_1/ρ_a) * exp(-q_2 / (T_sfc + C_to_K)) * (-q_2 / (T_sfc + C_to_K)^2)
    F_l[step]  = C_l * (Q_a[1] - Q_sf)
    dF_l[step] = -C_l * dQ_sf

    return nothing
end

# Performs one step in the iteration for set_atm_helper_values
@inline function set_atm_helper_values_step(c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, T_sfc, Q_sfc)

    # Update turbulent scales
    atm_u_star[1] = c_u[1] * max(U_dmin, (U_a[1]^2+U_a[2]^2+U_a[3]^2)^.5)
    Θ_star        = c_Θ[1] * (Θ_a[1] - T_sfc)
    Q_star        = c_q[1] * (Q_a[1] - Q_sfc)

    # Update Y and compute χ
    Y = (κ*g*z_deg)*(Θ_star / (Θ_a[1]*(1+0.606Q_a[1])) + Q_star / (1.0/0.606 + Q_a[1]))/(atm_u_star[1]^2)

    # Compute ψ_m and ψ_s TODO: add alternative if Y is unstable? (< 0)
    ψ_m = -(0.7Y + 0.75*(Y-14.3)*exp(-0.35Y) + 10.7)
    ψ_s = ψ_m

    if Y < 0.0
        χ   = (1 - 16Y)^0.25
        ψ_m = 2*log(0.5*(1+χ)) + log(0.5*(1+χ^2)) - 2*atan(χ) + π/2
        ψ_s = 2*log(0.5*(1+χ^2))
    end

    # Update exchange coefficients
    c_u[1] = c_u[1] / (1 + c_u[1]*(λ - ψ_m)/κ)
    c_Θ[1] = c_Θ[1] / (1 + c_Θ[1]*(λ - ψ_s)/κ)
    c_q[1] = c_Θ[1]

    return nothing
end
