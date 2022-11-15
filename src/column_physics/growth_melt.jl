# Written by Joseph Kump (josek97@utexas.edu)
# Determines thermodynamic-induced growth and melting of sea ice,
# and rebalances layers obeying conservation of energy.

# Runs a single time step of the ice growth/melt in bottom layer and rebalancing
@inline function step_growth_melt(N_i, N_s, S, T_frz, Δh, T_new, K, q_i, q_inew, z_old, z_new, Δt, u_star, T_w)

    # Updating ice enthalpy:
    generate_q_i_from_T(q_i, N_i, N_s, T_new, S)

    # Updating bottom ice thickness:
    bottom_ice_change(N_i, N_s, Δh, T_new, q_i, K, T_frz, Δt, u_star, T_w)

    # Restore the sea ice thicknesses so that they are uniform again, and calculate new
    # enthalpies to preserve conservation of energy:
    rebalance_ice_layers(Δh, q_i, N_i, N_s, q_inew, z_old, z_new)

    # Finally, recompute layer tempreatures to match new enthalpies:
    generate_T_from_q_i(T_new, N_i, N_s, q_i, S)

    return nothing
end

# Computes the ice enthalpy q_i from the ice temerature T
@inline function generate_q_i_from_T(q_i, N_i, N_s, T, S)

    for k in 1:N_i
        T_m  = -μ*S[k]
        q_i[k] = -ρ_i*(c_0*(T_m - T[k+N_s+1]) + L_0*(1 - T_m/T[k+N_s+1]) - c_w*T_m)
    end

    return nothing
end

# Computes the ice temperature T from the ice enthalpy q_i
@inline function generate_T_from_q_i(T, N_i, N_s, q_i, S)

    for k in 1:N_i
        T_m    = -μ*S[k]
        b      = (c_w-c_0)*T_m - (q_i[k]/ρ_i) - L_0
        T[k+N_s+1] = (-b - sqrt(b^2 - 4c_0*L_0*T_m))/(2c_0)
    end

    return nothing
end

# Computes the change in thickness of the bottom ice layer due to growth / melting
# This function is based on CICE documentation
@inline function bottom_ice_change(N_i, N_s, Δh, T_n, q_i, K, T_frz, Δt, u_star, T_w)

    # Get conductive heat flux at bottom surface:
    F_cb = (T_n[N_i+N_s+1] - T_frz) * K[N_i+N_s+1]/Δh[N_i+N_s+1]

    # Get net downward heat flux:
    F_bot = -ρ_w*c_w*c_h*u_star*(T_w - T_frz)

    # Now calculate the change in ice thickness in the bottom layer:
    δz = (F_cb - F_bot)*Δt / q_i[N_i]

    # And finally update the bottom later thickness:
    Δh[N_i+N_s+1] -= δz

    return nothing
end

# Redistributes the ice thicknesses so they are uniform, then updates the
# enthalpy q_i to match.
@inline function rebalance_ice_layers(Δh, q_i, N_i, N_s, q_new, z_old, z_new)

    # Get new total ice thickness
    H_inew = 0.0
    for k in (N_s+2):(N_i+N_s+1)
        H_inew += Δh[k]
    end

    # Split evenly
    Δh_inew = H_inew / N_i

    # Update the ice thicknesses
    z_old[1] = 0.0
    z_new[1] = 0.0
    for k in 2:(N_i+1)
        z_old[k] = Δh[k+N_s] + z_old[k-1]
        z_new[k] = Δh_inew * (k-1)
    end

    # Update enthalpies for each layer:
    for k in 1:N_i
        q_new[k] = 0.0
        for m in 1:N_i
        
            η_km      = min(z_old[m+1], z_new[k+1]) - max(z_old[m], z_new[k])
            η_km      = max(η_km, 0.0)
            q_new[k] += η_km*q_i[m]
        end
        q_new[k] /= Δh_inew
    end

    # Now that we have q_new and Δh_inew, we can replace current q_i and Δh with these:
    for k in 1:N_i
        q_i[k]  = q_new[k]
        Δh[k+N_s+1] = Δh_inew
    end

end