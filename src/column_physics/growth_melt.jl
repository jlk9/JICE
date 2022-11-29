# Written by Joseph Kump (josek97@utexas.edu)
# Determines thermodynamic-induced growth and melting of sea ice,
# and rebalances layers obeying conservation of energy.

# Runs a single time step of the ice growth/melt in bottom layer and rebalancing
@inline function step_growth_melt(N_i, N_s, S, T_frz, Δh, T_new, K, q, q_new, z_old, z_new, Δt, u_star, T_w)

    # Updating ice enthalpy:
    generate_q_from_T(q, N_i, N_s, T_new, S)

    # Updating bottom ice thickness:
    bottom_ice_change(N_i, N_s, Δh, T_new, q, K, T_frz, Δt, u_star, T_w)

    # Restore the sea ice thicknesses so that they are uniform again, and calculate new
    # enthalpies to preserve conservation of energy:
    rebalance_ice_layers(Δh, q, N_i, N_s, q_new, z_old, z_new)

    # Finally, recompute layer tempreatures to match new enthalpies:
    generate_T_from_q(T_new, N_i, N_s, q, S)

    return nothing
end

# Computes the ice enthalpy q from the ice temperature T
@inline function generate_q_from_T(q, N_i, N_s, T, S)

    # For the surface, using either ice or snow formulation
    q[1] = -ρ_s*(-c_0*T[1]+L_0) + min(N_s, 1)*(-ρ_i*(c_0*(-T[1]) + L_0) + ρ_s*(-c_0*T[1]+L_0))

    # For the snow layers
    for k in 1:N_s
        q[k+1] = -ρ_s*(-c_0*T[k+1]+L_0)
    end

    # For the ice layers
    for k in 1:N_i
        T_m  = -μ*S[k]
        q[k+N_s+1] = -ρ_i*(c_0*(T_m - T[k+N_s+1]) + L_0*(1 - T_m/T[k+N_s+1]) - c_w*T_m)
    end

    return nothing
end

# Computes the ice temperature T from the ice enthalpy q
@inline function generate_T_from_q(T, N_i, N_s, q, S)

    # For the snow layers
    for k in 1:N_s
        T[k+1] = (q[k+1]/ρ_s + L_0)/c_0
    end

    # For the ice layers
    for k in 1:N_i
        T_m    = -μ*S[k]
        b      = (c_w-c_0)*T_m - (q[k+N_s+1]/ρ_i) - L_0
        T[k+N_s+1] = (-b - sqrt(b^2 - 4c_0*L_0*T_m))/(2c_0)
    end

    return nothing
end

# Computes the melt at the surface ice/snow, if such melting occurs
@inline function surface_ice_snow_change(T_new, q, K̄)

    # Get the conditional surface flux
    F_ct = K̄[1] * (T_new[1] - T_new[2])

end

# Computes the change in thickness of the bottom ice layer due to growth / melting
# This function is based on CICE documentation
@inline function bottom_ice_change(N_i, N_s, Δh, T_n, q, K, T_frz, Δt, u_star, T_w)

    # Get conductive heat flux at bottom surface:
    F_cb = (T_n[N_i+N_s+1] - T_frz) * K[N_i+N_s+1]/Δh[N_i+N_s+1]

    # Get net downward heat flux:
    F_bot = -ρ_w*c_w*c_h*u_star*(T_w - T_frz)

    # Now calculate the change in ice thickness in the bottom layer:
    δz = (F_cb - F_bot)*Δt / q[N_i+N_s+1]

    # And finally update the bottom later thickness:
    Δh[N_i+N_s+1] -= δz

    return nothing
end

# Redistributes the ice thicknesses so they are uniform, then updates the
# enthalpy q to match.
@inline function rebalance_ice_layers(Δh, q, N_i, N_s, q_new, z_old, z_new)

    # TODO: rebalance snow layers here
    H_snew = 0.0
    for k in 1:N_s
        H_snew += Δh[k+1]
    end
    # Split evenly
    Δh_snew = H_snew / N_s

    # Update the snow thicknesses
    z_old[1] = 0.0
    z_new[1] = 0.0
    for k in 2:(N_s+1)
        z_old[k] = Δh[k] + z_old[k]
        z_new[k] = Δh_snew * (k-1)
    end

    # Update enthalpies for each layer:
    for k in 1:N_s
        q_new[k+1] = 0.0
        for m in 1:N_s
        
            η_km            = min(z_old[m+1], z_new[k+1]) - max(z_old[m], z_new[k])
            η_km            = max(η_km, 0.0)
            q_new[k+1] += η_km*q[m+1]
        end
        q_new[k+1] /= Δh_snew
    end


    # NEXT, rebalance ice layers

    # Get new total ice thickness
    H_inew = 0.0
    for k in 1:N_i
        H_inew += Δh[k+N_s+1]
    end

    # Split evenly
    Δh_inew = H_inew / N_i

    # Update the ice thicknesses
    z_old[N_s+2] = 0.0
    z_new[N_s+2] = 0.0
    for k in 2:(N_i+1)
        z_old[k+N_s+1] = Δh[k+N_s] + z_old[k+N_s]
        z_new[k+N_s+1] = Δh_inew * (k-1)
    end

    # Update enthalpies for each layer:
    for k in 1:N_i
        q_new[k+N_s+1] = 0.0
        for m in 1:N_i
        
            η_km            = min(z_old[m+N_s+2], z_new[k+N_s+2]) - max(z_old[m+N_s+1], z_new[k+N_s+1])
            η_km            = max(η_km, 0.0)
            q_new[k+N_s+1] += η_km*q[m+N_s+1]
        end
        q_new[k+N_s+1] /= Δh_inew
    end

    # Now that we have q_new and Δh_inew, we can replace current q and Δh with these:
    for k in 1:N_i
        q[k+N_s+1] = q_new[k+N_s+1]
        Δh[k+N_s+1]  = Δh_inew
    end

end