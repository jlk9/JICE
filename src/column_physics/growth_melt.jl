# Written by Joseph Kump (josek97@utexas.edu)
# Determines thermodynamic-induced growth and melting of sea ice,
# and rebalances layers obeying conservation of energy.

# Runs a single time step of the ice growth/melt in bottom layer and rebalancing
@inline function step_growth_melt(N_i, S, T_frz, Δh, T_old, T_new, K, q_i, q_inew, z_old, z_new, Δt, u_star, T_w)

    # Updating ice enthalpy:
    generate_q_i_from_T(q_i, N_i, T_new, S)

    # Updating bottom ice thickness:
    bottom_ice_change(N_i, Δh, T_new, q_i, K, T_frz, Δt, u_star, T_w)

    # Restore the sea ice thicknesses so that they are uniform again, and calculate new
    # enthalpies to preserve conservation of energy:
    #(Since K is no longer needed at this time step, we use it for intermediate steps)
    rebalance_ice_layers(Δh, q_i, N_i, q_inew, z_old, z_new)

    # Finally, recompute layer tempreatures to match new enthalpies:
    generate_T_from_q_i(T_new, N_i, q_i, S)

    return nothing
end

# Computes the enthalpy q from the temerature T
@inline function generate_q_i_from_T(q_i, N_i, T, S)

    for k in 1:N_i+1
        T_m  = -μ*S[k]
        q_i[k] = -ρ_i*(c_0*(T_m - T[k]) + L_0*(1 - T_m/T[k]) - c_w*T_m)
    end

    return nothing
end

# Computes the temperature T from the enthalpy q
@inline function generate_T_from_q_i(T, N_i, q_i, S)

    # Surface temperature is not changed
    for k in 1:(N_i+1)
        T_m  = -μ*S[k]
        b    = (c_w-c_0)*T_m - (q_i[k]/ρ_i) - L_0
        T[k] = (-b - sqrt(b^2 - 4c_0*L_0*T_m))/(2c_0)
    end

end

# Computes the change in thickness of the bottom ice layer due to growth / melting
# This function is based on CICE documentation
@inline function bottom_ice_change(N_i, Δh, T_n, q_i, K, T_frz, Δt, u_star, T_w)

    # Get conductive heat flux at bottom surface:
    F_cb = (T_n[N_i+1] - T_frz) * K[N_i+1]/Δh[N_i+1]

    # Get net downward heat flux:
    F_bot = -ρ_w*c_w*c_h*u_star*(T_w - T_frz)

    # Now calculate the change in ice thickness in the bottom layer:
    δz = (F_cb - F_bot)*Δt / q_i[N_i+1]

    # And finally update the bottom later thickness:
    Δh[N_i+1] -= δz

    return nothing
end

# Redistributes the ice thicknesses so they are uniform, then updates the
# enthalpy q_i to match.
@inline function rebalance_ice_layers(Δh, q_i, N_i, q_new, z_old, z_new)

    # Get new total thickness
    H_new = 0.0
    for k in 1:(N_i+1)
        H_new += Δh[k]
    end

    # Split evenly
    Δh_new = H_new / N_i

    # Update the ice thicknesses
    z_old[1] = Δh[1]
    z_new[1] = 0.0
    for k in 2:(N_i+1)
        z_old[k] = Δh[k] + z_old[k-1]
        z_new[k] = Δh_new * (k-1)
    end

    # Update enthalpies for each layer:
    for k in 2:(N_i+1)
        q_new[k] = 0.0
        for m in 2:(N_i+1)
        
            η_km      = min(z_old[m], z_new[k]) - max(z_old[m-1], z_new[k-1])
            η_km      = max(η_km, 0.0)
            q_new[k] += η_km*q_i[m]
        end
        q_new[k] /= Δh_new
    end

    # Now that we have q_new and Δh_new, we can replace current q_i and Δh with these:
    for k in 2:(N_i+1)
        q_i[k] = q_new[k]
        Δh[k]  = Δh_new
    end

end