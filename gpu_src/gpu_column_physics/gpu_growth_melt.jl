# Written by Joseph Kump (josek97@utexas.edu)
# Determines thermodynamic-induced growth and melting of sea ice,
# and rebalances layers obeying conservation of energy.

const m1 = 1.6e-6       # (m/s/deg^(-m2)) constant from Maykut & Perovich
const m2 = 1.36         # (unitless)      constant from Maykut & Perovich
const floeshape = 0.66  # (unitless)      constant
const floediam  = 300.0 # (m)             effective floe diameter for lateral melt

# Computre the preliminary bottom and lateral freeze/melt
@inline function frzmlt_bottom_lateral(Δt, N_i, N_s, H_i, H_s, T_w, T_frz, frzmlt, area, u_star_min, rside, fside, f_bot)

    rside[1] = 0.0
    fside[1] = 0.0
    T_bot    = T_frz
    f_bot[1] = 0.0
    w_lat    = 0.0

    if (area > puny) && (frzmlt < 0.0)

        # TODO: include ice-ocean stress to better compute u_star
        δT        = max((T_w - T_bot), 0.0)
        u_star    = u_star_min
        cpchr     = -0.006c_ocn*ρ_w
        f_bot[1]  = cpchr * u_star * δT
        f_bot[1]  = max(f_bot[1], frzmlt) # limiting how much bottom flux can take place

        # Next we compute rside
        w_lat    = m1 * δT^m2
        rside[1] = w_lat*Δt*pi / (floeshape*floediam)
        rside[1] = max(0.0, min(rside, 1.0))


        # Compute heat flux associated with this value of rside:
        e_tot_s = 0.0
        e_tot_i = 0.0

        for k in 2:(N_s+1)
            e_tot_s += q[k]
        end
        for k in (N_s+2):(N_i+N_s+1)
            e_tot_i += q[k]
        end
        e_tot_s *= H_s * area
        e_tot_i *= H_i * area

        # Different here if we trace fsd:
        fside[1] += rside[1]*(e_tot_s + e_tot_i) / Δt

        # Limit bottom and lateral heat fluxes if necessary:
        x_tmp     = frzmlt / (f_bot[1] + f_side[1] + puny)
        x_tmp     = min(x_tmp, 1.0)
        rside[1] *= x_tmp
        fside[1] *= x_tmp
        f_bot[1] *= x_tmp

    end

    return nothing
end

# Runs a single time step of the ice growth/melt in bottom layer and rebalancing
@inline function step_growth_melt(N_i, N_s, H_s, S, T_frz, Δh, T_new, K, K̄, q, q_new, z_old, z_new, Δt, u_star, T_w, F_0)

    # Updating ice enthalpy:
    generate_q_from_T(q, N_i, N_s, H_s, T_new, S)

    # Updating top ice/snow thickness:
    surface_ice_snow_change(N_s, H_s, T_new, q, K̄, F_0, Δh, Δt)

    # Updating bottom ice thickness:
    bottom_ice_change(N_i, N_s, Δh, T_new, q, K, T_frz, Δt, u_star, T_w)

    # Restore the sea ice thicknesses so that they are uniform again, and calculate new
    # enthalpies to preserve conservation of energy:
    rebalance_ice_layers(Δh, q, N_i, N_s, H_s, q_new, z_old, z_new)

    # Finally, recompute layer tempreatures to match new enthalpies:
    generate_T_from_q(T_new, N_i, N_s, H_s, q, S)

    return nothing
end

# Computes the ice enthalpy q from the ice temperature T
@inline function generate_q_from_T(q, N_i, N_s, H_s, T, S)

    # For the surface, using either ice or snow formulation
    w    = -c_0*T[1] + L_0 # intermediate variable
    q[1] = -ρ_i*w

    # For the snow layers
    if H_s >= puny
        q[1] = -ρ_s*w
        for k in 2:(N_s+1)
            q[k] = -ρ_s*(-c_0*T[k]+L_0)
        end
    end

    # For the ice layers
    for k in 1:N_i
        index    = k + N_s + 1
        T_m      = -μ*S[k]
        q[index] = -ρ_i*(c_0*(T_m - T[index]) + L_0*(1 - T_m/T[index]) - c_w*T_m)
    end

    return nothing
end

# Computes the ice temperature T from the ice enthalpy q
@inline function generate_T_from_q(T, N_i, N_s, H_s, q, S)

    # For the snow layers
    if H_s >= puny
        for k in 1:(N_s+1)
            T[k] = (q[k]/ρ_s + L_0)/c_0
        end
    else # need surface temp based on ice
        b    = - (q[1]/ρ_i + L_0)
        T[1] = (-b - sqrt(b^2))/(2c_0)
    end

    # For the ice layers
    for k in 1:N_i
        index    = k+N_s+1
        T_m      = -μ*S[k]
        b        = (c_w-c_0)*T_m - (q[index]/ρ_i) - L_0
        T[index] = (-b - sqrt(b^2 - 4c_0*L_0*T_m))/(2c_0)
    end

    return nothing
end

# Computes the melt at the surface ice/snow, if such melting occurs
@inline function surface_ice_snow_change(N_s, H_s, T_new, q, K̄, F_0, Δh, Δt)

    # Get the first subsurface layer, ice if there's no snow
    k = 2
    if H_s < puny
        k += N_s
    end

    # Get the conditional surface flux
    F_ct = K̄[k-1] * (T_new[1] - T_new[k])

    # Then the possible ice/snow melt
    δh     = max((Δt*(F_0 - F_ct) / q[1]), 0.0)
    Δh[k] -= δh

    return nothing
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
@inline function rebalance_ice_layers(Δh, q, N_i, N_s, H_s, q_new, z_old, z_new)

    # Rebalance snow layers here
    if H_s >= puny
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
            z_old[k] = Δh[k] + z_old[k-1]
            z_new[k] = Δh_snew * (k-1)
        end

        # Update enthalpies for each layer:
        for k in 1:N_s
            q_new[k+1] = 0.0
            for m in 1:N_s
                η_km        = min(z_old[m+1], z_new[k+1]) - max(z_old[m], z_new[k])
                η_km        = max(η_km, 0.0)
                q_new[k+1] += η_km*q[m+1]
            end
            q_new[k+1] /= Δh_snew
        end

        # Now that we have q_new and Δh_snew, we can replace current q and Δh with these:
        for k in 2:(N_s+1)
            q[k]  = q_new[k]
            Δh[k] = Δh_snew
        end
    end

    # NEXT, rebalance ice layers

    # Get new total ice thickness
    H_inew = 0.0
    for k in (N_s+2):(N_i+N_s+1)
        H_inew += Δh[k]
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
        k_index = k+N_s+1
        q_new[k_index] = 0.0
        for m in 1:N_i
        
            m_index = m+N_s+1
            η_km    = min(z_old[m_index+1], z_new[k_index+1]) - max(z_old[m_index], z_new[k_index])
            η_km    = max(η_km, 0.0)

            q_new[k_index] += η_km*q[m_index]
        end
        q_new[k_index] /= Δh_inew
    end

    # Now that we have q_new and Δh_inew, we can replace current q and Δh with these:
    for k in 1:N_i
        index     = k+N_s+1
        q[index]  = q_new[index]
        Δh[index] = Δh_inew
    end

end