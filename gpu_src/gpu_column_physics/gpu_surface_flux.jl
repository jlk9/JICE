# Written by Joseph Kump (josek97@utexas.edu)
# Contains surface heat flux functions for running the sea ice model, featuring
# the ATModel struct

using CUDA

# Computes the (constant) atmospheric flux affecting the model
@inline function step_surface_flux(N_c, N_i, N_layers, α_vdr_i, α_idr_i, α_vdf_i, α_idf_i, α_vdr_s, α_idr_s, α_vdf_s, α_idf_s, T_n, H_i, H_s, F_0, dF_0, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l,
                                    F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, I_pen, c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, ρ_a, c_p, Q_sfc, F_SWsfc, F_SWpen, onGPU)
    #=
    # Computes the current albedo
    if onGPU
        @cuda generate_α_i(N_c, N_layers, H_i, α_vdr_i, α_idr_i, α_vdf_i, α_idf_i, T_n)
        @cuda generate_α_s(N_c, N_layers, α_vdr_s, α_idr_s, α_vdf_s, α_idf_s, T_n)
    else
        generate_α_i(N_c, N_layers, H_i, α_vdr_i, α_idr_i, α_vdf_i, α_idf_i, T_n)
        generate_α_s(N_c, N_layers, α_vdr_s, α_idr_s, α_vdf_s, α_idf_s, T_n)
    end
    =#
    generate_α(N_layers, H_i, α_vdr_i, α_idr_i, α_vdf_i, α_idf_i, α_vdr_s, α_idr_s, α_vdf_s, α_idf_s, T_n)
    
    # Compute atmospheric fluxes dependent on ice:
    set_atm_flux_values(N_c, N_layers, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l, c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, ρ_a, c_p, Q_sfc, T_n, H_i, onGPU)
    
    # Get the fractional snow cover:
    F_SWsfc .= max.(H_s ./ (H_s .+ 0.02), puny / (0.02 + puny))
    F_SWpen .= 1.0 .- F_SWsfc
    
    # Reduce shortwave fluxes with albedo, getting absorbed shortwave fluxes.
    # Note that the albedo terms include ice and snow components:
    # F_SWabs = F_SWabsv .+ F_SWabsi
    F_SWsfc .= (F_SWvdr .* ((1.0.-α_vdr_i).*F_SWpen .+ (1.0.-α_vdr_s).*F_SWsfc) .+ F_SWvdf .* ((1.0.-α_vdf_i).*F_SWpen .+ (1.0.-α_vdf_s).*F_SWsfc)
            .+  F_SWidr .* ((1.0.-α_idr_i).*F_SWpen .+ (1.0.-α_idr_s).*F_SWsfc) .+ F_SWidf .* ((1.0.-α_idf_i).*F_SWpen .+ (1.0.-α_idf_s).*F_SWsfc))

    F_SWpen .*= (F_SWvdr .* (1.0.-α_vdr_i) .+ F_SWvdf .* (1.0.-α_vdf_i)) .* i0vis # F_SWpenvdr + F_SWpenvdr
    # This is the SW surface flux, which is added to F_0 (F_SWsfc = F_SWabs .- F_SWpen):
    F_SWsfc .-= F_SWpen

    # Transmittance at the top of the ice:
    # let trantop be c_u, tranbot be c_q:
    c_u .= 1.0
    c_q .= 0.0
    for k in 1:N_i
        c_q              .= exp.((-κ_i*k/N_i).*H_i)
        I_pen[k:N_i:end] .= F_SWpen .* (c_u .- c_q)
        c_u              .= c_q
    end

    # TODO: want to compute dF_SWsfc. Can we use enzyme for this?
    
    
    # Now compute total surface flux:
    #println(F_SWsfc)
    #println(F_Ld + F_Lu[step] + F_l[step] + F_s[step])
    F_0 .= F_SWsfc .+ F_Ld .+ F_Lu .+ F_l .+ F_s
    # And now compute derivative of flux:
    dF_0 .= dF_Lu .+ dF_s .+ dF_l

    return nothing
end

# Sets helper values needed to compute flux
@inline function set_atm_flux_values(N_c, N_layers, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l, c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, ρ_a, c_p, Q_sfc, T_n, H_i, onGPU)

    # Compute initial exchange coefficients:
    c_u .= κ ./ log.(z_ref ./ H_i)
    copyto!(c_Θ, c_u)
    copyto!(c_q, c_u)
    
    Q_sfc .= (q_1./ρ_a).*exp.(-q_2 ./ (T_n[begin:N_layers:end] .+ C_to_K))
    
    if onGPU
        for i in 1:5
            @cuda set_atm_helper_values_step(N_c, N_layers, c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, T_n, Q_sfc)
        end
    else
        for i in 1:5
            set_atm_helper_values_step(N_c, N_layers, c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, T_n, Q_sfc)
        end
    end
    
    # Compute heatflux coefficients. TODO: wind stress?
    dF_l .= (L_vap + L_ice) .* ρ_a .* atm_u_star .* c_q #C_l
    dF_s .= ρ_a .* c_p .* atm_u_star .* c_Θ .+ 1 #C_s

    # Outgoing longwave flux:
    F_Lu  .= (emissivity*sbc) .* (T_n[begin:N_layers:end] .+ C_to_K).^4
    dF_Lu .= (4.0emissivity*sbc) .* (T_n[begin:N_layers:end] .+ C_to_K).^3

    # Sensible heat flux:
    F_s   .= dF_s .* (Θ_a .- (T_n[begin:N_layers:end] .+ C_to_K))
    dF_s .*= -1.0

    
    # Latent heat flux:
    F_l  .= dF_l .* (Q_a .- Q_sfc)
    dF_l .= -dF_l .* Q_sfc .* (-q_2 ./ (T_n[begin:N_layers:end] .+ C_to_K).^2) #-C_l * dQ_sf
    
    return nothing
end

# Performs one step in the iteration for set_atm_helper_values
@inline function set_atm_helper_values_step(N_c, N_layers, c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, T_n, Q_sfc)

    for index in 1:N_c
        # Update turbulent scales
        U_a_index = (index - 1) * 3
        atm_u_star[index] = c_u[index] * max(U_dmin, (U_a[U_a_index+1]^2+U_a[U_a_index+2]^2+U_a[U_a_index+3]^2)^.5)

        Θ_star = c_Θ[index] * (Θ_a[index] - T_n[(index-1)*N_layers + 1])
        Q_star = c_q[index] * (Q_a[index] - Q_sfc[index])

        # Update Y and compute χ
        Y = (κ*g*z_deg)*(Θ_star / (Θ_a[index]*(1+0.606Q_a[index])) + Q_star / (1.0/0.606 + Q_a[index]))/(atm_u_star[index]^2)

        # Compute ψ_m and ψ_s TODO: add alternative if Y is unstable? (< 0)
        ψ_m = -(0.7Y + 0.75*(Y-14.3)*exp(-0.35Y) + 10.7)
        ψ_s = ψ_m

        if Y < 0.0
            # try to simplify this, might be an issue with derivative
            χ = (1 - 16.0Y)^0.25
            ψ_m = 2*log(0.5*(1+χ)) + log(0.5*(1+χ^2)) - 2*atan(χ) + π/2
            ψ_s = 2*log(0.5*(1+χ^2))
        end

        # Update exchange coefficients
        c_u[index] /= 1 + c_u[index]*(λ - ψ_m)/κ
        c_Θ[index] /= 1 + c_Θ[index]*(λ - ψ_s)/κ
        c_q[index]  = c_Θ[index]
    end

    return nothing
end

# Computes the albedo for this column's ice, assumed to be constant throughout
# model run (at least for now)
@inline function generate_α_i(N_c, N_layers, H_i, α_vdr_i, α_idr_i, α_vdf_i, α_idf_i, T_n)

    for index in 1:N_c

        # Get the albedo values for bare ice:
        fh       = min(atan(4.0H_i[index])/atan(4.0ahmax), 1.0)
        albo     = α_o*(1.0-fh)
        α_vdf_i[index] = α_icev*fh + albo
        α_idf_i[index] = α_icei*fh + albo

        # Temperature dependence component:
        dTs       = -T_n[(index-1)*N_layers + 1]
        fT        = min(dTs - 1.0, 0.0)
        α_vdf_i[index] -= dα_mlt*fT
        α_idf_i[index] -= dα_mlt*fT

        # Prevent negative albedos:
        α_vdf_i[index] = max(α_vdf_i[index], α_o)
        α_idf_i[index] = max(α_idf_i[index], α_o)

        # Direct albedos (same as diffuse albedos for now)
        α_vdr_i[index]  = α_vdf_i[index]
        α_idr_i[index]  = α_idf_i[index]
    end

    return nothing
end

# Computes the albedo for this column's snow, assumed to be constant throughout
# model run (at least for now)
@inline function generate_α_s(N_c, N_layers, α_vdr_s, α_idr_s, α_vdf_s, α_idf_s, T_n)
    
    for index in 1:N_c

        # Temperature dependence component:
        dTs       = -T_n[(index-1)*N_layers + 1]
        fT        = min(dTs - 1.0, 0.0)

        # Snow albedo:
        α_vdf_s[index] = α_snowv - dα_mltv*fT
        α_idf_s[index] = α_snowi - dα_mlti*fT

        # Direct albedos (same as diffuse albedos for now)
        α_vdr_s[index] = α_vdf_s[index]
        α_idr_s[index] = α_idf_s[index]
    end

    return nothing
end

# Computes all albedos using array programming
@inline function generate_α(N_layers, H_i, α_vdr_i, α_idr_i, α_vdf_i, α_idf_i, α_vdr_s, α_idr_s, α_vdf_s, α_idf_s, T_n)

    # Get the albedo values for bare ice:
    α_vdf_i .= min.(atan.((4.0/atan(4.0ahmax)) .* H_i), 1.0) #fh

    α_idf_i .= (α_icei - α_o).*α_vdf_i .+ α_o
    α_vdf_i .= (α_icev - α_o).*α_vdf_i .+ α_o

    # Temperature dependence component:
    α_vdf_s .= -T_n[begin:N_layers:end]  # dTs
    α_vdf_s .= min.(α_vdf_s .- 1.0, 0.0) # fT

    α_vdf_i .-= dα_mlt.*α_vdf_s
    α_idf_i .-= dα_mlt.*α_vdf_s

    # Snow albedos
    α_idf_s .= α_snowi .- dα_mlti.*α_vdf_s
    α_vdf_s .= α_snowv .- dα_mltv.*α_vdf_s

    # Prevent negative ice albedos:
    α_vdf_i .= max.(α_vdf_i, α_o)
    α_idf_i .= max.(α_idf_i, α_o)

    # Direct albedos (same as diffuse albedos for now)
    copyto!(α_vdr_i, α_vdf_i)
    copyto!(α_idr_i, α_idf_i)
    copyto!(α_vdr_s, α_vdf_s)
    copyto!(α_idr_s, α_idf_s)

end