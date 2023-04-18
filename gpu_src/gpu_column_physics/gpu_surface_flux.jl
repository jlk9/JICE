# Written by Joseph Kump (josek97@utexas.edu)
# Contains surface heat flux functions for running the sea ice model, featuring
# the ATModel struct

# Computes the (constant) atmospheric flux affecting the model
@inline function step_surface_flux(N_i, α_vdr, α_idr, α_vdf, α_idf, T_sfc, H_i, H_s, F_0, dF_0, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l,
                                    F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, I_pen, c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, ρ_a, c_p, step)

    # Compute atmospheric fluxes dependent on ice:
    set_atm_flux_values(F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l, c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, ρ_a, c_p, T_sfc, H_i, step)

    # Get the fractional snow cover:
    area_snow = 0.0
    if H_s > puny
        area_snow = H_s / (H_s+0.02)
    end
    area_no_snow = 1.0 - area_snow

    # Reduce shortwave fluxes with albedo, getting absorbed shortwave fluxes.
    # Note that the albedo terms include ice and snow components:
    F_SWabsv = F_SWvdr * ((1.0-α_vdr[1])*area_no_snow + (1.0-α_vdr[2])*area_snow) + F_SWvdf * ((1.0-α_vdf[1])*area_no_snow + (1.0-α_vdf[2])*area_snow)
    F_SWabsi = F_SWidr * ((1.0-α_idr[1])*area_no_snow + (1.0-α_idr[2])*area_snow) + F_SWidf * ((1.0-α_idf[1])*area_no_snow + (1.0-α_idf[2])*area_snow)
    F_SWabs  = F_SWabsv + F_SWabsi

    F_SWpenvdr = F_SWvdr * (1.0-α_vdr[1]) * area_no_snow * i0vis
    F_SWpenvdf = F_SWvdf * (1.0-α_vdf[1]) * area_no_snow * i0vis

    # The total penetrating SW radiation, only visible flux that enters the ice penetrates
    F_SWpen = F_SWpenvdr + F_SWpenvdf
    # This is the SW surface flux, which is added to F_0:
    F_SWsfc = F_SWabs - F_SWpen
    
    # Transmittance at the top of the ice:
    trantop = 1.0
    tranbot = 0.0
    for k in 1:N_i
        tranbot  = exp(-κ_i*(H_i/N_i)*k)
        I_pen[k] = F_SWpen * (trantop - tranbot)
        trantop  = tranbot
    end

    # TODO: want to compute dF_SWsfc. Can we use enzyme for this?
    

    # Now compute total surface flux:
    #println(F_SWsfc)
    #println(F_Ld + F_Lu[step] + F_l[step] + F_s[step])
    F_0[step] = F_SWsfc + F_Ld + F_Lu[step] + F_l[step] + F_s[step]
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
    for i in 1:5
        set_atm_helper_values_step(c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, T_sfc, Q_sfc)
    end

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
    Q_sf       = (q_1/ρ_a) * exp(-q_2 / (T_sfc + C_to_K))
    dQ_sf      = Q_sf * (-q_2 / (T_sfc + C_to_K)^2)
    F_l[step]  = C_l  * (Q_a[1] - Q_sf)
    dF_l[step] = -C_l * dQ_sf

    return nothing
end

#=
llpow(x::Float64, y::Float64) = Base.llvmcall(("""
declare double @llvm.pow.f64(double, double)
define double @f(double %x, double %y) {
entry:
  %z = call double @llvm.pow.f64(double %x, double %y)
  ret double %z
}
""","f"), Float64, Tuple{Float64, Float64}, x, y)
=#
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
        # try to simplify this, might be an issue with derivative
        #χ   = llpow(1 - 16Y, 0.25)
        χ = (1 - 16.0Y)^0.25
        ψ_m = 2*log(0.5*(1+χ)) + log(0.5*(1+χ^2)) - 2*atan(χ) + π/2
        ψ_s = 2*log(0.5*(1+χ^2))
    end

    # Update exchange coefficients
    c_u[1] = c_u[1] / (1 + c_u[1]*(λ - ψ_m)/κ)
    c_Θ[1] = c_Θ[1] / (1 + c_Θ[1]*(λ - ψ_s)/κ)
    c_q[1] = c_Θ[1]

    return nothing
end

# Computes the albedo for this column's ice and snow, assumed to be constant throughout
# model run (at least for now)
@inline function generate_α(H_i, α_vdr, α_idr, α_vdf, α_idf, T_sfc)

    # Get the albedo values for bare ice:
    fh       = min(atan(4.0H_i)/atan(4.0ahmax), 1.0)
    albo     = α_o*(1.0-fh)
    α_vdf[1] = α_icev*fh + albo
    α_idf[1] = α_icei*fh + albo

    # Temperature dependence component:
    dTs       = -T_sfc
    fT        = min(dTs - 1.0, 0.0)
    α_vdf[1] -= dα_mlt*fT
    α_idf[1] -= dα_mlt*fT

    # Prevent negative albedos:
    α_vdf[1] = max(α_vdf[1], α_o)
    α_idf[1] = max(α_idf[1], α_o)

    # Snow albedo:
    α_vdf[2] = α_snowv - dα_mltv*fT
    α_idf[2] = α_snowi - dα_mlti*fT

    # Direct albedos (same as diffuse albedos for now)
    α_vdr[1] = α_vdf[1]
    α_idr[1] = α_idf[1]
    α_vdr[2] = α_vdf[2]
    α_idr[2] = α_idf[2]

    return nothing
end
