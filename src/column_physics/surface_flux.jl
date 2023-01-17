# Written by Joseph Kump (josek97@utexas.edu)
# Contains surface heat flux functions for running the sea ice model, featuring
# the ATModel struct

# Computes the (constant) atmospheric flux affecting the model
@inline function step_surface_flux(jcolumn, atmodel, step)

    # Compute atmospheric fluxes dependent on ice:
    set_atm_flux_values(jcolumn, atmodel, jcolumn.H_i, step)

    # Get the fractional snow cover:
    area_snow = 0.0
    if jcolumn.H_s > puny
        area_snow = jcolumn.H_s / (jcolumn.H_s+0.02)
    end
    area_no_snow = 1.0 - area_snow
    
    # Reduce shortwave fluxes with albedo, getting absorbed shortwave fluxes.
    # Note that the albedo terms include ice and snow components:
    F_SWabsv = atmodel.F_SWvdr * ((1.0-jcolumn.α_vdr[1])*area_no_snow + (1.0-jcolumn.α_vdr[2])*area_snow) + atmodel.F_SWvdf * ((1.0-jcolumn.α_vdf[1])*area_no_snow + (1.0-jcolumn.α_vdf[2])*area_snow)
    F_SWabsi = atmodel.F_SWidr * ((1.0-jcolumn.α_idr[1])*area_no_snow + (1.0-jcolumn.α_idr[2])*area_snow) + atmodel.F_SWidf * ((1.0-jcolumn.α_idf[1])*area_no_snow + (1.0-jcolumn.α_idf[2])*area_snow)
    F_SWabs  = F_SWabsv + F_SWabsi

    F_SWpenvdr = atmodel.F_SWvdr * (1.0-jcolumn.α_vdr[1]) * area_no_snow * i0vis
    F_SWpenvdf = atmodel.F_SWvdf * (1.0-jcolumn.α_vdf[1]) * area_no_snow * i0vis

    # The total penetrating SW radiation, only visible flux that enters the ice penetrates
    F_SWpen = F_SWpenvdr + F_SWpenvdf
    # This is the SW surface flux, which is added to F_0:
    F_SWsfc = F_SWabs - F_SWpen
    
    # Transmittance at the top of the ice:
    trantop = 1.0
    tranbot = 0.0
    for k in 1:jcolumn.N_i
        tranbot          = exp(-κ_i*(jcolumn.H_i/jcolumn.N_i)*k)
        jcolumn.I_pen[k] = F_SWpen * (trantop - tranbot)
        trantop          = tranbot
    end

    # TODO: want to compute dF_SWsfc. Can we use enzyme for this?
    

    # Now compute total surface flux:
    #println(F_SWsfc)
    #println(F_Ld + F_Lu[step] + F_l[step] + F_s[step])
    jcolumn.F_0[step] = F_SWsfc + atmodel.F_Ld + jcolumn.F_Lu[step] + jcolumn.F_l[step] + jcolumn.F_s[step]
    # And now compute derivative of flux:
    jcolumn.dF_0[step] = jcolumn.dF_Lu[step] + jcolumn.dF_s[step] + jcolumn.dF_l[step]

    return nothing
end

# Sets helper values needed to compute flux
@inline function set_atm_flux_values(jcolumn, atmodel, z_ice, step)

    T_sfc = jcolumn.T_n[1]

    # Compute initial exchange coefficients:
    atmodel.c_u[1] = κ / log(z_ref / z_ice)
    atmodel.c_Θ[1] = atmodel.c_u[1]
    atmodel.c_q[1] = atmodel.c_u[1]

    Q_sfc = (q_1/atmodel.ρ_a)*exp(-q_2 / (T_sfc + C_to_K))

    # Iterate and update exchange coefficients:
    #=
    for k in 1:natmiter
        set_atm_helper_values_step(c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star, T_sfc, Q_sfc)
    end
    =#
    # Approach for work with AD
    set_atm_helper_values_step(atmodel, T_sfc, Q_sfc)
    set_atm_helper_values_step(atmodel, T_sfc, Q_sfc)
    set_atm_helper_values_step(atmodel, T_sfc, Q_sfc)
    set_atm_helper_values_step(atmodel, T_sfc, Q_sfc)
    set_atm_helper_values_step(atmodel, T_sfc, Q_sfc)

    # Compute heatflux coefficients. TODO: wind stress?
    C_l = atmodel.ρ_a * (L_vap + L_ice) * atmodel.atm_u_star[1] * atmodel.c_q[1]
    C_s = atmodel.ρ_a * atmodel.c_p     * atmodel.atm_u_star[1] * atmodel.c_Θ[1] + 1

    # Outgoing longwave flux:
    jcolumn.F_Lu[step]  = emissivity * sbc * (T_sfc + C_to_K)^4
    jcolumn.dF_Lu[step] = 4.0 * emissivity * sbc * (T_sfc + C_to_K)^3

    # Sensible heat flux:
    jcolumn.F_s[step]  = C_s * (atmodel.Θ_a[1] - (T_sfc + C_to_K))
    jcolumn.dF_s[step] = -C_s

    # Latent heat flux:
    Q_sf       = (q_1/atmodel.ρ_a) * exp(-q_2 / (T_sfc + C_to_K))
    dQ_sf      = Q_sf * (-q_2 / (T_sfc + C_to_K)^2)
    jcolumn.F_l[step]  = C_l  * (atmodel.Q_a[1] - Q_sf)
    jcolumn.dF_l[step] = -C_l * dQ_sf

    return nothing
end

llpow(x::Float64, y::Float64) = Base.llvmcall(("""
declare double @llvm.pow.f64(double, double)
define double @f(double %x, double %y) {
entry:
  %z = call double @llvm.pow.f64(double %x, double %y)
  ret double %z
}
""","f"), Float64, Tuple{Float64, Float64}, x, y)

# Performs one step in the iteration for set_atm_helper_values
@inline function set_atm_helper_values_step(atmodel, T_sfc, Q_sfc)

    # Update turbulent scales
    atmodel.atm_u_star[1] = atmodel.c_u[1] * max(U_dmin, (atmodel.U_a[1]^2+atmodel.U_a[2]^2+atmodel.U_a[3]^2)^.5)
    Θ_star        = atmodel.c_Θ[1] * (atmodel.Θ_a[1] - T_sfc)
    Q_star        = atmodel.c_q[1] * (atmodel.Q_a[1] - Q_sfc)

    # Update Y and compute χ
    Y = (κ*g*z_deg)*(Θ_star / (atmodel.Θ_a[1]*(1+0.606atmodel.Q_a[1])) + Q_star / (1.0/0.606 + atmodel.Q_a[1]))/(atmodel.atm_u_star[1]^2)

    # Compute ψ_m and ψ_s TODO: add alternative if Y is unstable? (< 0)
    ψ_m = -(0.7Y + 0.75*(Y-14.3)*exp(-0.35Y) + 10.7)
    ψ_s = ψ_m

    if Y < 0.0
        # try to simplify this, might be an issue with derivative
        χ   = llpow(1 - 16Y, 0.25)
        ψ_m = 2*log(0.5*(1+χ)) + log(0.5*(1+χ^2)) - 2*atan(χ) + π/2
        ψ_s = 2*log(0.5*(1+χ^2))
    end

    # Update exchange coefficients
    atmodel.c_u[1] = atmodel.c_u[1] / (1 + atmodel.c_u[1]*(λ - ψ_m)/κ)
    atmodel.c_Θ[1] = atmodel.c_Θ[1] / (1 + atmodel.c_Θ[1]*(λ - ψ_s)/κ)
    atmodel.c_q[1] = atmodel.c_Θ[1]

    return nothing
end

# Computes the albedo for this column's ice and snow, assumed to be constant throughout
# model run (at least for now)
@inline function generate_α(jcolumn)

    # Get the albedo values for bare ice:
    fh       = min(atan(4.0jcolumn.H_i)/atan(4.0ahmax), 1.0)
    albo     = α_o*(1.0-fh)

    jcolumn.α_vdf[1] = α_icev*fh + albo
    jcolumn.α_idf[1] = α_icei*fh + albo

    # Temperature dependence component:
    dTs               = -jcolumn.T_n[1]
    fT                = min(dTs - 1.0, 0.0)
    jcolumn.α_vdf[1] -= dα_mlt*fT
    jcolumn.α_idf[1] -= dα_mlt*fT

    # Prevent negative albedos:
    jcolumn.α_vdf[1] = max(jcolumn.α_vdf[1], α_o)
    jcolumn.α_idf[1] = max(jcolumn.α_idf[1], α_o)

    # Snow albedo:
    jcolumn.α_vdf[2] = α_snowv - dα_mltv*fT
    jcolumn.α_idf[2] = α_snowi - dα_mlti*fT

    # Direct albedos (same as diffuse albedos for now)
    jcolumn.α_vdr[1] = jcolumn.α_vdf[1]
    jcolumn.α_idr[1] = jcolumn.α_idf[1]
    jcolumn.α_vdr[2] = jcolumn.α_vdf[2]
    jcolumn.α_idr[2] = jcolumn.α_idf[2]

    return nothing
end
