# Written by Joseph Kump (josek97@utexas.edu)
# Basic code that is just ice temperature in column, with no formulation or melting.
# This version eliminates vectorized operations, enabling use of Enzyme.jl
# It also reduces the number of memory allocations substantially, increasing runtime
# and lowering memory footprint.

include("../gpu_atmosphere_model/gpu_atmodel_struct.jl")

include("./gpu_surface_flux.jl")
include("./gpu_temp_change.jl")
#include("./gpu_growth_melt.jl")

include("./gpu_jicecolumn_struct.jl")

#using CUDA

# Runs one step of ice process. 
@inline function run_column_step(jarrays::JICEColumnArrays, atmodels::ATModelArrays, step)

    
    # Computes the current albedo
    if jarrays.onGPU
        @cuda generate_α_i(jarrays.N_c, jarrays.N_layers, jarrays.H_i, jarrays.α_vdr_i, jarrays.α_idr_i, jarrays.α_vdf_i, jarrays.α_idf_i, jarrays.T_n)
        @cuda generate_α_s(jarrays.N_c, jarrays.N_layers, jarrays.α_vdr_s, jarrays.α_idr_s, jarrays.α_vdf_s, jarrays.α_idf_s, jarrays.T_n)
    else
        generate_α_i(jarrays.N_c, jarrays.N_layers, jarrays.H_i, jarrays.α_vdr_i, jarrays.α_idr_i, jarrays.α_vdf_i, jarrays.α_idf_i, jarrays.T_n)
        generate_α_s(jarrays.N_c, jarrays.N_layers, jarrays.α_vdr_s, jarrays.α_idr_s, jarrays.α_vdf_s, jarrays.α_idf_s, jarrays.T_n)
    end
    #=
    # Computes the surface fluxes at this time step
    step_surface_flux(N_i, α_vdr, α_idr, α_vdf, α_idf, T_n[1], H_i[1], H_s[1],
                         F_0, dF_0, F_Lu, F_s, F_l, dF_Lu, dF_s, dF_l,
                         F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, I_pen,
                         c_u, c_Θ, c_q, U_a, Θ_a, Q_a, atm_u_star,
                         ρ_a, c_p, step)
    =#
    # Computes the temperature changes at this step
    #=
    step_temp_change(jarrays.N_c, jarrays.N_i, jarrays.N_s, jarrays.N_layers, jarrays.H_s, jarrays.S, jarrays.T_frz, jarrays.Δh, jarrays.T_n,
                     jarrays.T_nplus, jarrays.c_i, jarrays.K, jarrays.K̄, jarrays.I_pen, jarrays.F_0, jarrays.dF_0,
                     jarrays.maindiag, jarrays.subdiag, jarrays.supdiag, jarrays.Δt, jarrays.onGPU, step)
    =#
    #=
    # Gets the growth/melt and rebalances
    step_growth_melt(N_i, N_s, H_s[1], S, T_frz, Δh, T_nplus, K, K̄,
                      q, q_new, z_old, z_new, Δt, u_star, T_w, F_0[step])

    # Add up the layer thicknesses to get the new total thickness
    readd_total_thickness(N_i, N_s, Δh, H_i, H_iold, H_s)
    =#
    #=
    # Update T_n
    for i in 1:(jarrays.N_c*jarrays.N_layers)
        jarrays.T_n[i] = jarrays.T_nplus[i]
    end
    =#
end

# Julia's sum() operation is memory inefficient for slices of arrays, this is oddly much faster
function readd_total_thickness(N_i, N_s, Δh, H_i, H_iold, H_s)

    H_iold[1] = H_i[1]
    H_i[1]    = 0.0
    H_s[1]    = 0.0

    for k in 1:(N_s+1)
        H_s[1] += Δh[k]
    end
    for k in (N_s+2):(N_i+N_s+1)
        H_i[1] += Δh[k]
    end
end