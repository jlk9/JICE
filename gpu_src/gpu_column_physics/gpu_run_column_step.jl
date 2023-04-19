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

    
    
    # Computes the surface fluxes at this time step
    step_surface_flux(jarrays.N_c, jarrays.N_i, jarrays.N_layers, jarrays.α_vdr_i, jarrays.α_idr_i, jarrays.α_vdf_i, jarrays.α_idf_i, jarrays.α_vdr_s, jarrays.α_idr_s, jarrays.α_vdf_s, jarrays.α_idf_s,
                      jarrays.T_n, jarrays.H_i, jarrays.H_s,
                      jarrays.F_0, jarrays.dF_0, jarrays.F_Lu, jarrays.F_s, jarrays.F_l, jarrays.dF_Lu, jarrays.dF_s, jarrays.dF_l,
                      atmodels.F_SWvdr, atmodels.F_SWidr, atmodels.F_SWvdf, atmodels.F_SWidf, atmodels.F_Ld, jarrays.I_pen,
                      atmodels.c_u, atmodels.c_Θ, atmodels.c_q, atmodels.U_a, atmodels.Θ_a, atmodels.Q_a, atmodels.atm_u_star,
                      atmodels.ρ_a, atmodels.c_p, atmodels.Q_sfc, atmodels.F_SWsfc, atmodels.F_SWpen, jarrays.onGPU)

    # Computes the temperature changes at this step
    
    step_temp_change(jarrays.N_c, jarrays.N_i, jarrays.N_s, jarrays.N_layers, jarrays.H_s, jarrays.S, jarrays.T_frz, jarrays.Δh, jarrays.T_n,
                     jarrays.T_nplus, jarrays.c_i, jarrays.K, jarrays.K̄, jarrays.I_pen, jarrays.F_0, jarrays.dF_0,
                     jarrays.maindiag, jarrays.subdiag, jarrays.supdiag, jarrays.Δt, jarrays.onGPU)
    
    #=
    # Gets the growth/melt and rebalances
    step_growth_melt(N_i, N_s, H_s[1], S, T_frz, Δh, T_nplus, K, K̄,
                      q, q_new, z_old, z_new, Δt, u_star, T_w, F_0[step])

    # Add up the layer thicknesses to get the new total thickness
    readd_total_thickness(N_i, N_s, Δh, H_i, H_iold, H_s)
    =#
    # Update T_n
    jarrays.T_n .= jarrays.T_nplus
    
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