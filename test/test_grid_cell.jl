# Written by Joseph Kump (josek97@utexas.edu)
# Tests the single grid cell functionality

include("../src/grid_cell/run_grid_cell_adjoint.jl")

using Printf

function test_cell_run(N_cat, N_t, Δt, T_frz, T_w, frzmlt, rside, N_i, N_s, H_i_cols, H_s_cols, u_star, T_0_cols, F_SWvdr, F_SWidr,
                        F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, thickness_bds, areas)

    jcell = initialize_JICECell(N_cat, N_t, Δt, T_frz, T_w, frzmlt, rside, N_i, N_s, H_i_cols, H_s_cols, u_star, T_0_cols,
                                F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, thickness_bds, areas)

    run_ice_cell(jcell)

    for k in 1:N_cat
        
        @printf("Initial and final temps for category %d are:\n", k)
        println(jcell.columns[k].T_array[:,1])
        println(jcell.columns[k].T_array[:,N_t+1])
        
        @printf("Initial and final thicknesses for category %d are:\n", k)
        println(jcell.columns[k].Δh_array[:,1])
        println(jcell.columns[k].Δh_array[:,N_t+1])
        
        @printf("Initial and final total thickness for category %d is:\n", k)
        println(sum(jcell.columns[k].Δh_array[:,1]))
        println(sum(jcell.columns[k].H_i_array[N_t+1] + jcell.columns[k].H_s_array[N_t+1]))
        
    end

    if jcell.i_vol_itd_change
        println("Oh no! The total change in ice volume from horizontal transport is above acceptable machine precision.")
    end
    if jcell.s_vol_itd_change
        println("Oh no! The total change in snow volume from horizontal transport is above acceptable machine precision.")
    end
    if jcell.i_energy_itd_change
        println("Oh no! The total change in ice energy from horizontal transport is above acceptable machine precision.")
    end
    if jcell.s_energy_itd_change
        println("Oh no! The total change in snow energy from horizontal transport is above acceptable machine precision.")
    end

    if jcell.i_vol_new_change
        println("Oh no! The total change in ice volume from adding new ice (when factoring in added component) is above acceptable machine precision.")
    end
    if jcell.i_energy_new_change
        println("Oh no! The total change in ice energy from adding new ice (when factoring in added component) is above acceptable machine precision.")
    end

    println("Category areas are:")
    println(jcell.areas)

    @time run_ice_cell(jcell)
end

function test_adjoint(N_cat, N_t, Δt, T_frz, T_w, frzmlt, rside, N_i, N_s, H_i_cols, H_s_cols, u_star, T_0_cols, F_SWvdr, F_SWidr,
                        F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, thickness_bds, areas)

    jcell = initialize_JICECell(N_cat, N_t, Δt, T_frz, T_w, frzmlt, rside, N_i, N_s, deepcopy(H_i_cols), deepcopy(H_s_cols), u_star, deepcopy(T_0_cols),
                                F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, deepcopy(thickness_bds), deepcopy(areas))

    ad_H_i_cols = [1.0, 0.0, 0.0, 0.0, 0.0]
    ad_T_cols = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]


    @time d_jcell = run_grid_cell_autodiff(jcell, ad_H_i_cols, ad_T_cols)

    println("Partial derivatives (∂H_i of first column) / (∂T of first column):")
    println(d_jcell.columns[1].T_n)
    println("In particular, the derivative of the first column thickness in terms of the second-from-bottom ice temperature is:")
    println(d_jcell.columns[1].T_n[N_i+N_s])

    println("To test this, we'll see the result we get from finite differences. The ϵ values are:")

    step_sizes  = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13]
    diffs       = []
    rel_errors  = []

    #println(T_0)

    #
    for ϵ in step_sizes

        # Needs to be the same as our original T_0
        T_ϵp              = deepcopy(T_0_cols)
        T_ϵp[1][N_i+N_s] += ϵ
        jcellp            = initialize_JICECell(N_cat, N_t, Δt, T_frz, T_w, frzmlt, rside, N_i, N_s, deepcopy(H_i_cols), deepcopy(H_s_cols), u_star, T_ϵp,
                                                F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, deepcopy(thickness_bds), deepcopy(areas))
        run_ice_cell(jcellp)
        
        H_ϵp_value  = jcellp.columns[1].H_i_array[N_t+1]

        T_ϵn              = deepcopy(T_0_cols)
        T_ϵn[1][N_i+N_s] -= ϵ
        jcelln            = initialize_JICECell(N_cat, N_t, Δt, T_frz, T_w, frzmlt, rside, N_i, N_s, deepcopy(H_i_cols), deepcopy(H_s_cols), u_star, T_ϵn,
                                                F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, deepcopy(thickness_bds), deepcopy(areas))
        run_ice_cell(jcelln)

        H_ϵn_value  = jcelln.columns[1].H_i_array[N_t+1]

        diff = (H_ϵp_value - H_ϵn_value) / (2*ϵ)

        error = abs(d_jcell.columns[1].T_n[N_i+N_s] - diff) / abs(d_jcell.columns[1].T_n[N_i+N_s])

        push!(diffs, diff)
        push!(rel_errors, error)

    end

    println(step_sizes)
    println("The finite diffs are:")
    println(diffs)
    println("and the errors are:")
    println(rel_errors)

    @time run_grid_cell_autodiff(jcell, ad_H_i_cols, ad_T_cols)
    
end

N_cat = 5
N_t   = 600
N_i   = 5
N_s   = 2
H_i_cols      = [0.45, 1.0, 1.85, 3.57, 5.1]
#H_i_cols      = [0.15, 0.45, 1.0, 1.6, 2.1] #WMO Case
#H_s_cols      = [0.0, 0.0, 0.01, 0.05, 0.09]
#H_s_cols      = [0.0, 0.0, 0.0, 0.0, 0.0]
H_s_cols      = [0.01, 0.02, 0.05, 0.07, 0.09]
T_frz         = 271.35 - 273.15
T_0_cols      = [0.0 .- [20.8, 20.5, 20.0, 19.0, 14.5, 10.0, 5.5, 1.0],
                 0.0 .- [20.9, 20.6, 20.1, 19.1, 14.6, 10.1, 5.6, 1.1],
                 0.0 .- [21.0, 20.7, 20.2, 19.2, 14.7, 10.2, 5.7, 1.2],
                 0.0 .- [21.1, 20.8, 20.3, 19.3, 14.8, 10.3, 5.8, 1.3],
                 0.0 .- [21.2, 20.9, 20.4, 19.4, 14.9, 10.4, 5.9, 1.4]]

Δt            = 1.0
u_star        = 0.0005 # recommended minimum value of u_star in CICE
T_w           = 274.47 - 273.15 # typical temp in C for sea surface in arctic
thickness_bds = [0.0, 0.64, 1.39, 2.47, 4.57, 999.9]
# thickness_bds = [0.0, 0.3, 0.7, 1.2, 2.0, 999.9] # WMO Case
areas         = [0.2, 0.1, 0.1, 0.1, 0.2] # Area of open water is 1 - sum(areas), in this case 0.3

# "Educated" guess for some normal atmospheric values
F_Ld    = 10.0
F_SWvdr = 120.0
F_SWidr = 0.0
F_SWvdf = 0.0
F_SWidf = 0.0
T_a     = -34.0
Θ_a     = T_a*(1000.0/1045.6)^0.286
ρ_a     = 1.4224
Q_a     = 0.005 #?
c_p     = 0.7171
U_a     = zeros(Float64, 3)

# Heat transfer coefficient of water is 598.5 W/ m^2 K.
frzmlt = 0.0 #598.5*(T_w - -1.0)
# fraction of ice that melts laterally
rside  = 0.000

println("Testing cell run...")
test_cell_run(N_cat, N_t, Δt, T_frz, T_w, frzmlt, rside, N_i, N_s, H_i_cols, H_s_cols, u_star, T_0_cols,
              F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, thickness_bds, areas)
#=
println("Testing adjoint...")
test_adjoint(N_cat, N_t, Δt, T_frz, T_w, frzmlt, rside, N_i, N_s, H_i_cols, H_s_cols, u_star, T_0_cols,
             F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, thickness_bds, areas)
=#