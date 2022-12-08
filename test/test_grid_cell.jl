# Written by Joseph Kump (josek97@utexas.edu)
# Tests the single grid cell functionality

include("../src/grid_cell/run_grid_cell.jl")

using Printf

function test_cell_run(N_cat, N_t, Δt, T_frz, T_w, N_i_cols, N_s_cols, H_i_cols, H_s_cols, u_star_cols, T_0_cols, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, thickness_bds, areas)

    jcell = initialize_JICECell(N_cat, N_t, Δt, T_frz, T_w, N_i_cols, N_s_cols, H_i_cols, H_s_cols, u_star_cols, T_0_cols, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, thickness_bds, areas)

    run_ice_cell(jcell)

    for k in 1:N_cat
        @printf("Initial and final temps for layer %d are:\n", k)
        println(jcell.columns[k].T_array[:,1])
        println(jcell.columns[k].T_array[:,N_t+1])

        @printf("Initial and final thicknesses for layer %d are:\n", k)
        println(jcell.columns[k].Δh_array[:,1])
        println(jcell.columns[k].Δh_array[:,N_t+1])

    end
    @time run_ice_cell(jcell)
end

N_cat         = 5
N_t           = 4
N_i_cols      = [1, 2, 3, 4, 5]
N_s_cols      = [0, 0, 0, 1, 2]
H_i_cols      = [0.15, 0.45, 1.0, 1.6, 2.1]
H_s_cols      = [0.0, 0.0, 0.0, 0.05, 0.09]
T_frz         = 271.35 - 273.15
T_0_cols      = [0 .- [7.0, 6.0],
                 0 .- [9.0, 8.0, 4.0],
                 0 .- [13.75, 12.5, 8.0, 3.5],
                 0 .- [20.0, 19.0, 14.5, 10.0, 5.5, 1.0],
                 0 .- [21.0, 20.5, 20.0, 19.0, 14.5, 10.0, 5.5, 1.0]]

Δt            = 1.0
u_star_cols   = 0.0005 .+ zeros(Float64, N_cat) # recommended minimum value of u_star in CICE
T_w           = 274.47 - 273.15 # typical temp in C for sea surface in arctic
thickness_bds = [0.0, 0.3, 0.7, 1.2, 2.0, 999.9]
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

println("Testing cell run...")
test_cell_run(N_cat, N_t, Δt, T_frz, T_w, N_i_cols, N_s_cols, H_i_cols, H_s_cols, u_star_cols, T_0_cols,
              F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, thickness_bds, areas)
