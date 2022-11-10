# Written by Joseph Kump (josek97@utexas.edu)
# Tests the single grid cell functionality

include("../src/grid_cell/run_grid_cell.jl")

function test_cell_run(N_t, Δt, T_frz, T_w, N_i_cols, N_s_cols, H_i_cols, H_s_cols, i_0_cols, κ_i_cols, u_star_cols, T_0_cols, F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, areas)

    jcell = initialize_JICECell(N_t, Δt, T_frz, T_w, N_i_cols, N_s_cols, H_i_cols, H_s_cols, i_0_cols, κ_i_cols, u_star_cols, T_0_cols, F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, areas)

    run_ice_cell(jcell)

    println("Resulting temps for 3-layer column are:")
    println(jcell.columns[1].T_array[:,N_t+1])
    println("Resulting temps for 5-layer column are:")
    println(jcell.columns[2].T_array[:,N_t+1])

end

N_t         = 4
N_i_cols    = [3, 5]
N_s_cols    = [0, 0]
H_i_cols    = [1.1, 2.0]
H_s_cols    = [0.0, 0.0]
T_frz       = 271.35 - 273.15
i_0_cols    = [0.7, 0.7]
κ_i_cols    = [1.4, 1.4]
T_0_cols    = [0 .- [14.75, 12.5, 8, 3.5], 0 .- [21.25, 19.0, 14.5, 10.0, 5.5, 1.0]]
Δt          = 1.0
u_star_cols = [0.0005, 0.0005] # recommended minimum value of u_star in CICE
T_w         = 274.47 - 273.15 # typical temp in C for sea surface in arctic
areas       = [0.2, 0.8]

# "Educated" guess for some normal atmospheric values
F_Ld  = 10.0
F_sw  = 120.0 
T_a   = -34.0
Θ_a   = T_a*(1000.0/1045.6)^0.286
ρ_a   = 1.4224
Q_a   = 0.005 #?
c_p   = 0.7171
U_a   = zeros(Float64, 3)

println("Testing cell run...")
test_cell_run(N_t, Δt, T_frz, T_w, N_i_cols, N_s_cols, H_i_cols, H_s_cols, i_0_cols, κ_i_cols, u_star_cols, T_0_cols, F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, areas)
