
include("../gpu_src/gpu_column_physics/gpu_jicecolumn_struct.jl")

N_t    = 100
N_c    = 10
N_i    = 5
N_s    = 2
H_i    = 2.0 .+ zeros(Float64, N_c)
H_s    = 0.10 .+ zeros(Float64, N_c)
T_frz  = 271.35 - 273.15
T_0    = repeat([-21.0, -20.5, -20.0, -19.0, -14.5, -10.0, -5.5, -1.0], N_c)
Δt     = 1.0
u_star = 0.0005 # recommended minimum value of u_star in CICE
T_w    = (274.47 - 273.15) .+ zeros(Float64, N_c) # typical temp in C for sea surface in arctic

jarrays = initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0)

println(jarrays)