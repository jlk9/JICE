
include("../gpu_src/gpu_atmosphere_model/gpu_atmodel_struct.jl")
include("../gpu_src/gpu_column_physics/gpu_run_ice_column.jl")

using BenchmarkTools, CUDA

N_t    = 10
N_c    = 2
N_i    = 5
N_s    = 2
H_i    = 2.0 .+ zeros(Float64, N_c)
H_s    = 0.10 .+ zeros(Float64, N_c)
T_frz  = 271.35 - 273.15
T_0    = repeat([-21.0, -20.5, -20.0, -19.0, -14.5, -10.0, -5.5, -1.0], N_c)
Δt     = 1.0
u_star = 0.0005 # recommended minimum value of u_star in CICE
T_w    = (274.47 - 273.15) .+ zeros(Float64, N_c) # typical temp in C for sea surface in arctic

# "Educated" guess for some normal atmospheric values
F_Ld    = 10.0 .+ zeros(Float64, N_c)
F_SWvdr = 120.0 .+ zeros(Float64, N_c)
F_SWidr = 0.0 .+ zeros(Float64, N_c)
F_SWvdf = 0.0 .+ zeros(Float64, N_c)
F_SWidf = 0.0 .+ zeros(Float64, N_c)
T_a     = -34.0 .+ zeros(Float64, N_c)
Θ_a     = (T_a*(1000.0/1045.6)^0.286) .+ zeros(Float64, N_c)
ρ_a     = 1.4224 .+ zeros(Float64, N_c)
Q_a     = 0.005 .+ zeros(Float64, N_c) #?
c_p     = 0.7171 .+ zeros(Float64, N_c)
U_a     = zeros(Float64, 3*N_c)

onDevice = CUDA.has_cuda()

atmodels = initialize_ATModelArrays(N_t, N_c, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onDevice)
jarrays  = initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onDevice)

#@btime run_ice_column(jarrays, atmodels)

run_ice_column(jarrays, atmodels)

#=
println(jarrays.S)
println(jarrays.Δh)
println(jarrays.K)
println(jarrays.c_i)
println(jarrays.maindiag)
println(jarrays.subdiag)
println(jarrays.supdiag)
=#
println(jarrays.T_nplus)

#=
println(jarrays.α_vdr_i)
println(jarrays.α_idr_i)
println(jarrays.α_vdf_i)
println(jarrays.α_idf_i)

println(jarrays.α_vdr_s)
println(jarrays.α_idr_s)
println(jarrays.α_vdf_s)
println(jarrays.α_idf_s)
=#
#=
println(atmodels.c_u)
println(atmodels.c_Θ)
println(atmodels.c_q)
println(atmodels.Q_sfc)
=#
#=
println(jarrays.F_Lu)
println(jarrays.F_s)
println(jarrays.F_l)
println(jarrays.dF_Lu)
println(jarrays.dF_s)
println(jarrays.dF_l)
println(jarrays.I_pen)
println(jarrays.F_0)
println(jarrays.dF_0)
=#