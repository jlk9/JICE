
include("../gpu_src/gpu_atmosphere_model/gpu_atmodel_struct.jl")
include("../gpu_src/gpu_column_physics/gpu_run_ice_column.jl")

using BenchmarkTools, CUDA

function test_model_run(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onDevice)

    @btime atmodels = initialize_ATModelArrays(N_t, N_c, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onDevice)
    @btime jarrays  = initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onDevice)

    @btime run_ice_column(jarrays, atmodels)

    return jarrays.T_nplus
end


N_t    = 1
N_c    = 100000
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
println("Performance for: initializing atmodel, initializing jarrays, running on 10^5 col for 1 time step:")
CUDA.@time atmodels = initialize_ATModelArrays(N_t, N_c, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onDevice)
CUDA.@time jarrays  = initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onDevice)
CUDA.@time run_ice_column(jarrays, atmodels)

println("Doing a second run as a sanity check for kernel compilation:")
CUDA.@time atmodels = initialize_ATModelArrays(N_t, N_c, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onDevice)
CUDA.@time jarrays  = initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onDevice)
CUDA.@time run_ice_column(jarrays, atmodels)

println("Now for the CPU comparison:")
onDevice = false
@time atmodels = initialize_ATModelArrays(N_t, N_c, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onDevice)
@time jarrays  = initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onDevice)
println(jarrays.onGPU)
@time run_ice_column(jarrays, atmodels)

println("Second run:")
onDevice = false
@time atmodels = initialize_ATModelArrays(N_t, N_c, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onDevice)
@time jarrays  = initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onDevice)
println(jarrays.onGPU)
@time run_ice_column(jarrays, atmodels)


# LONGER TIME
N_t    = 1000
N_c    = 100000
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
println("Performance for: initializing atmodel, initializing jarrays, running on 10^5 col for 10^3 time steps:")
CUDA.@time atmodels = initialize_ATModelArrays(N_t, N_c, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onDevice)
CUDA.@time jarrays  = initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onDevice)
CUDA.@time run_ice_column(jarrays, atmodels)

println("Second run to check effect of compiling:")
CUDA.@time atmodels = initialize_ATModelArrays(N_t, N_c, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onDevice)
CUDA.@time jarrays  = initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onDevice)
CUDA.@time run_ice_column(jarrays, atmodels)

println("Now for the CPU comparison:")
onDevice = false
@time atmodels = initialize_ATModelArrays(N_t, N_c, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onDevice)
@time jarrays  = initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onDevice)
println(jarrays.onGPU)
@time run_ice_column(jarrays, atmodels)

println("Second run:")
onDevice = false
@time atmodels = initialize_ATModelArrays(N_t, N_c, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onDevice)
@time jarrays  = initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onDevice)
println(jarrays.onGPU)
@time run_ice_column(jarrays, atmodels)


N_t    = 10000
N_c    = 100000
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
println("Performance for: initializing atmodel, initializing jarrays, running on 10^5 col for 10^4 time steps:")
CUDA.@time atmodels = initialize_ATModelArrays(N_t, N_c, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onDevice)
CUDA.@time jarrays  = initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onDevice)
CUDA.@time run_ice_column(jarrays, atmodels)

println("Second run to check effect of compiling:")
CUDA.@time atmodels = initialize_ATModelArrays(N_t, N_c, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onDevice)
CUDA.@time jarrays  = initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onDevice)
CUDA.@time run_ice_column(jarrays, atmodels)

println("Now for the CPU comparison:")
onDevice = false
@time atmodels = initialize_ATModelArrays(N_t, N_c, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onDevice)
@time jarrays  = initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onDevice)
println(jarrays.onGPU)
@time run_ice_column(jarrays, atmodels)

println("Second run:")
onDevice = false
@time atmodels = initialize_ATModelArrays(N_t, N_c, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, onDevice)
@time jarrays  = initialize_JICEColumnArrays(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, onDevice)
println(jarrays.onGPU)
@time run_ice_column(jarrays, atmodels)



#run_ice_column(jarrays, atmodels)

#println("First test correctness of CPU vs gpu model. On CPU:")
#=
cpu_temps = test_model_run(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, false)
gpu_temps = test_model_run(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, true)

println(cpu_temps)
println("Then on GPU:")
println(gpu_temps)
#println("Relative error is:")
#println(abs(cpu_temps - gpu_temps) / abs(cpu_temps))
println("For performance, in a problem with 1000 columns:")
N_c = 1000
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

println("CPU model time is:")
test_model_run(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, false)
println("GPU model time is:")
test_model_run(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, true)


println("For performance, in a problem with 100000 columns:")
N_c = 100000
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

println("CPU model time is:")
test_model_run(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, false)
println("GPU model time is:")

test_model_run(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, true)

#CUDA.@profile test_model_run(N_t, N_c, N_i, N_s, H_i, H_s, T_frz, Δt, u_star, T_w, T_0, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, true)
=#