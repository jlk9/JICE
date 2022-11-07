# Written by Joseph Kump (josek97@utexas.edu)
# Basic code that is just ice temperature in column, with no formulation or melting

include("../src/column_physics/run_ice_adjoint.jl")

using Enzyme, Test, Printf, LinearAlgebra

# Recall the run_ice_column function:
# Inputs:
#   H      (m)            total ice thickness, float
#   N_i      (dim'less)   number of layers, int (we consider the surface "skin layer")
#   S      (partsperthou) salinity, assumed to be constant throughout ice here, float
#   C_0    (J m^-3 k^-1)  initial volumetric heat capacity of ice (length K+1, including skin layer)
#   L      (J kg^-1)      latent heat of sublimation, float
#   T_frz  (K)            freezing point of water at salinity S, float
#   G_snow (W m^-2)       the heat flux from snow to ice, float
#   T_0    (K)            initial temps, array of floats (length K+1, including skin layer)
#   N_t    (dim'less)     number of time steps, int
#   Δt     (s)            size of each time step, float
#   u_star (m/s)          friction velocity
#   T_w    (K)            temperature of sea surface
#
# Output:
#   T_n    (K)            the sea ice layer temperatures at initial time + N_t*dt, array of K+1 floats

# Basic test with some common values for sea ice:
function test_temp_thickness(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w, T_0,
                             F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, U_a)

    atmodel = initialize_ATModel(N_t, F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, U_a)
    jcmodel = initialize_JICEColumn(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w, T_0)
    run_ice_column(jcmodel, atmodel)

    println("1. Basic temp test. Initial temps are:")
    println(T_0)
    println("Temps over time are: (should be close to initial)")
    println(jcmodel.T_array[:,N_t+1])

    println("Initial thicknesses are:")
    println(jcmodel.Δh_array[:,1])
    println("Thicknesses over time are: (should be close to initial in bottom)")
    println(jcmodel.Δh_array[:,N_t+1])

    atmodel = initialize_ATModel(N_t, F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, U_a)
    jcmodel = initialize_JICEColumn(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w, T_0)
    @time run_ice_column(jcmodel, atmodel)
    #@time compute_surface_flux(jcmodel, atmodel)

end

function test_tridiagonal_solve(T_0, N_i)
    println("4. Let's try autodiff on the tridiagonal solve, which is a little more complicated. First assume an identity matrix. Original ∂z_∂x is")
    maindiag = ones(Float64, N_i+1)
    subdiag  = zeros(Float64, N_i)
    supdiag  = zeros(Float64, N_i)

    x      = deepcopy(T_0)
    ∂z_∂x  = rand(N_i+1)
    ∂z_∂ma = zeros(N_i+1)
    ∂z_∂sb = zeros(N_i)
    ∂z_∂sp = zeros(N_i)

    println(∂z_∂x)
    autodiff(tridiagonal_solve, Const, Duplicated(x, ∂z_∂x), Const(N_i), Duplicated(maindiag, ∂z_∂ma), Duplicated(subdiag, ∂z_∂sb), Duplicated(supdiag, ∂z_∂sp))

    println("The new partial should be the same as the old one:")
    println(∂z_∂x)
    #println(∂z_∂ma)
    #println(∂z_∂sb)
    #println(∂z_∂sp)

    println("Now we'll assume -0.25 on sub and superdiagonals. Original ∂z_∂x is")
    maindiag = ones(Float64, N_i+1)
    subdiag  = -0.25 .* ones(Float64, N_i)
    supdiag  = -0.25 .* ones(Float64, N_i)

    x      = deepcopy(T_0)
    ∂z_∂x  = rand(N_i+1)
    ∂z_∂ma = zeros(N_i+1)
    ∂z_∂sb = zeros(N_i)
    ∂z_∂sp = zeros(N_i)

    mat            = Tridiagonal(subdiag, maindiag, supdiag)
    expected_∂z_∂x = mat\∂z_∂x

    println("Given ∂z_∂x is")
    println(∂z_∂x)
    println("What we expect from reverse-mode autodiff is the matrix inverse * original ∂z_∂x:")
    println(expected_∂z_∂x)

    autodiff(tridiagonal_solve, Const, Duplicated(x, ∂z_∂x), Const(N_i), Duplicated(maindiag, ∂z_∂ma), Duplicated(subdiag, ∂z_∂sb), Duplicated(supdiag, ∂z_∂sp))

    println("The actual partial ∂z_∂x after autodiff is:")
    println(∂z_∂x)
    #println(∂z_∂ma)
    #println(∂z_∂sb)
    #println(∂z_∂sp)
end

function test_run_ice_column_one_step(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w, T_0, F_0, dF_0)
    println("5. Now let's try sensitivity analysis on step_temp_change itself. In particular, the ice thickness to the water temperature:")

    # First we need to get our intermediate variables:
    Δh, Δh̄, S, c_i, K, K̄, I_pen, q_i, q_inew, z_old, z_new, maindiag, subdiag, supdiag, T_array, Δh_array = allocate_memory(N_i, N_t)
    for k in 1:N_i
        Δh[k+1] = H / N_i
    end
    T_n     = T_0
    T_nplus = zeros(Float64, N_i+1)
    generate_S(S, N_i)
    generate_I_pen(I_pen, i_0, κ_i, H, N_i)

    println("Let's begin by seeing what one step of the ice model yields. New thicknesses and temps:")
    step_temp_change(N_i, S, L, T_frz, κ_i, Δh, Δh̄, T_n, T_nplus, c_i, K, K̄, I_pen, q_i, q_inew,
                 z_old, z_new, F_0, dF_0, maindiag, subdiag, supdiag, Δt, u_star, T_w)

    println(Δh)
    println(T_nplus)

    # Now we'll reset for autodiff test:
    Δh, Δh̄, S, c_i, K, K̄, I_pen, q_i, q_inew, z_old, z_new, maindiag, subdiag, supdiag, T_array, Δh_array = allocate_memory(N_i, N_t)
    for k in 1:N_i
        Δh[k+1] = H / N_i
    end
    T_n    = T_0
    T_nplus = zeros(Float64, N_i+1)
    generate_S(S, N_i)
    generate_I_pen(I_pen, i_0, κ_i, H, N_i)

    ∂f_∂S        = ones(Float64, N_i+1)
    ∂f_∂h        = ones(Float64, N_i+1)
    ∂f_∂h̄        = ones(Float64, N_i)
    ∂f_∂T_n      = ones(Float64, N_i+1)
    ∂f_∂T_nplus  = ones(Float64, N_i+1)
    ∂f_∂c_i      = ones(Float64, N_i+1)
    ∂f_∂K        = ones(Float64, N_i+1)
    ∂f_∂K̄        = ones(Float64, N_i)
    ∂f_∂I_pen    = ones(Float64, N_i+1)
    ∂f_∂q_i      = ones(Float64, N_i+1)
    ∂f_∂q_inew   = ones(Float64, N_i+1)
    ∂f_∂z_old    = ones(Float64, N_i+1)
    ∂f_∂z_new    = ones(Float64, N_i+1)
    ∂f_∂maindiag = ones(Float64, N_i+1)
    ∂f_∂subdiag  = ones(Float64, N_i)
    ∂f_∂supdiag  = ones(Float64, N_i)

    # Running autodiff here
    ∂f_∂T_w = autodiff(step_temp_change, Const, Const(N_i), Duplicated(S, ∂f_∂S), Const(L), Const(T_frz),
    Const(κ_i), Duplicated(Δh, ∂f_∂h), Duplicated(Δh̄, ∂f_∂h̄), Duplicated(T_n, ∂f_∂T_n),
    Duplicated(T_nplus, ∂f_∂T_nplus), Duplicated(c_i, ∂f_∂c_i), Duplicated(K, ∂f_∂K), Duplicated(K̄, ∂f_∂K̄),
    Duplicated(I_pen, ∂f_∂I_pen), Duplicated(q_i, ∂f_∂q_i), Duplicated(q_inew, ∂f_∂q_inew), Duplicated(z_old, ∂f_∂z_old),
    Duplicated(z_new, ∂f_∂z_new), Const(F_0[1]), Const(dF_0[1]), Duplicated(maindiag, ∂f_∂maindiag),
    Duplicated(subdiag, ∂f_∂subdiag), Duplicated(supdiag, ∂f_∂supdiag), Const(Δt), Const(u_star), Active(T_w))

    println("Now the autodiff test. Thicknesses and temps should be the same:")
    println(Δh)
    println(T_nplus)

    println("∂f_∂T_w is")
    println(∂f_∂T_w)
end

function test_adjoint_temp(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w, T_0,
                            F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, U_a)

    println("6. A test of the adjoint method implementation as shown in tutorial by Sarah Williamson. First we'll run a forward loop:")

    atmodel = initialize_ATModel(N_t, F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, U_a)
    jcmodel = initialize_JICEColumn(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w, T_0)
    run_ice_column(jcmodel, atmodel)

    println("Now we'll go backwards.")
    # run one step of the adjoint function, Enzyme preferred that I did this one and then all of the others, I'm not sure why.......
    ad_h = [0.0;0.0;0.0;0.0;0.0;0.0]
    ad_T = [0.0;0.0;0.0;0.0;0.0;1.0]

    ∂h, ∂T_old = run_ice_column_adjoint_hT(jcmodel, atmodel, ad_h, ad_T)

    println("∂h:")
    println(∂h)
    println("∂T_old:")
    println(∂T_old)
    println("In particular, the diff of the second-from-bottom temp based on bottom temp is")
    init_layer = N_i
    println(∂T_old[init_layer])
    
    println("To test this, we'll see the result we get from finite differences. The ϵ values are:")

    step_sizes  = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13]
    diffs       = []
    rel_errors  = []

    for ϵ in step_sizes

        # Needs to be the same as our original T_0
        T_ϵp              = deepcopy(T_0)
        T_ϵp[init_layer] += ϵ
        jcmodelp          = initialize_JICEColumn(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w, T_ϵp)
        run_ice_column(jcmodelp, atmodel)
        
        T_ϵp_value  = jcmodelp.T_array[N_i+1, N_t+1]
        Δz_ϵp_value = jcmodelp.Δh_array[N_i+1, N_t+1]

        T_ϵn              = deepcopy(T_0)
        T_ϵn[init_layer] -= ϵ
        jcmodeln          = initialize_JICEColumn(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w, T_ϵn)
        run_ice_column(jcmodeln, atmodel)

        T_ϵn_value  = jcmodeln.T_array[N_i+1, N_t+1]
        Δz_ϵn_value = jcmodeln.Δh_array[N_i+1, N_t+1]

        diff = (T_ϵp_value - T_ϵn_value) / (2*ϵ)

        error = abs(∂T_old[init_layer] - diff) / abs(∂T_old[init_layer])

        push!(diffs, diff)
        push!(rel_errors, error)

    end

    println(step_sizes)
    println("The finite diffs are:")
    println(diffs)
    println("and the errors are:")
    println(rel_errors)
end

function test_adjoint_T_w(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w, T_0,
                            F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, U_a)

    println("7. A test of the adjoint method on the adjacent water temperature and its effect on thickness:")

    atmodel = initialize_ATModel(N_t, F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, U_a)
    jcmodel = initialize_JICEColumn(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w, T_0)
    run_ice_column(jcmodel, atmodel)
    
    println("Now we'll go backwards.")
    # run one step of the adjoint function, Enzyme preferred that I did this one and then all of the others, I'm not sure why.......
    ad_h = [0.0;0.0;0.0;0.0;0.0;1.0]
    ad_T = [0.0;0.0;0.0;0.0;0.0;0.0]

    ∂T_w = run_ice_column_adjoint_Tw(jcmodel, atmodel, ad_h, ad_T)

    println("∂T_w:")
    println(∂T_w)
    
    println("To test this, we'll see the result we get from finite differences. The ϵ values are:")

    step_sizes  = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13]
    diffs       = []
    rel_errors  = []

    for ϵ in step_sizes

        # Needs to be the same as our original T_0
        jcmodelp = initialize_JICEColumn(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w + ϵ, T_0)
        run_ice_column(jcmodelp, atmodel)
        
        T_ϵp_value  = jcmodelp.T_array[N_i+1, N_t+1]
        Δh_ϵp_value = jcmodelp.Δh_array[N_i+1, N_t+1]

        jcmodeln = initialize_JICEColumn(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w - ϵ, T_0)
        run_ice_column(jcmodeln, atmodel)

        T_ϵn_value  = jcmodeln.T_array[N_i+1, N_t+1]
        Δh_ϵn_value = jcmodeln.Δh_array[N_i+1, N_t+1]

        diff = (Δh_ϵp_value - Δh_ϵn_value) / (2*ϵ) #(T_ϵp_value - T_ϵn_value) / (2*ϵ)

        error = abs(∂T_w - diff) / abs(∂T_w)

        push!(diffs, diff)
        push!(rel_errors, error)

    end

    println(step_sizes)
    println("The finite diffs are:")
    println(diffs)
    println("and the errors are:")
    println(rel_errors)
end

N_i    = 5
N_t    = 4
H      = 2.0
#L      = 2260000.0
T_frz  = 271.35 - 273.15
i_0    = 0.7
κ_i    = 1.4
T_0    = 0 .- [20.0, 16.5, 13.0, 9.5, 6.0, 2.5]
Δt     = 1.0
u_star = 0.0005 # recommended minimum value of u_star in CICE
T_w    = 274.47 - 273.15 # typical temp in C for sea surface in arctic

# "Educated" guess for some normal atmospheric values
F_Ld  = 0.0
F_sw  = 120.0 
T_a   = -34.0
Θ_a   = T_a*(1000.0/1045.6)^0.286
ρ_a   = 1.4224
Q_a   = 0.005 #?
c_p   = 0.7171
U_a   = zeros(Float64, 3)

test_temp_thickness(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w, T_0,
                    F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, U_a)
println("")
#test_tridiagonal_solve(T_0, N_i)
println("")
#test_run_ice_column_one_step(N_i, N_t, H, L, T_frz, i_0, κ_i, Δt, u_star, T_w, T_0, 0.0, 0.0)
println("")

T_0  = 0 .- [20.0, 16.5, 13.0, 9.5, 6.0, 2.5]
N_t  = 1 # other variables are same as before

test_adjoint_temp(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w, T_0,
                    F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, U_a)
println("")
test_adjoint_T_w(N_i, N_t, H, T_frz, i_0, κ_i, Δt, u_star, T_w, T_0,
                    F_Ld, F_sw, T_a, Θ_a, ρ_a, Q_a, c_p, U_a)
