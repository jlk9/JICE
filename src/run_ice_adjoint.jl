# Written by Joseph Kump (josek97@utexas.edu)
# Adjoint method code for run_ice_inline.jl


include("run_ice_inline.jl")

using Enzyme

# Finds the adjoint of a completed JCModel (i.e. run_ice has already been called)
function run_ice_adjoint_hT(jcmodel, ad_h, ad_T)

    # run one step of the adjoint function:
    ∂h, ∂T_old = run_ice_adjoint_hT_step(jcmodel, jcmodel.N_t, ad_h, ad_T)

    # then the rest back to the original time step:
    for j in (jcmodel.N_t-1):-1:1

        ad_h[:] = ∂h
        ad_T[:] = ∂T_old

        ∂h, ∂T_old = run_ice_adjoint_hT_step(jcmodel, j, ad_h, ad_T)

    end

    return ∂h, ∂T_old

end

# Finds the adjoint of a completed JCModel (i.e. run_ice has already been called)
# This one gets the adjoint for T_w, which is tricky because T_w is applied
# separately at each step
function run_ice_adjoint_Tw(jcmodel, ad_h, ad_T)

    # run one step of the adjoint function, for both T_w and h, T_old:
    ∂h, ∂T_old = run_ice_adjoint_hT_step(jcmodel, jcmodel.N_t, ad_h, ad_T)
    ∂T_w = run_ice_adjoint_Tw_step(jcmodel, jcmodel.N_t, ad_h, ad_T)

    # then the rest back to the original time step:
    for j in (jcmodel.N_t-1):-1:1

        ad_h[:] = ∂h
        ad_T[:] = ∂T_old

        ∂h, ∂T_old = run_ice_adjoint_hT_step(jcmodel, j, ad_h, ad_T)

        ∂T_w += run_ice_adjoint_Tw_step(jcmodel, j, ad_h, ad_T)

    end

    return ∂T_w

end

# Runs the AD steps for computing how initial Δh and T_old affect changes in Δh and T_new
# TODO: keep track of time step for flux
function run_ice_adjoint_hT_step(jcmodel, step, ∂f_∂h, ∂f_∂T_new)

    # Producing derivatives for inputs / intermdiate variables:
    ∂f_∂S        = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂h̄        = zeros(Float64, jcmodel.N_i)
    ∂f_∂T_old    = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂c_i      = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂K        = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂K̄        = zeros(Float64, jcmodel.N_i)
    ∂f_∂I_pen    = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂q_i      = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂maindiag = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂subdiag  = zeros(Float64, jcmodel.N_i)
    ∂f_∂supdiag  = zeros(Float64, jcmodel.N_i)


    autodiff(run_ice_step, Const, Const(jcmodel.N_i), Duplicated(jcmodel.S, ∂f_∂S), Const(jcmodel.L), Const(jcmodel.T_frz),
    Const(jcmodel.κ_i), Duplicated(jcmodel.Δh_array[:,step], ∂f_∂h), Duplicated(jcmodel.Δh̄, ∂f_∂h̄), Duplicated(jcmodel.T_array[:,step], ∂f_∂T_old),
    Duplicated(jcmodel.T_array[:,step+1], ∂f_∂T_new), Duplicated(jcmodel.c_i, ∂f_∂c_i), Duplicated(jcmodel.K, ∂f_∂K), Duplicated(jcmodel.K̄, ∂f_∂K̄),
    Duplicated(jcmodel.I_pen, ∂f_∂I_pen), Duplicated(jcmodel.q_i, ∂f_∂q_i), Const(jcmodel.F_0[step]), Const(jcmodel.dF_0[step]), Duplicated(jcmodel.maindiag, ∂f_∂maindiag),
    Duplicated(jcmodel.subdiag, ∂f_∂subdiag), Duplicated(jcmodel.supdiag, ∂f_∂supdiag), Const(jcmodel.Δt), Const(jcmodel.u_star),
    Const(jcmodel.T_w))
    

    return ∂f_∂h, ∂f_∂T_old

end

# Runs the AD steps for computing how initial T_w affects changes in Δh and T_new
# TODO: keep track of time step for flux
function run_ice_adjoint_Tw_step(jcmodel, step, ∂f_∂h, ∂f_∂T_new)

    # Producing derivatives for inputs / intermdiate variables:
    ∂f_∂S        = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂h̄        = zeros(Float64, jcmodel.N_i)
    ∂f_∂T_old    = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂c_i      = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂K        = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂K̄        = zeros(Float64, jcmodel.N_i)
    ∂f_∂I_pen    = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂q_i      = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂maindiag = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂subdiag  = zeros(Float64, jcmodel.N_i)
    ∂f_∂supdiag  = zeros(Float64, jcmodel.N_i)


    ∂f_∂T_w = autodiff(run_ice_step, Const, Const(jcmodel.N_i), Duplicated(jcmodel.S, ∂f_∂S), Const(jcmodel.L), Const(jcmodel.T_frz),
    Const(jcmodel.κ_i), Duplicated(jcmodel.Δh_array[:,step], ∂f_∂h), Duplicated(jcmodel.Δh̄, ∂f_∂h̄), Duplicated(jcmodel.T_array[:,step], ∂f_∂T_old),
    Duplicated(jcmodel.T_array[:,step+1], ∂f_∂T_new), Duplicated(jcmodel.c_i, ∂f_∂c_i), Duplicated(jcmodel.K, ∂f_∂K), Duplicated(jcmodel.K̄, ∂f_∂K̄),
    Duplicated(jcmodel.I_pen, ∂f_∂I_pen), Duplicated(jcmodel.q_i, ∂f_∂q_i), Const(jcmodel.F_0[step]), Const(jcmodel.dF_0[step]), Duplicated(jcmodel.maindiag, ∂f_∂maindiag),
    Duplicated(jcmodel.subdiag, ∂f_∂subdiag), Duplicated(jcmodel.supdiag, ∂f_∂supdiag), Const(jcmodel.Δt), Const(jcmodel.u_star),
    Active(jcmodel.T_w))
    

    return ∂f_∂T_w[1]

end