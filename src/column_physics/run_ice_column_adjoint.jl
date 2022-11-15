# Written by Joseph Kump (josek97@utexas.edu)
# Adjoint method code for run_ice_inline.jl


include("run_ice_column.jl")

using Enzyme

# Finds the adjoint of a completed JCModel (i.e. run_ice has already been called)
function run_ice_column_adjoint_hT(jcmodel, atmodel, ad_h, ad_T)

    # run one step of the adjoint function:
    ∂h, ∂T_old = run_ice_adjoint_hT_step(jcmodel, atmodel, jcmodel.N_t, ad_h, ad_T)

    # then the rest back to the original time step:
    for j in (jcmodel.N_t-1):-1:1

        ad_h[:] = ∂h
        ad_T[:] = ∂T_old

        ∂h, ∂T_old = run_ice_adjoint_hT_step(jcmodel, atmodel, j, ad_h, ad_T)

    end

    return ∂h, ∂T_old

end

# Finds the adjoint of a completed JCModel (i.e. run_ice has already been called)
# This one gets the adjoint for T_w, which is tricky because T_w is applied
# separately at each step
function run_ice_column_adjoint_Tw(jcmodel, atmodel, ad_h, ad_T)

    # run one step of the adjoint function, for both T_w and h, T_old:
    ∂h, ∂T_old = run_ice_adjoint_hT_step(jcmodel, atmodel, jcmodel.N_t, ad_h, ad_T)
    ∂T_w = run_ice_adjoint_Tw_step(jcmodel, atmodel, jcmodel.N_t, ad_h, ad_T)

    # then the rest back to the original time step:
    for j in (jcmodel.N_t-1):-1:1

        ad_h[:] = ∂h
        ad_T[:] = ∂T_old

        ∂h, ∂T_old = run_ice_adjoint_hT_step(jcmodel, atmodel, j, ad_h, ad_T)

        ∂T_w += run_ice_adjoint_Tw_step(jcmodel, atmodel, j, ad_h, ad_T)

    end

    return ∂T_w

end

# Runs the AD steps for computing how initial Δh and T_old affect changes in Δh and T_new
# TODO: keep track of time step for flux
function run_ice_adjoint_hT_step(jcmodel, atmodel, step, ∂f_∂h, ∂f_∂T_new)

    # Producing derivatives for inputs / intermdiate variables:
    ∂f_∂F_0      = zeros(Float64, jcmodel.N_t)
    ∂f_∂dF_0     = zeros(Float64, jcmodel.N_t)
    ∂f_∂h̄        = zeros(Float64, jcmodel.N_i+jcmodel.N_s)
    ∂f_∂S        = zeros(Float64, jcmodel.N_i)
    ∂f_∂T_old    = zeros(Float64, jcmodel.N_i+jcmodel.N_s+1)
    ∂f_∂c_i      = zeros(Float64, jcmodel.N_i+jcmodel.N_s+1)
    ∂f_∂K        = zeros(Float64, jcmodel.N_i+jcmodel.N_s+1)
    ∂f_∂K̄        = zeros(Float64, jcmodel.N_i+jcmodel.N_s)
    ∂f_∂I_pen    = zeros(Float64, jcmodel.N_i)
    ∂f_∂q_i      = zeros(Float64, jcmodel.N_i)
    ∂f_∂q_inew   = zeros(Float64, jcmodel.N_i)
    ∂f_∂z_old    = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂z_new    = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂maindiag = zeros(Float64, jcmodel.N_i+jcmodel.N_s+1)
    ∂f_∂subdiag  = zeros(Float64, jcmodel.N_i+jcmodel.N_s)
    ∂f_∂supdiag  = zeros(Float64, jcmodel.N_i+jcmodel.N_s)
    ∂f_∂F_Lu     = zeros(Float64, jcmodel.N_t)
    ∂f_∂F_s      = zeros(Float64, jcmodel.N_t)
    ∂f_∂F_l      = zeros(Float64, jcmodel.N_t)
    ∂f_∂dF_Lu    = zeros(Float64, jcmodel.N_t)
    ∂f_∂dF_s     = zeros(Float64, jcmodel.N_t)
    ∂f_∂dF_l     = zeros(Float64, jcmodel.N_t)

    ∂f_∂α_vdr = zeros(Float64, 2)
    ∂f_∂α_idr = zeros(Float64, 2)
    ∂f_∂α_vdf = zeros(Float64, 2)
    ∂f_∂α_idf = zeros(Float64, 2)

    ∂f_∂c_u        = zeros(Float64, 1)
    ∂f_∂c_Θ        = zeros(Float64, 1)
    ∂f_∂c_q        = zeros(Float64, 1)
    ∂f_∂atm_u_star = zeros(Float64, 1)
    ∂f_∂U_a        = zeros(Float64, 3)


    autodiff(run_column_step, Const, Const(jcmodel.N_i), Const(jcmodel.N_s), Const(jcmodel.N_t),
            Const(sum(jcmodel.Δh_array[:,step])), Const(jcmodel.H_s), Const(jcmodel.T_frz),
            Const(jcmodel.Δt), Const(jcmodel.u_star), Const(jcmodel.T_w),
            Duplicated(jcmodel.α_vdr, ∂f_∂α_vdr), Duplicated(jcmodel.α_idr, ∂f_∂α_idr),
            Duplicated(jcmodel.α_vdf, ∂f_∂α_vdf), Duplicated(jcmodel.α_idf, ∂f_∂α_idf),
            Duplicated(jcmodel.F_0, ∂f_∂F_0), Duplicated(jcmodel.dF_0, ∂f_∂dF_0),
            Duplicated(jcmodel.Δh_array[:,step], ∂f_∂h), Duplicated(jcmodel.Δh̄, ∂f_∂h̄),
            Duplicated(jcmodel.S, ∂f_∂S), Duplicated(jcmodel.c_i, ∂f_∂c_i),
            Duplicated(jcmodel.K, ∂f_∂K), Duplicated(jcmodel.K̄, ∂f_∂K̄), Duplicated(jcmodel.I_pen, ∂f_∂I_pen),
            Duplicated(jcmodel.q_i, ∂f_∂q_i), Duplicated(jcmodel.q_inew, ∂f_∂q_inew),
            Duplicated(jcmodel.z_old, ∂f_∂z_old), Duplicated(jcmodel.z_new, ∂f_∂z_new),
            Duplicated(jcmodel.maindiag, ∂f_∂maindiag), Duplicated(jcmodel.subdiag, ∂f_∂subdiag), Duplicated(jcmodel.supdiag, ∂f_∂supdiag),
            Duplicated(jcmodel.F_Lu, ∂f_∂F_Lu), Duplicated(jcmodel.F_s, ∂f_∂F_s), Duplicated(jcmodel.F_l, ∂f_∂F_l),
            Duplicated(jcmodel.dF_Lu, ∂f_∂dF_Lu), Duplicated(jcmodel.dF_s, ∂f_∂dF_s), Duplicated(jcmodel.dF_l, ∂f_∂dF_l),
            Const(atmodel.F_SWvdr), Const(atmodel.F_SWidr), Const(atmodel.F_SWvdf), Const(atmodel.F_SWidf),
            Const(atmodel.F_Ld), Const(atmodel.T_a),
            Const(atmodel.Θ_a), Const(atmodel.ρ_a), Const(atmodel.Q_a), Const(atmodel.c_p),
            Duplicated(atmodel.c_u, ∂f_∂c_u), Duplicated(atmodel.c_Θ, ∂f_∂c_Θ), Duplicated(atmodel.c_q, ∂f_∂c_q),
            Duplicated(atmodel.atm_u_star, ∂f_∂atm_u_star), Duplicated(atmodel.U_a, ∂f_∂U_a),
            Duplicated(jcmodel.T_array[:,step], ∂f_∂T_old), Duplicated(jcmodel.T_array[:,step+1], ∂f_∂T_new),
            Const(step))

    return ∂f_∂h, ∂f_∂T_old
end

# Runs the AD steps for computing how initial T_w affects changes in Δh and T_new
# TODO: keep track of time step for flux
function run_ice_adjoint_Tw_step(jcmodel, atmodel, step, ∂f_∂h, ∂f_∂T_new)

    # Producing derivatives for inputs / intermdiate variables:
    ∂f_∂F_0      = zeros(Float64, jcmodel.N_t)
    ∂f_∂dF_0     = zeros(Float64, jcmodel.N_t)
    ∂f_∂h̄        = zeros(Float64, jcmodel.N_i+jcmodel.N_s)
    ∂f_∂S        = zeros(Float64, jcmodel.N_i)
    ∂f_∂T_old    = zeros(Float64, jcmodel.N_i+jcmodel.N_s+1)
    ∂f_∂c_i      = zeros(Float64, jcmodel.N_i+jcmodel.N_s+1)
    ∂f_∂K        = zeros(Float64, jcmodel.N_i+jcmodel.N_s+1)
    ∂f_∂K̄        = zeros(Float64, jcmodel.N_i+jcmodel.N_s)
    ∂f_∂I_pen    = zeros(Float64, jcmodel.N_i)
    ∂f_∂q_i      = zeros(Float64, jcmodel.N_i)
    ∂f_∂q_inew   = zeros(Float64, jcmodel.N_i)
    ∂f_∂z_old    = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂z_new    = zeros(Float64, jcmodel.N_i+1)
    ∂f_∂maindiag = zeros(Float64, jcmodel.N_i+jcmodel.N_s+1)
    ∂f_∂subdiag  = zeros(Float64, jcmodel.N_i+jcmodel.N_s)
    ∂f_∂supdiag  = zeros(Float64, jcmodel.N_i+jcmodel.N_s)
    ∂f_∂F_Lu     = zeros(Float64, jcmodel.N_t)
    ∂f_∂F_s      = zeros(Float64, jcmodel.N_t)
    ∂f_∂F_l      = zeros(Float64, jcmodel.N_t)
    ∂f_∂dF_Lu    = zeros(Float64, jcmodel.N_t)
    ∂f_∂dF_s     = zeros(Float64, jcmodel.N_t)
    ∂f_∂dF_l     = zeros(Float64, jcmodel.N_t)

    ∂f_∂α_vdr = zeros(Float64, 2)
    ∂f_∂α_idr = zeros(Float64, 2)
    ∂f_∂α_vdf = zeros(Float64, 2)
    ∂f_∂α_idf = zeros(Float64, 2)

    ∂f_∂c_u        = zeros(Float64, 1)
    ∂f_∂c_Θ        = zeros(Float64, 1)
    ∂f_∂c_q        = zeros(Float64, 1)
    ∂f_∂atm_u_star = zeros(Float64, 1)
    ∂f_∂U_a        = zeros(Float64, 3)


    ∂f_∂T_w = autodiff(run_column_step, Const, Const(jcmodel.N_i), Const(jcmodel.N_s), Const(jcmodel.N_t),
                Const(sum(jcmodel.Δh_array[:,step])), Const(jcmodel.H_s), Const(jcmodel.T_frz),
                Const(jcmodel.Δt), Const(jcmodel.u_star), Active(jcmodel.T_w),
                Duplicated(jcmodel.α_vdr, ∂f_∂α_vdr), Duplicated(jcmodel.α_idr, ∂f_∂α_idr),
                Duplicated(jcmodel.α_vdf, ∂f_∂α_vdf), Duplicated(jcmodel.α_idf, ∂f_∂α_idf),
                Duplicated(jcmodel.F_0, ∂f_∂F_0), Duplicated(jcmodel.dF_0, ∂f_∂dF_0),
                Duplicated(jcmodel.Δh_array[:,step], ∂f_∂h), Duplicated(jcmodel.Δh̄, ∂f_∂h̄),
                Duplicated(jcmodel.S, ∂f_∂S), Duplicated(jcmodel.c_i, ∂f_∂c_i),
                Duplicated(jcmodel.K, ∂f_∂K), Duplicated(jcmodel.K̄, ∂f_∂K̄), Duplicated(jcmodel.I_pen, ∂f_∂I_pen),
                Duplicated(jcmodel.q_i, ∂f_∂q_i), Duplicated(jcmodel.q_inew, ∂f_∂q_inew),
                Duplicated(jcmodel.z_old, ∂f_∂z_old), Duplicated(jcmodel.z_new, ∂f_∂z_new),
                Duplicated(jcmodel.maindiag, ∂f_∂maindiag), Duplicated(jcmodel.subdiag, ∂f_∂subdiag), Duplicated(jcmodel.supdiag, ∂f_∂supdiag),
                Duplicated(jcmodel.F_Lu, ∂f_∂F_Lu), Duplicated(jcmodel.F_s, ∂f_∂F_s), Duplicated(jcmodel.F_l, ∂f_∂F_l),
                Duplicated(jcmodel.dF_Lu, ∂f_∂dF_Lu), Duplicated(jcmodel.dF_s, ∂f_∂dF_s), Duplicated(jcmodel.dF_l, ∂f_∂dF_l),
                Const(atmodel.F_SWvdr), Const(atmodel.F_SWidr), Const(atmodel.F_SWvdf), Const(atmodel.F_SWidf),
                Const(atmodel.F_Ld), Const(atmodel.T_a),
                Const(atmodel.Θ_a), Const(atmodel.ρ_a), Const(atmodel.Q_a), Const(atmodel.c_p),
                Duplicated(atmodel.c_u, ∂f_∂c_u), Duplicated(atmodel.c_Θ, ∂f_∂c_Θ), Duplicated(atmodel.c_q, ∂f_∂c_q),
                Duplicated(atmodel.atm_u_star, ∂f_∂atm_u_star), Duplicated(atmodel.U_a, ∂f_∂U_a),
                Duplicated(jcmodel.T_array[:,step], ∂f_∂T_old), Duplicated(jcmodel.T_array[:,step+1], ∂f_∂T_new),
                Const(step))

    return ∂f_∂T_w[1]
end
