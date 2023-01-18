# Written by Joseph Kump (josek97@utexas.edu)
# Adjoint method code for run_ice_inline.jl


include("run_ice_column.jl")

using Enzyme

#= Runs the single column sea ice model and computes the derivatives for all of its fields. The various
    fields of d_jcolumn are the partial derivatives of those fields given the gradients for thickness and
    current temperature (ad_h and ad_T).

NOTE:   In cases where there is no snow and empty (thickness 0) snow layers, the partial derivative in terms of
        thickness can be NaN.
=#
function run_ice_column_autodiff(jcolumn, atmodel, ad_h, ad_T)

    # first we need to allocate d_jcolumn and d_atmodel:
    d_jcolumn = initialize_JICEColumn(jcolumn.N_t, jcolumn.N_i, jcolumn.N_s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, zeros(Float64, jcolumn.N_i+jcolumn.N_s+1))
    d_atmodel = initialize_ATModel(jcolumn.N_t, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, zeros(Float64, 3))

    # Want to reset some derivatives to 0:
    d_jcolumn.S[:] = zeros(Float64, jcolumn.N_i)
    d_jcolumn.N_t  = 0
    d_jcolumn.N_i  = 0
    d_jcolumn.N_s  = 0
    
    d_jcolumn.T_nplus = ad_T
    d_jcolumn.Δh      = ad_h

    autodiff(run_ice_column, Const, Duplicated(jcolumn, d_jcolumn), Duplicated(atmodel, d_atmodel))


    return d_jcolumn #d_jcolumn.Δh, d_jcolumn.T_n

end

#=
# Finds the adjoint of a completed JCModel (i.e. run_ice has already been called)
function run_ice_column_adjoint_hT(jcolumn, atmodel, ad_h, ad_T)

    # run one step of the adjoint function:
    d_jcol = run_ice_adjoint_step(jcolumn, atmodel, jcolumn.N_t, ad_h, ad_T)

    # then the rest back to the original time step:
    for j in (jcolumn.N_t-1):-1:1

        ad_h[:] = d_jcol.Δh
        ad_T[:] = d_jcol.T_n

        d_jcol = run_ice_adjoint_step(jcolumn, atmodel, j, ad_h, ad_T)

    end

    return d_jcol.Δh, d_jcol.T_n
end

# Finds the adjoint of a completed JCModel (i.e. run_ice has already been called)
# This one gets the adjoint for T_w, which is tricky because T_w is applied
# separately at each step
function run_ice_column_adjoint_Tw(jcolumn, atmodel, ad_h, ad_T)

    # run one step of the adjoint function, for both T_w and h, T_old:
    d_jcol = run_ice_adjoint_step(jcolumn, atmodel, jcolumn.N_t, ad_h, ad_T)
    ∂T_w   = d_jcol.T_w

    # then the rest back to the original time step:
    for j in (jcolumn.N_t-1):-1:1

        ad_h[:] = d_jcol.Δh
        ad_T[:] = d_jcol.T_n

        d_jcol     = run_ice_adjoint_step(jcolumn, atmodel, j, ad_h, ad_T)
        ∂T_w      += d_jcol.T_w

    end

    return ∂T_w

end

# Runs the AD steps for computing how initial variables affect changes in the outputs of our model, Δh and T_new.
# It applies enzyme directly to the struct, which requires us to create a second (empty)
# struct as the derivative
function run_ice_adjoint_step(jcolumn, atmodel, step, ∂f_∂h, ∂f_∂T_new)

    # TODO: autodiff includes forward loop and reverse mode AD pullback

    # first we need to allocate d_jcolumn and d_atmodel:
    d_jcolumn = initialize_JICEColumn(jcolumn.N_t, jcolumn.N_i, jcolumn.N_s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, zeros(Float64, jcolumn.N_i+jcolumn.N_s+1))
    d_atmodel = initialize_ATModel(jcolumn.N_t, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, zeros(Float64, 3))
    
    # Want to reset some derivatives to 0:
    d_jcolumn.S[:] = zeros(Float64, jcolumn.N_i)
    d_jcolumn.N_t  = 0
    d_jcolumn.N_i  = 0
    d_jcolumn.N_s  = 0
    
    d_jcolumn.T_nplus = ∂f_∂T_new
    d_jcolumn.Δh      = ∂f_∂h

    autodiff(run_column_step, Const, Duplicated(jcolumn, d_jcolumn), Duplicated(atmodel, d_atmodel), Const(step))

    return d_jcolumn
end


# Runs the AD steps for computing how initial Δh and T_old affect changes in Δh and T_new
# OLD VERSION that applies enzyme on all the struct's separate fields
function run_ice_adjoint_hT_step_old(jcolumn, atmodel, step, ∂f_∂h, ∂f_∂T_new)

    # Producing derivatives for inputs / intermdiate variables:
    ∂f_∂F_0      = zeros(Float64, jcolumn.N_t)
    ∂f_∂dF_0     = zeros(Float64, jcolumn.N_t)
    ∂f_∂h̄        = zeros(Float64, jcolumn.N_i+jcolumn.N_s)
    ∂f_∂S        = zeros(Float64, jcolumn.N_i)
    ∂f_∂T_old    = zeros(Float64, jcolumn.N_i+jcolumn.N_s+1)
    ∂f_∂c_i      = zeros(Float64, jcolumn.N_i+jcolumn.N_s+1)
    ∂f_∂K        = zeros(Float64, jcolumn.N_i+jcolumn.N_s+1)
    ∂f_∂K̄        = zeros(Float64, jcolumn.N_i+jcolumn.N_s)
    ∂f_∂I_pen    = zeros(Float64, jcolumn.N_i)
    ∂f_∂q_i      = zeros(Float64, jcolumn.N_i)
    ∂f_∂q_inew   = zeros(Float64, jcolumn.N_i)
    ∂f_∂z_old    = zeros(Float64, jcolumn.N_i+1)
    ∂f_∂z_new    = zeros(Float64, jcolumn.N_i+1)
    ∂f_∂maindiag = zeros(Float64, jcolumn.N_i+jcolumn.N_s+1)
    ∂f_∂subdiag  = zeros(Float64, jcolumn.N_i+jcolumn.N_s)
    ∂f_∂supdiag  = zeros(Float64, jcolumn.N_i+jcolumn.N_s)
    ∂f_∂F_Lu     = zeros(Float64, jcolumn.N_t)
    ∂f_∂F_s      = zeros(Float64, jcolumn.N_t)
    ∂f_∂F_l      = zeros(Float64, jcolumn.N_t)
    ∂f_∂dF_Lu    = zeros(Float64, jcolumn.N_t)
    ∂f_∂dF_s     = zeros(Float64, jcolumn.N_t)
    ∂f_∂dF_l     = zeros(Float64, jcolumn.N_t)

    ∂f_∂α_vdr = zeros(Float64, 2)
    ∂f_∂α_idr = zeros(Float64, 2)
    ∂f_∂α_vdf = zeros(Float64, 2)
    ∂f_∂α_idf = zeros(Float64, 2)

    ∂f_∂c_u        = zeros(Float64, 1)
    ∂f_∂c_Θ        = zeros(Float64, 1)
    ∂f_∂c_q        = zeros(Float64, 1)
    ∂f_∂atm_u_star = zeros(Float64, 1)
    ∂f_∂U_a        = zeros(Float64, 3)


    autodiff(run_column_step, Const, Const(jcolumn.N_i), Const(jcolumn.N_s), Const(jcolumn.N_t),
            Const(sum(jcolumn.Δh_array[:,step])), Const(jcolumn.H_s), Const(jcolumn.T_frz),
            Const(jcolumn.Δt), Const(jcolumn.u_star), Const(jcolumn.T_w),
            Duplicated(jcolumn.α_vdr, ∂f_∂α_vdr), Duplicated(jcolumn.α_idr, ∂f_∂α_idr),
            Duplicated(jcolumn.α_vdf, ∂f_∂α_vdf), Duplicated(jcolumn.α_idf, ∂f_∂α_idf),
            Duplicated(jcolumn.F_0, ∂f_∂F_0), Duplicated(jcolumn.dF_0, ∂f_∂dF_0),
            Duplicated(jcolumn.Δh_array[:,step], ∂f_∂h), Duplicated(jcolumn.Δh̄, ∂f_∂h̄),
            Duplicated(jcolumn.S, ∂f_∂S), Duplicated(jcolumn.c_i, ∂f_∂c_i),
            Duplicated(jcolumn.K, ∂f_∂K), Duplicated(jcolumn.K̄, ∂f_∂K̄), Duplicated(jcolumn.I_pen, ∂f_∂I_pen),
            Duplicated(jcolumn.q_i, ∂f_∂q_i), Duplicated(jcolumn.q_inew, ∂f_∂q_inew),
            Duplicated(jcolumn.z_old, ∂f_∂z_old), Duplicated(jcolumn.z_new, ∂f_∂z_new),
            Duplicated(jcolumn.maindiag, ∂f_∂maindiag), Duplicated(jcolumn.subdiag, ∂f_∂subdiag), Duplicated(jcolumn.supdiag, ∂f_∂supdiag),
            Duplicated(jcolumn.F_Lu, ∂f_∂F_Lu), Duplicated(jcolumn.F_s, ∂f_∂F_s), Duplicated(jcolumn.F_l, ∂f_∂F_l),
            Duplicated(jcolumn.dF_Lu, ∂f_∂dF_Lu), Duplicated(jcolumn.dF_s, ∂f_∂dF_s), Duplicated(jcolumn.dF_l, ∂f_∂dF_l),
            Const(atmodel.F_SWvdr), Const(atmodel.F_SWidr), Const(atmodel.F_SWvdf), Const(atmodel.F_SWidf),
            Const(atmodel.F_Ld), Const(atmodel.T_a),
            Const(atmodel.Θ_a), Const(atmodel.ρ_a), Const(atmodel.Q_a), Const(atmodel.c_p),
            Duplicated(atmodel.c_u, ∂f_∂c_u), Duplicated(atmodel.c_Θ, ∂f_∂c_Θ), Duplicated(atmodel.c_q, ∂f_∂c_q),
            Duplicated(atmodel.atm_u_star, ∂f_∂atm_u_star), Duplicated(atmodel.U_a, ∂f_∂U_a),
            Duplicated(jcolumn.T_array[:,step], ∂f_∂T_old), Duplicated(jcolumn.T_array[:,step+1], ∂f_∂T_new),
            Const(step))

    return ∂f_∂h, ∂f_∂T_old
end

# Runs the AD steps for computing how initial T_w affects changes in Δh and T_new
# OLD VERSION that applies enzyme on all the struct's separate fields
function run_ice_adjoint_Tw_step_old(jcolumn, atmodel, step, ∂f_∂h, ∂f_∂T_new)

    # Producing derivatives for inputs / intermdiate variables:
    ∂f_∂F_0      = zeros(Float64, jcolumn.N_t)
    ∂f_∂dF_0     = zeros(Float64, jcolumn.N_t)
    ∂f_∂h̄        = zeros(Float64, jcolumn.N_i+jcolumn.N_s)
    ∂f_∂S        = zeros(Float64, jcolumn.N_i)
    ∂f_∂T_old    = zeros(Float64, jcolumn.N_i+jcolumn.N_s+1)
    ∂f_∂c_i      = zeros(Float64, jcolumn.N_i+jcolumn.N_s+1)
    ∂f_∂K        = zeros(Float64, jcolumn.N_i+jcolumn.N_s+1)
    ∂f_∂K̄        = zeros(Float64, jcolumn.N_i+jcolumn.N_s)
    ∂f_∂I_pen    = zeros(Float64, jcolumn.N_i)
    ∂f_∂q_i      = zeros(Float64, jcolumn.N_i)
    ∂f_∂q_inew   = zeros(Float64, jcolumn.N_i)
    ∂f_∂z_old    = zeros(Float64, jcolumn.N_i+1)
    ∂f_∂z_new    = zeros(Float64, jcolumn.N_i+1)
    ∂f_∂maindiag = zeros(Float64, jcolumn.N_i+jcolumn.N_s+1)
    ∂f_∂subdiag  = zeros(Float64, jcolumn.N_i+jcolumn.N_s)
    ∂f_∂supdiag  = zeros(Float64, jcolumn.N_i+jcolumn.N_s)
    ∂f_∂F_Lu     = zeros(Float64, jcolumn.N_t)
    ∂f_∂F_s      = zeros(Float64, jcolumn.N_t)
    ∂f_∂F_l      = zeros(Float64, jcolumn.N_t)
    ∂f_∂dF_Lu    = zeros(Float64, jcolumn.N_t)
    ∂f_∂dF_s     = zeros(Float64, jcolumn.N_t)
    ∂f_∂dF_l     = zeros(Float64, jcolumn.N_t)

    ∂f_∂α_vdr = zeros(Float64, 2)
    ∂f_∂α_idr = zeros(Float64, 2)
    ∂f_∂α_vdf = zeros(Float64, 2)
    ∂f_∂α_idf = zeros(Float64, 2)

    ∂f_∂c_u        = zeros(Float64, 1)
    ∂f_∂c_Θ        = zeros(Float64, 1)
    ∂f_∂c_q        = zeros(Float64, 1)
    ∂f_∂atm_u_star = zeros(Float64, 1)
    ∂f_∂U_a        = zeros(Float64, 3)


    ∂f_∂T_w = autodiff(run_column_step, Const, Const(jcolumn.N_i), Const(jcolumn.N_s), Const(jcolumn.N_t),
                Const(sum(jcolumn.Δh_array[:,step])), Const(jcolumn.H_s), Const(jcolumn.T_frz),
                Const(jcolumn.Δt), Const(jcolumn.u_star), Active(jcolumn.T_w),
                Duplicated(jcolumn.α_vdr, ∂f_∂α_vdr), Duplicated(jcolumn.α_idr, ∂f_∂α_idr),
                Duplicated(jcolumn.α_vdf, ∂f_∂α_vdf), Duplicated(jcolumn.α_idf, ∂f_∂α_idf),
                Duplicated(jcolumn.F_0, ∂f_∂F_0), Duplicated(jcolumn.dF_0, ∂f_∂dF_0),
                Duplicated(jcolumn.Δh_array[:,step], ∂f_∂h), Duplicated(jcolumn.Δh̄, ∂f_∂h̄),
                Duplicated(jcolumn.S, ∂f_∂S), Duplicated(jcolumn.c_i, ∂f_∂c_i),
                Duplicated(jcolumn.K, ∂f_∂K), Duplicated(jcolumn.K̄, ∂f_∂K̄), Duplicated(jcolumn.I_pen, ∂f_∂I_pen),
                Duplicated(jcolumn.q_i, ∂f_∂q_i), Duplicated(jcolumn.q_inew, ∂f_∂q_inew),
                Duplicated(jcolumn.z_old, ∂f_∂z_old), Duplicated(jcolumn.z_new, ∂f_∂z_new),
                Duplicated(jcolumn.maindiag, ∂f_∂maindiag), Duplicated(jcolumn.subdiag, ∂f_∂subdiag), Duplicated(jcolumn.supdiag, ∂f_∂supdiag),
                Duplicated(jcolumn.F_Lu, ∂f_∂F_Lu), Duplicated(jcolumn.F_s, ∂f_∂F_s), Duplicated(jcolumn.F_l, ∂f_∂F_l),
                Duplicated(jcolumn.dF_Lu, ∂f_∂dF_Lu), Duplicated(jcolumn.dF_s, ∂f_∂dF_s), Duplicated(jcolumn.dF_l, ∂f_∂dF_l),
                Const(atmodel.F_SWvdr), Const(atmodel.F_SWidr), Const(atmodel.F_SWvdf), Const(atmodel.F_SWidf),
                Const(atmodel.F_Ld), Const(atmodel.T_a),
                Const(atmodel.Θ_a), Const(atmodel.ρ_a), Const(atmodel.Q_a), Const(atmodel.c_p),
                Duplicated(atmodel.c_u, ∂f_∂c_u), Duplicated(atmodel.c_Θ, ∂f_∂c_Θ), Duplicated(atmodel.c_q, ∂f_∂c_q),
                Duplicated(atmodel.atm_u_star, ∂f_∂atm_u_star), Duplicated(atmodel.U_a, ∂f_∂U_a),
                Duplicated(jcolumn.T_array[:,step], ∂f_∂T_old), Duplicated(jcolumn.T_array[:,step+1], ∂f_∂T_new),
                Const(step))

    return ∂f_∂T_w[1]
end
=#