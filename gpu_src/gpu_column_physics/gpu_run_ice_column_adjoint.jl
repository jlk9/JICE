# Written by Joseph Kump (josek97@utexas.edu)
# Adjoint method code for run_ice_inline.jl


include("./gpu_run_ice_column.jl")

using Enzyme

#= Runs the single column sea ice model and computes the derivatives for all of its fields. The various
    fields of d_jarrays are the partial derivatives of those fields given the gradients for thickness and
    current temperature (ad_h and ad_T).

NOTE:   In cases where there is no snow and empty (thickness 0) snow layers, the partial derivative in terms of
        thickness can be NaN.
=#
function run_ice_column_autodiff(jarrays, atmodels, ad_H_i, ad_T)

    # first we need to allocate d_jarrays and d_atmodels:
    d_jarrays  = initialize_JICEColumnArrays(jarrays.N_t, jarrays.N_c, jarrays.N_i, jarrays.N_s, zeros(Float64, jarrays.N_c), zeros(Float64, jarrays.N_c),
                                             0.0, 0.0, 0.0, zeros(Float64, jarrays.N_c), zeros(Float64, (jarrays.N_i+jarrays.N_s+1)*jarrays.N_c),
                                             jarrays.onGPU)

    d_atmodels = initialize_ATModelArrays(atmodels.N_t, atmodels.N_c, zeros(Float64, atmodels.N_c), zeros(Float64, atmodels.N_c), zeros(Float64, atmodels.N_c),
                                          zeros(Float64, atmodels.N_c), zeros(Float64, atmodels.N_c), zeros(Float64, atmodels.N_c), zeros(Float64, atmodels.N_c),
                                          zeros(Float64, atmodels.N_c), zeros(Float64, atmodels.N_c), zeros(Float64, atmodels.N_c), zeros(Float64, 3atmodels.N_c),
                                          jarrays.onGPU)

    # Want to reset some derivatives to 0:
    d_jarrays.S  .= 0
    d_jarrays.N_t = 0
    d_jarrays.N_c = 0
    d_jarrays.N_i = 0
    d_jarrays.N_s = 0
    
    d_jarrays.T_nplus = ad_T
    d_jarrays.H_i     = ad_H_i

    autodiff(Reverse, run_ice_column, Const, Duplicated(jarrays, d_jarrays), Duplicated(atmodels, d_atmodels))


    return d_jarrays #d_jarrays.Δh, d_jarrays.T_n

end

