# Written by Joseph Kump (josek97@utexas.edu)
# Code that runs the adjoint of the time step for a single grid cell object in
# our sea ice model. Ideally we want to apply enzyme to the entire process as
# we do for the single-column case.

include("./run_grid_cell.jl")

using Enzyme


function run_grid_cell_autodiff(jcell, ad_H_i_cols, ad_T_cols)

    d_H_i_cols      = zeros(Float64, jcell.N_cat)
    d_H_s_cols      = zeros(Float64, jcell.N_cat)
    d_u_star_cols   = zeros(Float64, jcell.N_cat)
    d_thickness_bds = zeros(Float64, jcell.N_cat+1)
    d_areas         = zeros(Float64, jcell.N_cat)

    d_T_0_cols    = Vector{Vector{Float64}}()
    for n in 1:jcell.N_cat
        push!(d_T_0_cols, zeros(Float64, jcell.N_i+jcell.N_s+1))
    end

    # First we need to allocate memory for the derivatives of jcell:
    d_jcell = initialize_JICECell(jcell.N_cat, jcell.N_t, 0.0, 0.0, 0.0, 0.0, 0.0, jcell.N_i, jcell.N_s,
                                    d_H_i_cols, d_H_s_cols, d_u_star_cols, d_T_0_cols,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, zeros(Float64,3), d_thickness_bds, d_areas)

    # Now we zero-out nonzero terms in d_jcell:
    d_jcell.N_cat = 0
    d_jcell.N_t   = 0
    d_jcell.N_i   = 0
    d_jcell.N_s   = 0

    # And add our gradient values of interest:
    for n in 1:jcell.N_cat
        d_jcell.columns[n].H_i     = ad_H_i_cols[n]
        d_jcell.columns[n].T_nplus = ad_T_cols[n]
    end

    # And run the model with autodiff:
    autodiff(run_ice_cell, Const, Duplicated(jcell, d_jcell))

    return d_jcell
end