# Written by Joseph Kump (josek97@utexas.edu)
# Outlines JICECell, a Julia structure for a single grid cell in sea ice physics consisting of multiple columns,
# their areas, and shared fields (including an atmosphere model).

# MAJOR ASSUMPTION: each column has the same number of layers

include("../column_physics/jicecolumn_struct.jl")
include("../atmosphere_model/atmodel_struct.jl")

#= JICECell struct

=#
struct JICECell

    N_cat::Int64
    N_t::Int64

    T_frz::Float64
    T_w::Float64

    frzmlt::Float64 # freezing / melting potential (W/m^2)

    H_bds::Vector{Float64}
    H_bnew::Vector{Float64}
    dH::Vector{Float64}
    areas::Vector{Float64}
    vol_i::Vector{Float64}
    vol_s::Vector{Float64}
    areas_old::Vector{Float64}
    vol_i_old::Vector{Float64}
    vol_s_old::Vector{Float64}
    columns::Vector{JICEColumn}

    g0::Vector{Float64}
    g1::Vector{Float64}
    hL::Vector{Float64}
    hR::Vector{Float64}

    dareas::Vector{Float64}
    dvol_i::Vector{Float64}
    dvol_s::Vector{Float64}
    donor::Vector{Int64}


    i_energy::Vector{Float64}
    s_energy::Vector{Float64}
    i_energy_old::Vector{Float64}
    s_energy_old::Vector{Float64}

    d_area_i_new::Vector{Float64}
    d_area_total::Vector{Float64}

    vi0_new_array::Vector{Float64}

    atm::ATModel

end

# Creates a JICECell object.
# Fields ending with '_cols' are lists of the given data, organized for each column in the JICECell object.
# It is assumed areas and all fields ending with '_cols' are the same length.
function initialize_JICECell(N_cat, N_t, Δt, T_frz, T_w, frzmlt, N_i_cols, N_s_cols, H_i_cols, H_s_cols, u_star_cols, T_0_cols, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, thickness_bds, areas)

    # TODO: error checking to make sure lengths are consistent
    # Want all variables ending with 'cols' and areas to be length N_cat, and H_bds to be length N_cat+1

    # Initialize the atmosphere object:
    atm = initialize_ATModel(N_t, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a)

    # Create empty array of columns and add each column model to it:
    columns = Vector{JICEColumn}()

    for k in 1:N_cat

        column = initialize_JICEColumn(N_t, N_i_cols[k], N_s_cols[k], H_i_cols[k], H_s_cols[k], T_frz, Δt, u_star_cols[k], T_w, T_0_cols[k])
        push!(columns, column)
    end

    H_bnew = zeros(Float64, N_cat+1)
    dH     = zeros(Float64, N_cat)
    vol_i  = zeros(Float64, N_cat)
    vol_s  = zeros(Float64, N_cat)

    areas_old = zeros(Float64, N_cat)
    vol_i_old = zeros(Float64, N_cat)
    vol_s_old = zeros(Float64, N_cat)

    g0 = zeros(Float64, N_cat+1)
    g1 = zeros(Float64, N_cat+1)
    hL = zeros(Float64, N_cat+1)
    hR = zeros(Float64, N_cat+1)

    dareas = zeros(Float64, N_cat)
    dvol_i = zeros(Float64, N_cat)
    dvol_s = zeros(Float64, N_cat)
    donor  = zeros(Int64, N_cat)

    i_energy     = zeros(Float64, N_cat)
    s_energy     = zeros(Float64, N_cat)
    i_energy_old = zeros(Float64, N_cat)
    s_energy_old = zeros(Float64, N_cat)

    d_area_i_new  = zeros(Float64, N_cat)
    d_area_total  = zeros(Float64, N_cat)
    vi0_new_array = zeros(Float64, N_cat)

    # Now construct grid cell object:
    grid_cell = JICECell(N_cat, N_t, T_frz, T_w, frzmlt, thickness_bds, H_bnew, dH, areas, vol_i, vol_s, areas_old, vol_i_old, vol_s_old, columns,
                         g0, g1, hL, hR, dareas, dvol_i, dvol_s, donor, i_energy, s_energy, i_energy_old, s_energy_old,
                         d_area_i_new, d_area_total, vi0_new_array, atm)
    return grid_cell
end