# Written by Joseph Kump (josek97@utexas.edu)
# Outlines JICECell, a Julia structure for a single grid cell in sea ice physics consisting of multiple columns,
# their areas, and shared fields (including an atmosphere model).

# MAJOR ASSUMPTION: each column has the same number of layers

include("../column_physics/jicecolumn_struct.jl")
include("../atmosphere_model/atmodel_struct.jl")

const sea_sal          = 35.0    # typical sea salt in ppt
const dSin0_frazil     = 3.0     # bulk salinity reduction of newly formed frazil
const ice_ref_salinity = 4.0     # reference salinity of sea ice, ppt

#= JICECell struct
Properties:
    N_cat (dim'less)   Number of thickness categories, int
    N_t (dim'less) Number of time steps in model, int

    T_frz      (C)          freezing point of water at salinity S, float
    T_w        (C)          Temperature of water at ocean surface, float
    frzmlt     (W / m^2)    freezing / melting potential of ice,   float
    rside      (dim'less)   fraction of ice that metls laterally,  float
    fresh_flux ()           freshwater flux to the ocean, float
    salt_flux  ()           saltwater flux to the ocean,  float
    heat_flux  ()           heat flux to the ocean,       float

    H_bds       (m)         upper and lower thickness bounds for each thickness category,                            Vector{Float64} of length N_cat+1
    H_bnew      (m)         upper and lower thickness bounds for temporary categories for linear remapping,          Vector{Float64} of length N_cat+1
    dH          (m)         amount the ice column height H_i changed in each category after vertical thermodynamics, Vector{Float64} of length N_cat
    areas       (dim'less)  percent area of the gridcell covered by each thickness category,                         Vector{Float64} of length N_cat
    vol_i       (m^3)       volume of ice in each thickness category,                                                Vector{Float64} of length N_cat
    vol_s       (m^3)       volume of snow in each thickness category,                                               Vector{Float64} of length N_cat
    areas_old   (dim'less)  old areas for conservation computations,                                                 Vector{Float64} of length N_cat
    vol_i_old   (m^3)       old ice volumes for conservation computations,                                           Vector{Float64} of length N_cat
    vol_s_old   (m^3)       old snow volumes for conservation computations,                                          Vector{Float64} of length N_cat
    columns     (a lot)     the ice columns that represent each thickness category,                                  Vector{JICEColumn} of length N_cat

    (NOTE: total area is in m^2 but it is fixed for the grid cell)

    # Variables for the fit_line helper in linear remapping, using a linear approximation of the thickness distribution function
    g0 () intercept of linear approximation, Vector{Float64}
    g1 () slope of linear approximation,     Vector{Float64}
    hL () new category lower bound,          Vector{Float64}
    hR () new category upper bound,          Vector{Float64}

    dareas (dim'less)   change in area after linear remapping, Vector{Float64}
    dvol_i (m^3)        change in ice volume after linear remapping, Vector{Float64}
    dvol_s (m^3)        change in snow volume after linear remapping, Vector{Float64}
    donor  (dim'less)   source of new ice/snow for each column after remapping, Vector{Int}


    i_energy     (E)    ice energy after grid cell process (linear remapping, adding new ice),   Vector{Float64}
    s_energy     (E)    snow energy after grid cell process (linear remapping),                  Vector{Float64}
    i_energy_old (E)    ice energy before grid cell process (linear remapping, adding new ice),  Vector{Float64}
    s_energy_old (E)    snow energy before grid cell process (linear remapping),                 Vector{Float64}

    d_area_i_new (dim'less)     ice area from new growth to add to each thickness category, Vector{Float64}
    d_area_total (dim'less)     total ice area to add to each thickness category, Vector{Float64}

    (NOTE: for now d_area_total = d_area_i_new, but if fsd is included then it will also include lateral growth)

    vi0_new_array (m^3)         volume of new ice to add into each category in add_new_ice,             Vector{Float64}
    rside_array   (dim'less)    fraction of ice melted from each thickness category at each time step,  Vector{Float64}

    atm (a lot)     accompanying atmosphere "model" with relevant atmosphere variables, ATModel
=#
mutable struct JICECell

    N_cat::Int64
    N_t::Int64

    T_frz::Float64
    T_w::Float64

    frzmlt::Float64 # freezing / melting potential (W/m^2)
    rside::Float64 # fraction of ice that metls laterally

    # freshwater, saltwater, and heat fluxes to the ocean
    fresh_flux::Float64
    salt_flux::Float64
    heat_flux::Float64

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
    rside_array::Vector{Float64}

    atm::ATModel

end

# Creates a JICECell object.
# Fields ending with '_cols' are lists of the given data, organized for each column in the JICECell object.
# It is assumed areas and all fields ending with '_cols' are the same length.
function initialize_JICECell(N_cat, N_t, Δt, T_frz, T_w, frzmlt, rside, N_i_cols, N_s_cols, H_i_cols, H_s_cols, u_star_cols, T_0_cols, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, thickness_bds, areas)

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
    rside_array   = zeros(Float64, N_cat)

    # Now construct grid cell object:
    grid_cell = JICECell(N_cat, N_t, T_frz, T_w, frzmlt, rside, 0.0, 0.0, 0.0, thickness_bds, H_bnew, dH, areas, vol_i, vol_s, areas_old, vol_i_old, vol_s_old, columns,
                         g0, g1, hL, hR, dareas, dvol_i, dvol_s, donor, i_energy, s_energy, i_energy_old, s_energy_old,
                         d_area_i_new, d_area_total, vi0_new_array, rside_array, atm)
    return grid_cell
end