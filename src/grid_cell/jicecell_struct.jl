# Written by Joseph Kump (josek97@utexas.edu)
# Outlines JICECell, a Julia structure for a single grid cell in sea ice physics consisting of multiple columns,
# their areas, and shared fields (including an atmosphere model).

include("../column_physics/jicecolumn_struct.jl")
include("../atmosphere_model/atmodel_struct.jl")

#= JICECell struct

=#
struct JICECell

    N_cat::Int64
    N_t::Int64

    T_frz::Float64
    T_w::Float64

    thickness_cats::Vector{Float64}
    areas::Vector{Float64}
    columns::Vector{JICEColumn}

    atm::ATModel

end

# Creates a JICECell object.
# Fields ending with '_cols' are lists of the given data, organized for each column in the JICECell object.
# It is assumed areas and all fields ending with '_cols' are the same length.
function initialize_JICECell(N_cat, N_t, Δt, T_frz, T_w, N_i_cols, N_s_cols, H_i_cols, H_s_cols, u_star_cols, T_0_cols, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a, thickness_cats, areas)

    # TODO: error checking to make sure lengths are consistent

    # Initialize the atmosphere object:
    atm = initialize_ATModel(N_t, F_SWvdr, F_SWidr, F_SWvdf, F_SWidf, F_Ld, T_a, Θ_a, ρ_a, Q_a, c_p, U_a)

    # Create empty array of columns and add each column model to it:
    columns = Vector{JICEColumn}()

    for k in 1:length(areas)

        column = initialize_JICEColumn(N_t, N_i_cols[k], N_s_cols[k], H_i_cols[k], H_s_cols[k], T_frz, Δt, u_star_cols[k], T_w, T_0_cols[k])
        push!(columns, column)
    end

    # Now construct grid cell object:
    grid_cell = JICECell(N_cat, N_t, T_frz, T_w, thickness_cats, areas, columns, atm)
    return grid_cell
end