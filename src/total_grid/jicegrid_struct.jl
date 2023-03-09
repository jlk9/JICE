# Written by Joseph Kump (josek97@utexas.edu)
# Outlines JICEGrid, a Julia structure for a grid of JICE grid cells, organized in a
# longitude/latitude storage format, longitude major.

include("../grid_cell/jicecell_struct.jl")

using Oceananigans

#= JICEGrid struct
Properties:

cell_count   total number of grid cells
cells        array of JICECell objects
latlong      LatitudeLongitudeGrid object from Oceananigans to track positional coordinates

=#
mutable struct JICEGrid

    cell_count::Int64

    cells::Vector{JICECell}

    latlong::LatitudeLongitudeGrid{}

end