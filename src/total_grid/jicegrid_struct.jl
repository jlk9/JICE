# Written by Joseph Kump (josek97@utexas.edu)
# Outlines JICEGrid, a Julia structure for a grid of JICE grid cells, organized in a
# longitude/latitude storage format

include("../grid_cell/jicecell_struct.jl")

#= JICEGrid struct
Properties:

=#
mutable struct JICEGrid

    cells::Vector{JICECell}

end