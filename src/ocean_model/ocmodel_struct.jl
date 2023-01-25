# Written by Joseph Kump (josek97@utexas.edu)
# Object that provides important adjacent variables for the sea ice model,
# particularly for getting ocean energy flux and ice growth/melt
# Initially, we'll assume values are constant on our time scale

#= Properties
    T_w     (K)         temperature of sea surface
    Frzmlt  (W/m^2)     freezing/melting potential for add new ice
    sea_sal (ppt)       sea surface salinity
=#
mutable struct OCModel

    T_w::Float64
    Frzmlt::Float64
    sea_sal::Float64

end