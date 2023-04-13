# Written by Joseph Kump (josek97@utexas.edu)
# Basic code that is just ice temperature in column, with no formulation or melting.
# This version eliminates vectorized operations, enabling use of Enzyme.jl
# It also reduces the number of memory allocations substantially, increasing runtime
# and lowering memory footprint.

include("./gpu_run_column_step.jl")

#= Model Function
Input:
    jcolumn         an initialized JICEColumn object
    atmodel         an initialized ATModel object for surface flux computations

Output:
    Technically nothing, but the function updates jcolumn's T_array and Δh_arrays to have a log
    of updated ice temperatures and thicknesses
=#
@time function run_ice_column(jarrays::JICEColumnArrays, atmodels::ATModelArrays)
    
    # Main loop of temperature modifications:
    for step in 1:jarrays.N_t
        
        run_column_step(jarrays, atmodels, step)
    end
    
    return nothing
end