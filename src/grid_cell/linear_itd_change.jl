# Written by Joseph Kump (josek97@utexas.edu)
# Code that computes the remapping of ice column areas
# between thickness categories within one grid cell

include("./jicecell_struct.jl")

@inline function linear_itd_change(jcell)

    # Compute thickness changes in each category
    thickness_changes = zeros(Float64, jcell.N_cat)
    for k in 1:jcell.N_cat

        thickness_changes[k] = jcell.columns[k].H_i - jcell.columns[k].H_iold
    end

    # Compute new category boundaries


    # Boundary check (each one lies between adjacen values of H_i)

    return nothing
end