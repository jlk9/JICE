# Written by Joseph Kump (josek97@utexas.edu)

include("jicedrv_init.jl")
include("jicedrv_run.jl")

# The main function for our sea ice column physics model:
function jicedrv_main()

    jicedrv_initialize()

    jicedrv_run()

    println("JICE COMPLETED SUCCESSFULLY")

end

jice_main()