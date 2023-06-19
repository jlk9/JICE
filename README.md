# JICE
Sea ice model written in Julia, ported from CICE/Icepack.

This model implements the vertical thermodynamics, atmospheric flux calculations, and growth/melt from Icepack (in ```src/column_physics```). It also implements linear remapping, lateral melting, and new ice formulation as gridcell level processes (in ```src/grid_cell```). The folder ```total_grid``` contains a skeleton for the overall grid structure. Since this model only represents thermodynamics and grid cell / column processes, individual grid cells don't interact with each other.

The folder ```gpu_src``` contains experiments for parallelizing the code with CUDA.jl, and the folder ```test``` contains some basic numeric tests for both the original model and the one modified to work with CUDA.
