# Written by Joseph Kump (josek97@utexas.edu)

# Runs and executes the model for ice column physics
function jicedrv_run()
    
    # ANY NECESSARY PRELIMINARIES HERE:
    # Integers
    stop_now = 0
    istep    = 0

    # Floats
    time = 0.0 # Total Elapsed Time (s)
    dt   = 0.0 # Thermodynamics timestep (s)


    # TIMESTEP LOOP:
    while true

        # Physics step:
        ice_step()

        istep = istep + 1

        cond(stop_now >= 1) || break

        # Here is calling query parameters, check warning flush, etc.

    end


end

# Runs the physics for a single step of the ice model
function ice_step()

    # Preliminary work setting variables, determine what is needed for thermodynamics
    #=
    use icedrv_calendar, only: dt, dt_dyn, ndtd, diagfreq, write_restart, istep
      use icedrv_diagnostics, only: runtime_diags, init_mass_diags
!     use icedrv_diagnostics, only: icedrv_diagnostics_debug
      use icedrv_diagnostics_bgc, only: hbrine_diags, zsal_diags, bgc_diags
      use icedrv_flux, only: init_history_therm, init_history_bgc, &
          daidtt, daidtd, dvidtt, dvidtd, dagedtt, dagedtd, init_history_dyn
      use icedrv_history, only: history_cdf, history_write
      use icedrv_restart, only: dumpfile, final_restart
      use icedrv_restart_bgc, only: write_restart_bgc
      use icedrv_step, only: prep_radiation, step_therm1, step_therm2, &
          update_state, step_dyn_ridge, step_snow, step_radiation, &
          biogeochemistry, step_dyn_wave

      integer (kind=int_kind) :: &
         k               ! dynamics supercycling index

      logical (kind=log_kind) :: &
         calc_Tsfc, skl_bgc, solve_zsal, z_tracers, tr_brine, &  ! from icepack
         tr_fsd, wave_spec, tr_snow

      real (kind=dbl_kind) :: &
         offset          ! d(age)/dt time offset

      character(len=*), parameter :: subname='(ice_step)'
    =#

    # Query icepack parameters (icepack_query_parameters) to get needed parameters



    # Initialize diagnostics if needed:
    #=
    call init_mass_diags   ! diagnostics per timestep
      call init_history_therm
      call init_history_bgc
    =#

    # Scale radiation fields
    if calc_Tsfc
        prep_radiation()
    end

    # thermodynamics and biogeochemistry
    step_therm1(dt)     # vertical thermodynamics
    # LATER: biogeochemistry between step_therm1 and step_therm2
    step_therm2(dt)     # ice thickness distribution thermo

    # LATER: dynamics, transport, ridging, snow redistribution and metamorphosis
    #offset = dt
    #update_state(dt, daidtt, dvidtt, dagedtt, offset)

    # albedo, shortwave radiation
    step_radiation(dt)


    # Get ready for coupling / next time step
    coupling_prep()


    # Write output data



end


# Coupling with ocean/atmospheric model, next time step
coupling_prep()



end