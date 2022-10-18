# Written by Joseph Kump (josek97@utexas.edu)

# Initializes parameters for ice column physics
function jicedrv_initialize()
    
    # Variables used for other packages, figure out what we need:
    #=
    use icedrv_arrays_column, only: hin_max, c_hi_range
      use icedrv_arrays_column, only: floe_rad_l, floe_rad_c, &
          floe_binwidth, c_fsd_range
      use icedrv_calendar, only: dt, time, istep, istep1, &
          init_calendar, calendar
      use icepack_intfc, only: icepack_init_itd, icepack_init_itd_hist
      use icepack_intfc, only: icepack_init_fsd_bounds
      use icepack_intfc, only: icepack_init_snow
      use icepack_intfc, only: icepack_warnings_flush
      use icedrv_domain_size, only: ncat, nfsd
!     use icedrv_diagnostics, only: icedrv_diagnostics_debug
      use icedrv_flux, only: init_coupler_flux, init_history_therm, &
          init_flux_atm_ocn
      use icedrv_forcing, only: init_forcing, get_forcing, get_wave_spec
      use icedrv_forcing_bgc, only: get_forcing_bgc, faero_default, fiso_default, init_forcing_bgc 
      use icedrv_restart_shared, only: restart
      use icedrv_init, only: input_data, init_state, init_grid2, init_fsd
      use icedrv_init_column, only: init_thermo_vertical, init_shortwave, init_zbgc
      use icepack_intfc, only: icepack_configure
      =#

    # Logical variables initialized here
    #=
    logical (kind=log_kind) :: &
         skl_bgc, &    ! from icepack
         z_tracers, &  ! from icepack
         tr_snow, &    ! from icepack
         tr_aero, &    ! from icepack
         tr_iso, &     ! from icepack
         tr_zaero, &   ! from icepack
         tr_fsd, wave_spec
    =#

    # Not sure what all we need to initialize for this stripped down model, but we'll put it here:

end