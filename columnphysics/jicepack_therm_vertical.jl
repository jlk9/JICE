# Written by Joseph Kump (josek97@utexas.edu)

#=
! Update ice and snow internal temperatures and compute
! thermodynamic growth rates and atmospheric fluxes.
!
! NOTE: The thermodynamic calculation is split in two for load balancing.
!       First icepack_therm_vertical computes vertical growth rates and coupler
!       fluxes.  Then icepack_therm_itd does thermodynamic calculations not
!       needed for coupling.
=#

#=
Driver for thermodynamic changes not needed for coupling:
transport in thickness space, lateral growth and melting.

TODO: figure out external variables as always, decide what to do with allocating optional arguments
    and optional parts like aerosols and melt ponds

Variables:

integer (kind=int_kind), intent(in) :: &
         ncat    , & ! number of thickness categories
         nilyr   , & ! number of ice layers
         nslyr       ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         dt          , & ! time step
         uvel        , & ! x-component of velocity (m/s)
         vvel        , & ! y-component of velocity (m/s)
         strax       , & ! wind stress components (N/m^2)
         stray       , & ! 
         yday            ! day of year

      logical (kind=log_kind), intent(in) :: &
         lmask_n     , & ! northern hemisphere mask
         lmask_s         ! southern hemisphere mask

      logical (kind=log_kind), intent(in), optional :: &
         prescribed_ice  ! if .true., use prescribed ice instead of computed

      real (kind=dbl_kind), intent(inout) :: &
         aice        , & ! sea ice concentration
         vice        , & ! volume per unit area of ice          (m)
         vsno        , & ! volume per unit area of snow         (m)
         zlvl        , & ! atm level height for momentum (and scalars if zlvs is not present) (m)
         uatm        , & ! wind velocity components (m/s)
         vatm        , &
         wind        , & ! wind speed (m/s)
         potT        , & ! air potential temperature  (K)
         Tair        , & ! air temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         rhoa        , & ! air density (kg/m^3)
         frain       , & ! rainfall rate (kg/m^2 s)
         fsnow       , & ! snowfall rate (kg/m^2 s)
         fpond       , & ! fresh water flux to ponds (kg/m^2/s)
         fresh       , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt       , & ! salt flux to ocean (kg/m^2/s)
         fhocn       , & ! net heat flux to ocean (W/m^2)
         fswthru     , & ! shortwave penetrating to ocean (W/m^2)
         fsurf       , & ! net surface heat flux (excluding fcondtop)(W/m^2)
         fcondtop    , & ! top surface conductive flux        (W/m^2)
         fcondbot    , & ! bottom surface conductive flux     (W/m^2)
         fsens       , & ! sensible heat flux (W/m^2)
         flat        , & ! latent heat flux   (W/m^2)
         fswabs      , & ! shortwave flux absorbed in ice and ocean (W/m^2)
         flw         , & ! incoming longwave radiation (W/m^2)
         flwout      , & ! outgoing longwave radiation (W/m^2)
         evap        , & ! evaporative water flux (kg/m^2/s)
         evaps       , & ! evaporative water flux over snow (kg/m^2/s)
         evapi       , & ! evaporative water flux over ice (kg/m^2/s)
         congel      , & ! basal ice growth         (m/step-->cm/day)
         snoice      , & ! snow-ice formation       (m/step-->cm/day)
         Tref        , & ! 2m atm reference temperature (K)
         Qref        , & ! 2m atm reference spec humidity (kg/kg)
         Uref        , & ! 10m atm reference wind speed (m/s)
         Cdn_atm     , & ! atm drag coefficient
         Cdn_ocn     , & ! ocn drag coefficient
         hfreebd     , & ! freeboard (m)
         hdraft      , & ! draft of ice + snow column (Stoessel1993)
         hridge      , & ! ridge height
         distrdg     , & ! distance between ridges
         hkeel       , & ! keel depth
         dkeel       , & ! distance between keels
         lfloe       , & ! floe length
         dfloe       , & ! distance between floes
         Cdn_atm_skin, & ! neutral skin drag coefficient
         Cdn_atm_floe, & ! neutral floe edge drag coefficient
         Cdn_atm_pond, & ! neutral pond edge drag coefficient
         Cdn_atm_rdg , & ! neutral ridge drag coefficient
         Cdn_ocn_skin, & ! skin drag coefficient
         Cdn_ocn_floe, & ! floe edge drag coefficient
         Cdn_ocn_keel, & ! keel drag coefficient
         Cdn_atm_ratio,& ! ratio drag atm / neutral drag atm
         strairxT    , & ! stress on ice by air, x-direction
         strairyT    , & ! stress on ice by air, y-direction
         strocnxT    , & ! ice-ocean stress, x-direction
         strocnyT    , & ! ice-ocean stress, y-direction
         fbot        , & ! ice-ocean heat flux at bottom surface (W/m^2)
         frzmlt      , & ! freezing/melting potential (W/m^2)
         rside       , & ! fraction of ice that melts laterally
         fside       , & ! lateral heat flux (W/m^2)
         sst         , & ! sea surface temperature (C)
         Tf          , & ! freezing temperature (C)
         Tbot        , & ! ice bottom surface temperature (deg C)
         Tsnice      , & ! snow ice interface temperature (deg C)
         sss         , & ! sea surface salinity (ppt)
         meltt       , & ! top ice melt             (m/step-->cm/day)
         melts       , & ! snow melt                (m/step-->cm/day)
         meltb       , & ! basal ice melt           (m/step-->cm/day)
         mlt_onset   , & ! day of year that sfc melting begins
         frz_onset       ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), intent(inout), optional :: &
         fswthru_vdr , & ! vis dir shortwave penetrating to ocean (W/m^2)
         fswthru_vdf , & ! vis dif shortwave penetrating to ocean (W/m^2)
         fswthru_idr , & ! nir dir shortwave penetrating to ocean (W/m^2)
         fswthru_idf , & ! nir dif shortwave penetrating to ocean (W/m^2)
         dsnow       , & ! change in snow depth     (m/step-->cm/day)
         meltsliq    , & ! mass of snow melt (kg/m^2)
         fsloss          ! rate of snow loss to leads (kg/m^2/s)

      real (kind=dbl_kind), dimension(:), optional, intent(inout) :: &
         Qa_iso      , & ! isotope specific humidity (kg/kg)
         Qref_iso    , & ! isotope 2m atm reference spec humidity (kg/kg)
         fiso_atm    , & ! isotope deposition rate (kg/m^2 s)
         fiso_ocn    , & ! isotope flux to ocean  (kg/m^2/s)
         fiso_evap   , & ! isotope evaporation (kg/m^2/s)
         meltsliqn       ! mass of snow melt (kg/m^2)

      real (kind=dbl_kind), dimension(:,:), optional, intent(inout) :: &
         rsnwn       , & ! snow grain radius (10^-6 m)
         smicen      , & ! tracer for mass of ice in snow (kg/m^3)
         smliqn          ! tracer for mass of liquid in snow (kg/m^3)

      real (kind=dbl_kind), optional, intent(in) :: &
         HDO_ocn     , & ! ocean concentration of HDO (kg/kg)
         H2_16O_ocn  , & ! ocean concentration of H2_16O (kg/kg)
         H2_18O_ocn  , & ! ocean concentration of H2_18O (kg/kg)
         zlvs            ! atm level height for scalars (if different than zlvl) (m)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         aicen_init  , & ! fractional area of ice
         vicen_init  , & ! volume per unit area of ice (m)
         vsnon_init  , & ! volume per unit area of snow (m)
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         vsnon       , & ! volume per unit area of snow         (m)
         Tsfc        , & ! ice/snow surface temperature, Tsfcn
         alvl        , & ! level ice area fraction
         vlvl        , & ! level ice volume fraction
         apnd        , & ! melt pond area fraction
         hpnd        , & ! melt pond depth (m)
         ipnd        , & ! melt pond refrozen lid thickness (m)
         iage        , & ! volume-weighted ice age
         FY          , & ! area-weighted first-year ice area
         fsurfn      , & ! net flux to top surface, excluding fcondtop
         fcondtopn   , & ! downward cond flux at top surface (W m-2)
         fcondbotn   , & ! downward cond flux at bottom surface (W m-2)
         flatn       , & ! latent heat flux (W m-2)
         fsensn      , & ! sensible heat flux (W m-2)
         fsurfn_f    , & ! net flux to top surface, excluding fcondtop
         fcondtopn_f , & ! downward cond flux at top surface (W m-2)
         flatn_f     , & ! latent heat flux (W m-2)
         fsensn_f    , & ! sensible heat flux (W m-2)
         fswsfcn     , & ! SW absorbed at ice/snow surface (W m-2)
         fswthrun    , & ! SW through ice to ocean            (W/m^2)
         fswintn     , & ! SW absorbed in ice interior, below surface (W m-2)
         faero_atm   , & ! aerosol deposition rate (kg/m^2 s)
         faero_ocn   , & ! aerosol flux to ocean  (kg/m^2/s)
         dhsn        , & ! depth difference for snow on sea ice and pond ice
         ffracn      , & ! fraction of fsurfn used to melt ipond
         meltsn      , & ! snow melt                       (m)
         melttn      , & ! top ice melt                    (m)
         meltbn      , & ! bottom ice melt                 (m)
         congeln     , & ! congelation ice growth          (m)
         snoicen     , & ! snow-ice growth                 (m)
         dsnown          ! change in snow thickness (m/step-->cm/day)

      real (kind=dbl_kind), optional, dimension(:), intent(inout) :: &
         fswthrun_vdr , & ! vis dir SW through ice to ocean            (W/m^2)
         fswthrun_vdf , & ! vis dif SW through ice to ocean            (W/m^2)
         fswthrun_idr , & ! nir dir SW through ice to ocean            (W/m^2)
         fswthrun_idf     ! nir dif SW through ice to ocean            (W/m^2)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         zqsn        , & ! snow layer enthalpy (J m-3)
         zqin        , & ! ice layer enthalpy (J m-3)
         zSin        , & ! internal ice layer salinities
         Sswabsn     , & ! SW radiation absorbed in snow layers (W m-2)
         Iswabsn         ! SW radiation absorbed in ice layers (W m-2)

      real (kind=dbl_kind), dimension(:,:,:), intent(inout) :: &
         aerosno    , &  ! snow aerosol tracer (kg/m^2)
         aeroice         ! ice aerosol tracer (kg/m^2)

      real (kind=dbl_kind), dimension(:,:), optional, intent(inout) :: &
         isosno     , &  ! snow isotope tracer (kg/m^2)
         isoice          ! ice isotope tracer (kg/m^2)
=#
function icepack_step_therm1(dt, ncat, nilyr, nslyr, aicen_init, vicen_init, vsnon_init,
    aice, aicen,
    vice, vicen,
    vsno, vsnon,
    uvel, vvel,
    Tsfc, zqsn,
    zqin, zSin,
    alvl, vlvl,
    apnd, hpnd,
    ipnd,
    iage, FY,
    aerosno, aeroice,
    isosno, isoice,
    uatm, vatm,
    wind, zlvl,
    Qa, rhoa,
    Qa_iso,
    Tair, Tref,
    Qref, Uref,
    Qref_iso,
    Cdn_atm_ratio,
    Cdn_ocn, Cdn_ocn_skin,
    Cdn_ocn_floe, Cdn_ocn_keel,
    Cdn_atm, Cdn_atm_skin,
    Cdn_atm_floe, Cdn_atm_pond,
    Cdn_atm_rdg, hfreebd,
    hdraft, hridge,
    distrdg, hkeel,
    dkeel, lfloe,
    dfloe,
    strax, stray,
    strairxT, strairyT,
    potT, sst,
    sss, Tf,
    strocnxT, strocnyT,
    fbot,
    Tbot, Tsnice,
    frzmlt, rside,
    fside,
    fsnow, frain,
    fpond, fsloss,
    fsurf, fsurfn,
    fcondtop, fcondtopn,
    fcondbot, fcondbotn,
    fswsfcn, fswintn,
    fswthrun,
    fswthrun_vdr,
    fswthrun_vdf,
    fswthrun_idr,
    fswthrun_idf,
    fswabs,
    flwout,
    Sswabsn, Iswabsn,
    flw,
    fsens, fsensn,
    flat, flatn,
    evap,
    evaps, evapi,
    fresh, fsalt,
    fhocn,
    fswthru,
    fswthru_vdr,
    fswthru_vdf,
    fswthru_idr,
    fswthru_idf,
    flatn_f, fsensn_f,
    fsurfn_f, fcondtopn_f,
    faero_atm, faero_ocn,
    fiso_atm, fiso_ocn,
    fiso_evap,
    HDO_ocn, H2_16O_ocn,
    H2_18O_ocn,
    dhsn, ffracn,
    meltt, melttn,
    meltb, meltbn,
    melts, meltsn,
    congel, congeln,
    snoice, snoicen,
    dsnow, dsnown,
    meltsliq, meltsliqn,
    rsnwn,
    smicen, smliqn,
    lmask_n, lmask_s,
    mlt_onset, frz_onset,
    yday, prescribed_ice,
    zlvs)

    # local variables
    rnslyr = 1 / nslyr
    worka  = 0.0
    workb  = 0.0

    # Float 2D coupler variables (computed for each category, then aggregated)
    fswabsn         = 0.0  # shortwave absorbed by ice          (W/m^2)
    flwoutn         = 0.0  # upward LW at surface               (W/m^2)
    evapn           = 0.0  # flux of vapor, atmos to ice   (kg m-2 s-1)
    evapsn          = 0.0  # flux of vapor, atmos to ice over snow  (kg m-2 s-1)
    evapin          = 0.0  # flux of vapor, atmos to ice over ice  (kg m-2 s-1)
    freshn          = 0.0  # flux of water, ice to ocean     (kg/m^2/s)
    fsaltn          = 0.0  # flux of salt, ice to ocean      (kg/m^2/s)
    fhocnn          = 0.0  # fbot corrected for leftover energy (W/m^2)
    strairxn        = 0.0  # air/ice zonal  stress,             (N/m^2)
    strairyn        = 0.0  # air/ice meridional stress,         (N/m^2)
    Cdn_atm_ratio_n = 0.0  # drag coefficient ratio
    Trefn           = 0.0  # air tmp reference level                (K)
    Urefn           = 0.0  # air speed reference level            (m/s)
    Qrefn           = 0.0  # air sp hum reference level         (kg/kg)
    shcoef          = 0.0  # transfer coefficient for sensible heat
    lhcoef          = 0.0  # transfer coefficient for latent heat
    rfrac           = 0.0  # water fraction retained for melt ponds

    # Array of floats, dimension(nslyr,ncat)
    massicen = zeros(Float64, nslyr, ncat) # mass of ice in snow (kg/m^2)
    massliqn = zeros(Float64, nslyr, ncat) # mass of liquid in snow (kg/m^2)

    # Array of floats, dimension(n_iso)
    Qrefn_iso  = zeros(Float64, n_iso) # isotope air sp hum reference level (kg/kg)
    fiso_ocnn  = zeros(Float64, n_iso) # isotope flux to ocean  (kg/m^2/s)
    fiso_evapn = zeros(Float64, n_iso) # isotope evaporation (kg/m^2/s)

    # Allocatable, might be optional
    #=
    real (kind=dbl_kind), allocatable, dimension(:,:) :: &
         l_isosno   , &  ! local snow isotope tracer (kg/m^2)
         l_isoice        ! local ice isotope tracer (kg/m^2)

      real (kind=dbl_kind), allocatable, dimension(:) :: &
         l_Qa_iso    , & ! local isotope specific humidity (kg/kg)
         l_Qref_iso  , & ! local isotope 2m atm reference spec humidity (kg/kg)
         l_fiso_atm  , & ! local isotope deposition rate (kg/m^2 s)
         l_fiso_ocn  , & ! local isotope flux to ocean  (kg/m^2/s)
         l_fiso_evap , & ! local isotope evaporation (kg/m^2/s)
         l_meltsliqn     ! mass of snow melt (kg/m^2)

      real (kind=dbl_kind), allocatable, dimension(:,:) :: &
         l_rsnw      , & ! snow grain radius (10^-6 m)
         l_smice     , & ! tracer for mass of ice in snow (kg/m^3)
         l_smliq         ! tracer for mass of liquid in snow (kg/m^3)
    =#

    # More floats
    l_fsloss     = 0.0  # rate of snow loss to leads (kg/m^2/s)
    l_meltsliq   = 0.0  # mass of snow melt (kg/m^2)
    l_HDO_ocn    = 0.0  # local ocean concentration of HDO (kg/kg)
    l_H2_16O_ocn = 0.0  # local ocean concentration of H2_16O (kg/kg)
    l_H2_18O_ocn = 0.0  # local ocean concentration of H2_18O (kg/kg)

    # Even more floats
    l_fswthru_vdr = 0.0  # vis dir SW through ice to ocean            (W/m^2)
    l_fswthru_vdf = 0.0  # vis dif SW through ice to ocean            (W/m^2)
    l_fswthru_idr = 0.0  # nir dir SW through ice to ocean            (W/m^2)
    l_fswthru_idf = 0.0  # nir dif SW through ice to ocean            (W/m^2)

    # Another optional (?) allocatable
    #=
    real (kind=dbl_kind), dimension(:), allocatable :: &
         l_fswthrun_vdr , & ! vis dir SW through ice to ocean            (W/m^2)
         l_fswthrun_vdf , & ! vis dif SW through ice to ocean            (W/m^2)
         l_fswthrun_idr , & ! nir dir SW through ice to ocean            (W/m^2)
         l_fswthrun_idf     ! nir dif SW through ice to ocean            (W/m^2)
    =#

    # water retained in ponds (m)
    pond = 0.0

    # Allocate local optional arguments
    #=
          if (present(isosno)    ) then
         allocate(l_isosno(size(isosno,dim=1),size(isosno,dim=2)))
         l_isosno     = isosno
      else
         allocate(l_isosno(1,1))
         l_isosno     = c0
      endif

      if (present(isoice)    ) then
         allocate(l_isoice(size(isoice,dim=1),size(isoice,dim=2)))
         l_isoice     = isoice
      else
         allocate(l_isoice(1,1))
         l_isoice     = c0
      endif

      if (present(Qa_iso)    ) then
         allocate(l_Qa_iso(size(Qa_iso)))
         l_Qa_iso     = Qa_iso
      else
         allocate(l_Qa_iso(1))
         l_Qa_iso     = c0
      endif

      if (present(Qref_iso)    ) then
         allocate(l_Qref_iso(size(Qref_iso)))
         l_Qref_iso     = Qref_iso
      else
         allocate(l_Qref_iso(1))
         l_Qref_iso     = c0
      endif

      if (present(fiso_atm)  ) then
         allocate(l_fiso_atm(size(fiso_atm)))
         l_fiso_atm = fiso_atm
      else
         allocate(l_fiso_atm(1))
         l_fiso_atm   = c0
      endif

      if (present(fiso_ocn)  ) then
         allocate(l_fiso_ocn(size(fiso_ocn)))
         l_fiso_ocn = fiso_ocn
      else
         allocate(l_fiso_ocn(1))
         l_fiso_ocn   = c0
      endif

      if (present(fiso_evap)  ) then
         allocate(l_fiso_evap(size(fiso_evap)))
         l_fiso_evap = fiso_evap
      else
         allocate(l_fiso_evap(1))
         l_fiso_evap   = c0
      endif

      l_fsloss     = c0
      if (present(fsloss)    ) l_fsloss     = fsloss

      l_HDO_ocn    = c0
      if (present(HDO_ocn)   ) l_HDO_ocn    = HDO_ocn

      l_H2_16O_ocn = c0
      if (present(H2_16O_ocn)) l_H2_16O_ocn = H2_16O_ocn

      l_H2_18O_ocn = c0
      if (present(H2_18O_ocn)) l_H2_18O_ocn = H2_18O_ocn

      l_fswthru_vdr = c0
      if (present(fswthru_vdr)) l_fswthru_vdr = fswthru_vdr

      l_fswthru_vdf = c0
      if (present(fswthru_vdf)) l_fswthru_vdf = fswthru_vdf

      l_fswthru_idr = c0
      if (present(fswthru_idr)) l_fswthru_idr = fswthru_idr

      l_fswthru_idf = c0
      if (present(fswthru_idf)) l_fswthru_idf = fswthru_idf

      allocate(l_fswthrun_vdr(ncat))
      l_fswthrun_vdr = c0
      if (present(fswthrun_vdr)) l_fswthrun_vdr = fswthrun_vdr

      allocate(l_fswthrun_vdf(ncat))
      l_fswthrun_vdf = c0
      if (present(fswthrun_vdf)) l_fswthrun_vdf = fswthrun_vdf

      allocate(l_fswthrun_idr(ncat))
      l_fswthrun_idr = c0
      if (present(fswthrun_idr)) l_fswthrun_idr = fswthrun_idr

      allocate(l_fswthrun_idf(ncat))
      l_fswthrun_idf = c0
      if (present(fswthrun_idf)) l_fswthrun_idf = fswthrun_idf

      allocate(l_meltsliqn(ncat))
      l_meltsliqn = c0
      if (present(meltsliqn)) l_meltsliqn = meltsliqn
      l_meltsliq  = c0
      if (present(meltsliq )) l_meltsliq  = meltsliq

      allocate(l_rsnw(nslyr,ncat))
      l_rsnw = rsnw_fall
      if (present(rsnwn)) l_rsnw = rsnwn

      allocate(l_smice(nslyr,ncat))
      l_smice = c0
      if (present(smicen)) l_smice = smicen

      allocate(l_smliq(nslyr,ncat))
      l_smliq = c0
      if (present(smliqn)) l_smliq = smliqn
    =#

    #=
    !-----------------------------------------------------------------
    ! Done with allocating, onto model
    !-----------------------------------------------------------------
    =#

    # Initialize rate of snow loss to leads
    l_fsloss = fsnow * (c1 - aice)

    # snow redistribution using snwlvlfac:  precip factor
    if (snwredist == 'bulk')
        worka = c0
        if (aice > puny)
            for n = 1:ncat
                worka = worka + alvl[n]*aicen[n]
            end
            worka  = worka * (snwlvlfac/(c1+snwlvlfac)) / aice
        end
        l_fsloss = l_fsloss + fsnow*worka
        fsnow    = fsnow*(c1-worka)
    end # snwredist

    # solid and liquid components of snow mass
    massicen[:,:] = c0
    massliqn[:,:] = c0
    if (snwgrain)
        rnslyr = c1 / nslyr
        for n = 1:ncat
            for k = 1:nslyr
                massicen[k,n] = l_smice[k,n] * vsnon[n] * rnslyr # kg/m^2
                massliqn[k,n] = l_smliq[k,n] * vsnon[n] * rnslyr
            end
        end
    end

    #=
    Update the neutral drag coefficients to account for form drag
    Oceanic and atmospheric drag coefficients
    =#
    if (formdrag)
        neutral_drag_coeffs(apnd,
                            hpnd, ipnd,
                            alvl, vlvl,
                            aice, vice,
                            vsno, aicen,
                            vicen,
                            Cdn_ocn, Cdn_ocn_skin,
                            Cdn_ocn_floe, Cdn_ocn_keel,
                            Cdn_atm, Cdn_atm_skin,
                            Cdn_atm_floe, Cdn_atm_pond,
                            Cdn_atm_rdg, hfreebd,
                            hdraft, hridge,
                            distrdg, hkeel,
                            dkeel, lfloe,
                            dfloe, ncat)
        # if (icepack_warnings_aborted(subname)) return
    end

    #=
    Adjust frzmlt to account for ice-ocean heat fluxes since last
      !  call to coupler.
      ! Compute lateral and bottom heat fluxes.
    =#
    frzmlt_bottom_lateral(dt, ncat,
                            nilyr, nslyr,
                            aice, frzmlt,
                            vicen, vsnon,
                            zqin, zqsn,
                            sst, Tf,
                            ustar_min,
                            fbot_xfer_type,
                            strocnxT, strocnyT,
                            Tbot, fbot,
                            rside, Cdn_ocn,
                            fside)

    # if (icepack_warnings_aborted(subname)) return

    for n = 1:ncat
        meltsn[n]    = c0
        melttn[n]    = c0
        meltbn[n]    = c0
        congeln[n]   = c0
        snoicen[n]   = c0
        dsnown[n]    = c0
        meltsliqn[n] = c0

        Trefn  = c0
        Qrefn  = c0
        Qrefn_iso[:]  = c0
        fiso_ocnn[:]  = c0
        fiso_evapn[:] = c0
        Urefn  = c0
        lhcoef = c0
        shcoef = c0
        worka  = c0
        workb  = c0

    if (aicen_init[n] > puny)
        if (calc_Tsfc || calc_strair)

        #=
        ! Atmosphere boundary layer calculation; compute coefficients
        ! for sensible and latent heat fluxes.
        !
        ! NOTE: The wind stress is computed here for later use if 
        !       calc_strair = .true.   Otherwise, the wind stress
        !       components are set to the data values.
        =#
 
            jicepack_atm_boundary('ice',
                                Tsfc[n],  potT,
                                uatm,     vatm,
                                wind,     zlvl,
                                Qa,       rhoa,
                                strairxn, strairyn,
                                Trefn,    Qrefn,
                                worka,    workb,
                                lhcoef,   shcoef,
                                Cdn_atm,
                                Cdn_atm_ratio_n,
                                l_Qa_iso,
                                Qrefn_iso,
                                uvel, vvel,
                                Urefn, zlvs)

            # if (icepack_warnings_aborted(subname)) return

        end   # calc_Tsfc or calc_strair

        if !(calc_strair)
            #ifndef CICE_IN_NEMO
            # Set to data values (on T points)
            strairxn = strax
            strairyn = stray
            #=
    #else
          ! NEMO wind stress is supplied on u grid, multipied 
          ! by ice concentration and set directly in evp, so
          ! strairxT/yT = 0. Zero u-components here for safety.
          strairxn = c0
          strairyn = c0
    #endif
    =#
        end

        #=
        Update ice age
        This is further adjusted for freezing in the thermodynamics.
        Melting does not alter the ice age.
        =#

        if (tr_iage)
            increment_age(dt, iage[n])
        end
        # if (icepack_warnings_aborted(subname)) return
        if (tr_FY)
            update_FYarea(dt, lmask_n, lmask_s, yday, FY[n])
        end
        # if (icepack_warnings_aborted(subname)) return

        # Vertical thermodynamics: Heat conduction, growth and melting.
        if !(calc_Tsfc)

            #=
          ! If not calculating surface temperature and fluxes, set 
          ! surface fluxes (flatn, fsurfn, and fcondtopn) to be used 
          ! in thickness_changes
            =#

            # hadgem routine sets fluxes to default values in ice-only mode
            set_sfcflux(aicen[n], flatn_f[n], fsensn_f[n], fcondtopn_f[n],                 
                        fsurfn_f[n], flatn[n], fsensn[n], fsurfn[n], fcondtopn[n])
            # if (icepack_warnings_aborted(subname)) return
        end

        thermo_vertical(nilyr, nslyr, 
                            dt,               aicen[n],     
                            vicen[n], vsnon[n],     
                            Tsfc[n], zSin[:,n],      
                            zqin[:,n], zqsn[:,n],      
                            apnd[n], hpnd[n],     
                            flw,             potT,              
                            Qa,               rhoa,              
                            fsnow,         fpond,            
                            fbot,           Tbot,              
                            Tsnice,       sss,                
                            l_rsnw[:,n], 
                            lhcoef,       shcoef,          
                            fswsfcn[n], fswintn[n],    
                            Sswabsn[:,n], Iswabsn[:,n],    
                            fsurfn[n], fcondtopn[n], 
                            fcondbotn[n],                      
                            fsensn[n], flatn[n],     
                            flwoutn,     evapn,            
                            evapsn,       evapin,          
                            freshn,       fsaltn,          
                            fhocnn,       frain,            
                            melttn[n], meltsn[n],     
                            meltbn[n], l_meltsliqn[n],
                            l_smice[:,n], massicen[:,n],  
                            l_smliq[:,n], massliqn[:,n],  
                            congeln[n], snoicen[n],    
                            mlt_onset, frz_onset,    
                            yday,           dsnown[n],     
                            prescribed_ice=prescribed_ice)
        #=
        if (icepack_warnings_aborted(subname)) then
            write(warnstr,*) subname, ' ice: Vertical thermo error, cat ', n
            call icepack_warnings_add(warnstr)
            return
        end
        =#

        # Total absorbed shortwave radiation
        fswabsn = fswsfcn[n] + fswintn[n] + fswthrun[n]
    #=
 !-----------------------------------------------------------------
 ! Aerosol update
 !-----------------------------------------------------------------

       if (tr_aero) then
          call update_aerosol (dt,                             &
                               nilyr, nslyr, n_aero,           &
                               melttn     (n), meltsn     (n), &
                               meltbn     (n), congeln    (n), &
                               snoicen    (n), fsnow,          &
                               aerosno(:,:,n), aeroice(:,:,n), &
                               aicen_init (n), vicen_init (n), &
                               vsnon_init (n),                 &
                               vicen      (n), vsnon      (n), &
                               aicen      (n),                 &
                               faero_atm     ,  faero_ocn)
          if (icepack_warnings_aborted(subname)) return
       endif

       if (tr_iso) then
          call update_isotope (dt = dt, &
                               nilyr = nilyr, nslyr = nslyr, &
                               meltt = melttn(n),melts = meltsn(n),     &
                               meltb = meltbn(n),congel=congeln(n),    &
                               snoice=snoicen(n),evap=evapn,         & 
                               fsnow=fsnow,      Tsfc=Tsfc(n),       &
                               Qref_iso=Qrefn_iso(:),                 &
                               isosno=l_isosno(:,n),isoice=l_isoice(:,n), &
                               aice_old=aicen_init(n),vice_old=vicen_init(n), &
                               vsno_old=vsnon_init(n),                &
                               vicen=vicen(n),vsnon=vsnon(n),      &
                               aicen=aicen(n),                     &
                               fiso_atm=l_fiso_atm(:),                  &
                               fiso_evapn=fiso_evapn(:),                &
                               fiso_ocnn=fiso_ocnn(:),                 &
                               HDO_ocn=l_HDO_ocn,H2_16O_ocn=l_H2_16O_ocn,    &
                               H2_18O_ocn=l_H2_18O_ocn)
          if (icepack_warnings_aborted(subname)) return
       endif
            =#
    end   # aicen_init

    if (snwgrain && use_smliq_pnd)
        drain_snow(nslyr, vsnon[n], aicen[n], massicen[:,n], massliqn[:,n], l_meltsliqn[n])
        # if (icepack_warnings_aborted(subname)) return
    end
    #=
 !-----------------------------------------------------------------
 ! Melt ponds
 ! If using tr_pond_cesm, the full calculation is performed here.
 ! If using tr_pond_topo, the rest of the calculation is done after
 ! the surface fluxes are merged, below.
 !-----------------------------------------------------------------

    !call ice_timer_start(timer_ponds)
    if (tr_pond) then
          
       if (tr_pond_cesm) then
          rfrac = rfracmin + (rfracmax-rfracmin) * aicen(n) 
          call compute_ponds_cesm(dt=dt,           &
                                  hi_min=hi_min,   &
                                  rfrac=rfrac,     &
                                  meltt=melttn(n), &
                                  melts=meltsn(n), &
                                  frain=frain,     &
                                  aicen=aicen (n), &
                                  vicen=vicen (n), &
                                  Tsfcn=Tsfc  (n), &
                                  apnd=apnd   (n), &
                                  hpnd=hpnd   (n), &
                                  meltsliqn=l_meltsliqn(n))
          if (icepack_warnings_aborted(subname)) return
             
       elseif (tr_pond_lvl) then
          rfrac = rfracmin + (rfracmax-rfracmin) * aicen(n)
          call compute_ponds_lvl (dt=dt,            &
                                  nilyr=nilyr,      &
                                  ktherm=ktherm,    &
                                  hi_min=hi_min,    &
                                  dpscale=dpscale,  &
                                  frzpnd=frzpnd,    &
                                  rfrac=rfrac,      &
                                  meltt=melttn (n), &
                                  melts=meltsn (n), &
                                  frain=frain,      &
                                  Tair=Tair,        &
                                  fsurfn=fsurfn(n), &
                                  dhs=dhsn     (n), &
                                  ffrac=ffracn (n), &
                                  aicen=aicen  (n), &
                                  vicen=vicen  (n), &
                                  vsnon=vsnon  (n), &
                                  qicen=zqin (:,n), &
                                  sicen=zSin (:,n), &
                                  Tsfcn=Tsfc   (n), &
                                  alvl=alvl    (n), &
                                  apnd=apnd    (n), &
                                  hpnd=hpnd    (n), &
                                  ipnd=ipnd    (n), &
                                  meltsliqn=l_meltsliqn(n))
          if (icepack_warnings_aborted(subname)) return
             
       elseif (tr_pond_topo) then
          if (aicen_init(n) > puny) then
                
             ! collect liquid water in ponds
             ! assume salt still runs off
             rfrac = rfracmin + (rfracmax-rfracmin) * aicen(n)
             if (snwgrain .and. use_smliq_pnd) then
                pond = rfrac/rhofresh * (melttn(n)*rhoi &
                     +                 l_meltsliqn(n))
             else
                pond = rfrac/rhofresh * (melttn(n)*rhoi &
                     +                   meltsn(n)*rhos &
                     +                   frain *dt)
             endif

             ! if pond does not exist, create new pond over full ice area
             ! otherwise increase pond depth without changing pond area
             if (apnd(n) < puny) then
                hpnd(n) = c0
                apnd(n) = c1
             endif
             hpnd(n) = (pond + hpnd(n)*apnd(n)) / apnd(n)
             fpond = fpond + pond * aicen(n) ! m
          endif ! aicen_init
       endif

    endif ! tr_pond
    !call ice_timer_stop(timer_ponds)
    =#
 !-----------------------------------------------------------------
 ! Increment area-weighted fluxes.
 !-----------------------------------------------------------------

    if (aicen_init(n) > puny) &
       call merge_fluxes (aicen=aicen_init(n),            &
                          flw=flw, & 
                          strairxn=strairxn, strairyn=strairyn,&
                          Cdn_atm_ratio_n=Cdn_atm_ratio_n,     &
                          fsurfn=fsurfn(n), fcondtopn=fcondtopn(n),&
                          fcondbotn=fcondbotn(n),              &
                          fsensn=fsensn(n),  flatn=flatn(n),   &
                          fswabsn=fswabsn,   flwoutn=flwoutn,  &
                          evapn=evapn,                         &
                          evapsn=evapsn,     evapin=evapin,    &
                          Trefn=Trefn,       Qrefn=Qrefn,      &
                          freshn=freshn,     fsaltn=fsaltn,    &
                          fhocnn=fhocnn,                       &
                          fswthrun=fswthrun(n),                &
                          fswthrun_vdr=l_fswthrun_vdr(n),      &
                          fswthrun_vdf=l_fswthrun_vdf(n),      &
                          fswthrun_idr=l_fswthrun_idr(n),      &
                          fswthrun_idf=l_fswthrun_idf(n),      &
                          strairxT=strairxT, strairyT=strairyT,&
                          Cdn_atm_ratio=Cdn_atm_ratio,         &
                          fsurf=fsurf,       fcondtop=fcondtop,&
                          fcondbot=fcondbot,                   &
                          fsens=fsens,       flat=flat,        &
                          fswabs=fswabs,     flwout=flwout,    &
                          evap=evap,                           &
                          evaps=evaps,       evapi=evapi,      &
                          Tref=Tref,         Qref=Qref,        &
                          fresh=fresh,       fsalt=fsalt,      &
                          fhocn=fhocn,                         &
                          fswthru=fswthru,                     &
                          fswthru_vdr=l_fswthru_vdr,           &
                          fswthru_vdf=l_fswthru_vdf,           &
                          fswthru_idr=l_fswthru_idr,           &
                          fswthru_idf=l_fswthru_idf,           &
                          melttn=melttn (n), meltsn=meltsn(n), &
                          meltbn=meltbn (n), congeln=congeln(n),&
                          meltt=meltt,       melts=melts,      &
                          meltb=meltb,       snoicen=snoicen(n),&
                          dsnow=dsnow,       dsnown=dsnown(n), &
                          congel=congel,     snoice=snoice,    &
                          meltsliq=l_meltsliq,      &
                          meltsliqn=l_meltsliqn(n), &
                          Uref=Uref,  Urefn=Urefn,  &
                          Qref_iso=l_Qref_iso,      &
                          Qrefn_iso=Qrefn_iso,      &
                          fiso_ocn=l_fiso_ocn,      &
                          fiso_ocnn=fiso_ocnn,      &
                          fiso_evap=l_fiso_evap,    &
                          fiso_evapn=fiso_evapn)

    # if (icepack_warnings_aborted(subname)) return

    end # ncat

    #=
    !-----------------------------------------------------------------
      ! reload snow mass tracers
      !-----------------------------------------------------------------

      if (snwgrain) then
         do n = 1, ncat
            if (vsnon(n) > puny) then
               do k = 1, nslyr
                  l_smice(k,n) = massicen(k,n) / (vsnon(n) * rnslyr)
                  l_smliq(k,n) = massliqn(k,n) / (vsnon(n) * rnslyr)
                  worka = l_smice(k,n) + l_smliq(k,n)
                  if (worka > puny) then
                     l_smice(k,n) = rhos * l_smice(k,n) / worka
                     l_smliq(k,n) = rhos * l_smliq(k,n) / worka
                  endif
               enddo
            else ! reset to default values
               do k = 1, nslyr
                  l_smice(k,n) = rhos
                  l_smliq(k,n) = c0
               enddo
            endif
         enddo
      endif

      if (present(isosno      )) isosno       = l_isosno
      if (present(isoice      )) isoice       = l_isoice
      if (present(Qa_iso      )) Qa_iso       = l_Qa_iso
      if (present(Qref_iso    )) Qref_iso     = l_Qref_iso
      if (present(fiso_atm    )) fiso_atm     = l_fiso_atm
      if (present(fiso_ocn    )) fiso_ocn     = l_fiso_ocn
      if (present(fiso_evap   )) fiso_evap    = l_fiso_evap
      if (present(fswthrun_vdr)) fswthrun_vdr = l_fswthrun_vdr
      if (present(fswthrun_vdf)) fswthrun_vdf = l_fswthrun_vdf
      if (present(fswthrun_idr)) fswthrun_idr = l_fswthrun_idr
      if (present(fswthrun_idf)) fswthrun_idf = l_fswthrun_idf
      if (present(fswthru_vdr )) fswthru_vdr  = l_fswthru_vdr
      if (present(fswthru_vdf )) fswthru_vdf  = l_fswthru_vdf
      if (present(fswthru_idr )) fswthru_idr  = l_fswthru_idr
      if (present(fswthru_idf )) fswthru_idf  = l_fswthru_idf
      if (present(fsloss      )) fsloss       = l_fsloss
      if (present(meltsliqn   )) meltsliqn    = l_meltsliqn
      if (present(meltsliq    )) meltsliq     = l_meltsliq
      if (present(rsnwn       )) rsnwn        = l_rsnw
      if (present(smicen      )) smicen       = l_smice
      if (present(smliqn      )) smliqn       = l_smliq
      deallocate(l_isosno)
      deallocate(l_isoice)
      deallocate(l_Qa_iso)
      deallocate(l_Qref_iso)
      deallocate(l_fiso_atm)
      deallocate(l_fiso_ocn)
      deallocate(l_fiso_evap)
      deallocate(l_fswthrun_vdr)
      deallocate(l_fswthrun_vdf)
      deallocate(l_fswthrun_idr)
      deallocate(l_fswthrun_idf)
      deallocate(l_meltsliqn)
      deallocate(l_rsnw)
      deallocate(l_smice)
      deallocate(l_smliq)

      !-----------------------------------------------------------------
      ! Calculate ponds from the topographic scheme
      !-----------------------------------------------------------------
      !call ice_timer_start(timer_ponds)
      if (tr_pond_topo) then
         call compute_ponds_topo(dt,       ncat,      nilyr,     &
                                 ktherm,   heat_capacity,        &
                                 aice,     aicen,                &
                                 vice,     vicen,                &
                                 vsno,     vsnon,                &
                                 meltt,                &
                                 fsurf,    fpond,                &
                                 Tsfc,     Tf,                   &
                                 zqin,     zSin,                 &
                                 apnd,     hpnd,      ipnd       )
         if (icepack_warnings_aborted(subname)) return
      endif
      !call ice_timer_stop(timer_ponds)
    =#

end