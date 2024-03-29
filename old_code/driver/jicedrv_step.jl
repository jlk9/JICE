# Written by Joseph Kump (josek97@utexas.edu)

# Does prep and scale radiation
# TODO: warning check, optional variables from jicepack_prep_radiation
function prep_radiation()

    #=
    use icedrv_domain_size, only: ncat, nilyr, nslyr, nx
    use icedrv_flux, only: scale_factor, swvdr, swvdf, swidr, swidf
    use icedrv_flux, only: alvdr_ai, alvdf_ai, alidr_ai, alidf_ai
    use icedrv_flux, only: alvdr_init, alvdf_init, alidr_init, alidf_init
    use icedrv_arrays_column, only: fswsfcn, fswintn
    use icedrv_arrays_column, only: fswthrun, fswthrun_vdr, fswthrun_vdf, fswthrun_idr, fswthrun_idf
    use icedrv_arrays_column, only: fswpenln, Sswabsn, Iswabsn
    use icedrv_state, only: aice, aicen
    =#

    # column package includes
    #use icepack_intfc, only: icepack_prep_radiation

    # Compute netsw scaling factor (new netsw / old netsw)
    for i = 1:nx

            alvdr_init[i] = alvdr_ai[i]
            alvdf_init[i] = alvdf_ai[i]
            alidr_init[i] = alidr_ai[i]
            alidf_init[i] = alidf_ai[i]

            icepack_prep_radiation(ncat, nilyr, nslyr, aice[i], aicen[i,:], swvdr[i], swvdf[i],
                         swidr[i], swidf[i], alvdr_ai[i], alvdf_ai[i], alidr_ai[i], alidf_ai[i],
                         scale_factor[i], fswsfcn[i,:], fswintn[i,:], fswthrun[i,:],
                         fswthrun_vdr[i,:], fswthrun_vdf[i,:], fswthrun_idr[i,:],
                         fswthrun_idf[i,:], fswpenln[i,:,:], Sswabsn[i,:,:], Iswabsn[i,:,:])

    end

end

# Driver for updating ice and snow internal temperatures and
# computing thermodynamic growth rates and coupler fluxes.
# dt is time step
function step_therm_1(dt)

    #=
    use icedrv_arrays_column, only: ffracn, dhsn
      use icedrv_arrays_column, only: Cdn_ocn, Cdn_ocn_skin, Cdn_ocn_floe
      use icedrv_arrays_column, only: Cdn_ocn_keel, Cdn_atm_ratio
      use icedrv_arrays_column, only: Cdn_atm, Cdn_atm_skin, Cdn_atm_floe
      use icedrv_arrays_column, only: Cdn_atm_rdg, Cdn_atm_pond
      use icedrv_arrays_column, only: hfreebd, hdraft, hridge, distrdg
      use icedrv_arrays_column, only: hkeel, dkeel, lfloe, dfloe
      use icedrv_arrays_column, only: fswsfcn, fswintn, Sswabsn, Iswabsn
      use icedrv_arrays_column, only: fswthrun, fswthrun_vdr, fswthrun_vdf, fswthrun_idr, fswthrun_idf
      use icedrv_arrays_column, only: meltsliqn, meltsliq
      use icedrv_calendar, only: yday
      use icedrv_domain_size, only: ncat, nilyr, nslyr, n_aero, n_iso, nx
      use icedrv_flux, only: frzmlt, sst, Tf, strocnxT, strocnyT, rside, fside, &
                             fbot, Tbot, Tsnice
      use icedrv_flux, only: meltsn, melttn, meltbn, congeln, snoicen, uatm, vatm
      use icedrv_flux, only: wind, rhoa, potT, Qa, Qa_iso, zlvl, strax, stray, flatn
      use icedrv_flux, only: fsensn, fsurfn, fcondtopn, fcondbotn
      use icedrv_flux, only: flw, fsnow, fpond, sss, mlt_onset, frz_onset, fsloss
      use icedrv_flux, only: frain, Tair, strairxT, strairyT, fsurf
      use icedrv_flux, only: fcondtop, fcondbot, fsens, fresh, fsalt, fhocn
      use icedrv_flux, only: flat, fswabs, flwout, evap, evaps, evapi
      use icedrv_flux, only: Tref, Qref, Qref_iso, Uref
      use icedrv_flux, only: meltt, melts, meltb, congel, snoice
      use icedrv_flux, only: fswthru, fswthru_vdr, fswthru_vdf, fswthru_idr, fswthru_idf
      use icedrv_flux, only: flatn_f, fsensn_f, fsurfn_f, fcondtopn_f
      use icedrv_flux, only: dsnow, dsnown, faero_atm, faero_ocn
      use icedrv_flux, only: fiso_atm, fiso_ocn, fiso_evap
      use icedrv_flux, only: HDO_ocn, H2_16O_ocn, H2_18O_ocn
      use icedrv_init, only: lmask_n, lmask_s
      use icedrv_state, only: aice, aicen, aice_init, aicen_init, vicen_init
      use icedrv_state, only: vice, vicen, vsno, vsnon, trcrn, uvel, vvel, vsnon_init

      ! column packge includes
      use icepack_intfc, only: icepack_step_therm1

      logical (kind=log_kind) :: & 
         prescribed_ice ! if .true., use prescribed ice instead of computed
    =#

    # Local Variables
    # Integers:
    ntrcr     = 0
    nt_apnd   = 0
    nt_hpnd   = 0
    nt_ipnd   = 0
    nt_alvl   = 0
    nt_vlvl   = 0
    nt_Tsfc   = 0
    nt_iage   = 0
    nt_FY     = 0
    nt_qice   = 0
    nt_sice   = 0
    nt_qsno   = 0
    nt_aero   = 0
    nt_isosno = 0
    nt_isoice = 0
    nt_rsnw   = 0
    nt_smice  = 0
    nt_smliq  = 0

    # logical:
    tr_iage   = false
    tr_FY     = false
    tr_aero   = false
    tr_iso    = false
    calc_Tsfc = false
    tr_snow   = false

    # float
    puny = 0.0

    # float arrays of dimension(n_aero,2,ncat), kg/m^2
    aerosno = zeros(Float64, n_aero, 2, ncat)
    aeroice = zeros(Float64, n_aero, 2, ncat)

    # float arrays of dimension(n_iso,ncat), kg/m^2
    isosno = zeros(Float64, n_iso, ncat)
    isoice = zeros(Float64, n_iso, ncat)

    # float arrays of dimension(nslyr,ncat)
    rsnwn  = zeros(Float64, nslyr, ncat)
    smicen = zeros(Float64, nslyr, ncat)
    smliqn = zeros(Float64, nslyr, ncat)

    # query icepack values

    #=
    call icepack_query_parameters(puny_out=puny)
      call icepack_query_parameters(calc_Tsfc_out=calc_Tsfc)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call icepack_query_tracer_sizes( &
           ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call icepack_query_tracer_flags( &
         tr_iage_out=tr_iage, tr_FY_out=tr_FY, &
         tr_aero_out=tr_aero, tr_iso_out=tr_iso, tr_snow_out=tr_snow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      call icepack_query_tracer_indices( &
         nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd, &
         nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, nt_Tsfc_out=nt_Tsfc, &
         nt_iage_out=nt_iage, nt_FY_out=nt_FY, &
         nt_qice_out=nt_qice, nt_sice_out=nt_sice, &
         nt_aero_out=nt_aero, nt_qsno_out=nt_qsno, &
         nt_rsnw_out=nt_rsnw, nt_smice_out=nt_smice, nt_smliq_out=nt_smliq, &
         nt_isosno_out=nt_isosno, nt_isoice_out=nt_isoice)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)
    =#

    prescribed_ice = false
    aerosno[:,:,:] = c0
    aeroice[:,:,:] = c0
    isosno[:,:]    = c0
    isoice[:,:]    = c0
    rsnwn [:,:]    = c0
    smicen[:,:]    = c0
    smliqn[:,:]    = c0

    for i = 1:nx

        #=
        Save the ice area passed to the coupler (so that history fields
        can be made consistent with coupler fields).
        Save the initial ice area and volume in each category.
        =#

        aice_init[i] = aice[i]

        for n = 1:ncat
            aicen_init[i,n] = aicen[i,n]
            vicen_init[i,n] = vicen[i,n]
            vsnon_init[i,n] = vsnon[i,n]
        end
    end

    for i = 1:nx
        # echmod    do i = 3,3
        if (tr_aero)
            # trcrn(nt_aero) has units kg/m^3
            for n = 1:ncat
                for k = 1:n_aero
                    aerosno[k,:,n] = trcrn[i,nt_aero+(k-1)*4:nt_aero+(k-1)*4+1,n]   * vsnon_init[i,n]
                    aeroice[k,:,n] = trcrn[i,nt_aero+(k-1)*4+2:nt_aero+(k-1)*4+3,n] * vicen_init[i,n]
                end
            end
        end # tr_aero
        
        if (tr_iso)
            # trcrn(nt_isosno/ice) has units kg/m^3
            for n = 1:ncat
                for k = 1:n_iso
                    isosno[k,n] = trcrn[i,nt_isosno+k-1,n] * vsnon_init[i,n]
                    isoice[k,n] = trcrn[i,nt_isoice+k-1,n] * vicen_init[i,n]
                end
            end
        end # tr_iso

        if (tr_snow)
            for n = 1:ncat
                for k = 1:nslyr
                    rsnwn[k,n]  = trcrn[i,nt_rsnw +k-1,n]
                    smicen[k,n] = trcrn[i,nt_smice+k-1,n]
                    smliqn[k,n] = trcrn[i,nt_smliq+k-1,n]
                end
            end
        end # tr_snow

        # Clean up :(
        jicepack_step_therm1(dt, ncat, nilyr, nslyr,
            aicen_init[i,:],
            vicen_init[i,:],
            vsnon_init[i,:],
            aice[i],   aicen[i,:],
            vice[i],   vicen[i,:],
            vsno[i],   vsnon[i,:],
            uvel[i],   vvel[i],   
            trcrn[i,nt_Tsfc,:],
            trcrn[i,nt_qsno:nt_qsno+nslyr-1,:],
            trcrn[i,nt_qice:nt_qice+nilyr-1,:],
            trcrn[i,nt_sice:nt_sice+nilyr-1,:],
            trcrn[i,nt_alvl,:],
            trcrn[i,nt_vlvl,:],
            trcrn[i,nt_apnd,:],
            trcrn[i,nt_hpnd,:],
            trcrn[i,nt_ipnd,:],
            trcrn[i,nt_iage,:],
            trcrn[i,nt_FY,:],
            rsnwn[:,:],
            smicen[:,:],
            smliqn[:,:],
            aerosno[:,:,:],
            aeroice[:,:,:],
            isosno[:,:],
            isoice[:,:],
            uatm[i], vatm[i],
            wind[i], zlvl[i],
            Qa[i],   rhoa[i],
            Qa_iso[i,:],
            Tair[i], Tref[i],
            Qref[i], Uref[i],
            Qref_iso[i,:],
            Cdn_atm_ratio[i],
            Cdn_ocn[i],
            Cdn_ocn_skin[i],
            Cdn_ocn_floe[i],
            Cdn_ocn_keel[i],
            Cdn_atm[i],
            Cdn_atm_skin[i],
            Cdn_atm_floe[i],
            Cdn_atm_pond[i],
            Cdn_atm_rdg[i], 
            hfreebd[i],    hkeel[i],
            hdraft[i],     hridge[i],
            distrdg[i],    dkeel[i],
            lfloe[i],      dfloe[i],
            strax[i],      stray[i],
            strairxT[i],   strairyT[i],
            potT[i],       sst[i],
            sss[i],        Tf[i],
            strocnxT[i],   strocnyT[i],
            fbot[i],       frzmlt[i],
            Tbot[i],       Tsnice[i],
            rside[i],      fside[i],
            fsnow[i],      frain[i],
            fpond[i],      fsloss[i],
            fsurf[i],      fsurfn[i,:],
            fcondtop[i],   fcondtopn[i,:],
            fcondbot[i],   fcondbotn[i,:],
            fswsfcn[i,:],  fswintn[i,:],
            fswthrun[i,:],
            fswthrun_vdr[i,:],
            fswthrun_vdf[i,:],
            fswthrun_idr[i,:],
            fswthrun_idf[i,:],
            fswabs[i],
            flwout[i],     flw[i],
            fsens[i],      fsensn[i,:],
            flat[i],       flatn[i,:],
            fresh[i],      fsalt[i],
            fhocn[i],
            fswthru[i],
            fswthru_vdr[i],
            fswthru_vdf[i],
            fswthru_idr[i],
            fswthru_idf[i],
            flatn_f[i,:],  fsensn_f[i,:],
            fsurfn_f[i,:],
            fcondtopn_f[i,:],
            faero_atm[i,1:n_aero],
            faero_ocn[i,1:n_aero],
            fiso_atm[i,:],
            fiso_ocn[i,:],
            fiso_evap[i,:],
            HDO_ocn[i],
            H2_16O_ocn[i],
            H2_18O_ocn[i],
            Sswabsn[i,:,:], Iswabsn[i,:,:],
            evap[i], evaps[i], evapi[i],
            dhsn[i,:],     ffracn[i,:],
            meltt[i],      melttn[i,:],
            meltb[i],      meltbn[i,:],
            melts[i],      meltsn[i,:],
            congel[i],     congeln[i,:],
            snoice[i],     snoicen[i,:],
            dsnow[i],      dsnown[i,:],
            meltsliqn[i,:],
            lmask_n[i],    lmask_s[i],
            mlt_onset[i],   frz_onset[i]
            yday, prescribed_ice)

        if (tr_aero)
            for n = 1:ncat
                if (vicen[i,n] > puny)
                    aeroice[:,:,n] = aeroice[:,:,n]/vicen[i,n]
                end
                if (vsnon[i,n] > puny)
                    aerosno[:,:,n] = aerosno[:,:,n]/vsnon[i,n]
                end

                for k = 1:n_aero
                    for kk = 1:2
                        trcrn[i,nt_aero+(k-1)*4+kk-1,n] = aerosno[k,kk,n]
                        trcrn[i,nt_aero+(k-1)*4+kk+1,n] = aeroice[k,kk,n]
                    end
                end
            end
        end # tr_aero
        
        if (tr_iso)
            for n = 1:ncat
                if (vicen(i,n) > puny)
                    isoice[:,n] = isoice[:,n]/vicen[i,n]
                end
                if (vsnon(i,n) > puny)
                    isosno[:,n] = isosno[:,n]/vsnon[i,n]
                end
                for k = 1:n_iso
                    trcrn[i,nt_isosno+k-1,n] = isosno[k,n]
                    trcrn[i,nt_isoice+k-1,n] = isoice[k,n]
                end
            end
        end # tr_iso
        
        if (tr_snow)
            for n = 1:ncat
                for k = 1:nslyr
                    trcrn[i,nt_rsnw +k-1,n] = rsnwn[k,n]
                    trcrn[i,nt_smice+k-1,n] = smicen[k,n]
                    trcrn[i,nt_smliq+k-1,n] = smliqn[k,n]
                end
            end
        end # tr_snow

    end # i

    #= Warning check
    call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__, line=__LINE__)
    =#

end

# Driver for thermodynamic changes not needed for coupling: transport in thickness space,
# lateral growth and melting.
function step_therm_2(dt)



end

# Computes radiation fields
function step_radiation()



end