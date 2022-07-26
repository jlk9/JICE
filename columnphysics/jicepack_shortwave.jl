# Written by Joseph Kump (josek97@utexas.edu)

# The albedo and absorbed/transmitted flux parameterizations for
# snow over ice, bare ice and ponded ice.

#= Scales radiation fields computed on the previous time step. Variables:

TODO:   determine variables that need to be saved as external constants
        figure out what to do with optional variables

      integer (kind=int_kind), intent(in) :: &
         ncat    , & ! number of ice thickness categories
         nilyr   , & ! number of ice layers
         nslyr       ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         aice        , & ! ice area fraction
         swvdr       , & ! sw down, visible, direct  (W/m^2)
         swvdf       , & ! sw down, visible, diffuse (W/m^2)
         swidr       , & ! sw down, near IR, direct  (W/m^2)
         swidf       , & ! sw down, near IR, diffuse (W/m^2)
         ! grid-box-mean albedos aggregated over categories (if calc_Tsfc)
         alvdr_ai    , & ! visible, direct   (fraction)
         alidr_ai    , & ! near-ir, direct   (fraction)
         alvdf_ai    , & ! visible, diffuse  (fraction)
         alidf_ai        ! near-ir, diffuse  (fraction)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         aicen           ! ice area fraction in each category

      real (kind=dbl_kind), intent(inout) :: &
         scale_factor    ! shortwave scaling factor, ratio new:old

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         fswsfcn     , & ! SW absorbed at ice/snow surface (W m-2)
         fswintn     , & ! SW absorbed in ice interior, below surface (W m-2)
         fswthrun        ! SW through ice to ocean (W/m^2)

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         fswthrun_vdr , & ! vis dir SW through ice to ocean (W/m^2)
         fswthrun_vdf , & ! vis dif SW through ice to ocean (W/m^2)
         fswthrun_idr , & ! nir dir SW through ice to ocean (W/m^2)
         fswthrun_idf     ! nir dif SW through ice to ocean (W/m^2)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         fswpenln    , & ! visible SW entering ice layers (W m-2)
         Iswabsn     , & ! SW radiation absorbed in ice layers (W m-2)
         Sswabsn         ! SW radiation absorbed in snow layers (W m-2)
=#
function jicepack_prep_radiation(ncat, nilyr, nslyr, aice, aicen, swvdr, swvdf, swidr, swidf, alvdr_ai, alvdf_ai,
                                 alidr_ai, alidf_ai, scale_factor, fswsfcn, fswthrun, fswthrun_vdr, fswthrun_vdf,
                                 fswthrun_idr, fswthrun_idf, fswpenln, Sswabsn, Iswabsn)

    # Local variables:
    k     = 0 # vertical index
    n     = 0 # thickness category index
    netsw = 0.0

    #character(len=*),parameter :: subname='(icepack_prep_radiation)'

    # Compute netsw scaling factor (new netsw / old netsw):
    if (aice > c0 && scale_factor > puny)
        netsw        = swvdr*(c1 - alvdr_ai) + swvdf*(c1 - alvdf_ai)
                     + swidr*(c1 - alidr_ai) + swidf*(c1 - alidf_ai)
        scale_factor = netsw / scale_factor
    else
        scale_factor = c1
    end

    for n = 1:ncat
        if (aicen(n) > puny)

            # Scale absorbed solar radiation for change in net shortwave

            fswsfcn[n]  = scale_factor*fswsfcn[n]
            fswintn[n]  = scale_factor*fswintn[n]
            fswthrun[n] = scale_factor*fswthrun[n]

            #=
            if (present(fswthrun_vdr)) fswthrun_vdr[n] = scale_factor*fswthrun_vdr[n]
            if (present(fswthrun_vdf)) fswthrun_vdf[n] = scale_factor*fswthrun_vdf[n]
            if (present(fswthrun_idr)) fswthrun_idr[n] = scale_factor*fswthrun_idr[n]
            if (present(fswthrun_idf)) fswthrun_idf[n] = scale_factor*fswthrun_idf[n]
            =#

            for k = 1:(nilyr+1)
                fswpenln[k,n] = scale_factor*fswpenln[k,n]
            end
            for k=1:nslyr
                Sswabsn[k,n] = scale_factor*Sswabsn[k,n]
            end
            for k=1:nilyr
                Iswabsn[k,n] = scale_factor*Iswabsn[k,n]
            end

        end
    end

end