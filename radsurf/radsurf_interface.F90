! radsurf_interface.F90 - Perform surface radiative transfer calculation
!
! Copyright (C) 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radsurf_interface

contains

  subroutine radsurf(config, canopy_props, facet_props, volume_props, &
       &             bc_out, &
       &             istartcol, iendcol, &
       &             sw_norm_dir, sw_norm_diff, &
       &             lw_internal, lw_norm)

    use parkind1,                   only : jpim, jprb
    use yomhook,                    only : lhook, dr_hook
    use radiation_io,               only : nulout
    use radsurf_config,             only : config_type
    use radsurf_canopy_properties,  only : ITileFlat,  ITileForest, &
         &                                 ITileUrban, ITileVegetatedUrban, &
         &                                 canopy_properties_type
    use radsurf_facet_properties,   only : facet_properties_type
    use radsurf_volume_properties,  only : volume_properties_type
    use radsurf_boundary_conds_out, only : boundary_conds_out_type
    use radsurf_canopy_flux,        only : canopy_flux_type
    use radsurf_forest_sw,          only : spartacus_forest_sw
    use radsurf_urban_sw,           only : spartacus_urban_sw
    use radsurf_urban_lw,           only : spartacus_urban_lw
    
    implicit none

    type(config_type),             intent(in)  :: config
    type(canopy_properties_type),  intent(in)  :: canopy_props
    type(facet_properties_type),   intent(in), target  :: facet_props
    type(volume_properties_type),  intent(in)  :: volume_props
    type(boundary_conds_out_type), intent(inout) :: bc_out
    type(canopy_flux_type),        intent(inout), optional &
         &  :: sw_norm_dir, &  ! SW fluxes normalized by top-of-canopy direct
         &     sw_norm_diff, & ! SW fluxes normalized by top-of-canopy diffuse
         &     lw_internal, &  ! LW fluxes from internal emission
         &     lw_norm         ! LW fluxes normalized by top-of-canopy down

    integer(kind=jpim), optional, intent(in) :: istartcol, iendcol

    ! Actual start and end columns
    integer(kind=jpim) :: icol1, icol2

    ! Start and end layers for current column
    integer(kind=jpim) :: ilay1, ilay2

    ! Loop index for column
    integer(kind=jpim) :: jcol

    ! Representation for current column (flat, vegetated, urban etc)
    integer(kind=jpim) :: irep

    ! Is the current column layered?
    logical :: do_canopy

    ! Direct ground shortwave albedo, in case not provided by user
    real(kind=jprb), pointer :: ground_sw_albedo_direct(:,:) ! (nswinterval,ncol)

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_interface:radsurf',0,hook_handle)

    ! Find which range of columns to process
    if (present(istartcol)) then
      icol1 = istartcol
    else
      icol1 = 1
    end if

    if (present(iendcol)) then
      icol2 = iendcol
    else
      icol2 = canopy_props%ncol
    end if

    if (config%use_sw_direct_albedo) then
      ground_sw_albedo_direct => facet_props%ground_sw_albedo_direct
    else
      ground_sw_albedo_direct => facet_props%ground_sw_albedo
    end if

    ! Loop through columns calculating radiative transfer on each
    do jcol = icol1,icol2

      irep = canopy_props%i_representation(jcol)
      
      do_canopy = (.not. irep == ITileFlat)

      if (do_canopy) then
        ilay1 = canopy_props%istartlay(jcol)
        ilay2 = ilay1 + canopy_props%nlay(jcol)-1
      else
        ! If this is ever used then should be detected with Fortran
        ! bounds checking on
        ilay1 = 0
        ilay2 = 0
      end if

      select case (irep)
      case (ITileFlat)
        ! Flat tiles simply involve copying the ground facet
        ! properties to the outward boundary conditions, and copying
        ! the incoming fluxes down to the ground facet
        if (config%iverbose >= 4) then
          write(nulout,'(i5,a)') jcol, ': Flat'
        end if

        if (config%do_sw) then
          ! Copy shortwave albedo
          bc_out%sw_albedo(:,jcol)        = facet_props%ground_sw_albedo(:,jcol)
          bc_out%sw_albedo_direct(:,jcol) = ground_sw_albedo_direct(:,jcol)
          ! Rate of change of ground and top-of-canopy fluxes with
          ! respect to direct flux at top-of-canopy
          sw_norm_dir%ground_dn_direct(:,jcol) = 1.0_jprb
          sw_norm_dir%ground_dn(:,jcol)  = 1.0_jprb
          sw_norm_dir%ground_net(:,jcol) = 1.0_jprb - ground_sw_albedo_direct(:,jcol)
          sw_norm_dir%top_dn_direct(:,jcol) = 1.0_jprb
          sw_norm_dir%top_dn(:,jcol)  = 1.0_jprb
          sw_norm_dir%top_net(:,jcol) = 1.0_jprb - ground_sw_albedo_direct(:,jcol)
          ! Rate of change of ground and top-of-canopy fluxes with
          ! respect to diffuse downward flux at top-of-canopy
          sw_norm_diff%ground_dn_direct(:,jcol) = 0.0_jprb
          sw_norm_diff%ground_dn(:,jcol)  = 1.0_jprb
          sw_norm_diff%ground_net(:,jcol) = 1.0_jprb - facet_props%ground_sw_albedo(:,jcol)
          sw_norm_diff%top_dn_direct(:,jcol) = 0.0_jprb
          sw_norm_diff%top_dn(:,jcol)  = 1.0_jprb
          sw_norm_diff%top_net(:,jcol) = 1.0_jprb - facet_props%ground_sw_albedo(:,jcol)
        end if

        if (config%do_lw) then
          ! Copy longwave albedo and upward emission
          bc_out%lw_emissivity(:,jcol)   = facet_props%ground_lw_emissivity(:,jcol)
          bc_out%lw_emission(:,jcol)     = facet_props%ground_lw_emission(:,jcol)
          ! Longwave fluxes due to surface emission
          lw_internal%ground_dn(:,jcol)  = 0.0_jprb
          lw_internal%ground_net(:,jcol) = -facet_props%ground_lw_emission(:,jcol)
          lw_internal%top_dn(:,jcol)  = 0.0_jprb
          lw_internal%top_net(:,jcol) = -facet_props%ground_lw_emission(:,jcol)
          ! Rate of change of ground and top-of-canopy fluxes with
          ! respect to diffuse downward flux at top-of-canopy
          lw_norm%ground_dn(:,jcol) = 1.0_jprb
          lw_norm%ground_net(:,jcol) = facet_props%ground_lw_emissivity(:,jcol)
          lw_norm%top_dn(:,jcol) = 1.0_jprb
          lw_norm%top_net(:,jcol) = facet_props%ground_lw_emissivity(:,jcol)
        end if

      case (ITileForest)
        if (config%iverbose >= 4) then
          write(nulout,'(i5,a,i0,a,i0,a,i0,a)') jcol, ': Forest,            ', &
               &  canopy_props%nlay(jcol), ' layers, ', &
               &  config%lg_forest%nstream, ' diffuse streams per hemisphere, ', &
               &  config%n_vegetation_region_forest+1, ' regions'
        end if
        if (config%do_sw) then
          call spartacus_forest_sw(config, config%nswinternal, &
               &  config%lg_forest%nstream, config%n_vegetation_region_forest+1, &
               &  canopy_props%nlay(jcol), jcol, ilay1, ilay2, &
               &  config%lg_forest, canopy_props%cos_sza(jcol), &
               &  canopy_props, volume_props, &
               &  facet_props%ground_sw_albedo(:,jcol), &
               &              ground_sw_albedo_direct(:,jcol), &
               &  bc_out%sw_albedo(:,jcol), bc_out%sw_albedo_direct(:,jcol), &
               &  sw_norm_dir, sw_norm_diff)
        end if

      case (ITileUrban)
        if (config%iverbose >= 4) then
          write(nulout,'(i5,a,i0,a,i0,a)') jcol, ': Unvegetated urban, ', &
               &  canopy_props%nlay(jcol), ' layers, ', &
               &  config%lg_urban%nstream, ' diffuse streams per hemisphere, 1 region'
        end if

        if (config%do_sw) then
          call sw_norm_dir%zero(jcol,ilay1,ilay2)
          call sw_norm_diff%zero(jcol,ilay1,ilay2)
          call spartacus_urban_sw(config, config%nswinternal, &
               &  config%lg_urban%nstream, 1, &
               &  canopy_props%nlay(jcol), jcol, ilay1, ilay2, &
               &  config%lg_urban, canopy_props%cos_sza(jcol), &
               &  canopy_props, volume_props, facet_props, &
               &  facet_props%ground_sw_albedo(:,jcol), &
               &              ground_sw_albedo_direct(:,jcol), &
               &  bc_out%sw_albedo(:,jcol), bc_out%sw_albedo_direct(:,jcol), &
               &  sw_norm_dir, sw_norm_diff)
        end if
        if (config%do_lw) then
          call lw_internal%zero(jcol,ilay1,ilay2)
          call lw_norm%zero(jcol,ilay1,ilay2)
          call spartacus_urban_lw(config, config%nlwinternal, &
               &  config%lg_urban%nstream, 1, &
               &  canopy_props%nlay(jcol), jcol, ilay1, ilay2, &
               &  config%lg_urban, canopy_props, volume_props, facet_props, &
               &  bc_out%lw_emissivity(:,jcol), bc_out%lw_emission(:,jcol), &
               &  lw_internal, lw_norm)
        end if
        
      case (ITileVegetatedUrban)
        if (config%iverbose >= 4) then
          write(nulout,'(i5,a,i0,a,i0,a,i0,a)') jcol, ': Vegetated urban,   ', &
               &  canopy_props%nlay(jcol), ' layers, ', &
               &  config%lg_urban%nstream, ' diffuse streams per hemisphere, ', &
               &  config%n_vegetation_region_urban+1, ' regions'
        end if
        if (config%do_sw) then
          call spartacus_urban_sw(config, config%nswinternal, &
               &  config%lg_urban%nstream, config%n_vegetation_region_urban+1, &
               &  canopy_props%nlay(jcol), jcol, ilay1, ilay2, &
               &  config%lg_urban, canopy_props%cos_sza(jcol), &
               &  canopy_props, volume_props, facet_props, &
               &  facet_props%ground_sw_albedo(:,jcol), &
               &              ground_sw_albedo_direct(:,jcol), &
               &  bc_out%sw_albedo(:,jcol), bc_out%sw_albedo_direct(:,jcol), &
               &  sw_norm_dir, sw_norm_diff)
        end if
        if (config%do_lw) then
          call spartacus_urban_lw(config, config%nlwinternal, &
               &  config%lg_urban%nstream, config%n_vegetation_region_urban+1, &
               &  canopy_props%nlay(jcol), jcol, ilay1, ilay2, &
               &  config%lg_urban, canopy_props, volume_props, facet_props, &
               &  bc_out%lw_emissivity(:,jcol), bc_out%lw_emission(:,jcol), &
               &  lw_internal, lw_norm)
        end if

      end select

    end do

    if (lhook) call dr_hook('radiation_interface:radsurf',0,hook_handle)

  end subroutine radsurf

end module radsurf_interface

