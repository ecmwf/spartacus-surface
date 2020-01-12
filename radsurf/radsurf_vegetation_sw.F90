! radsurf_vegetation_sw.F90 - SPARTACUS shortwave vegetation solver
!
! Copyright (C) 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radsurf_vegetation_sw

contains

  subroutine spartacus_vegetation_sw(config, &
       &  nsw, ns, nreg, nlay, ilay1, ilay2, &
       &  lg, cos_sza, &
       &  canopy_props, volume_props, &
       &  ground_sw_albedo, ground_sw_albedo_direct, &
       &  bc_out, &
       &  sw_norm_dir, sw_norm_diff, lw_internal, lw_norm)
    
    use parkind1,                   only : jpim, jprb
    use yomhook,                    only : lhook, dr_hook
    use radsurf_config,             only : config_type
    use radtool_legendre_gauss,     only : legendre_gauss_type
    use radsurf_canopy_properties,  only : canopy_properties_type
    use radsurf_volume_properties,  only : volume_properties_type
    use radsurf_boundary_conds_out, only : boundary_conds_out_type
    use radsurf_canopy_flux,        only : canopy_flux_type

    implicit none

    type(config_type),             intent(in)  :: config
    ! Number of spectral intervals, number of layers
    integer(kind=jpim),            intent(in)  :: nsw, nlay
    ! Index of first and last layer of the current column
    integer(kind=jpim),            intent(in)  :: ilay1, ilay2
    ! Number of regions, number of diffuse streams in each hemisphere
    integer(kind=jpim),            intent(in)  :: nreg, ns
    ! Legendre-Gauss coefficients
    type(legendre_gauss_type),     intent(in)  :: lg
    ! Cosine of the solar zenith angle
    real(kind=jprb),               intent(in)  :: cos_sza
    ! Properties of the canopy
    type(canopy_properties_type),  intent(in)  :: canopy_props
    type(volume_properties_type),  intent(in)  :: volume_props
    real(kind=jprb), dimension(nsw),intent(in) :: ground_sw_albedo, ground_sw_albedo_direct
    type(boundary_conds_out_type), intent(out) :: bc_out
    type(canopy_flux_type),        intent(out), optional &
         &  :: sw_norm_dir, &  ! SW fluxes normalized by top-of-canopy direct
         &     sw_norm_diff, & ! SW fluxes normalized by top-of-canopy diffuse
         &     lw_internal, &  ! LW fluxes from internal emission
         &     lw_norm         ! LW fluxes normalized by top-of-canopy down

    ! Transmittance and reflectance of a layer to diffuse radiation
    real(kind=jprb), dimension(nsw,nreg*ns,nreg*ns,nlay) &
         &  :: trans_dif, ref_dif
    ! Reflectance of a layer to direct radiation, and fraction of
    ! direct radiation that is scattered and passes out through the
    ! base of the layer
    real(kind=jprb), dimension(nsw,nreg*ns,nreg,nlay) &
         &   :: ref_dir, trans_dir_dif
    ! Direct (unscattered) transmittance
    real(kind=jprb), dimension(nsw,nreg,nreg,nlay) :: trans_dir_dir

    ! Area fraction of each region in each layer
    real(kind=jprb), dimension(nreg,nlay) :: frac
    
    ! Components of the Gamma matrices
    real(kind=jprb), dimension(nsw,nreg,nreg,nlay)       :: gamma0
    real(kind=jprb), dimension(nsw,nreg*ns,nreg*ns,nlay) :: gamma1, gamma2
    real(kind=jprb), dimension(nsw,nreg*ns,nreg,nlay)    :: gamma3

    integer(kind=jpim) :: jlay, jreg, jsw
    
    associate( &
         &  dz           => canopy_props%dz(ilay1:ilay2), &
         &  veg_fraction => canopy_props%veg_fraction(ilay1:ilay2), &
         &  veg_scale    => canopy_props%veg_scale(ilay1:ilay2), &
         &  veg_fsd      => canopy_props%veg_fsd(ilay1:ilay2), &
         &  air_sw_ext   => volume_props%air_sw_ext(:,ilay1:ilay2), &
         &  air_sw_ssa   => volume_props%air_sw_ssa(:,ilay1:ilay2), &
         &  veg_sw_ext   => volume_props%veg_sw_ext(:,ilay1:ilay2), &
         &  veg_sw_ssa   => volume_props%veg_sw_ssa(:,ilay1:ilay2))

    end associate

  end subroutine spartacus_vegetation_sw


end module radsurf_vegetation_sw
