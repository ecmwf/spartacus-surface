! radsurf_forest_lw.F90 - SPARTACUS longwave solver for forests
!
! Copyright (C) 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radsurf_forest_lw

contains

  ! ------------------------------------------------------------------
  ! This routine implements the SPARTACUS longwave radiative transfer
  ! algorithm for computing the propagation of thermal-infrared
  ! radiation in a forest canopy
  ! 
  ! Sections:
  !   1. Declare variables and arrays
  !   2. Prepare general variables and arrays
  !   3. First loop over layers
  !     3a. Prepare the properties of the current layer
  !     3b. Compute Gamma matrices
  !     3c. Compute reflection/transmission matrices for this layer
  !   4. Albedo/emission of scene at each layer interface
  !   5. Compute normalized flux profile
  !
  subroutine spartacus_forest_lw(config, &
       &  nlw, ns, nreg, nlay, icol, ilay1, ilay2, &
       &  lg, canopy_props, lw_spectral_props, &
       &  top_emissivity, top_emission, &
       &  lw_internal, lw_norm)
    
    use parkind1,                   only : jpim, jprb
    use yomhook,                    only : lhook, dr_hook
    use radsurf_config,             only : config_type
    use radtool_legendre_gauss,     only : legendre_gauss_type
    use radsurf_canopy_properties,  only : canopy_properties_type
    use radsurf_lw_spectral_properties,  only : lw_spectral_properties_type
    use radsurf_canopy_flux,        only : canopy_flux_type
    use radtool_calc_matrices_lw_eig,only: calc_matrices_lw_eig
    use radiation_constants,        only : Pi
    use radtool_matrix,             only : identity_minus_mat_x_mat, &
         &  mat_x_mat, singlemat_x_vec, mat_x_vec, rect_mat_x_vec, &
         &  solve_mat, rect_mat_x_mat, rect_expandedmat_x_mat, &
         &  rect_mat_x_expandedmat, rect_expandedmat_x_vec, solve_vec, &
         &  solve_rect_mat, rect_mat_x_singlemat, rect_singlemat_x_vec
    use radsurf_overlap,            only : calc_overlap_matrices

!#define PRINT_ARRAYS 1

#ifdef PRINT_ARRAYS
    use print_matrix_mod
#endif
    
    implicit none

    ! --------------------------------------------------------
    ! Section 1: Declare variables and arrays
    ! --------------------------------------------------------

    ! Inputs

    ! Algorithm configuration
    type(config_type),             intent(in)  :: config
    ! Number of spectral intervals, number of layers
    integer(kind=jpim),            intent(in)  :: nlw, nlay
    ! Index of current column and first and last layer of the current column
    integer(kind=jpim),            intent(in)  :: icol, ilay1, ilay2
    ! Number of regions, number of diffuse streams in each hemisphere
    integer(kind=jpim),            intent(in)  :: nreg, ns
    ! Legendre-Gauss coefficients
    type(legendre_gauss_type),     intent(in)  :: lg
    ! Geometric and other spectrally independent properties of the canopy
    type(canopy_properties_type),  intent(in)  :: canopy_props
    ! Spectral properties of the air and vegetation
    type(lw_spectral_properties_type),  intent(in)  :: lw_spectral_props

    ! Outputs

    ! Top-of-canopy spectral emissivity and emission (W m-2)
    real(kind=jprb), dimension(nlw),intent(out):: top_emissivity, &
         &                                        top_emission
    ! Flux outputs
    type(canopy_flux_type), intent(inout), optional &
         &  :: lw_internal, & ! LW fluxes from internal emission
         &     lw_norm        ! LW fluxes normalized by top-of-canopy down

    ! Local variables and arrays

    ! Transmittance and reflectance of a layer
    real(kind=jprb), dimension(nlw,nreg*ns,nreg*ns,nlay) &
         &  :: trans, ref

    ! Radiation emerging up from the top or down from the base of a
    ! particular layer (in a given region and stream) due to emission
    ! within the layer
    real(kind=jprb), dimension(nlw,nreg*ns,nlay) :: source_lay

    ! Area fraction of each region in each layer, plus a pseudo-layer
    ! at the end representing the free atmosphere above
    real(kind=jprb), dimension(nreg,nlay+1) :: frac
    
    ! Components of the Gamma matrix
    real(kind=jprb), dimension(nlw,nreg*ns,nreg*ns) :: gamma1, gamma2

    ! Emission rate per unit height, "b" in Eq. 32 (W m-3)
    real(kind=jprb), dimension(nlw,nreg*ns) :: emiss_rate

    ! Emission rate per unit height for walls, air and vegetation (the
    ! last two in regions) considering all emission directions (W m-3)
    real(kind=jprb) :: emiss_wall(nlw,nlay), emiss_reg(nlw,nreg,nlay), &
         &             emiss_air(nlw,2:nreg,nlay), emiss_veg(nlw,2:nreg,nlay)

    ! Geometric factor needed in computing the emiss_* arrays above
    real(kind=jprb) :: emiss_factor
    
    ! Normalized vegetation perimeter length (perimeter length divided
    ! by domain area), m-1.  If nreg=2 then there is a clear-sky and a
    ! vegetation region, and norm_perim(1) is the normalized length
    ! between the two regions, while norm_perim(2) is unused.  If
    ! nreg=3 then region 1 is clear-sky, region 2 is low optical depth
    ! vegetation and region 3 is high optical depth
    ! vegetation. norm_perim(1) is the normalized length between
    ! regions 1 and 2, norm_perim(2) is that between regions 2 and 3,
    ! and norm_perim(3) is that between regions 3 and 1.
    real(kind=jprb) :: norm_perim(nreg)

    ! Rate of exchange between regions, excluding the tangent term,
    ! where the dimensions are in the sense of
    ! f_exchange(region_to,region_from)
    real(kind=jprb) :: f_exchange(nreg,nreg)

    ! Spectral emission from the air/leaves in a region
    real(kind=jprb) :: volume_emiss(nlw)
    
    ! Extinction (m-1), single scattering albedo and Planck function
    ! of each region
    real(kind=jprb) :: ext_reg(nlw,nreg), ssa_reg(nlw,nreg), planck_reg(nlw,nreg)

    ! Optical depth scaling of vegetation optical depth to represent
    ! inhomogeneity
    real(kind=jprb) :: od_scaling(2:nreg,nlay)

    ! Albedo matrix just above and below each interface. The "above"
    ! albedo includes surface (interface=1) and top-of-canopy
    ! (interface=nlay+1) so has one element more than the number of
    ! layers. The "below" albedo uses the same indexing of interfaces,
    ! but the surface is never used so this dimension starts at 2.
    real(kind=jprb) :: a_above(nlw,nreg*ns,nreg*ns,nlay+1)
    real(kind=jprb) :: a_below(nlw,nreg*ns,nreg*ns,2:nlay+1)

    ! Upward radiation just above and below a layer interface due
    ! purely to emission below that interface. The dimensioning is as
    ! for a_above and a_below but with one fewer dimension.
    real(kind=jprb) :: source_above(nlw,nreg*ns,nlay+1)
    real(kind=jprb) :: source_below(nlw,nreg*ns,1:nlay+1)
    
    ! Denominator matrix used in the matrix adding method
    real(kind=jprb) :: denominator(nlw,nreg*ns,nreg*ns,nlay)

    ! Directional overlap matrices expressing how the fluxes in each
    ! region of one layer pass into the layer above (u_overlap) or
    ! below (v_overlap)
    real(kind=jprb), dimension(nreg,nreg,nlay+1) :: u_overlap
    real(kind=jprb), dimension(nreg,nreg,nlay+1) :: v_overlap

    ! Fluxes just above or just below a layer interface, either due to
    ! canopy sources or normalized by downwelling longwave flux at
    ! canopy top
    real(kind=jprb), dimension(nlw,nreg*ns)     :: flux_dn_above, flux_up_above
    real(kind=jprb), dimension(nlw,nreg*ns) :: flux_dn_below, flux_up_below
    
    ! The following matrices are defined such that the integrated flux
    ! (u_hat+v_hat in Eq. 29) is
    !   int_flux_mat * (u_conv+v_conv) + int_source
    ! where u_conv+v_conv is the convergence of fluxes into the layer,
    ! i.e. the sum of the fluxes entering the layer from the top or
    ! base, minus the fluxes leaving the layer from top or
    ! base. "int_source" is the integrated fluxes (u_hat+v_hat) due to
    ! emission within the layer.
    real(kind=jprb) &
         &  :: int_flux_mat(nlw,nreg*ns,nreg*ns,nlay), &
         &     int_source(nlw,nreg*ns,nlay)

    ! Integrated flux across a layer
    real(kind=jprb) :: int_flux(nlw,nreg*ns)

    ! Loop counters
    integer(kind=jpim) :: jlay, jreg, jreg_fr, jreg_to, js, js_fr, js_to

    ! Index to layer inside canopy_flux object
    integer(kind=jpim) :: ilay

    ! Index to the matrix dimension expressing regions "from" and "to"
    ! in combination with a particular stream
    integer(kind=jpim) :: ifr, ito
    
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radsurf_forest_lw:spartacus_forest_lw',0,hook_handle)

    associate( &
         &  dz                    => canopy_props%dz(ilay1:ilay2), &
         &  veg_fraction          => canopy_props%veg_fraction(ilay1:ilay2), &
         &  veg_scale             => canopy_props%veg_scale(ilay1:ilay2), &
         &  veg_ext               => canopy_props%veg_ext(ilay1:ilay2), &
         &  veg_fsd               => canopy_props%veg_fsd(ilay1:ilay2), &
         &  veg_ssa               => lw_spectral_props%veg_ssa(:,ilay1:ilay2), &
         &  veg_planck            => lw_spectral_props%veg_planck(:,ilay1:ilay2), &
         &  air_ext               => lw_spectral_props%air_ext(:,ilay1:ilay2), &
         &  air_ssa               => lw_spectral_props%air_ssa(:,ilay1:ilay2), &
         &  clear_air_planck      => lw_spectral_props%clear_air_planck(:,ilay1:ilay2), &
         &  veg_air_planck        => lw_spectral_props%veg_air_planck(:,ilay1:ilay2), &
         &  ground_emissivity     => lw_spectral_props%ground_emissivity(:,icol), &
         &  ground_emission       => lw_spectral_props%ground_emission(:,icol))


      ! --------------------------------------------------------
      ! Section 2: Prepare general variables and arrays
      ! --------------------------------------------------------
      
      ! Set the area fraction of each region
      frac(1,nlay+1)  = 1.0_jprb ! Free atmosphere
      frac(1,1:nlay)  = 1.0_jprb - canopy_props%veg_fraction(ilay1:ilay2)
      frac(2:,1:nlay) = spread(1.0_jprb - frac(1,1:nlay),1,nreg-1) &
           &          / real(nreg-1,jprb)
      frac(2:,nlay+1) = 0.0_jprb

#ifdef PRINT_ARRAYS
      call print_vector('veg_fraction',veg_fraction)
      call print_vector('veg_scale', veg_scale);
      call print_matrix('frac', frac);
#endif

      ! Compute overlap matrices
      call calc_overlap_matrices(nlay,nreg,frac,u_overlap,v_overlap, &
           &  config%min_vegetation_fraction);

#ifdef PRINT_ARRAYS
      call print_array3('u_overlap',u_overlap)
      call print_array3('v_overlap',v_overlap)
#endif      

      ! --------------------------------------------------------
      ! Section 3: First loop over layers
      ! --------------------------------------------------------
      ! Loop up through the canopy computing the Gamma matrices, and
      ! from those the transmission and reflection matrices
      do jlay = 1,nlay

        ! --------------------------------------------------------
        ! Section 3a: Prepare the properties of the current layer
        ! --------------------------------------------------------
        ! Compute the extinction coefficient and single-scattering
        ! albedo of each region
        ext_reg(:,1)    = air_ext(:,jlay)
        ssa_reg(:,1)    = air_ssa(:,jlay)
        planck_reg(:,1) = clear_air_planck(:,jlay)
        if (nreg == 2) then
          ext_reg(:,2) = air_ext(:,jlay) + veg_ext(jlay)
          ssa_reg(:,2) = (ext_reg(:,1)*ssa_reg(:,1) + veg_ext(jlay)*veg_ssa(:,jlay)) &
               &       / max(ext_reg(:,2), 1.0e-8_jprb)
          planck_reg(:,2) = (ext_reg(:,1)*(1.0_jprb-ssa_reg(:,1))*veg_air_planck(:,jlay) &
               &   + veg_ext(jlay)*(1.0_jprb-veg_ssa(:,jlay))*veg_planck(:,jlay)) &
               &       / max(ext_reg(:,2)*(1.0_jprb-ssa_reg(:,2)), 1.0e-8_jprb)
          od_scaling(2,jlay) = 1.0_jprb
        else if (nreg == 3) then
          ! Approximate method to approximate a Gamma distribution
          od_scaling(2,jlay) = exp(-veg_fsd(jlay)*(1.0_jprb + 0.5_jprb*veg_fsd(jlay) &
               &                            *(1.0_jprb + 0.5_jprb*veg_fsd(jlay))))
          od_scaling(3,jlay) = 2.0_jprb - od_scaling(2,jlay)
          ext_reg(:,2) = air_ext(:,jlay) + od_scaling(2,jlay)*veg_ext(jlay)
          ext_reg(:,3) = air_ext(:,jlay) + od_scaling(3,jlay)*veg_ext(jlay)
          ssa_reg(:,2) = (ext_reg(:,1)*ssa_reg(:,1) &
               &          + od_scaling(2,jlay)*veg_ext(jlay)*veg_ssa(:,jlay)) &
               &       / max(ext_reg(:,2), 1.0e-8_jprb)
          ssa_reg(:,3) = (ext_reg(:,1)*ssa_reg(:,1) &
               &          + od_scaling(3,jlay)*veg_ext(jlay)*veg_ssa(:,jlay)) &
               &       / max(ext_reg(:,3), 1.0e-8_jprb)
          planck_reg(:,2) = (ext_reg(:,1)*(1.0_jprb-ssa_reg(:,1))*veg_air_planck(:,jlay) &
               &   + od_scaling(2,jlay) * veg_ext(jlay)*(1.0_jprb-veg_ssa(:,jlay))*veg_planck(:,jlay)) &
               &       / max(ext_reg(:,2)*(1.0_jprb-ssa_reg(:,2)), 1.0e-8_jprb)
          planck_reg(:,3) = (ext_reg(:,1)*(1.0_jprb-ssa_reg(:,1))*veg_air_planck(:,jlay) &
               &   + od_scaling(3,jlay) * veg_ext(jlay)*(1.0_jprb-veg_ssa(:,jlay))*veg_planck(:,jlay)) &
               &       / max(ext_reg(:,3)*(1.0_jprb-ssa_reg(:,3)), 1.0e-8_jprb)
        end if

        norm_perim = 0.0_jprb
        if (nreg > 1) then
          ! Compute the normalized vegetation perimeter length
          if (config%use_symmetric_vegetation_scale_forest) then
            norm_perim(1) = 4.0_jprb * veg_fraction(jlay) * (1.0_jprb - veg_fraction(jlay)) &
                 &        / veg_scale(jlay)
          else
            norm_perim(1) = 4.0_jprb * veg_fraction(jlay) / veg_scale(jlay)
          end if

          if (nreg > 2) then
            ! Share the clear-air/vegetation perimeter between the two
            ! vegetated regions
            norm_perim(nreg) = config%vegetation_isolation_factor_forest * norm_perim(1)
            norm_perim(1) = (1.0_jprb - config%vegetation_isolation_factor_forest) &
                 &          * norm_perim(1) 
            ! We assume that the horizontal scale of the vegetation
            ! inhomogeneities is the same as the scale of the tree
            ! crowns themselves. Therefore, to compute the interface
            ! between the two vegetated regions, we use the same
            ! formula as before but with the fraction associated with
            ! one of the two vegetated regions, which is half the
            ! total vegetation fraction.
            if (config%use_symmetric_vegetation_scale_forest) then
              norm_perim(2) = (1.0_jprb - config%vegetation_isolation_factor_forest) &
                   &  * 4.0_jprb * (0.5_jprb*veg_fraction(jlay)) &
                   &  * (1.0_jprb - (0.5_jprb*veg_fraction(jlay))) &
                   &  / veg_scale(jlay)
            else
              !            norm_perim(2) = (1.0_jprb - config%vegetation_isolation_factor_forest) &
              !                 &  * 4.0_jprb * (0.5_jprb*veg_fraction(jlay)) / veg_scale(jlay)
              ! Lollipop model - see Hogan, Quaife and Braghiere (2018) explaining sqrt(2)
              norm_perim(2) = (1.0_jprb - config%vegetation_isolation_factor_forest) &
                   &  * 4.0_jprb * veg_fraction(jlay) / (sqrt(2.0_jprb)*veg_scale(jlay))
            end if
          else
            ! Only one vegetated region so the other column of norm_perim
            ! is unused
            norm_perim(2:) = 0.0_jprb
          end if
        end if
        
        ! Compute the rates of exchange between regions, excluding the
        ! tangent term
        f_exchange = 0.0_jprb
        do jreg = 1,nreg-1
          if (frac(jreg,jlay) <= config%min_vegetation_fraction) then
            f_exchange(jreg+1,jreg) = 0.0_jprb
          else
            f_exchange(jreg+1,jreg) = norm_perim(jreg) / (Pi * frac(jreg,jlay))
          end if
          if (frac(jreg+1,jlay) <= config%min_vegetation_fraction) then
            f_exchange(jreg,jreg+1) = 0.0_jprb
          else
            f_exchange(jreg,jreg+1) = norm_perim(jreg) / (Pi * frac(jreg+1,jlay))
          end if
        end do
        if (nreg > 2 .and. norm_perim(nreg) > 0.0_jprb) then
          if (frac(3,jlay) <= config%min_vegetation_fraction) then
            f_exchange(1,3) = 0.0_jprb
          else
            f_exchange(1,3) = norm_perim(jreg) / (Pi * frac(3,jlay))
          end if
          if (frac(1,jlay) <= config%min_vegetation_fraction) then
            f_exchange(3,1) = 0.0_jprb
          else
            f_exchange(3,1) = norm_perim(jreg) / (Pi * frac(1,jlay))
          end if
        end if

        ! --------------------------------------------------------
        ! Section 3b: Compute Gamma matrices
        ! --------------------------------------------------------
        ! Compute the Gamma matrix representing exchange of
        ! radiation between regions
        gamma1 = 0.0_jprb
        do jreg_fr = 1,nreg
          do jreg_to = 1,nreg
            if (jreg_fr /= jreg_to) then
              do js = 1,ns
                ifr = js + (jreg_fr-1)*ns
                ito = js + (jreg_to-1)*ns
                gamma1(:,ifr,ifr) = gamma1(:,ifr,ifr) &
                     &  - lg%tan_ang(js) * f_exchange(jreg_to,jreg_fr)
                gamma1(:,ito,ifr) = + lg%tan_ang(js) * f_exchange(jreg_to,jreg_fr)
              end do
            end if
          end do
        end do

        ! Diagonal elements of gamma1 also has loss term due to
        ! extinction through the regions
        do jreg = 1,nreg
          do js = 1,ns
            ifr = js + (jreg-1)*ns
            gamma1(:,ifr,ifr) = gamma1(:,ifr,ifr) - ext_reg(:,jreg)/lg%mu(js)
          end do
        end do

        ! Compute Gamma2, representing the rate of scattering of
        ! radiation from one diffuse stream to another
        gamma2 = 0.0_jprb
        do js_fr = 1,ns
          do js_to = 1,ns
            do jreg = 1,nreg
              ifr = js_fr + (jreg-1)*ns
              ito = js_to + (jreg-1)*ns
              ! Factor of 0.5 is because half of energy goes into
              ! other hemisphere
              gamma2(:,ito,ifr) = (0.5_jprb * lg%weight(js_to) / lg%mu(js_fr)) &
                   &            * ext_reg(:,jreg) * ssa_reg(:,jreg)
            end do
          end do
        end do

        ! Gamma1 also contains gain of diffuse radiation by
        ! scattering, so add Gamma2 to it
        gamma1 = gamma1 + gamma2

        ! Compute emission rates
        emiss_factor = 2.0_jprb * sum(lg%hweight / lg%mu)
        do jreg = 1,nreg
          volume_emiss = frac(jreg,jlay) * (ext_reg(:,jreg) * (1.0_jprb - ssa_reg(:,jreg)) &
               &                            * planck_reg(:,jreg))
          do js = 1,ns
            ifr = js + (jreg-1)*ns
            emiss_rate(:,ifr) = (lg%hweight(js) / lg%mu(js)) * volume_emiss
          end do
          emiss_reg(:,jreg,jlay) = emiss_factor * volume_emiss

          if (jreg > 1) then
            ! Note that in the following we use the extinction and
            ! single-scattering albedo of the clear air region
            emiss_air(:,jreg,jlay) = emiss_factor * frac(jreg,jlay) * ext_reg(:,1) &
                 &              * (1.0_jprb-ssa_reg(:,1)) * veg_air_planck(:,jlay)
            emiss_veg(:,jreg,jlay) = emiss_factor * frac(jreg,jlay) * veg_ext(jlay) &
                 &   * (1.0_jprb-veg_ssa(:,jlay)) * veg_planck(:,jlay) * od_scaling(jreg,jlay)
          end if
        end do
        
        
#ifdef PRINT_ARRAYS
        print *, 'PROPERTIES OF LAYER ', jlay
        call print_vector('ext_reg',ext_reg(1,:))
        call print_vector('ssa_reg',ssa_reg(1,:))
        call print_matrix('f_exchange',f_exchange)
        call print_vector('norm_perim', norm_perim)
        call print_matrix('gamma1',gamma1(1,:,:))
        call print_matrix('gamma2',gamma2(1,:,:))
        call print_vector('emiss_rate',emiss_rate(1,:))
        call print_vector('clear_air_planck', clear_air_planck(1,:))
        call print_vector('veg_air_planck', veg_air_planck(1,:))
        call print_vector('veg_planck', veg_planck(1,:))
#endif        

        ! --------------------------------------------------------
        ! Section 3c: Compute reflection/transmission matrices for this layer
        ! --------------------------------------------------------
        call calc_matrices_lw_eig(nlw, nreg*ns, dz(jlay), &
             &  gamma1, gamma2, emiss_rate, &
             &  ref(:,:,:,jlay), trans(:,:,:,jlay), source_lay(:,:,jlay), &
             &  int_flux_mat(:,:,:,jlay), int_source(:,:,jlay))

      end do ! Loop over layers to compute reflectance/transmittance matrices

      ! --------------------------------------------------------
      ! Section 4: Albedo of scene at each layer interface
      ! --------------------------------------------------------
      ! Set the surface albedo and upward emission
      a_above = 0.0_jprb
      source_above = 0.0_jprb
      do jreg = 1,nreg
        do js_to = 1,ns
          do js_fr = 1,ns
            a_above(:,js_to+(jreg-1)*ns,js_fr+(jreg-1)*ns,1) &
                 &  = (1.0_jprb - ground_emissivity) * lg%hweight(js_to)
          end do
        end do
        do js = 1,ns
          source_above(:,js+(jreg-1)*ns,1) = (lg%hweight(js) * frac(jreg,1)) &
               &                           * ground_emission
        end do
      end do

      ! Loop up through the half levels / interfaces computing albedos
      a_below = 0.0_jprb
      source_below = 0.0_jprb
      do jlay = 1,nlay
        ! Adding method
        denominator(:,:,:,jlay) = identity_minus_mat_x_mat(nlw,nlw,nreg*ns, &
             &  a_above(:,:,:,jlay), ref(:,:,:,jlay))
        a_below(:,:,:,jlay+1) = ref(:,:,:,jlay) &
             &  + mat_x_mat(nlw,nlw,nreg*ns,trans(:,:,:,jlay), &
             &  solve_mat(nlw,nlw,nreg*ns,denominator(:,:,:,jlay), &
             &  mat_x_mat(nlw,nlw,nreg*ns,a_above(:,:,:,jlay), &
             &  trans(:,:,:,jlay))))
        ! Eq. 34
        source_below(:,:,jlay+1) = source_lay(:,:,jlay) &
             &  + mat_x_vec(nlw,nlw,nreg*ns,trans(:,:,:,jlay), &
             &  solve_vec(nlw,nlw,nreg*ns,denominator(:,:,:,jlay), &
             &  source_above(:,:,jlay) &
             &  + mat_x_vec(nlw,nlw,nreg*ns,a_above(:,:,:,jlay),source_lay(:,:,jlay))))

        ! Overlap: Hogan (2019), equations 22 and 36
        a_above(:,:,:,jlay+1) = rect_expandedmat_x_mat(nlw,nreg,nreg,ns,nreg*ns, &
             &  u_overlap(:,:,jlay+1), &
             &  rect_mat_x_expandedmat(nlw,nreg,nreg,ns,nreg*ns, a_below(:,:,:,jlay+1), &
             &                          v_overlap(:,:,jlay+1)))
        source_above(:,:,jlay+1) = rect_expandedmat_x_vec(nlw,nreg,nreg,ns, &
             &  u_overlap(:,:,jlay+1), source_below(:,:,jlay+1))
      end do ! Loop over layers for upward pass to compute albedo

#ifdef PRINT_ARRAYS
      call print_array3('a_above', a_above(1,:,:,:))
      call print_matrix('source_above', source_above(1,:,:))
      call print_matrix('source_below', source_below(1,:,:))
      call print_array3('a_below', a_below(1,:,:,:))
      call print_array3('T', trans(1,:,:,:))
      call print_array3('R', ref(1,:,:,:))
      call print_matrix('source_lay', source_lay(1,:,:))
      call print_array3('int_flux_mat', int_flux_mat(1,:,:,:))
      call print_matrix('int_source', int_source(1,:,:))
      call print_matrix('emiss_reg', emiss_reg(1,:,:))
#endif

      ! Store top-of-canopy boundary conditions.  Isotropic emissivity
      ! is weighted average over the ns streams, noting that just
      ! above the "nlay+1" interface we are above the canopy so only
      ! need to consider the clear-sky region (indexed 1:ns).
      top_emissivity = 1.0_jprb - sum(mat_x_vec(nlw,nlw,ns,a_above(:,1:ns,1:ns,nlay+1), &
           &                             spread(lg%hweight,1,nlw)),2)
      ! Top-of-canopy emission is the sum over the streams in the
      ! clear-sky region
      top_emission = sum(source_above(:,1:ns,nlay+1),2)
      
      ! --------------------------------------------------------
      ! Section 5: Compute normalized flux profile
      ! --------------------------------------------------------

      ! Set to the flux components to zero initially
      call lw_internal%zero(icol, ilay1, ilay2)
      call lw_norm%zero( icol, ilay1, ilay2)

      ! First the fluxes associated with emission within the canopy

      ! Initial normalized downward flux is zero
      flux_dn_above = 0.0_jprb
      
      lw_internal%top_dn(:,icol)  = 0.0_jprb
      lw_internal%top_net(:,icol) = -top_emission

      ! Loop down through layers
      do jlay = nlay,1,-1
        ! Find index into output arrays
        ilay = ilay1 + jlay - 1

        ! Translate the downwelling flux component across the
        ! interface at the top of the layer
        flux_dn_below = rect_expandedmat_x_vec(nlw,nreg,nreg,ns, & 
             &  v_overlap(:,:,jlay+1), flux_dn_above)
        flux_up_below = mat_x_vec(nlw,nlw,nreg*ns,a_below(:,:,:,jlay+1),flux_dn_below) &
             &  + source_below(:,:,jlay+1)
        
        ! Compute fluxes at base of layer
        flux_dn_above = solve_vec(nlw,nlw,nreg*ns,denominator(:,:,:,jlay), &
             &  (mat_x_vec(nlw,nlw,nreg*ns,trans(:,:,:,jlay),flux_dn_below) &
             &   + mat_x_vec(nlw,nlw,nreg*ns,ref(:,:,:,jlay),source_above(:,:,jlay)) &
             &   + source_lay(:,:,jlay)))
        flux_up_above = mat_x_vec(nlw,nlw,nreg*ns,a_above(:,:,:,jlay),flux_dn_above) &
             &  + source_above(:,:,jlay)

        ! Compute integrated flux vectors, recalling that _above means
        ! above the just above the *base* of the layer, and _below
        ! means just below the *top* of the layer
        int_flux = mat_x_vec(nlw,nlw,nreg*ns,int_flux_mat(:,:,:,jlay), &
             &               flux_dn_below + flux_up_above) &
             &   + int_source(:,:,jlay)
       
        ! Absorption by clear-air region - see Eqs. 29 and 30
        lw_internal%clear_air_abs(:,ilay) = lw_internal%clear_air_abs(:,ilay) &
             &  + air_ext(:,jlay)*(1.0_jprb-air_ssa(:,jlay)) &
             &    * sum(int_flux(:,1:ns) * spread(1.0_jprb/lg%mu,nlw,1), 2) &
             &  - emiss_reg(:,1,jlay)*dz(jlay)
        do jreg = 2,nreg
          ! Absorption by clear-air in the vegetated regions
          lw_internal%veg_air_abs(:,ilay) = lw_internal%veg_air_abs(:,ilay) &
               &  + air_ext(:,jlay)*(1.0_jprb-air_ssa(:,jlay)) & ! Use clear-air properties
               &    * sum(int_flux(:,(jreg-1)*ns+1:jreg*ns) &
               &             * spread(1.0_jprb/lg%mu,nlw,1), 2) &
               &  - emiss_air(:,jreg,jlay)*dz(jlay)
          lw_internal%veg_abs(:,ilay) = lw_internal%veg_abs(:,ilay) &
               &  + veg_ext(jlay)*(1.0_jprb-veg_ssa(:,jlay)) & ! Use vegetation properties
               &    * sum(int_flux(:,(jreg-1)*ns+1:jreg*ns) &
               &             * spread(1.0_jprb/lg%mu,nlw,1), 2) * od_scaling(jreg,jlay) &
               &  - emiss_veg(:,jreg,jlay)*dz(jlay)
        end do

#ifdef PRINT_ARRAYS
        print *, 'ABSOLUTE FLUXES DUE TO INTERNAL EMISSION AT LAYER ', jlay
        call print_vector('  flux_dn_below ', flux_dn_below(1,:))
        call print_vector('  flux_dn_above ', flux_dn_above(1,:))
        call print_vector('  flux_up_above ', flux_up_above(1,:))
#endif

      end do
      lw_internal%ground_dn(:,icol) = sum(flux_dn_above,2)
      lw_internal%ground_net(:,icol) = lw_internal%ground_dn(:,icol) &
           &  - sum(flux_up_above,2)

      ! Second the fluxes due to and normalized by the downwelling
      ! flux at canopy top. This is identical to the shortwave case.

      ! Initial normalized diffuse flux sums to 1 in each stream, in
      ! each spectral interval, so use the Legendre-Gauss horizontal
      ! weights
      flux_dn_above           = 0.0_jprb
      flux_dn_above(:,1:ns)   = spread(lg%hweight,nlw,1)

      lw_norm%top_dn(:,icol)  = 1.0_jprb
      lw_norm%top_net(:,icol) = top_emissivity
      
      ! Loop down through layers
      do jlay = nlay,1,-1
        ! Find index into output arrays
        ilay = ilay1 + jlay - 1

        ! Translate the downwelling flux component across the
        ! interface at the top of the layer
        flux_dn_below = rect_expandedmat_x_vec(nlw,nreg,nreg,ns, &
             &  v_overlap(:,:,jlay+1), flux_dn_above)
        flux_up_below = mat_x_vec(nlw,nlw,nreg*ns,a_below(:,:,:,jlay+1),flux_dn_below)

        ! Compute fluxes at base of layer
        flux_dn_above = solve_vec(nlw,nlw,nreg*ns,denominator(:,:,:,jlay), &
             &  mat_x_vec(nlw,nlw,nreg*ns,trans(:,:,:,jlay),flux_dn_below))
        flux_up_above = mat_x_vec(nlw,nlw,nreg*ns,a_above(:,:,:,jlay), &
             &  flux_dn_above)

        ! Compute integrated flux vectors, recalling that _above means
        ! above the just above the *base* of the layer, and _below
        ! means just below the *top* of the layer
        !int_flux = mat_x_vec(nlw,nlw,nreg*ns,int_flux_mat(:,:,:,jlay), flux_dn_below(:,1:nreg*ns) &
        !     &                   - flux_dn_above - flux_up_below(:,1:nreg*ns) + flux_up_above)
        int_flux = mat_x_vec(nlw,nlw,nreg*ns,int_flux_mat(:,:,:,jlay), &
             &               flux_dn_below + flux_up_above)

        ! Absorption by clear-air region - see Eqs. 29 and 30
        lw_norm%clear_air_abs(:,ilay) = lw_norm%clear_air_abs(:,ilay) &
             &  + air_ext(:,jlay)*(1.0_jprb-air_ssa(:,jlay)) &
             &    * sum(int_flux(:,1:ns) * spread(1.0_jprb/lg%mu,nlw,1), 2)
        do jreg = 2,nreg
          ! Absorption by clear-air in the vegetated regions
          lw_norm%veg_air_abs(:,ilay) = lw_norm%veg_air_abs(:,ilay) &
               &  + air_ext(:,jlay)*(1.0_jprb-air_ssa(:,jlay)) & ! Use clear-air properties
               &    * sum(int_flux(:,(jreg-1)*ns+1:jreg*ns) &
               &             * spread(1.0_jprb/lg%mu,nlw,1), 2)
          lw_norm%veg_abs(:,ilay) = lw_norm%veg_abs(:,ilay) &
               &  + veg_ext(jlay)*(1.0_jprb-veg_ssa(:,jlay)) & ! Use vegetation properties
               &    * sum(int_flux(:,(jreg-1)*ns+1:jreg*ns) &
               &             * spread(1.0_jprb/lg%mu,nlw,1), 2) * od_scaling(jreg,jlay)
        end do

#ifdef PRINT_ARRAYS
        print *, 'NORMALIZED FLUXES W.R.T. DIFFUSE INCOMING RADIATION AT LAYER ', jlay
        call print_vector('  flux_dn_below ', flux_dn_below(1,:))
        call print_vector('  flux_dn_above ', flux_dn_above(1,:))
        call print_vector('  flux_up_above ', flux_up_above(1,:))
#endif

      end do
      lw_norm%ground_dn(:,icol) = sum(flux_dn_above,2)
      lw_norm%ground_net(:,icol) = lw_norm%ground_dn(:,icol) &
           &  - sum(flux_up_above,2)
    
#ifdef PRINT_ARRAYS
      print *, 'ABSOLUTE FLUXES DUE TO INTERNAL EMISSION'
      call print_vector('  clear_air_abs ', lw_internal%clear_air_abs(1,ilay1:ilay2))
      call print_vector('  veg_air_abs ', lw_internal%veg_air_abs(1,ilay1:ilay2))
      call print_vector('  veg_abs ', lw_internal%veg_abs(1,ilay1:ilay2))
      print *, '  ground_dn  = ', lw_internal%ground_dn(1,icol)
      print *, '  ground_net = ', lw_internal%ground_net(1,icol)
      call print_vector('  wall_in', lw_internal%wall_in(1,ilay1:ilay2))
      call print_vector('  wall_net', lw_internal%wall_net(1,ilay1:ilay2))
      print *, 'NORMALIZED FLUXES W.R.T. DIFFUSE INCOMING RADIATION'
      call print_vector('  clear_air_abs ', lw_norm%clear_air_abs(1,ilay1:ilay2))
      call print_vector('  veg_air_abs ', lw_norm%veg_air_abs(1,ilay1:ilay2))
      call print_vector('  veg_abs ', lw_norm%veg_abs(1,ilay1:ilay2))
      print *, '  ground_dn  = ', lw_norm%ground_dn(1,icol)
      print *, '  ground_net = ', lw_norm%ground_net(1,icol)
      call print_vector('  wall_in', lw_norm%wall_in(1,ilay1:ilay2))
      call print_vector('  wall_net', lw_norm%wall_net(1,ilay1:ilay2))
#endif

    end associate

    if (lhook) call dr_hook('radsurf_forest_lw:spartacus_forest_lw',1,hook_handle)

  end subroutine spartacus_forest_lw

end module radsurf_forest_lw
