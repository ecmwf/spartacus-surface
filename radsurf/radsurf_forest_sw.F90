! radsurf_forest_sw.F90 - SPARTACUS shortwave solver for forests
!
! Copyright (C) 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radsurf_forest_sw

contains

  subroutine spartacus_forest_sw(config, &
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
    use radtool_calc_matrices_sw_eig,only: calc_matrices_sw_eig
    use radiation_constants,        only : Pi

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
         &  :: trans_diff, ref_diff
    ! Reflectance of a layer to direct radiation, and fraction of
    ! direct radiation that is scattered and passes out through the
    ! base of the layer
    real(kind=jprb), dimension(nsw,nreg*ns,nreg,nlay) &
         &   :: ref_dir, trans_dir_diff
    ! Direct (unscattered) transmittance
    real(kind=jprb), dimension(nsw,nreg,nreg,nlay) :: trans_dir_dir

    ! Area fraction of each region in each layer, plus a pseudo-layer
    ! at the end representing the free atmosphere above
    real(kind=jprb), dimension(nreg,nlay+1) :: frac
    
    ! Components of the Gamma matrices
    real(kind=jprb), dimension(nsw,nreg,nreg)       :: gamma0
    real(kind=jprb), dimension(nsw,nreg*ns,nreg*ns) :: gamma1, gamma2
    real(kind=jprb), dimension(nsw,nreg*ns,nreg)    :: gamma3

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

    ! Tangent of solar zenith angle
    real(kind=jprb) :: tan0

    ! Rate of exchange between regions, excluding the tangent term,
    ! where the dimensions are in the sense of
    ! f_exchange(region_to,region_from)
    real(kind=jprb) :: f_exchange(nreg,nreg)

    ! Extinction (m-1) and single scattering albedo of each region
    real(kind=jprb) :: ext_reg(nsw,nreg), ssa_reg(nsw,nreg)

    ! Optical depth scaling of vegetation optical depth to represent
    ! inhomogeneity
    real(kind=jprb) :: od_scaling(2:nreg)

    ! Diffuse and direct albedo matrices just above and below each
    ! interface. The "above" albedos include surface (interface=1) and
    ! top-of-canopy (interface=nlay+1) so have one element more than
    ! the number of layers. The "below" albedos use the same indexing
    ! of interfaces, but the surface is never used so this dimension
    ! starts at 2.
    real(kind=jprb) :: a_above(nsw,nreg*ns,nreg*ns,nlay+1) ! Diffuse
    real(kind=jprb) :: d_above(nsw,nreg*ns,nreg,nlay+1)    ! Direct
    real(kind=jprb) :: a_below(nsw,nreg*ns,nreg*ns,2:nlay+1)
    real(kind=jprb) :: d_below(nsw,nreg*ns,nreg,2:nlay+1)

    ! Loop counters
    integer(kind=jpim) :: jlay, jreg, jreg_fr, jreg_to, jsw, js, js_fr, js_to

    ! Index to the matrix dimension expressing regions "from" and "to"
    ! in combination with a particular stream
    integer(kind=jpim) :: ifr, ito
    
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radsurf_forest_sw:spartacus_forest_sw',0,hook_handle)

    associate( &
         &  dz           => canopy_props%dz(ilay1:ilay2), &
         &  veg_fraction => canopy_props%veg_fraction(ilay1:ilay2), &
         &  veg_scale    => canopy_props%veg_scale(ilay1:ilay2), &
         &  veg_fsd      => canopy_props%veg_fsd(ilay1:ilay2), &
         &  air_sw_ext   => volume_props%air_sw_ext(:,ilay1:ilay2), &
         &  air_sw_ssa   => volume_props%air_sw_ssa(:,ilay1:ilay2), &
         &  veg_sw_ext   => volume_props%veg_sw_ext(:,ilay1:ilay2), &
         &  veg_sw_ssa   => volume_props%veg_sw_ssa(:,ilay1:ilay2))

      ! Tangent of solar zenith angle
      tan0 = sqrt(1.0_jprb - cos_sza*cos_sza) / max(cos_sza,1.0e-6_jprb)

      ! Set the area fraction of each region
      frac(1,1:nlay)  = 1.0_jprb - veg_fraction
      frac(2:,1:nlay) = spread(veg_fraction,1,nreg-1) / real(nreg-1,jprb)
      frac(1,nlay+1)  = 1.0_jprb
      frac(2:,nlay+1) = 0.0_jprb

      ! Loop up through the canopy computing the Gamma matrices, and
      ! from those the transmission and reflection matrices
      do jlay = 1,nlay

        ! Compute the extinction coefficient and single-scattering
        ! albedo of each region
        ext_reg(:,1) = air_sw_ext(:,jlay)
        ssa_reg(:,1) = veg_sw_ext(:,jlay)
        if (nreg == 2) then
          ext_reg(:,2) = air_sw_ext(:,jlay) + veg_sw_ext(:,jlay)
          ssa_reg(:,2) = (ext_reg(:,1)*ssa_reg(:,1) + veg_sw_ext(:,jlay)*veg_sw_ssa(:,jlay)) &
               &       / ext_reg(:,2)
        else
          ! Approximate method to approximate a Gamma distribution
          od_scaling(2) = exp(-veg_fsd(jlay)*(1.0_jprb + 0.5_jprb*veg_fsd(jlay) &
               &                            *(1.0_jprb + 0.5_jprb*veg_fsd(jlay))))
          od_scaling(3) = 2.0_jprb - od_scaling(2)
          ext_reg(:,2) = air_sw_ext(:,jlay) + od_scaling(2)*veg_sw_ext(:,jlay)
          ext_reg(:,3) = air_sw_ext(:,jlay) + od_scaling(3)*veg_sw_ext(:,jlay)
          ssa_reg(:,2) = (ext_reg(:,1)*ssa_reg(:,1) &
               &          + od_scaling(2)*veg_sw_ext(:,jlay)*veg_sw_ssa(:,jlay)) &
               &       / ext_reg(:,2)
          ssa_reg(:,3) = (ext_reg(:,1)*ssa_reg(:,1) &
               &          + od_scaling(3)*veg_sw_ext(:,jlay)*veg_sw_ssa(:,jlay)) &
               &       / ext_reg(:,3)
        end if

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
          norm_perim(3) = config%vegetation_isolation_factor_forest * norm_perim(3)
          norm_perim(1) = (1.0_jprb - config%vegetation_isolation_factor_forest) &
               &          * norm_perim(1) 
          ! We assume that the horizontal scale of the vegetation
          ! inhomogeneities is the same as the scale of the tree crowns
          ! themselves. Therefore, to compute the interface between the
          ! two vegetated regions, we use the same formula as before but
          ! with the fraction associated with one of the two vegetated
          ! regions, which is half the total vegetation fraction.
          if (config%use_symmetric_vegetation_scale_forest) then
            norm_perim(2) = (1.0_jprb - config%vegetation_isolation_factor_forest) &
                 &  * 4.0_jprb * (0.5_jprb*veg_fraction(jlay)) * (1.0_jprb - (0.5_jprb*veg_fraction(jlay))) &
                 &  / veg_scale(jlay)
          else
            norm_perim(2) = (1.0_jprb - config%vegetation_isolation_factor_forest) &
                 &  * 4.0_jprb * (0.5_jprb*veg_fraction(jlay)) / veg_scale(jlay)
          end if
        else
          ! Only one vegetated region so the other column of norm_perim
          ! is unused
          norm_perim(2:) = 0.0_jprb
        end if

        ! Compute the rates of exchange between regions, excluding the
        ! tangent term
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
        if (nreg > 2 .and. norm_perim(3) > 0.0_jprb) then
          if (frac(3,jlay) <= config%min_vegetation_fraction) then
            f_exchange(1,3) = 0.0_jprb
          else
            f_exchange(1,3) = norm_perim(jreg) / (Pi * frac(3,jlay))
          end if
          if (frac(1,jlay) <= config%min_vegetation_fraction) then
            f_exchange(3,1) = 0.0_jprb
          else
            f_exchange(3,1) = norm_perim(jreg) * tan0 / (Pi * frac(1,jlay))
          end if
        end if

        ! Compute the Gamma matrices representing exchange of direct
        ! and diffuse radiation between regions (gamma0 and gamma1
        ! respectively)
        gamma0 = 0.0_jprb
        gamma1 = 0.0_jprb
        do jreg_fr = 1,nreg
          do jreg_to = 1,nreg
            if (jreg_fr /= jreg_to) then
              gamma0(:,jreg_fr,jreg_fr) = - tan0 * f_exchange(jreg_to,jreg_fr)
              gamma0(:,jreg_to,jreg_fr) = + tan0 * f_exchange(jreg_to,jreg_fr)
              do js = 1,ns
                ifr = js + (jreg_fr-1)*ns
                ito = js + (jreg_to-1)*ns
                gamma1(:,ifr,ifr) = - lg%tan_ang(js) * f_exchange(jreg_to,jreg_fr)
                gamma1(:,ito,ifr) = + lg%tan_ang(js) * f_exchange(jreg_to,jreg_fr)
              end do
            end if
          end do
        end do
        
        ! Diagonal elements of gamma0 and gamma1 also have loss term
        ! due to extinction through the regions
        do jreg = 1,nreg
          gamma0(:,jreg,jreg) = gamma0(:,jreg,jreg) - ext_reg(:,jreg)/cos_sza
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
              gamma2(:,ito,ifr) = 0.5_jprb * lg%hweight(js_fr) &
                   &  * ext_reg(:,jreg) * ssa_reg(:,jreg) / lg%mu(js_fr)
            end do
          end do
        end do

        ! Gamma1 also contains gain of diffuse radiation by
        ! scattering, so add Gamma2 to it
        gamma1 = gamma1 + gamma2

        ! Compute Gamma3, representing scattering from the direct to
        ! the diffuse streams
        gamma3 = 0.0_jprb
        do jreg = 1,nreg
          do js = 1,ns
            ito = js + (jreg-1)*ns
            gamma3(:,ito,jreg) = 0.5_jprb * lg%hweight(js) &
                 &             * ext_reg(:,jreg) * ssa_reg(:,jreg) / cos_sza
          end do
        end do

        ! Compute reflection/transmission matrices for this layer
        call calc_matrices_sw_eig(nsw, nreg*ns, nreg, dz(jlay), cos_sza, &
             &  gamma0, gamma1, gamma2, gamma3, &
             &  ref_diff(:,:,:,jlay), trans_diff(:,:,:,jlay), &
             &  ref_dir(:,:,:,jlay), trans_dir_diff(:,:,:,jlay), &
             &  trans_dir_dir(:,:,:,jlay))

      end do ! Loop over layers to compute reflectance/transmittance matrices

      ! Set the surface albedo
      a_above = 0.0_jprb
      d_above = 0.0_jprb
      do jreg = 1,nreg
        do js_to = 1,ns
          d_above(:,js_to+(jreg-1)*ns,jreg,1) = cos_sza * ground_sw_albedo_direct * lg%hweight(js)
          do js_fr = 1,ns
            a_above(:,js_to+(jreg-1)*ns,jreg_fr+(jreg-1)*ns,1) = ground_sw_albedo * lg%hweight(js)
          end do
        end do
      end do

      ! Loop up through the half levels / interfaces computing albedos
      a_below = 0.0_jprb
      d_below = 0.0_jprb
      do jlay = 1,nlay
        


      end do ! Loop over layers for upward pass to compute albedos

    end associate

    if (lhook) call dr_hook('radsurf_forest_sw:spartacus_forest_sw',1,hook_handle)

  end subroutine spartacus_forest_sw

end module radsurf_forest_sw
