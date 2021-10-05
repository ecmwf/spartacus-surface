! radsurf_norm_perim.F90 - Compute normalized perimeter length between regions
!
! (C) Copyright 2020- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!

module radsurf_norm_perim

contains

  !---------------------------------------------------------------------
  ! This routine computes the normalized perimeter length (perimeter
  ! length divided by domain area), in m-1, between regions in
  ! forests. "nreg" is the number of permeable regions per layer,
  ! which can be 1 (clear), 2 (clear and vegetated) or 3 (clear and
  ! two vegetated). "norm_perim" has nreg elements per layer: if
  ! nreg=1 then norm_perim is unused, if nreg=2 then norm_perim(1,:)
  ! is the length between the two regions and if nreg=3 then
  ! norm_perim(i,:) is the length between regions i and i+1 except for
  ! norm_perim(3,:) which is the length between regions 1 and 3.
  subroutine calc_norm_perim_forest(config, nlay, nreg, &
       &  veg_fraction, veg_scale, norm_perim)
    
    use parkind1,                   only : jpim, jprb
    use yomhook,                    only : lhook, dr_hook
    use radsurf_config,             only : config_type

    implicit none

    ! Algorithm configuration
    type(config_type), intent(in)  :: config

    ! Number of layers and permeable regions
    integer(jpim), intent(in)  :: nlay, nreg
    
    ! Vegetation area fraction, and the horizontal scale of the
    ! vegetation (m), which are converted to edge length using Eqs. 19
    ! or 20 of Hogan et al. (GMD 2018)
    real(jprb), dimension(nlay), intent(in) :: veg_fraction, veg_scale

    ! Normalized perimeter length between regions (m-1)
    real(jprb), dimension(nreg,nlay), intent(out) :: norm_perim

    ! Layer index
    integer(jpim) :: jlay
    
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radsurf_norm_perim:calc_norm_perim_forest',0,hook_handle)

    norm_perim = 0.0_jprb

    do jlay = 1,nlay
      if (nreg > 1) then
        if (veg_fraction(jlay) > config%min_vegetation_fraction) then
          ! Compute the normalized vegetation perimeter length
          if (config%use_symmetric_vegetation_scale_forest) then
            ! If S=veg scale, v=veg fraction and c=clear fraction then
            ! normal formula for normalized perimeter length is
            ! L=4*v*c/S 
            norm_perim(1,jlay) = 4.0_jprb * veg_fraction(jlay) &
                 &  * max(0.0_jprb, 1.0_jprb - veg_fraction(jlay)) &
                 &  / veg_scale(jlay)
          else
            ! The Jensen et al. (JClim 2008) effective diameter "D" is
            ! used in a simpler formula for normalized perimeter
            ! length L=4v/D
            norm_perim(1,jlay) = 4.0_jprb * veg_fraction(jlay) / veg_scale(jlay)
          end if

          if (nreg > 2) then
            ! Share the clear-air/vegetation perimeter between the two
            ! vegetated regions
            norm_perim(nreg,jlay) = config%vegetation_isolation_factor_forest * norm_perim(1,jlay)
            norm_perim(1,jlay) = (1.0_jprb - config%vegetation_isolation_factor_forest) &
                 &             * norm_perim(1,jlay) 
            ! We assume that the horizontal scale of the vegetation
            ! inhomogeneities is the same as the scale of the tree
            ! crowns themselves. Therefore, to compute the interface
            ! between the two vegetated regions, we use the same
            ! formula as before but with the fraction associated with
            ! one of the two vegetated regions, which is half the
            ! total vegetation fraction.
            if (config%use_symmetric_vegetation_scale_forest) then
              norm_perim(2,jlay) = (1.0_jprb - config%vegetation_isolation_factor_forest) &
                   &  * 4.0_jprb * (0.5_jprb*veg_fraction(jlay)) &
                   &  * (1.0_jprb - (0.5_jprb*veg_fraction(jlay))) &
                   &  /  veg_scale(jlay)
            else
              !            norm_perim(2) = (1.0_jprb - config%vegetation_isolation_factor_forest) &
              !                 &  * 4.0_jprb * (0.5_jprb*veg_fraction(jlay)) / veg_scale(jlay)
              ! Lollipop model - see Hogan, Quaife and Braghiere (2018) explaining sqrt(2)
              norm_perim(2,jlay) = (1.0_jprb - config%vegetation_isolation_factor_forest) &
                   &  * 4.0_jprb * veg_fraction(jlay) / (sqrt(2.0_jprb)*veg_scale(jlay))
            end if
          else
            ! Only one vegetated region so the other column of
            ! norm_perim is unused
            norm_perim(2:,jlay) = 0.0_jprb
          end if
        end if
      end if
    end do
      
    if (lhook) call dr_hook('radsurf_norm_perim:calc_norm_perim_forest',1,hook_handle)

  end subroutine calc_norm_perim_forest

  
  !---------------------------------------------------------------------
  ! This routine computes the normalized perimeter lengths (perimeter
  ! length divided by domain area), in m-1, between regions and walls
  ! in an urban environment. "nreg" is the number of permeable regions
  ! per layer, which can be 1 (clear), 2 (clear and vegetated) or 3
  ! (clear and two vegetated). "norm_perim" has nreg elements per
  ! layer: if nreg=1 then norm_perim is unused, if nreg=2 then
  ! norm_perim(1,:) is the length between the two regions and if
  ! nreg=3 then norm_perim(i,:) is the length between regions i and
  ! i+1 except for norm_perim(3,:) which is the length between regions
  ! 1 and 3. "norm_perim_wall(i,:)" is the length of the interface
  ! between a wall and region i.
  subroutine calc_norm_perim_urban(config, nlay, nreg, building_fraction, &
       &  building_scale, veg_fraction, veg_scale, veg_contact_fraction, &
       &  norm_perim, norm_perim_wall)
    
    use parkind1,                   only : jpim, jprb
    use yomhook,                    only : lhook, dr_hook
    use radsurf_config,             only : config_type

    implicit none

    ! Algorithm configuration
    type(config_type), intent(in)  :: config

    ! Number of layers and permeable regions
    integer(jpim), intent(in)  :: nlay, nreg
    
    ! Building area fraction and the effective diameter of the
    ! buildings (m), and vegetation area fraction and the horizontal
    ! scale of the vegetation (m). The latter two are converted to
    ! vegetation edge length using Eqs. 19 or 20 of Hogan et al. (GMD
    ! 2018). The former two use Eq. 19 of this paper.
    real(jprb), dimension(nlay), intent(in) &
         &  :: building_fraction, building_scale, &
         &     veg_fraction, veg_scale, veg_contact_fraction

    ! Normalized perimeter length between regions, and between regions
    ! and walls (m-1)
    real(jprb), dimension(nreg,nlay), intent(out) &
         &  :: norm_perim, norm_perim_wall

    ! Layer index
    integer(jpim) :: jlay
    
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radsurf_norm_perim:calc_norm_perim_urban',0,hook_handle)

    norm_perim = 0.0_jprb
    norm_perim_wall = 0.0_jprb

    do jlay = 1,nlay
      ! First compute the perimeter length between permeable regions
      if (nreg > 1) then
        if (veg_fraction(jlay) > config%min_vegetation_fraction) then
          ! Compute the normalized vegetation perimeter length
          if (config%use_symmetric_vegetation_scale_urban) then
            ! If S=veg scale, v=veg fraction and c=clear fraction then
            ! normal formula for normalized perimeter length is
            ! L=4*v*c/S. But there are buildings present too with
            ! fraction b, and v+c+b=1. So we are dealing only with the
            ! veg+clear fraction, and need to normalize quantities by
            ! v+c: L/(v+c)=4*(v/(v+c))*(c/(v+c))/S leading to
            ! L=4*v*c/(S*(v+c)). Note that this is the length of the
            ! vegetation-clear interface and cannot be reassigned to
            ! the veg-building or clear-building interface.
            norm_perim(1,jlay) = 4.0_jprb * veg_fraction(jlay) &
                 &  * max(0.0_jprb, 1.0_jprb - veg_fraction(jlay) - building_fraction(jlay)) &
                 &  / (max(config%min_building_fraction, &
                 &         1.0_jprb - building_fraction(jlay)) * veg_scale(jlay))
          else
            ! The Jensen et al. (JClim 2008) effective diameter "D" is
            ! used in a simpler formula for normalized perimeter
            ! length L=4v/D. In this case if we normalize L and v by
            ! (v+c) they cancel leaving the same formula:
            norm_perim(1,jlay) = 4.0_jprb * veg_fraction(jlay) / veg_scale(jlay)
          end if

          if (nreg > 2) then
            ! Share the clear-air/vegetation perimeter between the two
            ! vegetated regions
            norm_perim(nreg,jlay) = config%vegetation_isolation_factor_urban * norm_perim(1,jlay)
            norm_perim(1,jlay) = (1.0_jprb - config%vegetation_isolation_factor_urban) &
                 &             * norm_perim(1,jlay) 
            ! We assume that the horizontal scale of the vegetation
            ! inhomogeneities is the same as the scale of the tree
            ! crowns themselves. Therefore, to compute the interface
            ! between the two vegetated regions, we use the same
            ! formula as before but with the fraction associated with
            ! one of the two vegetated regions, which is half the
            ! total vegetation fraction.
            if (config%use_symmetric_vegetation_scale_urban) then
              norm_perim(2,jlay) = (1.0_jprb - config%vegetation_isolation_factor_urban) &
                   &  * 4.0_jprb * (0.5_jprb*veg_fraction(jlay)) &
                   &  * (1.0_jprb - (0.5_jprb*veg_fraction(jlay)) - building_fraction(jlay)) &
                   &  / (max(config%min_building_fraction, &
                 &         1.0_jprb - building_fraction(jlay)) * veg_scale(jlay))
              else
                !            norm_perim(2) = (1.0_jprb - config%vegetation_isolation_factor_urban) &
                !                 &  * 4.0_jprb * (0.5_jprb*veg_fraction(jlay)) / veg_scale(jlay)
                ! Lollipop model - see Hogan, Quaife and Braghiere (2018) explaining sqrt(2)
                norm_perim(2,jlay) = (1.0_jprb - config%vegetation_isolation_factor_urban) &
                     &  * 4.0_jprb * veg_fraction(jlay) / (sqrt(2.0_jprb)*veg_scale(jlay))
              end if
            else
              ! Only one vegetated region so the other column of
              ! norm_perim is unused
              norm_perim(2:,jlay) = 0.0_jprb
            end if
          end if
        end if

        ! Second, compute the normalized length of the interface
        ! between each region and a building wall
        if (building_fraction(jlay) > config%min_building_fraction) then

          ! By default the only region in contact with walls is the
          ! clear region
          norm_perim_wall(1,jlay) = 4.0_jprb * building_fraction(jlay) / building_scale(jlay)

          if (nreg > 1) then
            if (1.0_jprb - veg_fraction(jlay) - building_fraction(jlay) &
                 &  <= config%min_vegetation_fraction) then
              ! There is no clear region (region 1), so all walls must
              ! be in contact with vegetation (region 2 and possibly
              ! 3)
              if (nreg == 2) then
                norm_perim_wall(2,jlay) = norm_perim_wall(1,jlay)
              else
                norm_perim_wall(2,jlay) = norm_perim_wall(1,jlay) &
                     &  * (1.0_jprb - config%vegetation_isolation_factor_urban)
                norm_perim_wall(3,jlay) = norm_perim_wall(1,jlay) &
                     &  * config%vegetation_isolation_factor_urban
              end if
              norm_perim_wall(1,jlay) = 0.0_jprb
            else if (veg_fraction(jlay) > config%min_vegetation_fraction) then
              ! Nominal case: clear and vegetated regions are present
              if (veg_contact_fraction(jlay) > 0.0_jprb) then
                ! Compute normalized length of interface between wall
                ! and any vegetation
                if (nreg == 2) then
                  norm_perim_wall(2,jlay) = norm_perim_wall(1,jlay)*veg_contact_fraction(jlay)
                else
                  norm_perim_wall(2,jlay) = norm_perim_wall(1,jlay)*veg_contact_fraction(jlay) &
                       &  * (1.0_jprb - config%vegetation_isolation_factor_urban)
                  norm_perim_wall(3,jlay) = norm_perim_wall(1,jlay)*veg_contact_fraction(jlay) &
                       &  * config%vegetation_isolation_factor_urban
                end if
                ! Reduce length of interface between wall and clear-air
                norm_perim_wall(1,jlay) = norm_perim_wall(1,jlay) &
                     &                  * (1.0_jprb - veg_contact_fraction(jlay))
              end if
            end if
          ! else no significant vegetation so norm_perim_wall(2:,jlay) = 0
          end if
        end if

    end do

    if (lhook) call dr_hook('radsurf_norm_perim:calc_norm_perim_urban',1,hook_handle)

  end subroutine calc_norm_perim_urban

end module radsurf_norm_perim
