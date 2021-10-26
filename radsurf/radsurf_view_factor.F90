! radsurf_view_factor.F90 - Compute view factors for single-layer urban canopy
!
! (C) Copyright 2021- ECMWF.
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

! Equation numbers in this file refer to Hogan (BLM 2019) "An
! exponential model of urban geometry for use in radiative transfer
! applications"

module radsurf_view_factor

contains

  !---------------------------------------------------------------------
  ! Calculate view factors making the single-layer infinite street
  ! assumption
  !  elemental 
  subroutine calc_view_factors_inf(height_width_ratio, &
       &  view_ground_sky, view_wall_wall, cos_sza, view_dir_ground)
    
    use parkind1, only : jpim, jprb
    use radiation_constants, only : Pi

    implicit none

    ! Ratio of building height to street width in the single-layer
    ! infinite street canyon assumption (H/W)
    real(kind=jprb), intent(in)  :: height_width_ratio
    ! Fraction of diffuse radiation emanating from the ground that
    ! reaches the sky and fraction emanating from a wall that reaches
    ! another wall
    real(kind=jprb), intent(out) :: view_ground_sky, view_wall_wall
    ! Cosine of solar zenith angle
    real(kind=jprb), intent(in),  optional :: cos_sza
    ! Fraction of direct radiation at the top-of-canopy that
    ! penetrates to the ground unscattered
    real(kind=jprb), intent(out), optional :: view_dir_ground

    ! Direct penetration factor, defined as (H/W)*tan(sza) - see Eq. 2
    real(kind=jprb) :: norm_x0
    ! Ratio of two variables in Eq. 16
    real(kind=jprb) :: y_over_w

    view_ground_sky = sqrt(height_width_ratio * height_width_ratio + 1.0_jprb) &
         &          - height_width_ratio
    view_wall_wall = sqrt(1.0_jprb / (height_width_ratio * height_width_ratio) + 1.0_jprb) &
         &         - 1.0_jprb / height_width_ratio

    if (present(cos_sza) .and. present(view_dir_ground)) then
      ! The sqrt term is tan(sza)
      norm_x0 = (Pi*0.5_jprb) * height_width_ratio * sqrt(1.0_jprb / (cos_sza*cos_sza) - 1.0_jprb)
      y_over_w = sqrt(max(norm_x0*norm_x0 - 1.0_jprb, 0.0_jprb))
      if (y_over_w > 0.0_jprb) then
        view_dir_ground = (2.0_jprb/Pi) * (y_over_w - norm_x0 + atan(1.0_jprb / y_over_w))
      else
        view_dir_ground = 1.0_jprb - 2.0_jprb*norm_x0/Pi
      end if
    end if

  end subroutine calc_view_factors_inf


  !---------------------------------------------------------------------
  ! Calculate view factors assuming the exponential model of urban
  ! geometry
  elemental subroutine calc_view_factors_exp(height_separation_ratio, &
       &  view_ground_sky, view_wall_wall, cos_sza, view_dir_ground)
    
    use parkind1, only : jpim, jprb
    use radiation_constants, only : Pi

    implicit none

    ! Parameters
    integer, parameter :: nw = 8 ! Number of weights
    real(kind=jprb), parameter :: weights(nw) &
         &  = [0.0506142681451884_jprb, 0.111190517226687_jprb, &
         &     0.156853322938944_jprb,  0.181341891689181_jprb, &
         &     0.181341891689181_jprb,  0.156853322938944_jprb, &
         &     0.111190517226687_jprb,  0.0506142681451884_jprb ]
    real(kind=jprb), parameter :: nodes(nw) &
         &  = [0.0198550717512319_jprb, 0.101666761293187_jprb, &
         &     0.237233795041836_jprb,  0.408282678752175_jprb, &
         &     0.591717321247825_jprb,  0.762766204958164_jprb, &
         &     0.898333238706813_jprb,  0.980144928248768_jprb ]

    ! Ratio of building height to the building horizontal separation
    ! scale in the exponential-urban-geometry assumption (H/X)
    real(kind=jprb), intent(in)  :: height_separation_ratio
    ! Fraction of diffuse radiation emanating from the ground that
    ! reaches the sky and fraction emanating from a wall that reaches
    ! another wall
    real(kind=jprb), intent(out) :: view_ground_sky, view_wall_wall
    ! Cosine of solar zenith angle
    real(kind=jprb), intent(in),  optional :: cos_sza
    ! Fraction of direct radiation at the top-of-canopy that
    ! penetrates to the ground unscattered
    real(kind=jprb), intent(out), optional :: view_dir_ground

    ! Direct penetration factor, defined as (H/X)*tan(sza) - see Eq. 2
    real(kind=jprb) :: norm_x0

    real(kind=jprb) :: tk(nw), exp_tk(nw)

    real(kind=jprb) :: hweight(nw), vweight(nw)

    ! Eqs. 12 and 17
    hweight = weights * nodes / sum(weights * nodes)
    vweight = weights * sqrt(1.0_jprb - nodes*nodes)
    vweight = vweight / sum(vweight)

    ! The sqrt term is tan(zenith angle) since node is the cosine of
    ! the zenith angle
    tk = height_separation_ratio * sqrt(1.0_jprb / (nodes*nodes) - 1.0_jprb)
    exp_tk = exp(-tk)

    ! Eqs. 41 and 42 (note that there is a mistake in the latter
    ! equation in the paper)
    view_ground_sky = sum(hweight * exp_tk)
    view_wall_wall  = 1.0_jprb - sum(vweight * (1.0_jprb-exp_tk) / tk) 

    if (present(cos_sza) .and. present(view_dir_ground)) then
      ! The sqrt term is tan(sza)
      norm_x0 = height_separation_ratio * sqrt(1.0_jprb / (cos_sza*cos_sza) - 1.0_jprb)
      view_dir_ground = exp(-norm_x0)
    end if

  end subroutine calc_view_factors_exp

end module radsurf_view_factor
