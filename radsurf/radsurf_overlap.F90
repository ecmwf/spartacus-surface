! radsurf_overlap.F90 - Compute overlap matrices for vegetation
!
! Copyright (C) 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radsurf_overlap

contains

  pure function calc_overlap_matrix_max_ran(nreg, &
       &  f_upper, f_lower) result(overlap_matrix)

    use parkind1, only : jprb

    implicit none
    
    integer, intent(in) :: nreg ! Number of regions

    ! Vegetation fraction in the upper and lower layers
    real(jprb), intent(in) :: f_upper, f_lower

    ! Output overlap matrix
    real(jprb) :: overlap_matrix(nreg,nreg)

    ! Combined vegetation cover of pair of layers
    real(jprb) :: pair_cover

    real(jprb) :: veg_unit

    pair_cover = max(f_upper,f_lower)

    ! Clear in both layers
    overlap_matrix(1,1) = 1.0_jprb - pair_cover
    if (nreg == 2) then
      ! Clear in upper layer, vegetation in lower layer
      overlap_matrix(1,2) = pair_cover - f_upper
      ! Clear in lower layer, vegetation in upper layer
      overlap_matrix(2,1) = pair_cover - f_lower
      ! Vegetation in both layers
      overlap_matrix(2,2) = f_upper + f_lower - pair_cover
    else
      ! The following assumes that the two vegetated regions are of
      ! equal area.
      ! Clear in upper layer, vegetated in lower layer
      overlap_matrix(1,2) = 0.5_jprb * (pair_cover - f_upper)
      overlap_matrix(1,3) = overlap_matrix(1,2)
      ! Clear in lower layer, vegetated in upper layer
      overlap_matrix(2,1) = 0.5_jprb * (pair_cover - f_lower)
      overlap_matrix(3,1) = overlap_matrix(2,1)
      ! Vegetated in both layers
      overlap_matrix(2,2) = 0.5_jprb * (f_upper + f_lower - pair_cover)
      overlap_matrix(2,3) = overlap_matrix(2,2)
      overlap_matrix(3,3) = overlap_matrix(2,2)
      overlap_matrix(3,2) = overlap_matrix(2,2)
    end if

  end function calc_overlap_matrix_max_ran


  !---------------------------------------------------------------------
  ! Compute the upward and downward overlap matrices u_overlap and
  ! v_overlap, respectively, where u_overlap is defined such that
  ! y=u_overlap*x, where x is a vector of upwelling fluxes in each
  ! region just below an interface, and y is a vector of upwelling
  ! fluxes in each region just above that interface. For nlay model
  ! layers there are nlay+1 interfaces including the ground and
  ! top-of-atmosphere, and so that is one of the dimensions of
  ! u_overlap and v_overlap.
  subroutine calc_overlap_matrices(nlay,nreg, &
       &     region_fracs, u_overlap, v_overlap, &
       &     frac_threshold)

    use parkind1,     only : jprb
    use yomhook,      only : lhook, dr_hook

    implicit none
    
    ! Number of layers and regions
    integer,  intent(in) :: nlay, nreg

    ! Area fraction of each region: region 1 is clear sky, and 2+ are
    ! the vegetated regions (only one or two vegetated regions are
    ! supported)
    real(jprb), intent(in), dimension(1:nreg,nlay)  :: region_fracs

    ! Output overlap matrices, one for each layer interface
    real(jprb), intent(out), dimension(nreg,nreg,nlay+1) &
         &  :: u_overlap, v_overlap

    ! Regions smaller than this are ignored
    real(jprb), intent(in), optional :: frac_threshold

    ! Loop indices for column, layer, region and the regions in the
    ! upper and lower layers for an interface
    integer  :: jcol, jlay, jupper, jlower

    ! Overlap matrix (non-directional)
    real(jprb) :: overlap_matrix(nreg,nreg)

    ! Fraction of the gridbox occupied by each region in the upper and
    ! lower layers for an interface
    real(jprb) :: frac_upper(nreg), frac_lower(nreg)

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radsurf_overlap:calc_overlap_matrices',0,hook_handle)

    ! The ground is treated as one clear-sky region, so the fractions
    ! are assigned as such
    frac_lower(1)      = 1.0_jprb
    frac_lower(2:nreg) = 0.0_jprb

    ! Loop up through the canopy, where jlay indexes each layer
    ! interface starting at 1 for the ground, as well as indexing each
    ! layer starting at 1 for the lowest layer.
    do jlay = 1,nlay+1
      ! Fraction of each region just below the interface
      if (jlay > nlay) then
        ! We are at the top of the canopy: treat as a single clear-sky
        ! region
        frac_upper(1)      = 1.0_jprb
        frac_upper(2:nreg) = 0.0_jprb
      else
        frac_upper = region_fracs(1:nreg,jlay)
      end if
   
      ! Simpler scheme assuming the two cloudy regions have the
      ! same fraction

      overlap_matrix = calc_overlap_matrix_max_ran(nreg, &
           1.0_jprb - frac_upper(1), 1.0_jprb - frac_lower(1))

      ! Convert to directional overlap matrices
      do jupper = 1,nreg
        do jlower = 1,nreg
          if (frac_lower(jlower) >= frac_threshold) then
            u_overlap(jupper,jlower,jlay) = overlap_matrix(jupper,jlower) &
                 &  / frac_lower(jlower)
          else
            u_overlap(jupper,jlower,jlay) = 0.0_jprb
          end if
          if (frac_upper(jupper) >= frac_threshold) then
            v_overlap(jlower,jupper,jlay) = overlap_matrix(jupper,jlower) &
                 &  / frac_upper(jupper)
          else
            v_overlap(jlower,jupper,jlay) = 0.0_jprb
          end if
        end do
      end do
      frac_lower = frac_upper
        
    end do ! layers

    if (lhook) call dr_hook('radsurf_overlap:calc_overlap_matrices',1,hook_handle)

  end subroutine calc_overlap_matrices
  
  
end module radsurf_overlap
