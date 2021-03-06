! radsurf_boundary_conds_out.F90 - Boundary conditions presented to atmosphere above
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

module radsurf_boundary_conds_out

  use parkind1, only : jprb
  
  implicit none

  ! Spectral properties presented to the atmosphere above in each
  ! spectral band
  type boundary_conds_out_type

    ! Shortwave albedo to diffuse and direct radiation
    real(kind=jprb), allocatable :: sw_albedo(:,:)     ! (nsw,ncol)
    real(kind=jprb), allocatable :: sw_albedo_dir(:,:) ! (nsw,ncol)

    ! Longwave emissivity, and the upware emission in W m-2
    real(kind=jprb), allocatable :: lw_emissivity(:,:)    ! (nlw,ncol)
    real(kind=jprb), allocatable :: lw_emission(:,:)      ! (nlw,ncol)

  contains

    procedure :: allocate   => allocate_boundary_conds_out
    procedure :: deallocate => deallocate_boundary_conds_out
    
  end type boundary_conds_out_type

contains

  subroutine allocate_boundary_conds_out(this, ncol, nsw, nlw)

    class(boundary_conds_out_type), intent(inout) :: this
    integer, intent(in) :: nsw, nlw, ncol

    call this%deallocate()

    if (nsw > 0) then
      allocate(this%sw_albedo(nsw,ncol))
      allocate(this%sw_albedo_dir(nsw,ncol))
    end if
    if (nlw > 0) then
      allocate(this%lw_emissivity(nlw,ncol))
      allocate(this%lw_emission(nlw,ncol))
    end if

  end subroutine allocate_boundary_conds_out

  subroutine deallocate_boundary_conds_out(this)

    class(boundary_conds_out_type), intent(inout) :: this

    if (allocated(this%sw_albedo))     deallocate(this%sw_albedo)
    if (allocated(this%sw_albedo_dir)) deallocate(this%sw_albedo_dir)
    if (allocated(this%lw_emissivity)) deallocate(this%lw_emissivity)
    if (allocated(this%lw_emission))   deallocate(this%lw_emission)

  end subroutine deallocate_boundary_conds_out
    
end module radsurf_boundary_conds_out
