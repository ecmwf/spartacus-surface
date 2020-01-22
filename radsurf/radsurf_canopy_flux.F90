! radsurf_canopy_flux.f90 - Derived type for fluxes within the vegetated/urban canopy
!
! Copyright (C) 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radsurf_canopy_flux

  use parkind1, only : jpim, jprb

  implicit none

  !---------------------------------------------------------------------
  ! Derived type storing spectral fluxes within the canopy in either
  ! the shortwave or longwave.  All quantities are in W m-2,
  ! which is radiative power per unit horizontal area of the entire
  ! domain. Net fluxes are in-minus-out of the facet.
  type canopy_flux_type

    ! Flux components are dn/in for fluxes into a surface, and net for
    ! fluxes into minus out of a surface
    real(kind=jprb), allocatable :: ground_dn(:,:)        ! (nspec,ncol)
    real(kind=jprb), allocatable :: ground_dn_direct(:,:) ! (nspec,ncol)
    real(kind=jprb), allocatable :: ground_net(:,:)       ! (nspec,ncol)
    real(kind=jprb), allocatable :: roof_in(:,:)          ! (nspec,ntotlay)
    real(kind=jprb), allocatable :: roof_net(:,:)         ! (nspec,ntotlay)
    real(kind=jprb), allocatable :: wall_in(:,:)          ! (nspec,ntotlay)
    real(kind=jprb), allocatable :: wall_net(:,:)         ! (nspec,ntotlay)

    ! Absorption by the clear-air region, the vegetation in the
    ! vegetated region, and the air in the vegetated region
    real(kind=jprb), allocatable :: clear_air_abs(:,:)    ! (nspec,ntotlay)
    real(kind=jprb), allocatable :: veg_abs(:,:)          ! (nspec,ntotlay)
    real(kind=jprb), allocatable :: veg_air_abs(:,:)      ! (nspec,ntotlay)

    ! Number of spectral intervals, number of columns total number of
    ! layers across all columns
    integer(kind=jpim) :: nspec, ncol, ntotlay

  contains

    procedure :: allocate   => allocate_canopy_flux
    procedure :: deallocate => deallocate_canopy_flux
    procedure :: scale      => scale_canopy_flux
    procedure :: zero       => zero_canopy_flux

  end type canopy_flux_type

contains

  subroutine allocate_canopy_flux(this, ncol, ntotlay, nspec, use_direct)

    class(canopy_flux_type), intent(inout) :: this
    integer(kind=jpim),      intent(in)    :: ncol, ntotlay, nspec
    logical, optional,       intent(in)    :: use_direct

    logical :: use_direct_local

    if (present(use_direct)) then
      use_direct_local = use_direct
    else
      use_direct_local = .true.
    end if
    
    call this%deallocate()

    allocate(this%ground_dn(nspec,ncol))
    if (use_direct_local) then
      allocate(this%ground_dn_direct(nspec,ncol))
    end if
    allocate(this%ground_net(nspec,ncol))
    allocate(this%roof_in(nspec,ntotlay))
    allocate(this%roof_net(nspec,ntotlay))
    allocate(this%wall_in(nspec,ntotlay))
    allocate(this%wall_net(nspec,ntotlay))
    allocate(this%clear_air_abs(nspec,ntotlay))
    allocate(this%veg_abs(nspec,ntotlay))
    allocate(this%veg_air_abs(nspec,ntotlay))

  end subroutine allocate_canopy_flux

  subroutine deallocate_canopy_flux(this)

    class(canopy_flux_type), intent(inout) :: this

    if (allocated(this%ground_dn))        deallocate(this%ground_dn)
    if (allocated(this%ground_dn_direct)) deallocate(this%ground_dn_direct)
    if (allocated(this%ground_net))       deallocate(this%ground_net)
    if (allocated(this%roof_in))          deallocate(this%roof_in)
    if (allocated(this%roof_net))         deallocate(this%roof_net)
    if (allocated(this%wall_in))          deallocate(this%wall_in)
    if (allocated(this%wall_net))         deallocate(this%wall_net)
    if (allocated(this%clear_air_abs))    deallocate(this%clear_air_abs)
    if (allocated(this%veg_abs))          deallocate(this%veg_abs)
    if (allocated(this%veg_air_abs))      deallocate(this%veg_air_abs)

  end subroutine deallocate_canopy_flux
  
  ! Typically the canopy_flux object initially contains normalized
  ! fluxes, e.g. for a top-of-canopy downwelling flux of unity.  Once
  ! the top-of-canopy downwelling flux is known, the canopy fluxes can
  ! be scaled.
  subroutine scale_canopy_flux(this, factor)
   
    class(canopy_flux_type), intent(inout) :: this
    real(kind=jprb),         intent(in)    :: factor

    this%ground_dn        = factor * this%ground_dn
    this%ground_net       = factor * this%ground_net
    if (allocated(this%ground_dn_direct)) then
      this%ground_dn_direct = factor * this%ground_dn_direct
    end if
    if (allocated(this%roof_in)) then
      this%roof_in  = factor * this%roof_in
      this%roof_net = factor * this%roof_net
      this%wall_in  = factor * this%wall_in
      this%wall_net = factor * this%wall_net
    end if
    if (allocated(this%clear_air_abs)) then
      this%clear_air_abs = factor * this%clear_air_abs
    end if
    if (allocated(this%veg_abs)) then
      this%veg_abs     = factor * this%veg_abs
      this%veg_air_abs = factor * this%veg_air_abs
    end if

  end subroutine scale_canopy_flux

  ! Set the fluxes to zero for a particular column
  subroutine zero_canopy_flux(this, icol, ilay1, ilay2)
   
    class(canopy_flux_type),      intent(inout) :: this
    integer(kind=jpim),           intent(in)    :: icol
    integer(kind=jpim), optional, intent(in)    :: ilay1, ilay2

    this%ground_dn(:,icol)        = 0.0_jprb
    this%ground_net(:,icol)       = 0.0_jprb
    if (allocated(this%ground_dn_direct)) then
      this%ground_dn_direct(:,icol) = 0.0_jprb
    end if
    if (present(ilay1) .and. present(ilay2)) then
      if (ilay2 >= ilay1) then
        if (allocated(this%roof_in)) then
          this%roof_in(:,ilay1:ilay2)  = 0.0_jprb
          this%roof_net(:,ilay1:ilay2) = 0.0_jprb
          this%wall_in(:,ilay1:ilay2)  = 0.0_jprb
          this%wall_net(:,ilay1:ilay2) = 0.0_jprb
        end if
        if (allocated(this%clear_air_abs)) then
          this%clear_air_abs(:,ilay1:ilay2) = 0.0_jprb
        end if
        if (allocated(this%veg_abs)) then
          this%veg_abs(:,ilay1:ilay2)     = 0.0_jprb
          this%veg_air_abs(:,ilay1:ilay2) = 0.0_jprb
        end if
      end if
    end if

  end subroutine zero_canopy_flux

end module radsurf_canopy_flux
