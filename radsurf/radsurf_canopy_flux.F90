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
    procedure :: sum        => sum_canopy_flux

  end type canopy_flux_type

contains

  !---------------------------------------------------------------------
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

    this%nspec   = nspec
    this%ncol    = ncol
    this%ntotlay = ntotlay

  end subroutine allocate_canopy_flux

  !---------------------------------------------------------------------
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
  
  !---------------------------------------------------------------------
  ! Typically the canopy_flux object initially contains normalized
  ! fluxes, e.g. for a top-of-canopy downwelling flux of unity.  Once
  ! the top-of-canopy downwelling flux is known, the canopy fluxes can
  ! be scaled.
  subroutine scale_canopy_flux(this, nlay, factor)

    use radiation_io, only : nulerr, radiation_abort

    class(canopy_flux_type), intent(inout) :: this
    integer(kind=jpim),      intent(in)    :: nlay(:)     ! (ncol)
    real(kind=jprb),         intent(in)    :: factor(:,:) ! (nspec,ncol)

    integer(kind=jpim) :: ncol, istartcol, jcol

    ! Index from layer to column
    integer(kind=jpim) :: indcol(this%ntotlay)

    if (ubound(factor,1) /= this%nspec) then
      write(nulerr,'(a)') '*** Error: spectral resolution mismatch when scaling canopy fluxes'
      print *, ubound(factor,1), this%nspec
      call radiation_abort()
    end if

    ncol = this%ncol

    ! Assign indices from layer to column
    istartcol = 1
    do jcol = 1,ncol
      if (nlay(jcol) > 0) then
        indcol(istartcol:istartcol-1+nlay(jcol)) = jcol
        istartcol = istartcol + nlay(jcol)
      end if
    end do
    
    this%ground_dn        = factor * this%ground_dn
    this%ground_net       = factor * this%ground_net
    if (allocated(this%ground_dn_direct)) then
      this%ground_dn_direct = factor * this%ground_dn_direct
    end if
    if (allocated(this%roof_in)) then
      this%roof_in  = factor(:,indcol) * this%roof_in
      this%roof_net = factor(:,indcol) * this%roof_net
      this%wall_in  = factor(:,indcol) * this%wall_in
      this%wall_net = factor(:,indcol) * this%wall_net
    end if
    if (allocated(this%clear_air_abs)) then
      this%clear_air_abs = factor(:,indcol) * this%clear_air_abs
    end if
    if (allocated(this%veg_abs)) then
      this%veg_abs     = factor(:,indcol) * this%veg_abs
      this%veg_air_abs = factor(:,indcol) * this%veg_air_abs
    end if

  end subroutine scale_canopy_flux

  !---------------------------------------------------------------------
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

  !---------------------------------------------------------------------
  ! this = flux1 + flux2
  subroutine sum_canopy_flux(this, flux1, flux2)

    class(canopy_flux_type), intent(inout) :: this
    type(canopy_flux_type),  intent(in)    :: flux1, flux2

    logical :: use_direct

    use_direct = allocated(flux1%ground_dn_direct)

    if (.not. allocated(this%ground_dn)) then
      call this%allocate(flux1%ncol, flux1%ntotlay, flux1%nspec, use_direct=use_direct)
    end if

    this%ground_dn = flux1%ground_dn + flux2%ground_dn
    if (use_direct) then
      this%ground_dn_direct = flux1%ground_dn_direct + flux2%ground_dn_direct
    end if
    this%ground_net = flux1%ground_net + flux2%ground_net
    this%roof_in = flux1%roof_in + flux2%roof_in
    this%roof_net = flux1%roof_net + flux2%roof_net
    this%wall_in = flux1%wall_in + flux2%wall_in
    this%wall_net = flux1%wall_net + flux2%wall_net
    this%clear_air_abs = flux1%clear_air_abs + flux2%clear_air_abs
    this%veg_abs = flux1%veg_abs + flux2%veg_abs
    this%veg_air_abs = flux1%veg_air_abs + flux2%veg_air_abs

  end subroutine sum_canopy_flux

end module radsurf_canopy_flux
