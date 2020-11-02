! radsurf_lw_spectral_properties.F90 - Derived type for longwave spectral properties
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

module radsurf_lw_spectral_properties

  use parkind1, only : jpim, jprb

  implicit none

  !---------------------------------------------------------------------
  ! Derived type storing longwave spectral properties of the canopy
  type lw_spectral_properties_type

    ! Volume properties
    real(kind=jprb), allocatable :: air_ext(:,:) ! m-1 (nspec,ntotlay)
    real(kind=jprb), allocatable :: air_ssa(:,:) !     (nspec,ntotlay)
    real(kind=jprb), allocatable :: veg_ssa(:,:) !     (nspec,ntotlay)

    ! Planck function at the temperature of the clear-air, leaves and
    ! the air within vegetation (W m-2)
    real(kind=jprb), allocatable :: clear_air_planck(:,:) !  (nspec,ntotlay)
    real(kind=jprb), allocatable :: veg_planck(:,:)       !  (nspec,ntotlay)
    real(kind=jprb), allocatable :: veg_air_planck(:,:)   !  (nspec,ntotlay)

    ! Facet properties
    real(kind=jprb), allocatable :: ground_emissivity(:,:) ! (nspec,ncol)
    real(kind=jprb), allocatable :: roof_emissivity(:,:)   ! (nspec,ntotlay)
    real(kind=jprb), allocatable :: wall_emissivity(:,:)   ! (nspec,ntotlay)

    ! Emission 
    real(kind=jprb), allocatable :: ground_emission(:,:) ! (nspec,ncol)
    real(kind=jprb), allocatable :: roof_emission(:,:)   ! (nspec,ntotlay)
    real(kind=jprb), allocatable :: wall_emission(:,:)   ! (nspec,ntotlay)

    ! Number of spectral intervals, number of columns and total number
    ! of layers
    integer(kind=jpim) :: nspec, ncol, ntotlay

  contains
    
    procedure :: allocate   => allocate_spectral
    procedure :: deallocate => deallocate_spectral
    procedure :: calc_monochromatic_emission

  end type lw_spectral_properties_type

contains

  !---------------------------------------------------------------------
  subroutine allocate_spectral(this, config, nspec, ncol, ntotlay, &
       &                       i_representation)

    use yomhook,        only : lhook, dr_hook
    use radsurf_config, only : config_type

    use radsurf_canopy_properties, only : ITileFlat,  ITileForest, &
         &                                ITileUrban, ITileVegetatedUrban

    class(lw_spectral_properties_type), intent(inout) :: this
    type(config_type),                  intent(in)    :: config
    integer(kind=jpim),                 intent(in)    :: nspec, ncol, ntotlay
    integer(kind=jpim), optional,       intent(in)    :: i_representation(:)

    logical :: do_vegetation, do_urban, do_canopy

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_lw_spectral_properties:allocate',0,hook_handle)

    call this%deallocate()

    this%ncol    = ncol
    this%ntotlay = ntotlay
    this%nspec   = nspec
    do_vegetation = config%do_vegetation
    do_urban      = config%do_urban

    if (present(i_representation)) then
      if (.not. any(i_representation == ITileForest &
           &        .or. i_representation == ITileVegetatedUrban)) then
        do_vegetation = .false.
      end if
      if (.not. any(i_representation == ITileUrban &
           &        .or. i_representation == ITileVegetatedUrban)) then
        do_urban = .false.
      end if
    end if

    do_canopy = (do_vegetation .or. do_urban)

    if (do_canopy) then
      allocate(this%air_ext(nspec,ntotlay))
      allocate(this%air_ssa(nspec,ntotlay))
      allocate(this%clear_air_planck(nspec,ntotlay))
      !
      ! FIXME: Allocate veg_ssa, veg_air_planck, veg_planck
      ! in both cases since they are a lot used in urban routines
      !
      allocate(this%veg_ssa(nspec,ntotlay))
      allocate(this%veg_air_planck(nspec,ntotlay))
      allocate(this%veg_planck(nspec,ntotlay))
      !
    end if
    !
    allocate(this%ground_emissivity(nspec,ncol))
    allocate(this%ground_emission(nspec,ncol))
    if (do_urban) then
      allocate(this%roof_emissivity(nspec,ntotlay))
      allocate(this%wall_emissivity(nspec,ntotlay))
      allocate(this%roof_emission(nspec,ntotlay))
      allocate(this%wall_emission(nspec,ntotlay))
    end if
    !
    if (lhook) call dr_hook('radiation_lw_spectral_properties:allocate',1,hook_handle)
    !
  end subroutine allocate_spectral


  !---------------------------------------------------------------------
  subroutine deallocate_spectral(this)

    use yomhook, only : lhook, dr_hook

    class(lw_spectral_properties_type), intent(inout) :: this

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_lw_spectral_properties:deallocate',0,hook_handle)

    if (allocated(this%air_ext))           deallocate(this%air_ext)
    if (allocated(this%air_ssa))           deallocate(this%air_ssa)
    if (allocated(this%veg_ssa))           deallocate(this%veg_ssa)
    if (allocated(this%clear_air_planck))  deallocate(this%clear_air_planck)
    if (allocated(this%veg_air_planck))    deallocate(this%veg_air_planck)
    if (allocated(this%veg_planck))        deallocate(this%veg_planck)
    if (allocated(this%ground_emissivity)) deallocate(this%ground_emissivity)
    if (allocated(this%roof_emissivity))   deallocate(this%roof_emissivity)
    if (allocated(this%wall_emissivity))   deallocate(this%wall_emissivity)
    if (allocated(this%ground_emission))   deallocate(this%ground_emission)
    if (allocated(this%roof_emission))     deallocate(this%roof_emission)
    if (allocated(this%wall_emission))     deallocate(this%wall_emission)

    if (lhook) call dr_hook('radiation_lw_spectral_properties:deallocate',1,hook_handle)

  end subroutine deallocate_spectral

  
  !---------------------------------------------------------------------
  subroutine calc_monochromatic_emission(this, canopy_props)

    use radiation_constants,       only : StefanBoltzmann
    use radsurf_canopy_properties, only : canopy_properties_type
    
    class(lw_spectral_properties_type), intent(inout) :: this
    type(canopy_properties_type),       intent(in)    :: canopy_props
    
    integer :: ntotlay, ncol

    ntotlay = canopy_props%ntotlay
    ncol    = canopy_props%ncol

    if (allocated(canopy_props%ground_temperature) &
         &  .and. allocated(this%ground_emissivity)) then
      if (.not. allocated(this%ground_emission)) then
        allocate(this%ground_emission(1,ncol))
      end if
      this%ground_emission(1,:) = StefanBoltzmann &
           &  * this%ground_emissivity(1,:) &
           &  * canopy_props%ground_temperature ** 4
    end if
    if (allocated(canopy_props%roof_temperature) &
         &  .and. allocated(this%roof_emissivity)) then
      if (.not. allocated(this%roof_emission)) then
        allocate(this%roof_emission(1,ntotlay))
      end if
      this%roof_emission(1,:) = StefanBoltzmann &
           &  * this%roof_emissivity(1,:) &
           &  * canopy_props%roof_temperature ** 4
      if (.not. allocated(this%wall_emission)) then
        allocate(this%wall_emission(1,ntotlay))
      end if
      this%wall_emission(1,:) = StefanBoltzmann &
           &  * this%wall_emissivity(1,:) &
           &  * canopy_props%wall_temperature ** 4
    end if

  end subroutine calc_monochromatic_emission
  
end module radsurf_lw_spectral_properties
  
