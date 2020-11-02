! radsurf_sw_spectral_properties.F90 - Derived type for shortwave spectral properties
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

module radsurf_sw_spectral_properties

  use parkind1, only : jpim, jprb

  implicit none

  !---------------------------------------------------------------------
  ! Derived type storing shortwave spectral properties of the canopy
  type sw_spectral_properties_type

    ! Volume properties
    real(kind=jprb), allocatable :: air_ext(:,:) ! m-1 (nspec,ntotlay)
    real(kind=jprb), allocatable :: air_ssa(:,:) !     (nspec,ntotlay)
    real(kind=jprb), allocatable :: veg_ssa(:,:) !     (nspec,ntotlay)

    ! Facet properties
    real(kind=jprb), allocatable :: ground_albedo(:,:) ! (nspec,ncol)
    real(kind=jprb), allocatable :: roof_albedo(:,:)   ! (nspec,ntotlay)
    real(kind=jprb), allocatable :: wall_albedo(:,:)   ! (nspec,ntotlay)

    real(kind=jprb), allocatable :: ground_albedo_dir(:,:)  ! (nspec,ncol)
    real(kind=jprb), allocatable :: roof_albedo_dir(:,:)    ! (nspec,ntotlay)
    real(kind=jprb), allocatable :: wall_specular_frac(:,:) ! (nspec,ntotlay)

    ! Number of spectral intervals, number of columns and total number
    ! of layers
    integer(kind=jpim) :: nspec, ncol, ntotlay

  contains

    procedure :: allocate   => allocate_spectral
    procedure :: deallocate => deallocate_spectral
    
  end type sw_spectral_properties_type

contains

  
  !---------------------------------------------------------------------
  subroutine allocate_spectral(this, config, ncol, ntotlay, nspec, &
       &                       i_representation)

    use yomhook,        only : lhook, dr_hook
    use radsurf_config, only : config_type

    use radsurf_canopy_properties, only : ITileFlat,  ITileForest, &
         &                                ITileUrban, ITileVegetatedUrban

    class(sw_spectral_properties_type), intent(inout) :: this
    type(config_type),                  intent(in)    :: config
    integer(kind=jpim),                 intent(in)    :: nspec, ncol, ntotlay
    integer(kind=jpim), optional,       intent(in)    :: i_representation(:)

    logical :: do_vegetation, do_urban, do_canopy

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_sw_spectral_properties:allocate',0,hook_handle)

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
      !
      ! FIXME: allocate veg_ssa in do_canopy case
      !
      allocate(this%veg_ssa(nspec,ntotlay))
      !
    end if
    !
    !if (do_vegetation) then
    !  allocate(this%veg_ssa(nspec,ntotlay))
    !end if
    !
    allocate(this%ground_albedo(nspec,ncol))
    if (config%use_sw_direct_albedo) then
      allocate(this%ground_albedo_dir(nspec,ncol))
    end if
    if (do_urban) then
      allocate(this%roof_albedo(nspec,ntotlay))
      allocate(this%wall_albedo(nspec,ntotlay))
      allocate(this%wall_specular_frac(nspec,ntotlay))
      if (config%use_sw_direct_albedo) then
        allocate(this%roof_albedo_dir(nspec,ntotlay))
      end if
    end if
    
    if (lhook) call dr_hook('radiation_sw_spectral_properties:allocate',1,hook_handle)
    
  end subroutine allocate_spectral


  

  !---------------------------------------------------------------------
  subroutine deallocate_spectral(this)

    use yomhook, only : lhook, dr_hook

    class(sw_spectral_properties_type), intent(inout) :: this

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_sw_spectral_properties:deallocate',0,hook_handle)

    if (allocated(this%air_ext))           deallocate(this%air_ext)
    if (allocated(this%air_ssa))           deallocate(this%air_ssa)
    if (allocated(this%veg_ssa))           deallocate(this%veg_ssa)
    if (allocated(this%ground_albedo))     deallocate(this%ground_albedo)
    if (allocated(this%roof_albedo))       deallocate(this%roof_albedo)
    if (allocated(this%wall_albedo))       deallocate(this%wall_albedo)
    if (allocated(this%ground_albedo_dir)) deallocate(this%ground_albedo_dir)
    if (allocated(this%roof_albedo_dir))   deallocate(this%roof_albedo_dir)
    if (allocated(this%wall_specular_frac))deallocate(this%wall_specular_Frac)

    if (lhook) call dr_hook('radiation_sw_spectral_properties:deallocate',1,hook_handle)

  end subroutine deallocate_spectral

  
end module radsurf_sw_spectral_properties
  
