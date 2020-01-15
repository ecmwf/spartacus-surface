! radsurf_facet_properties.f90 - Derived type for facet spectral properties
!
! Copyright (C) 2019 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radsurf_facet_properties

  use parkind1, only : jpim, jprb

  implicit none

  !---------------------------------------------------------------------
  ! Derived type storing spectral properties of the facets of the surface
  type facet_properties_type

    ! Shortwave spectral properties
    real(kind=jprb), allocatable :: ground_sw_albedo(:,:) ! (nsw,ncol)
    real(kind=jprb), allocatable :: ground_sw_albedo_direct(:,:) ! (nsw,ncol)
    real(kind=jprb), allocatable :: roof_sw_albedo(:,:) ! (nsw,ntotlay)
    real(kind=jprb), allocatable :: roof_sw_albedo_direct(:,:) ! (nsw,ntotlay)
    real(kind=jprb), allocatable :: wall_sw_albedo(:,:) ! (nsw,ntotlay)
    real(kind=jprb), allocatable :: wall_specular_fraction(:,:) ! (nsw,ntotlay)

    ! Longwave spectral properties: emissivity...
    real(kind=jprb), allocatable :: ground_lw_emissivity(:,:) ! (nlw,ncol)
    real(kind=jprb), allocatable :: roof_lw_emissivity(:,:) ! (nlw,ntotlay)
    real(kind=jprb), allocatable :: wall_lw_emissivity(:,:) ! (nlw,ntotlay)
    ! ...and outward emission (W m-2)
    real(kind=jprb), allocatable :: ground_lw_emission(:,:) ! (nlw,ncol)
    real(kind=jprb), allocatable :: roof_lw_emission(:,:) ! (nlw,ntotlay)
    real(kind=jprb), allocatable :: wall_lw_emission(:,:) ! (nlw,ntotlay)

    ! Number of columns and maximum number of layers
    integer(kind=jpim) :: ncol, ntotlay

    ! Number of shortwave albedo and longwave emissivity bands
    integer(kind=jpim) :: nsw, nlw

  contains
    procedure :: allocate   => allocate_facet
    procedure :: deallocate => deallocate_facet
    procedure :: calc_monochromatic_emission
    !procedure :: read => read_from_netcdf
    
  end type facet_properties_type


contains

  !---------------------------------------------------------------------
  subroutine allocate_facet(this, config, ncol, ntotlay, &
       &                    i_representation)

    use yomhook,        only : lhook, dr_hook
    use radsurf_config, only : config_type

    use radsurf_canopy_properties, only : ITileFlat,  ITileForest, &
         &                                ITileUrban, ITileVegetatedUrban

    class(facet_properties_type), intent(inout) :: this
    type(config_type),   intent(in)    :: config
    integer(kind=jpim),  intent(in)    :: ncol, ntotlay
    integer(kind=jpim),  intent(in), optional :: i_representation(:)

    logical :: do_vegetation, do_urban, do_canopy

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_facet_properties:allocate',0,hook_handle)

    call this%deallocate()

    this%ncol    = ncol
    this%ntotlay = ntotlay
    do_vegetation = config%do_vegetation
    do_urban      = config%do_urban
    this%nsw      = config%nsw
    this%nlw      = config%nlw

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

    if (config%do_sw) then
      allocate(this%ground_sw_albedo(this%nsw,ncol))
      if (do_urban) then
        allocate(this%roof_sw_albedo(this%nsw,ntotlay))
        allocate(this%wall_sw_albedo(this%nsw,ntotlay))
        allocate(this%wall_specular_fraction(this%nsw,ntotlay))
        this%wall_specular_fraction = 0.0_jprb
      end if
      if (config%use_sw_direct_albedo) then
        allocate(this%ground_sw_albedo_direct(this%nsw,ncol))
        if (do_urban) then
          allocate(this%roof_sw_albedo_direct(this%nsw,ntotlay))
        end if
      end if
    end if

    if (config%do_lw) then
      allocate(this%ground_lw_emissivity(this%nlw,ncol))
      allocate(this%ground_lw_emission(this%nlw,ncol))
      if (do_urban) then
        allocate(this%roof_lw_emissivity(this%nlw,ntotlay))
        allocate(this%wall_lw_emissivity(this%nlw,ntotlay))
        allocate(this%roof_lw_emission(this%nlw,ntotlay))
        allocate(this%wall_lw_emission(this%nlw,ntotlay))
      end if
    end if

    if (lhook) call dr_hook('radiation_facet_properties:allocate',1,hook_handle)

  end subroutine allocate_facet


  !---------------------------------------------------------------------
  subroutine deallocate_facet(this)

    use yomhook, only : lhook, dr_hook

    class(facet_properties_type), intent(inout) :: this

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_facet_properties:deallocate',0,hook_handle)

    if (allocated(this%ground_sw_albedo))        deallocate(this%ground_sw_albedo)
    if (allocated(this%ground_sw_albedo_direct)) deallocate(this%ground_sw_albedo_direct)
    if (allocated(this%roof_sw_albedo))          deallocate(this%roof_sw_albedo)
    if (allocated(this%roof_sw_albedo_direct))   deallocate(this%roof_sw_albedo_direct)
    if (allocated(this%wall_sw_albedo))          deallocate(this%wall_sw_albedo)
    if (allocated(this%wall_specular_fraction))  deallocate(this%wall_specular_fraction)
    if (allocated(this%ground_lw_emissivity))    deallocate(this%ground_lw_emissivity)
    if (allocated(this%roof_lw_emissivity))      deallocate(this%roof_lw_emissivity)
    if (allocated(this%wall_lw_emissivity))      deallocate(this%wall_lw_emissivity)
    if (allocated(this%ground_lw_emission))      deallocate(this%ground_lw_emission)
    if (allocated(this%roof_lw_emission))        deallocate(this%roof_lw_emission)
    if (allocated(this%wall_lw_emission))        deallocate(this%wall_lw_emission)

    if (lhook) call dr_hook('radiation_facet_properties:deallocate',1,hook_handle)

  end subroutine deallocate_facet

  
  !---------------------------------------------------------------------
  subroutine calc_monochromatic_emission(this, canopy_props)

    use radiation_constants, only : StefanBoltzmann
    use radsurf_canopy_properties, only : canopy_properties_type
    
    class(facet_properties_type), intent(inout) :: this
    type(canopy_properties_type), intent(in)    :: canopy_props
    
    integer :: jlay, jcol, ntotlay, ncol

    ntotlay = canopy_props%ntotlay
    ncol    = canopy_props%ncol

    if (allocated(canopy_props%ground_temperature) &
         &  .and. allocated(this%ground_lw_emissivity)) then
      if (.not. allocated(this%ground_lw_emission)) then
        allocate(this%ground_lw_emission(1,ncol))
      end if
      this%ground_lw_emission(1,:) = StefanBoltzmann &
           &  * this%ground_lw_emissivity(1,:) &
           &  * canopy_props%ground_temperature ** 4
    end if
    if (allocated(canopy_props%roof_temperature) &
         &  .and. allocated(this%roof_lw_emissivity)) then
      if (.not. allocated(this%roof_lw_emission)) then
        allocate(this%roof_lw_emission(1,ntotlay))
      end if
      this%roof_lw_emission(1,:) = StefanBoltzmann &
           &  * this%roof_lw_emissivity(1,:) &
           &  * canopy_props%roof_temperature ** 4
      if (.not. allocated(this%wall_lw_emission)) then
        allocate(this%wall_lw_emission(1,ntotlay))
      end if
      this%wall_lw_emission(1,:) = StefanBoltzmann &
           &  * this%wall_lw_emissivity(1,:) &
           &  * canopy_props%wall_temperature ** 4
    end if

  end subroutine calc_monochromatic_emission
  
end module radsurf_facet_properties
