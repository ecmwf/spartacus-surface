! radsurf_unpacked_canopy_properties.f90 - Derived type for general canopy properties
!
! Copyright (C) 2019 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radsurf_unpacked_canopy_properties

  use parkind1, only : jpim, jprb

  implicit none

  !---------------------------------------------------------------------
  ! Derived type storing a physical description of the properties of
  ! the surface tiles
  type unpacked_canopy_properties_type

    ! Number of layers for each column
    integer(kind=jpim), allocatable :: nlay(:) ! (ncol)

    ! Layer thickness (m)
    real(kind=jprb), allocatable :: dz(:,:) ! (nmaxlay,ncol)

    ! Cosine of solar zenith angle
    real(kind=jprb), allocatable :: cos_sza(:) ! (ncol)

    ! Skin temperature of various surfaces (K)
    real(kind=jprb), allocatable :: ground_temperature(:) ! (ncol)
    real(kind=jprb), allocatable :: roof_temperature(:,:) ! (nmaxlay,ncol)
    real(kind=jprb), allocatable :: wall_temperature(:,:) ! (nmaxlay,ncol)

    ! Air temperature in canopy (K)
    real(kind=jprb), allocatable :: air_temperature(:,:) ! (nmaxlay,ncol)
    real(kind=jprb), allocatable :: vegetation_temperature(:,:) ! (nmaxlay,ncol)

    ! Fractional coverage of buildings and vegetation
    real(kind=jprb), allocatable :: building_fraction(:,:) ! (nmaxlay,ncol)
    real(kind=jprb), allocatable :: vegetation_fraction(:,:) ! (nmaxlay,ncol)

    ! Horizontal scale of buildings and vegetation (m)
    real(kind=jprb), allocatable :: building_scale(:,:) ! (nmaxlay,ncol)
    real(kind=jprb), allocatable :: vegetation_scale(:,:) ! (nmaxlay,ncol)

    ! Fractional standard deviation of vegetation optical depth
    real(kind=jprb), allocatable :: vegetation_fsd(:,:) ! (nmaxlay,ncol)

    ! Fraction of vegetation edge in contact with buildings rather
    ! than air
    real(kind=jprb), allocatable :: vegetation_contact_fraction(:,:) ! (nmaxlay,ncol)

    ! Representation codes (ITileFlat etc) for each tile: dimensioning
    ! is (ncol)
    integer(kind=jpim), allocatable :: i_representation(:)

    ! Number of columns and maximum number of layers
    integer(kind=jpim) :: ncol, nmaxlay

  contains
    procedure :: allocate   => allocate_canopy
    procedure :: deallocate => deallocate_canopy
!    procedure :: read => read_from_netcdf

  end type unpacked_canopy_properties_type


contains

  !---------------------------------------------------------------------
  subroutine allocate_canopy(this, config, ncol, nmaxlay, &
       &                    i_representation)

    use yomhook,        only : lhook, dr_hook
    use radsurf_config, only : config_type
    use radsurf_canopy_properties,  only : ITileFlat,  ITileVegetation, &
         &                                 ITileUrban, ITileVegetatedUrban, &
         &                                 canopy_properties_type

    class(canopy_properties_type), intent(inout) :: this
    type(config_type),   intent(in)    :: config
    integer(kind=jpim),  intent(in)    :: ncol, nmaxlay
    integer(kind=jpim),  intent(in), optional :: i_representation(:)

    logical :: do_vegetation, do_urban, do_canopy

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_canopy_properties:allocate',0,hook_handle)

    call this%deallocate()

    this%ncol     = ncol
    this%nmaxlay  = nmaxlay
    do_vegetation = config%do_vegetation
    do_urban      = config%do_urban

    if (present(i_representation)) then
      if (.not. any(i_representation == ITileVegetation &
           &        .or. i_representation == ITileVegetatedUrban)) then
        do_vegetation = .false.
      end if
      if (.not. any(i_representation == ITileUrban &
           &        .or. i_representation == ITileVegetatedUrban)) then
        do_urban = .false.
      end if
    end if

    do_canopy = (do_vegetation .or. do_urban)

    allocate(this%cos_sza(ncol))
    allocate(this%nlay(ncol))
    if (do_canopy) then
      allocate(this%dz(nmaxlay,ncol))
    end if
    allocate(this%ground_temperature(ncol))
    if (do_urban) then
      allocate(this%roof_temperature(nmaxlay,ncol))
      allocate(this%wall_temperature(nmaxlay,ncol))
      allocate(this%building_fraction(nmaxlay,ncol))
      allocate(this%building_scale(nmaxlay,ncol))
    end if
    if (do_canopy) then
      allocate(this%air_temperature(nmaxlay,ncol))
    end if
 
    if (do_vegetation) then
      allocate(this%vegetation_temperature(nmaxlay,ncol))
      allocate(this%vegetation_fraction(nmaxlay,ncol))
      allocate(this%vegetation_scale(nmaxlay,ncol))
    end if

    if (config%nvegregion > 1) then
      allocate(this%vegetation_fsd(nmaxlay,ncol))
    end if

    if (do_urban .or. do_vegetation) then
      allocate(this%vegetation_contact_fraction(nmaxlay,ncol))
    end if

    if (lhook) call dr_hook('radiation_canopy_properties:allocate',1,hook_handle)

  end subroutine allocate_canopy


  !---------------------------------------------------------------------
  subroutine deallocate_canopy(this)

    use yomhook, only : lhook, dr_hook

    class(canopy_properties_type), intent(inout) :: this

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_canopy_properties:deallocate',0,hook_handle)

    if (allocated(this%i_representation))       deallocate(this%i_representation)
    if (allocated(this%cos_sza))                deallocate(this%cos_sza)
    if (allocated(this%ground_temperature))     deallocate(this%ground_temperature)
    if (allocated(this%roof_temperature))       deallocate(this%roof_temperature)
    if (allocated(this%wall_temperature))       deallocate(this%wall_temperature)
    if (allocated(this%air_temperature))        deallocate(this%air_temperature)
    if (allocated(this%vegetation_temperature)) deallocate(this%vegetation_temperature)
    if (allocated(this%building_fraction))      deallocate(this%building_fraction)
    if (allocated(this%vegetation_fraction))    deallocate(this%vegetation_fraction)
    if (allocated(this%building_scale))         deallocate(this%building_scale)
    if (allocated(this%vegetation_scale))       deallocate(this%vegetation_scale)
    if (allocated(this%vegetation_fsd))         deallocate(this%vegetation_fsd)
    if (allocated(this%vegetation_contact_fraction)) deallocate(this%vegetation_contact_fraction)

    if (lhook) call dr_hook('radiation_canopy_properties:deallocate',1,hook_handle)

  end subroutine deallocate_canopy


  !---------------------------------------------------------------------
  ! Print a description of the surface tile types
  subroutine print_surface_representation(i_representation)

    use radiation_io, only : nulout

    integer(kind=jpim), dimension(:), allocatable, intent(in) :: i_representation

    integer :: ntile, jtile

    write(nulout,'(a)') 'Surface tile representation:'
    if (.not. allocated(i_representation)) then
      write(nulout,'(a)') '  Simple (one flat tile)'
    else
      ntile = size(i_representation,1)
      do jtile = 1,ntile
        write(nulout,'(a,i0,a,a)') '  Tile ', jtile, ': ', trim(TileRepresentationName(i_representation(jtile)))
      end do
    end if
    
  end subroutine print_surface_representation

#ifdef PANTS

  !---------------------------------------------------------------------
  subroutine read_from_netcdf(this, file)

    use parkind1,           only : jprb, jpim
    use easy_netcdf,        only : netcdf_file
    
    implicit none

    type(netcdf_file),  intent(in)     :: file
    class(surface_type), intent(inout) :: this

    real(kind=jprb), allocatable       :: data_1d(:)

    integer(kind=jpim), parameter :: ipermute(3) = [2,3,1]

    call this%deallocate

    this%is_simple = .false.

    call file%get('skin_temperature', this%skin_temperature, do_transp=.true.)
    call file%get('canopy_temperature', this%canopy_temperature, do_transp=.true.)
    call file%get('sw_albedo', this%sw_albedo, ipermute=ipermute)
    if (file%exists('sw_albedo_direct')) then
      call file%get('sw_albedo', this%sw_albedo_direct, ipermute=ipermute)
    end if
    call file%get('lw_emissivity', this%lw_emissivity, ipermute=ipermute)
    call file%get('tile_representation', data_1d)
    this%ntile = size(data_1d)
    allocate(this%i_representation(this%ntile))
    this%i_representation = int(data_1d)
    call file%get('tile_fraction', this%tile_fraction, do_transp=.true.)
    call file%get('canopy_depth', this%canopy_depth, do_transp=.true.)
    call file%get('building_fraction', this%building_fraction, do_transp=.true.)
    if (file%exists('building_normalized_perimeter')) then
      call file%get('building_normalized_perimeter', &
           &  this%building_normalized_perimeter, do_transp=.true.)
    else
      ! Convert building scale to normalized perimeter
      call file%get('building_scale', this%building_normalized_perimeter, do_transp=.true.)
      this%building_normalized_perimeter = 4.0_jprb * this%building_fraction &
           &  * (1.0_jprb-this%building_fraction) / max(1.0e-8_jprb,this%building_normalized_perimeter)
    end if
    call file%get('vegetation_optical_depth', this%vegetation_optical_depth, do_transp=.true.)
    call file%get('vegetation_sw_albedo', this%vegetation_sw_albedo, ipermute=ipermute)
    call file%get('vegetation_lw_emissivity', this%vegetation_lw_emissivity, ipermute=ipermute)

    this%nregion = sum(NTileRegions(this%i_representation))
    this%ncol    = size(this%skin_temperature,1)
    this%nfacet  = size(this%skin_temperature,2)

    this%nalbedobands = size(this%sw_albedo,2)
    this%nemissbands = size(this%lw_emissivity,2)

    call this%set_facet_indices

  end subroutine read_from_netcdf

#endif

end module radsurf_canopy_properties
