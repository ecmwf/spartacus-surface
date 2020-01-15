! radsurf_volume_properties.F90 - Derived type for volume spectral properties
!
! Copyright (C) 2019 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radsurf_volume_properties

  use parkind1, only : jpim, jprb

  implicit none

  !---------------------------------------------------------------------
  ! Derived type storing spectral properties of the volumes of the surface
  type volume_properties_type

    ! Shortwave spectral properties
    real(kind=jprb), allocatable :: air_sw_ext(:,:) ! m-1 (nsw,ntotlay,ncol)
    real(kind=jprb), allocatable :: air_sw_ssa(:,:) !     (nsw,ntotlay)
    real(kind=jprb), allocatable :: veg_sw_ext(:,:) ! m-1 (nsw,ntotlay)
    real(kind=jprb), allocatable :: veg_sw_ssa(:,:) !     (nsw,ntotlay)

    ! Longwave spectral properties
    real(kind=jprb), allocatable :: air_lw_ext(:,:) ! m-1 (nlw,ntotlay)
    real(kind=jprb), allocatable :: air_lw_ssa(:,:) !     (nlw,ntotlay)
    real(kind=jprb), allocatable :: veg_lw_ext(:,:) ! m-1 (nlw,ntotlay)
    real(kind=jprb), allocatable :: veg_lw_ssa(:,:) !     (nlw,ntotlay)

    ! Number of columns and maximum number of layers
    integer(kind=jpim) :: ncol, ntotlay

    ! Number of shortwave albedo and longwave emissivity bands
    integer(kind=jpim) :: nsw, nlw

  contains
    procedure :: allocate   => allocate_volume
    procedure :: deallocate => deallocate_volume
    !procedure :: read => read_from_netcdf
    
  end type volume_properties_type

contains


  !---------------------------------------------------------------------
  subroutine allocate_volume(this, config, ncol, ntotlay, &
       &                    i_representation, is_internal)

    use yomhook,        only : lhook, dr_hook
    use radsurf_config, only : config_type
    use radsurf_canopy_properties, only : ITileFlat,  ITileForest, &
         &                                ITileUrban, ITileVegetatedUrban

    class(volume_properties_type), intent(inout) :: this
    type(config_type),   intent(in)    :: config
    integer(kind=jpim),  intent(in)    :: ncol, ntotlay
    integer(kind=jpim),  intent(in), optional :: i_representation(:)
    logical, intent(in), optional :: is_internal

    logical :: do_vegetation, do_urban, do_canopy

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_volume_properties:allocate',0,hook_handle)

    call this%deallocate()

    this%ncol     = ncol
    this%ntotlay  = ntotlay
    do_vegetation = config%do_vegetation
    do_urban      = config%do_urban
    this%nsw      = config%nsw
    this%nlw      = config%nlw
    if (present(is_internal)) then
      if (is_internal) then
        this%nsw  = config%nswinternal
        this%nlw  = config%nlwinternal
      end if
    end if

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
      if (config%do_sw) then      
        allocate(this%air_sw_ext(this%nsw,ntotlay))
        allocate(this%air_sw_ssa(this%nsw,ntotlay))
        if (do_vegetation) then
          allocate(this%veg_sw_ext(this%nsw,ntotlay))
          allocate(this%veg_sw_ssa(this%nsw,ntotlay))
        end if
      end if
      if (config%do_lw) then      
        allocate(this%air_lw_ext(this%nlw,ntotlay))
        allocate(this%air_lw_ssa(this%nlw,ntotlay))
        if (do_vegetation) then
          allocate(this%veg_lw_ext(this%nlw,ntotlay))
          allocate(this%veg_lw_ssa(this%nlw,ntotlay))
        end if
      end if
    end if

    if (lhook) call dr_hook('radiation_volume_properties:allocate',1,hook_handle)

  end subroutine allocate_volume

  !---------------------------------------------------------------------
  subroutine deallocate_volume(this)

    use yomhook, only : lhook, dr_hook

    class(volume_properties_type), intent(inout) :: this

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_volume_properties:deallocate',0,hook_handle)

    if (allocated(this%air_sw_ext)) deallocate(this%air_sw_ext)
    if (allocated(this%air_sw_ssa)) deallocate(this%air_sw_ssa)
    if (allocated(this%veg_sw_ext)) deallocate(this%veg_sw_ext)
    if (allocated(this%veg_sw_ssa)) deallocate(this%veg_sw_ssa)
    if (allocated(this%air_lw_ext)) deallocate(this%air_lw_ext)
    if (allocated(this%air_lw_ssa)) deallocate(this%air_lw_ssa)
    if (allocated(this%veg_lw_ext)) deallocate(this%veg_lw_ext)
    if (allocated(this%veg_lw_ssa)) deallocate(this%veg_lw_ssa)

    if (lhook) call dr_hook('radiation_volume_properties:deallocate',1,hook_handle)

  end subroutine deallocate_volume


end module radsurf_volume_properties
