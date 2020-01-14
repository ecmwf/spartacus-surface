! spartacus_surface_config.F90 - Configure driver for SPARTACUS-Surface radiation scheme
!
! Copyright (C) 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details

module spartacus_surface_config

  use parkind1,                      only : jprb, jpim

  implicit none

  type driver_config_type
    
     ! Parallel settings
     logical :: do_parallel = .true.
     integer(kind=jpim) :: nblocksize = 16 ! Number of columns processed at once

     ! Do we repeat the core call, for benchmarking?
     integer(kind=jpim) :: nrepeat = 1

     ! Process a limited number of columns (iendcol=0 indicates to
     ! process from istartcol up to the end)
     integer(kind=jpim) :: istartcol = 1
     integer(kind=jpim) :: iendcol = 0

     ! Control verbosity in driver routine: 0=none (no output to
     ! standard output; write to standard error only if an error
     ! occurs), 1=warning, 2=info, 3=progress, 4=detailed, 5=debug
     integer(kind=jpim) :: iverbose = 3

     ! Override values
     real(kind=jprb) :: cos_sza_override     = -1.0
     real(kind=jprb) :: ground_sw_albedo     = -1.0
     real(kind=jprb) :: roof_sw_albedo       = -1.0
     real(kind=jprb) :: wall_sw_albedo       = -1.0
     real(kind=jprb) :: ground_lw_emissivity = -1.0
     real(kind=jprb) :: roof_lw_emissivity   = -1.0
     real(kind=jprb) :: wall_lw_emissivity   = -1.0
     real(kind=jprb) :: veg_fsd              = -1.0
     real(kind=jprb) :: veg_sw_extinction    = -1.0
     real(kind=jprb) :: veg_sw_ssa           = -1.0
     real(kind=jprb) :: veg_lw_extinction    = -1.0
     real(kind=jprb) :: veg_lw_ssa           = -1.0
     
   contains
     procedure :: read => read_config_from_namelist

  end type driver_config_type

contains

  !---------------------------------------------------------------------
  ! This subroutine reads configuration data from a namelist file, and
  ! anything that is not in the namelists will be set to default
  ! values. If optional output argument "is_success" is present, then on
  ! error (e.g. missing file) it will be set to .false.; if this
  ! argument is missing then on error the program will be aborted.
  subroutine read_config_from_namelist(this, file_name, is_success)

    use yomhook,      only : lhook, dr_hook
    use radiation_io, only : nulerr, radiation_abort

    class(driver_config_type), intent(inout), target :: this
    character(*), intent(in)          :: file_name
    logical, intent(out), optional    :: is_success

    integer :: iosopen ! Status after calling open

    logical,            pointer :: do_parallel
    integer(kind=jpim), pointer :: nblocksize, istartcol, iendcol, iverbose, nrepeat
    real(kind=jprb),    pointer :: cos_solar_zenith_angle, veg_fsd
    real(kind=jprb),    pointer :: ground_sw_albedo, roof_sw_albedo, wall_sw_albedo
    real(kind=jprb),    pointer :: ground_lw_emissivity, roof_lw_emissivity, wall_lw_emissivity

    namelist /radsurf_driver/ do_parallel, nblocksize, nrepeat, istartcol, iendcol, &
         &  iverbose, cos_solar_zenith_angle, veg_fsd, &
         &  ground_sw_albedo, roof_sw_albedo, wall_sw_albedo, &
         &  ground_lw_emissivity, roof_lw_emissivity, wall_lw_emissivity

    do_parallel            => this%do_parallel
    nblocksize             => this%nblocksize
    nrepeat                => this%nrepeat
    istartcol              => this%istartcol
    iendcol                => this%iendcol
    iverbose               => this%iverbose
    cos_solar_zenith_angle => this%cos_sza_override
    ground_sw_albedo       => this%ground_sw_albedo
    roof_sw_albedo         => this%roof_sw_albedo
    wall_sw_albedo         => this%wall_sw_albedo
    ground_lw_emissivity   => this%ground_lw_emissivity
    roof_lw_emissivity     => this%roof_lw_emissivity
    wall_lw_emissivity     => this%wall_lw_emissivity
    veg_fsd                => this%veg_fsd

    ! Open the namelist file and read the radiation_driver namelist
    open(unit=10, iostat=iosopen, file=trim(file_name))
    if (iosopen /= 0) then
      ! An error occurred
      if (present(is_success)) then
        is_success = .false.
        ! We now continue the subroutine so that the default values
        ! are placed in the config structure
      else
        write(nulerr,'(a,a,a)') '*** Error: namelist file "', &
             &                trim(file_name), '" not found'
        call radiation_abort('Driver configuration error')
      end if
    else
      ! Read the radiation_driver namelist, noting that it is not an
      ! error if this namelist is not present, provided all the required
      ! variables are present in the NetCDF data file instead
      read(unit=10, nml=radsurf_driver)
      close(unit=10)
    end if

  end subroutine read_config_from_namelist

end module spartacus_surface_config