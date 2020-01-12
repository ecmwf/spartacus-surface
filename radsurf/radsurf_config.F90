! radsurf_config.f90 - Derived type for configuring the radiative treament of surfaces
!
! Copyright (C) 2019-2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
module radsurf_config

  use parkind1, only : jpim, jprb
  use radtool_legendre_gauss, only : legendre_gauss_type

  implicit none

  !---------------------------------------------------------------------
  ! Derived type containing all the configuration information needed
  ! to run the radiation scheme.  The intention is that this is fixed
  ! for a given model run.  The parameters are to list first those
  ! quantities that can be set directly by the user, for example using
  ! a namelist, and second those quantities that are computed
  ! afterwards from the user-supplied numbers.
  type config_type

    ! USER-CONFIGURABLE PARAMETERS

    ! Do we do shortwave and longwave calculations, and include
    ! shortwave direct albedo?
    logical :: do_sw = .true.
    logical :: do_lw = .true.
    logical :: use_sw_direct_albedo = .false.

    ! Do we include the capability to do vegetation and urban areas?
    ! If not, some arrays will not be allocated
    logical :: do_vegetation = .true.
    logical :: do_urban = .true.

    ! Number of regions to represent vegetation (1 or 2)
    integer :: nvegregion = 2

    ! Number of shortwave and longwave bands for user to supply
    ! properties of facets
    integer :: nsw = 1
    integer :: nlw = 1

    ! Number of diffuse streams to use in a single hemisphere for a
    ! pure vegetation tile and an urban tile
    integer :: n_stream_vegetation = 4
    integer :: n_stream_urban = 4

    integer :: iverbose = 3

    ! COMPUTED PARAMETERS

    ! Number of shortwave and longwave bands used internally to
    ! represent the spectral variation of atmospheric absorption
    integer :: nswinternal, nlwinternal

    type(legendre_gauss_type) :: lg_urban, lg_vegetation

   contains
     procedure :: read => read_config_from_namelist

  end type config_type

contains
  
  !---------------------------------------------------------------------
  ! This subroutine reads configuration data from a namelist file, and
  ! anything that is not in the namelists will be set to default
  ! values. If optional output argument "is_success" is present, then
  ! on error (e.g. missing file) it will be set to .false.; if this
  ! argument is missing then on error the program will be aborted. You
  ! may either specify the file_name or the unit of an open file to
  ! read, but not both.
  subroutine read_config_from_namelist(this, file_name, unit, is_success)

    use yomhook,      only : lhook, dr_hook
    use radiation_io, only : nulout, nulerr, nulrad, radiation_abort

    class(config_type), intent(inout), target :: this
    character(*),       intent(in),  optional :: file_name
    integer,            intent(in),  optional :: unit
    logical,            intent(out), optional :: is_success

    integer :: iosopen, iosread ! Status after calling open and read
    integer :: iunit ! Unit number of namelist file

    ! Namelist variables mirroring values in config_type
    
    logical, pointer :: do_sw, do_lw, use_sw_direct_albedo, do_vegetation, do_urban
    integer, pointer :: nvegregion, nsw, nlw, n_stream_vegetation, &
         &  n_stream_urban, iverbose

    namelist /radsurf/ do_sw, do_lw, use_sw_direct_albedo, do_vegetation, do_urban, &
         &  nvegregion, nsw, nlw, n_stream_vegetation, n_stream_urban, iverbose

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radsurf_config:read',0,hook_handle)

    do_sw                => this%do_sw
    do_lw                => this%do_lw
    use_sw_direct_albedo => this%use_sw_direct_albedo
    do_vegetation        => this%do_vegetation
    do_urban             => this%do_urban
    nvegregion           => this%nvegregion
    nsw                  => this%nsw
    nlw                  => this%nlw
    n_stream_vegetation  => this%n_stream_vegetation
    n_stream_urban       => this%n_stream_urban
    iverbose             => this%iverbose

    if (present(file_name) .and. present(unit)) then
      write(nulerr,'(a)') '*** Error: cannot specify both file_name and unit in call to config_type%read'
      call radiation_abort('Radiation configuration error')
    else if (.not. present(file_name) .and. .not. present(unit)) then
      write(nulerr,'(a)') '*** Error: neither file_name nor unit specified in call to config_type%read'
      call radiation_abort('Radiation configuration error')
    end if

    if (present(file_name)) then
      ! Open the namelist file
      iunit = nulrad
      open(unit=iunit, iostat=iosopen, file=trim(file_name))
    else
      ! Assume that iunit represents and open file
      iosopen = 0
      iunit = unit
    end if

    if (iosopen /= 0) then
      ! An error occurred opening the file
      if (present(is_success)) then
        is_success = .false.
        ! We now continue the subroutine so that the default values
        ! are placed in the config structure
      else
        write(nulerr,'(a,a,a)') '*** Error: namelist file "', &
             &                trim(file_name), '" not found'
        call radiation_abort('Radiation configuration error')
      end if
    else
      read(unit=iunit, iostat=iosread, nml=radsurf)
      if (iosread /= 0) then
        ! An error occurred reading the file
        if (present(is_success)) then
          is_success = .false.
          ! We now continue the subroutine so that the default values
          ! are placed in the config structure
        else if (present(file_name)) then
          write(nulerr,'(a,a,a)') '*** Error reading namelist "radsurf" from file "', &
               &      trim(file_name), '"'
          close(unit=iunit)
          call radiation_abort('Radiation configuration error')
        else
          write(nulerr,'(a,i0)') '*** Error reading namelist "radsurf" from unit ', &
               &      iunit
          call radiation_abort('Radiation configuration error')
        end if
      end if

      if (present(file_name)) then
        close(unit=iunit)
      end if
    end if

    if (lhook) call dr_hook('radsurf_config:read',1,hook_handle)

  end subroutine read_config_from_namelist

end module radsurf_config
