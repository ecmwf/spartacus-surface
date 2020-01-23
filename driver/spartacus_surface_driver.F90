! spartacus_surface_driver.F90 - Driver for SPARTACUS-Surface radiation algorithm
!
! Copyright (C) 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! This program takes three arguments:
! 1) Namelist file to configure the radiation calculation
! 2) Name of a NetCDF file containing one or more columns
! 3) Name of output NetCDF file

program spartacus_surface_driver

  ! --------------------------------------------------------
  ! Section 1: Declarations
  ! --------------------------------------------------------

  use parkind1,                     only : jprb ! Working precision

  use radiation_io,                 only : nulout
  use radsurf_config,               only : config_type
  use spartacus_surface_config,     only : driver_config_type
  use radsurf_canopy_properties,    only : canopy_properties_type
  use radsurf_facet_properties,     only : facet_properties_type
  use radsurf_volume_properties,    only : volume_properties_type
  use radsurf_boundary_conds_out,   only : boundary_conds_out_type
  use spartacus_surface_read_input, only : read_input
  use radsurf_canopy_flux,          only : canopy_flux_type
  use radsurf_interface,            only : radsurf
  use radsurf_save,                 only : save_canopy_fluxes
  use easy_netcdf

  implicit none

  ! Derived types for the inputs to the radiation scheme
  type(config_type)            :: config
  type(canopy_properties_type) :: canopy_props
  type(facet_properties_type)  :: facet_props
  type(volume_properties_type) :: volume_props
  type(boundary_conds_out_type):: bc_out
  ! Canopy flux components, the last three normalized by top-of-canopy
  ! downwelling
  type(canopy_flux_type)       :: lw_internal, lw_norm, sw_norm_diff, sw_norm_dir
  ! Total canopy fluxes
  type(canopy_flux_type)       :: lw_flux, sw_flux

  ! Top-of-canopy downward radiation, all dimensioned (nspec,ncol)
  real(kind=jprb), allocatable &
       &  :: top_flux_dn_sw(:,:), &        ! Total shortwave (direct+diffuse)
       &     top_flux_dn_direct_sw(:,:), & ! ...diffuse only
       &     top_flux_dn_lw(:,:)           ! longwave

  ! Configuration specific to this driver
  type(driver_config_type)     :: driver_config

  ! The NetCDF file containing the input profiles
  type(netcdf_file)  :: file

  ! Name of file names specified on command line
  character(len=512) :: file_name
  integer            :: istatus ! Result of command_argument_count

  ! For parallel processing of multiple blocks
  integer :: jblock, nblock ! Block loop index and number
  integer, external :: omp_get_thread_num

  ! Loop index for repeats (for benchmarking)
  integer :: jrepeat

  integer :: ncol, ntotlay
  integer :: istartcol, iendcol ! Range of columns to process

  ! --------------------------------------------------------
  ! Section 2: Configure
  ! --------------------------------------------------------

  ! Check program called with correct number of arguments
  if (command_argument_count() < 3) then
    stop 'Usage: spartacus_surface config.nam input_file.nc output_file.nc'
  end if

  ! Use namelist to configure the radiation calculation
  call get_command_argument(1, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of namelist file as string of length < 512'
  end if

  ! Read "radsurf" namelist into radiation configuration type
  call config%read(file_name=file_name)
  call config%consolidate()

  ! Read "radsurf_driver" namelist into radiation driver config type
  call driver_config%read(file_name)
  
  ! --------------------------------------------------------
  ! Section 3: Read input data file
  ! --------------------------------------------------------

  ! Get NetCDF input file name
  call get_command_argument(2, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of input NetCDF file as string of length < 512'
  end if

  ! Open the file and configure the way it is read
  call file%open(trim(file_name), iverbose=config%iverbose)

  ! Get NetCDF output file name
  call get_command_argument(3, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of output NetCDF file as string of length < 512'
  end if

  ! Read input variables from NetCDF file
  call read_input(file, config, driver_config, ncol, ntotlay, &
       &  canopy_props, facet_props, volume_props, &
       &  top_flux_dn_sw, top_flux_dn_direct_sw, top_flux_dn_lw)

  ! Set first and last columns to process
  if (driver_config%iendcol < 1 .or. driver_config%iendcol > ncol) then
    driver_config%iendcol = ncol
  end if

  if (driver_config%istartcol > driver_config%iendcol) then
    write(nulout,'(a,i0,a,i0,a,i0,a)') '*** Error: requested column range (', &
         &  driver_config%istartcol, &
         &  ' to ', driver_config%iendcol, ') is out of the range in the data (1 to ', &
         &  ncol, ')'
    stop 1
  end if

  call bc_out%allocate(ncol, config%nsw, config%nlw)
  if (config%do_sw) then
    call sw_norm_dir%allocate(ncol, ntotlay, config%nsw, use_direct=.true.)
    call sw_norm_diff%allocate(ncol, ntotlay, config%nsw, use_direct=.true.)
    call sw_flux%allocate(ncol, ntotlay, config%nsw, use_direct=.true.)
  end if
  if (config%do_lw) then
    call lw_internal%allocate(ncol, ntotlay, config%nlw, use_direct=.false.)
    call lw_norm%allocate(ncol, ntotlay, config%nlw, use_direct=.false.)
    call lw_flux%allocate(ncol, ntotlay, config%nlw, use_direct=.false.)
  end if
  
  ! --------------------------------------------------------
  ! Section 4: Call radiation scheme
  ! --------------------------------------------------------

  call facet_props%calc_monochromatic_emission(canopy_props)

  ! Option of repeating calculation multiple time for more accurate
  ! profiling
  do jrepeat = 1,driver_config%nrepeat

    if (driver_config%do_parallel) then
      ! Run radiation scheme over blocks of columns in parallel
        
      ! Compute number of blocks to process
      nblock = (driver_config%iendcol - driver_config%istartcol &
           &  + driver_config%nblocksize) / driver_config%nblocksize
        
      !$OMP PARALLEL DO PRIVATE(istartcol, iendcol) SCHEDULE(RUNTIME)
      do jblock = 1, nblock
        ! Specify the range of columns to process.
        istartcol = (jblock-1) * driver_config%nblocksize &
             &    + driver_config%istartcol
        iendcol = min(istartcol + driver_config%nblocksize - 1, &
             &        driver_config%iendcol)
        
        if (driver_config%iverbose >= 3) then
          write(nulout,'(a,i0,a,i0,a,i0)')  'Thread ', omp_get_thread_num(), &
               &  ' processing columns ', istartcol, '-', iendcol
        end if

        ! Call the SPARTACUS-Surface radiation scheme
        call radsurf(config, canopy_props, facet_props, volume_props, bc_out, &
             &       istartcol, iendcol, sw_norm_dir, sw_norm_diff, &
             &       lw_internal, lw_norm)
        
      end do
      !$OMP END PARALLEL DO
        
    else
      ! Run radiation scheme serially
      if (driver_config%iverbose >= 3) then
        write(nulout,'(a,i0,a)')  'Processing ', ncol, ' columns'
      end if
      
      ! Call the SPARTACUS-Surface radiation scheme
      call radsurf(config, canopy_props, facet_props, volume_props, bc_out, &
           &       istartcol, iendcol, sw_norm_dir, sw_norm_diff, &
           &       lw_internal, lw_norm)
      
    end if

    ! Scale the normalized fluxes
    if (config%do_sw) then
      call sw_norm_dir%scale(canopy_props%nlay, top_flux_dn_direct_sw)
      call sw_norm_diff%scale(canopy_props%nlay, &
           &  top_flux_dn_sw - top_flux_dn_direct_sw)
      call sw_flux%sum(sw_norm_dir, sw_norm_diff)
    end if

    if (config%do_lw) then
      call lw_norm%scale(canopy_props%nlay, top_flux_dn_lw)
      call lw_flux%sum(lw_internal, lw_norm)
    end if

  end do

  ! --------------------------------------------------------
  ! Section 5: Check and save output
  ! --------------------------------------------------------
  call save_canopy_fluxes(trim(file_name), config, canopy_props, &
       &  sw_flux, lw_flux, iverbose=driver_config%iverbose)

end program spartacus_surface_driver
