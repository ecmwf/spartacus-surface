! spartacus_surface_driver.F90 - Driver for SPARTACUS-Surface radiation algorithm
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
! This program takes three arguments:
! 1) Namelist file to configure the radiation calculation
! 2) Name of a NetCDF file containing one or more columns
! 3) Name of output NetCDF file

program spartacus_surface_driver

  ! --------------------------------------------------------
  ! Section 1: Declarations
  ! --------------------------------------------------------

  use parkind1,                     only : jprb, jprd ! Working precision

  use radiation_io,                 only : nulout
  use radsurf_config,               only : config_type
  use spartacus_surface_config,     only : driver_config_type
  use radsurf_canopy_properties,    only : canopy_properties_type
  use radsurf_sw_spectral_properties,only: sw_spectral_properties_type
  use radsurf_lw_spectral_properties,only: lw_spectral_properties_type
  use radsurf_boundary_conds_out,   only : boundary_conds_out_type
  use spartacus_surface_read_input, only : read_input
  use radsurf_canopy_flux,          only : canopy_flux_type
  use radsurf_interface,            only : radsurf
  use radsurf_save,                 only : save_canopy_fluxes
  use radsurf_simple_spectrum,      only : calc_simple_spectrum_lw
  use easy_netcdf
  use print_matrix_mod, only : print_vector, print_matrix
  
  implicit none

  ! Derived types for the inputs to the radiation scheme
  type(config_type)                 :: config
  type(canopy_properties_type)      :: canopy_props
  type(sw_spectral_properties_type) :: sw_spectral_props
  type(lw_spectral_properties_type) :: lw_spectral_props
  type(boundary_conds_out_type)     :: bc_out

  ! Canopy flux components, the last three normalized by top-of-canopy
  ! downwelling
  type(canopy_flux_type)       :: lw_internal,  lw_norm,&
       &                          sw_norm_diff, sw_norm_dir
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

#ifndef NO_OPENMP
  integer, external :: omp_get_thread_num
  double precision, external :: omp_get_wtime
  ! Start/stop time in seconds
  real(kind=jprd) :: tstart, tstop
#endif

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
  
  ! Read "radsurf_driver" namelist into radiation driver config type
  call driver_config%read(file_name)

  ! Print out configuration information
  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)') '------------------ OFFLINE SPARTACUS-SURFACE RADIATION SCHEME ------------------'
    write(nulout,'(a)') 'Copyright (C) 2019- ECMWF'
    write(nulout,'(a)') 'Contact: Robin Hogan (r.j.hogan@ecmwf.int)'
#ifdef SINGLE_PRECISION
    write(nulout,'(a)') 'Floating-point precision: single'
#else
    write(nulout,'(a)') 'Floating-point precision: double'
#endif
    call config%print(driver_config%iverbose)
  end if

  ! Act on any information in the configuration, e.g. loading config
  ! files
  call config%consolidate()
  
  ! --------------------------------------------------------
  ! Section 3: Read input data file
  ! --------------------------------------------------------

  ! Get NetCDF input file name
  call get_command_argument(2, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of input NetCDF file as string of length < 512'
  end if

  ! Open the file and configure the way it is read
  call file%open(trim(file_name), iverbose=driver_config%iverbose)

  ! Get NetCDF output file name
  call get_command_argument(3, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of output NetCDF file as string of length < 512'
  end if

  ! Read input variables from NetCDF file
  call read_input(file, config, driver_config, ncol, ntotlay, &
       &  canopy_props, sw_spectral_props, lw_spectral_props, &
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
    call sw_norm_dir%allocate(config, ncol, ntotlay, config%nsw, use_direct=.true.)
    call sw_norm_diff%allocate(config, ncol, ntotlay, config%nsw, use_direct=.true.)
    
    call sw_norm_dir%zero_all()
    call sw_norm_diff%zero_all()

    call sw_flux%allocate(config, ncol, ntotlay, config%nsw, use_direct=.true.)
  end if
  if (config%do_lw) then
    call lw_internal%allocate(config, ncol, ntotlay, config%nlw, use_direct=.false.)
    call lw_norm%allocate(config, ncol, ntotlay, config%nlw, use_direct=.false.)

    call lw_internal%zero_all()
    call lw_norm%zero_all()

    call lw_flux%allocate(config, ncol, ntotlay, config%nlw, use_direct=.false.)
  end if
  
  ! --------------------------------------------------------
  ! Section 4: Call radiation scheme
  ! --------------------------------------------------------

  call lw_spectral_props%calc_monochromatic_emission(canopy_props)

  ! Option of repeating calculation multiple time for more accurate
  ! profiling
#ifndef NO_OPENMP
  tstart = omp_get_wtime()
#endif
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
        
#ifndef NO_OPENMP
        if (driver_config%iverbose >= 3) then
          write(nulout,'(a,i0,a,i0,a,i0)')  'Thread ', omp_get_thread_num(), &
               &  ' processing columns ', istartcol, '-', iendcol
        end if
#endif

        if (config%do_lw) then
          ! Gas optics and spectral emission
          call calc_simple_spectrum_lw(config, canopy_props, lw_spectral_props, &
               &                       istartcol, iendcol)
        end if
        
        ! Call the SPARTACUS-Surface radiation scheme
        call radsurf(config, canopy_props, &
             &       sw_spectral_props, lw_spectral_props, bc_out, &
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
      call radsurf(config, canopy_props, &
           &       sw_spectral_props, lw_spectral_props, bc_out, &
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
#ifndef NO_OPENMP
  tstop = omp_get_wtime()
  if (driver_config%iverbose >= 2) then
    write(nulout, '(a,g11.5,a)') 'Time elapsed in radiative transfer: ', tstop-tstart, ' seconds'
  endif
#endif

  ! --------------------------------------------------------
  ! Section 5: Check and save output
  ! --------------------------------------------------------

  if (driver_config%do_conservation_check) then
    istartcol = driver_config%istartcol
    iendcol   = driver_config%iendcol
    if (iendcol <= 0) then
      iendcol = ncol
    end if
    if (config%do_sw) then
      write(nulout, '(a)') 'Direct shortwave budget: radiation originating from direct solar at canopy top'
      call sw_norm_dir%check(canopy_props, istartcol, iendcol)
      write(nulout, '(a)') 'Diffuse shortwave budget: radiation originating from downward diffuse solar at canopy top'
      call sw_norm_diff%check(canopy_props, istartcol, iendcol)
    end if
    if (config%do_lw) then
      write(nulout, '(a)') 'Internal longwave budget: radiation originating from emission within canopy'
      call lw_internal%check(canopy_props, istartcol, iendcol)
      write(nulout, '(a)') 'Incoming longwave budget: radiation originating from downward longwave at canopy top'
      call lw_norm%check(canopy_props, istartcol, iendcol)
    end if
  end if

  call save_canopy_fluxes(trim(file_name), config, canopy_props, &
       &  sw_flux, lw_flux, iverbose=driver_config%iverbose)

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)') '--------------------------------------------------------------------------------'
  end if

end program spartacus_surface_driver
