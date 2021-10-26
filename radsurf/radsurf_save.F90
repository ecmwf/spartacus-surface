! radsurf_save.f90 - Save results of surface radiation calculation
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

module radsurf_save

  use parkind1, only : jprb, jpim

  implicit none

  real(kind=jprb), parameter :: FillValueFlux = -9999.0_jprb

contains

  subroutine save_canopy_fluxes(file_name, config, canopy_props, &
       &                        flux_sw, flux_lw, &
       &                        iverbose, is_hdf5_file)

    use easy_netcdf
    use radiation_io,             only : nulout
    use radsurf_config,           only : config_type
    use radsurf_canopy_properties,only : canopy_properties_type
    use radsurf_canopy_flux,      only : canopy_flux_type

    character(len=*), intent(in)        :: file_name
    type(config_type), intent(in)       :: config
    type(canopy_properties_type), intent(in) :: canopy_props
    type(canopy_flux_type), intent(in)  :: flux_sw, flux_lw
    logical, optional, intent(in)       :: is_hdf5_file
    integer(kind=jpim), optional, intent(in) :: iverbose

    integer(kind=jpim)     :: i_local_verbose
    type(netcdf_file)      :: out_file
    integer(kind=jpim)     :: jcol, jlay, ncol, nlay, ninterface, itotlay
    logical                :: do_spectral_sw, do_spectral_lw
    logical                :: do_broadband_sw, do_broadband_lw

    ! Temporary variable at layer interfaces
    real(kind=jprb), allocatable :: tmp(:,:)

    if (present(iverbose)) then
      i_local_verbose = iverbose
    else
      i_local_verbose = config%iverbose
    end if

    if (config%do_sw) then
      do_spectral_sw  = config%do_save_spectral_flux
      ! .and. flux_sw%nspec > 1
      do_broadband_sw = config%do_save_broadband_flux
    else
      do_spectral_sw  = .false.
      do_broadband_sw = .false.
    end if

    if (config%do_lw) then
      do_spectral_lw  = config%do_save_spectral_flux
      ! .and. flux_lw%nspec > 1
      do_broadband_sw = config%do_save_broadband_flux
    else
      do_spectral_lw  = .false.
      do_broadband_sw = .false.
    end if

    ncol = canopy_props%ncol
    nlay = maxval(canopy_props%nlay)
    ninterface = nlay+1

    allocate(tmp(ninterface,ncol))

    ! Open the file
    call out_file%create(trim(file_name), iverbose=i_local_verbose, &
         &               is_hdf5_file=is_hdf5_file)

    ! Define dimensions
    call out_file%define_dimension("column", ncol)
    call out_file%define_dimension("layer", nlay)
    call out_file%define_dimension("layer_interface", ninterface)
    if (do_spectral_sw) then
      call out_file%define_dimension("band_sw", flux_sw%nspec)
    end if
    if (do_spectral_lw) then
      call out_file%define_dimension("band_lw", flux_lw%nspec)
    end if
    
    ! Put global attributes
    call out_file%put_global_attributes( &
         &  title_str="Radiative fluxes from the SPARTACUS-Surface radiation model", &
         &  references_str="Hogan, R. J., T. Quaife and R. Braghiere, 2018: Fast matrix treatment of 3-D radiative" &
         &  //" transfer in vegetation canopies: SPARTACUS-Vegetation 1.1. Geosci. Model Dev., 11, 339-350." &
         &  //NEW_LINE('A')//"Hogan, R. J., 2019: Flexible treatment of radiative transfer in complex urban" &
         &  // " canopies for use in weather and climate models. Boundary-Layer Meteorol., 173, 53-78.", &
         &  source_str="SPARTACUS-Surface offline radiation model", &
         &  comment_str="All fluxes and absorption rates are in terms of power per unit horizontal area of the domain. " &
         &  //"Net fluxes are downwelling (or incoming) minus upwelling (or outgoing).")

    ! Define general variables
    call out_file%define_variable("height", &
         &  dim2_name="column", dim1_name="layer_interface", fill_value=-1.0_jprb, &
         &  units_str="m", long_name="Height of layer interfaces above ground", &
         &  standard_name="height")
    call out_file%define_variable("surface_type", data_type_name="short", &
         &  dim1_name="column", long_name="Surface type")
    call out_file%put_attribute("surface_type", "definition", &
         &    "0: Flat"//NEW_LINE('A') &
         &  //"1: Forest"//NEW_LINE('A') &
         &  //"2: Unvegetated urban"//NEW_LINE('A') &
         &  //"3: Vegetated urban"//NEW_LINE('A') &
         &  //"4: Simple urban"//NEW_LINE('A') &
         &  //"5: Infinite street")
    call out_file%define_variable("nlayer", data_type_name="short", &
         &  dim1_name="column", long_name="Number of active layers")

    ! Define shortwave variables
    if (config%do_sw) then
      call define_canopy_flux_variables(out_file, "sw", "shortwave", flux_sw, &
           &  do_broadband_sw, do_spectral_sw)
    end if

    ! Define longwave variables
    if (config%do_lw) then
      call define_canopy_flux_variables(out_file, "lw", "longwave", flux_lw, &
           &  do_broadband_lw, do_spectral_lw)
     end if

    ! Write general variables
    tmp = -1.0_jprb
    tmp(1,:) = 0.0_jprb ! Surface height is zero m
    itotlay = 1
    do jcol = 1,ncol
      do jlay = 1,canopy_props%nlay(jcol)
        tmp(jlay+1,jcol) = tmp(jlay,jcol) + canopy_props%dz(itotlay)
        itotlay = itotlay + 1
      end do
    end do
    call out_file%put("height", tmp)
    call out_file%put("surface_type", canopy_props%i_representation)
    call out_file%put("nlayer", canopy_props%nlay)

    ! Write shortwave variables
    if (config%do_sw) then
      call write_canopy_flux_variables(out_file, "sw", nlay, &
       &  canopy_props%nlay, flux_sw, do_broadband_sw, do_spectral_sw)
    end if

    ! Write longwave variables
    if (config%do_lw) then
      call write_canopy_flux_variables(out_file, "lw", nlay, &
       &  canopy_props%nlay, flux_lw, do_broadband_lw, do_spectral_lw)
    end if

    ! Close file
    call out_file%close()

  end subroutine save_canopy_fluxes

  subroutine define_canopy_flux_variables(out_file, band_name, band_long_name, &
       &  flux, do_broadband, do_spectral)

    use easy_netcdf
    use radsurf_canopy_flux,      only : canopy_flux_type

    type(netcdf_file),      intent(inout) :: out_file
    character(len=*),       intent(in)    :: band_name, band_long_name
    type(canopy_flux_type), intent(in)    :: flux
    logical,                intent(in)    :: do_broadband, do_spectral

    ! Define wavelength-independent variables
    if (allocated(flux%ground_dn_dir)) then
      call out_file%define_variable("ground_sunlit_fraction", units_str="1", &
           &  long_name="Fraction of ground in direct sunlight", dim1_name="column")
    end if
    if (allocated(flux%roof_in_dir)) then
      call out_file%define_variable("roof_sunlit_fraction", units_str="1", &
           &  long_name="Fraction of roof in direct sunlight", &
           &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
    end if
    if (allocated(flux%wall_in_dir)) then
      call out_file%define_variable("wall_sunlit_fraction", units_str="1", &
           &  long_name="Fraction of wall in direct sunlight", &
           &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
    end if
    if (allocated(flux%veg_abs_dir)) then
      call out_file%define_variable("veg_sunlit_fraction", units_str="1", &
           &  long_name="Fraction of vegetation in direct sunlight", &
           &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
    end if

    ! Define broadband variables
    if (do_broadband) then
      call out_file%define_variable("ground_flux_dn_"//band_name, units_str="W m-2", &
           &  long_name="Downwelling "//band_long_name//" flux at ground", &
           &  dim1_name="column")
      call out_file%define_variable("ground_flux_net_"//band_name, units_str="W m-2", &
           &  long_name="Net "//band_long_name//" flux at ground", &
           &  dim1_name="column")
      if (allocated(flux%ground_dn_dir)) then
        call out_file%define_variable("ground_flux_dn_direct_"//band_name, units_str="W m-2", &
             &  long_name="Downwelling direct "//band_long_name//" flux at ground", &
             &  dim1_name="column")
        call out_file%define_variable("ground_flux_vertical_diffuse_"//band_name, units_str="W m-2", &
             &  long_name="Diffuse "//band_long_name//" flux into a vertical surface at ground level", &
             &  dim1_name="column")
      else
        call out_file%define_variable("ground_flux_vertical_"//band_name, units_str="W m-2", &
             &  long_name="Flux in "//band_long_name//" into a vertical surface at ground level", &
             &  dim1_name="column")
      end if
      call out_file%define_variable("top_flux_dn_"//band_name, units_str="W m-2", &
           &  long_name="Downwelling "//band_long_name//" flux at top of canopy", &
           &  dim1_name="column")
      call out_file%define_variable("top_flux_net_"//band_name, units_str="W m-2", &
           &  long_name="Net "//band_long_name//" flux at top of canopy", &
           &  dim1_name="column")
      if (allocated(flux%top_dn_dir)) then
        call out_file%define_variable("top_flux_dn_direct_"//band_name, units_str="W m-2", &
             &  long_name="Downwelling direct "//band_long_name//" flux at top of canopy", &
             &  dim1_name="column")
      end if
      if (allocated(flux%roof_in)) then
        call out_file%define_variable("roof_flux_in_"//band_name, units_str="W m-2", &
             &  long_name="Incoming "//band_long_name//" flux at roofs", &
             &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
        if (allocated(flux%roof_in_dir)) then
          call out_file%define_variable("roof_flux_in_direct_"//band_name, units_str="W m-2", &
               &  long_name="Direct incoming "//band_long_name//" flux at roofs", &
               &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
        end if
        call out_file%define_variable("roof_flux_net_"//band_name, units_str="W m-2", &
             &  long_name="Net "//band_long_name//" flux at roofs", &
             &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
        call out_file%define_variable("wall_flux_in_"//band_name, units_str="W m-2", &
             &  long_name="Incoming "//band_long_name//" flux at walls", &
             &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
        if (allocated(flux%wall_in_dir)) then
          call out_file%define_variable("wall_flux_in_direct_"//band_name, units_str="W m-2", &
               &  long_name="Direct incoming "//band_long_name//" flux at walls", &
               &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
        end if
        call out_file%define_variable("wall_flux_net_"//band_name, units_str="W m-2", &
             &  long_name="Net "//band_long_name//" flux at walls", &
             &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
      end if
      if (allocated(flux%clear_air_abs)) then
        call out_file%define_variable("clear_air_absorption_"//band_name, units_str="W m-2", &
             &  long_name="Absorbed "//band_long_name//" in clear air", &
             &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
      end if

      if (allocated(flux%veg_abs)) then
        call out_file%define_variable("veg_absorption_"//band_name, units_str="W m-2", &
             &  long_name="Absorbed "//band_long_name//" by vegetation", &
             &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
        call out_file%define_variable("veg_air_absorption_"//band_name, units_str="W m-2", &
             &  long_name="Absorbed "//band_long_name//" by air in vegetated regions", &
             &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
      end if
      if (allocated(flux%veg_abs_dir)) then
        call out_file%define_variable("veg_absorption_direct_"//band_name, units_str="W m-2", &
             &  long_name="Absorbed direct "//band_long_name//" by vegetation", &
             &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
      end if

      if (allocated(flux%flux_dn_layer_top)) then
        call out_file%define_variable("flux_dn_layer_top_"//band_name, units_str="W m-2", &
             &  long_name="Downwelling "//band_long_name//" flux at top of layer", &
             &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
        if (allocated(flux%flux_dn_dir_layer_top)) then
          call out_file%define_variable("flux_dn_direct_layer_top_"//band_name, units_str="W m-2", &
               &  long_name="Downwelling direct "//band_long_name//" flux at top of layer", &
               &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
        end if
        call out_file%define_variable("flux_up_layer_top_"//band_name, units_str="W m-2", &
             &  long_name="Upwelling "//band_long_name//" flux at top of layer", &
             &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
        call out_file%define_variable("flux_dn_layer_base_"//band_name, units_str="W m-2", &
             &  long_name="Downwelling "//band_long_name//" flux at base of layer", &
             &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
        if (allocated(flux%flux_dn_dir_layer_base)) then
          call out_file%define_variable("flux_dn_direct_layer_base_"//band_name, units_str="W m-2", &
               &  long_name="Downwelling direct "//band_long_name//" flux at base of layer", &
               &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
        end if
        call out_file%define_variable("flux_up_layer_base_"//band_name, units_str="W m-2", &
             &  long_name="Upwelling "//band_long_name//" flux at base of layer", &
             &  dim2_name="column", dim1_name="layer", fill_value=FillValueFlux)
      end if
    end if

    ! Define spectral fluxes
    if (do_spectral) then
      call out_file%define_variable("ground_spectral_flux_dn_"//band_name, units_str="W m-2", &
           &  long_name="Downwelling "//band_long_name//" spectral flux at ground", &
           &  dim2_name="column", dim1_name="band_"//band_name)
      call out_file%define_variable("ground_spectral_flux_net_"//band_name, units_str="W m-2", &
           &  long_name="Net "//band_long_name//" spectral flux at ground", &
           &  dim2_name="column", dim1_name="band_"//band_name)
      if (allocated(flux%ground_dn_dir)) then
        call out_file%define_variable("ground_spectral_flux_dn_direct_"//band_name, units_str="W m-2", &
             &  long_name="Downwelling direct "//band_long_name//" spectral flux at ground", &
             &  dim2_name="column", dim1_name="band_"//band_name)
        call out_file%define_variable("ground_spectral_flux_vertical_diffuse_"//band_name, units_str="W m-2", &
             &  long_name="Diffuse "//band_long_name//" spectral flux into a vertical surface at ground level", &
             &  dim2_name="column", dim1_name="band_"//band_name)
      else
        call out_file%define_variable("ground_spectral_flux_vertical_"//band_name, units_str="W m-2", &
             &  long_name="Flux in "//band_long_name//" into a vertical surface at ground level", &
             &  dim2_name="column", dim1_name="band_"//band_name)
      end if
      call out_file%define_variable("top_spectral_flux_dn_"//band_name, units_str="W m-2", &
           &  long_name="Downwelling "//band_long_name//" spectral flux at top of canopy", &
           &  dim2_name="column", dim1_name="band_"//band_name)
      call out_file%define_variable("top_spectral_flux_net_"//band_name, units_str="W m-2", &
           &  long_name="Net "//band_long_name//" spectral flux at top of canopy", &
           &  dim2_name="column", dim1_name="band_"//band_name)
      if (allocated(flux%top_dn_dir)) then
        call out_file%define_variable("top_spectral_flux_dn_direct_"//band_name, units_str="W m-2", &
             &  long_name="Downwelling direct "//band_long_name//" spectral flux at top of canopy", &
             &  dim2_name="column", dim1_name="band_"//band_name)
      end if
      if (allocated(flux%roof_in)) then
        call out_file%define_variable("roof_spectral_flux_in_"//band_name, units_str="W m-2", &
             &  long_name="Incoming "//band_long_name//" spectral flux at roofs", &
             &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
             &  fill_value=FillValueFlux)
        if (allocated(flux%roof_in_dir)) then
          call out_file%define_variable("roof_spectral_flux_in_direct_"//band_name, units_str="W m-2", &
               &  long_name="Direct incoming "//band_long_name//" spectral flux at roofs", &
               &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
               &  fill_value=FillValueFlux)
        end if
        call out_file%define_variable("roof_spectral_flux_net_"//band_name, units_str="W m-2", &
             &  long_name="Net "//band_long_name//" spectral flux at roofs", &
             &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
             &  fill_value=FillValueFlux)
        call out_file%define_variable("wall_spectral_flux_in_"//band_name, units_str="W m-2", &
             &  long_name="Incoming "//band_long_name//" spectral flux at walls", &
             &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
             &  fill_value=FillValueFlux)
        if (allocated(flux%wall_in_dir)) then
          call out_file%define_variable("wall_spectral_flux_in_direct_"//band_name, units_str="W m-2", &
               &  long_name="Direct incoming "//band_long_name//" spectral flux at walls", &
               &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
               &  fill_value=FillValueFlux)
        end if
        call out_file%define_variable("wall_spectral_flux_net_"//band_name, units_str="W m-2", &
             &  long_name="Net "//band_long_name//" spectral flux at walls", &
             &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
             &  fill_value=FillValueFlux)
      end if
      if (allocated(flux%clear_air_abs)) then
        call out_file%define_variable("clear_air_spectral_absorption_"//band_name, units_str="W m-2", &
             &  long_name="Absorbed "//band_long_name//" in clear air", &
             &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
             &  fill_value=FillValueFlux)
      end if

      if (allocated(flux%veg_abs)) then
        call out_file%define_variable("veg_spectral_absorption_"//band_name, units_str="W m-2", &
             &  long_name="Absorbed "//band_long_name//" by vegetation", &
             &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
             &  fill_value=FillValueFlux)
        call out_file%define_variable("veg_air_spectral_absorption_"//band_name, units_str="W m-2", &
             &  long_name="Absorbed "//band_long_name//" by air in vegetated regions", &
             &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
             &  fill_value=FillValueFlux)
      end if
      if (allocated(flux%veg_abs_dir)) then
        call out_file%define_variable("veg_spectral_absorption_direct_"//band_name, units_str="W m-2", &
             &  long_name="Absorbed direct "//band_long_name//" by vegetation", &
             &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
             &  fill_value=FillValueFlux)
      end if

      if (allocated(flux%flux_dn_layer_top)) then
        call out_file%define_variable("spectral_flux_dn_layer_top_"//band_name, units_str="W m-2", &
             &  long_name="Downwelling "//band_long_name//" spectral flux at top of layer", &
             &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
             &  fill_value=FillValueFlux)
        if (allocated(flux%flux_dn_dir_layer_top)) then
          call out_file%define_variable("spectral_flux_dn_direct_layer_top_"//band_name, units_str="W m-2", &
               &  long_name="Downwelling direct "//band_long_name//" spectral flux at top of layer", &
               &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
               &  fill_value=FillValueFlux)
        end if
        call out_file%define_variable("spectral_flux_up_layer_top_"//band_name, units_str="W m-2", &
             &  long_name="Upwelling "//band_long_name//" spectral flux at top of layer", &
             &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
             &  fill_value=FillValueFlux)
        call out_file%define_variable("spectral_flux_dn_layer_base_"//band_name, units_str="W m-2", &
             &  long_name="Downwelling "//band_long_name//" spectral flux at base of layer", &
             &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
             &  fill_value=FillValueFlux)
        if (allocated(flux%flux_dn_dir_layer_base)) then
          call out_file%define_variable("spectral_flux_dn_direct_layer_base_"//band_name, units_str="W m-2", &
               &  long_name="Downwelling direct "//band_long_name//" spectral flux at base of layer", &
               &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
               &  fill_value=FillValueFlux)
        end if
        call out_file%define_variable("spectral_flux_up_layer_base_"//band_name, units_str="W m-2", &
             &  long_name="Upwelling "//band_long_name//" spectral flux at base of layer", &
             &  dim3_name="column", dim2_name="layer", dim1_name="band_"//band_name, &
             &  fill_value=FillValueFlux)
      end if
    end if

  end subroutine define_canopy_flux_variables


  subroutine write_canopy_flux_variables(out_file, band_name, nmaxlay, &
       &  nlay, flux, do_broadband, do_spectral)

    use easy_netcdf
    use radsurf_canopy_flux,      only : canopy_flux_type

    type(netcdf_file),      intent(inout) :: out_file
    character(len=*),       intent(in)    :: band_name
    integer(kind=jpim),     intent(in)    :: nmaxlay
    integer(kind=jpim),     intent(in)    :: nlay(:)
    type(canopy_flux_type), intent(in)    :: flux
    logical,                intent(in)    :: do_broadband, do_spectral

    ! Temporary variable at layer interfaces
    real(kind=jprb) :: tmp(nmaxlay, flux%ncol)
    real(kind=jprb) :: tmp_spec(flux%nspec, nmaxlay, flux%ncol)
    
    ! Wavelength-independent quantities
    if (allocated(flux%ground_dn_dir)) then
      call out_file%put("ground_sunlit_fraction", flux%ground_sunlit_frac)
    end if
    if (allocated(flux%roof_in_dir)) then
      call unpack_variable(flux%ncol, nmaxlay, nlay, FillValueFlux, &
           &  flux%roof_sunlit_frac, tmp)
      call out_file%put("roof_sunlit_fraction", tmp)
      call unpack_variable(flux%ncol, nmaxlay, nlay, FillValueFlux, &
           &  flux%wall_sunlit_frac, tmp)
      call out_file%put("wall_sunlit_fraction", tmp)
    end if
    if (allocated(flux%veg_abs_dir)) then
      call unpack_variable(flux%ncol, nmaxlay, nlay, FillValueFlux, &
           &  flux%veg_sunlit_frac, tmp)
      call out_file%put("veg_sunlit_fraction", tmp)
    end if

    ! Broadband fluxes
    if (do_broadband) then
      call out_file%put("ground_flux_dn_"//band_name, sum(flux%ground_dn,1))
      call out_file%put("ground_flux_net_"//band_name, sum(flux%ground_net,1))
      if (allocated(flux%ground_dn_dir)) then
        call out_file%put("ground_flux_dn_direct_"//band_name, &
             &  sum(flux%ground_dn_dir,1))
        call out_file%put("ground_flux_vertical_diffuse_"//band_name, &
             &  sum(flux%ground_vertical_diff,1))
      else
        call out_file%put("ground_flux_vertical_"//band_name, &
             &  sum(flux%ground_vertical_diff,1))
      end if
      call out_file%put("top_flux_dn_"//band_name, sum(flux%top_dn,1))
      call out_file%put("top_flux_net_"//band_name, sum(flux%top_net,1))
      if (allocated(flux%top_dn_dir)) then
        call out_file%put("top_flux_dn_direct_"//band_name, &
             &  sum(flux%top_dn_dir,1))
      end if
      if (allocated(flux%roof_in)) then
        call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
             &  flux%roof_in, tmp)
        call out_file%put("roof_flux_in_"//band_name, tmp)
        call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
             &  flux%roof_net, tmp)
        call out_file%put("roof_flux_net_"//band_name, tmp)
        call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
             &  flux%wall_in, tmp)
        call out_file%put("wall_flux_in_"//band_name, tmp)
        call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
             &  flux%wall_net, tmp)
        call out_file%put("wall_flux_net_"//band_name, tmp)
        if (allocated(flux%roof_in_dir)) then
          call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
               &  flux%roof_in_dir, tmp)
          call out_file%put("roof_flux_in_direct_"//band_name, tmp)
          call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
               &  flux%wall_in_dir, tmp)
          call out_file%put("wall_flux_in_direct_"//band_name, tmp)
        end if
      end if
      if (allocated(flux%clear_air_abs)) then
        call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
             &  flux%clear_air_abs, tmp)
        call out_file%put("clear_air_absorption_"//band_name, tmp)    
      end if
      if (allocated(flux%veg_abs)) then
        call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
             &  flux%veg_abs, tmp)
        call out_file%put("veg_absorption_"//band_name, tmp)
        call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
             &  flux%veg_air_abs, tmp)
        call out_file%put("veg_air_absorption_"//band_name, tmp)
      end if
      if (allocated(flux%veg_abs_dir)) then
        call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
             &  flux%veg_abs_dir, tmp)
        call out_file%put("veg_absorption_direct_"//band_name, tmp)
      end if
      if (allocated(flux%flux_dn_layer_top)) then
        call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
             &  flux%flux_dn_layer_top, tmp)
        call out_file%put("flux_dn_layer_top_"//band_name, tmp)
        if (allocated(flux%flux_dn_dir_layer_top)) then
          call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
               &  flux%flux_dn_dir_layer_top, tmp)
          call out_file%put("flux_dn_direct_layer_top_"//band_name, tmp)
        end if
        call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
             &  flux%flux_up_layer_top, tmp)
        call out_file%put("flux_up_layer_top_"//band_name, tmp)
        call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
             &  flux%flux_dn_layer_base, tmp)
        call out_file%put("flux_dn_layer_base_"//band_name, tmp)
        if (allocated(flux%flux_dn_dir_layer_base)) then
          call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
               &  flux%flux_dn_dir_layer_base, tmp)
          call out_file%put("flux_dn_direct_layer_base_"//band_name, tmp)
        end if
        call unpack_variable_broadband(flux%ncol, nmaxlay, nlay, FillValueFlux, &
             &  flux%flux_up_layer_base, tmp)
        call out_file%put("flux_up_layer_base_"//band_name, tmp)
      end if
    end if

    ! Spectral fluxes
    if (do_spectral) then
      call out_file%put("ground_spectral_flux_dn_"//band_name, flux%ground_dn)
      call out_file%put("ground_spectral_flux_net_"//band_name, flux%ground_net)
      if (allocated(flux%ground_dn_dir)) then
        call out_file%put("ground_spectral_flux_dn_direct_"//band_name, &
             &  flux%ground_dn_dir)
        call out_file%put("ground_spectral_flux_vertical_diffuse_"//band_name, &
             &  flux%ground_vertical_diff)
      else
        call out_file%put("ground_spectral_flux_vertical_"//band_name, &
             &  flux%ground_vertical_diff)
      end if
      call out_file%put("top_spectral_flux_dn_"//band_name, flux%top_dn)
      call out_file%put("top_spectral_flux_net_"//band_name, flux%top_net)
      if (allocated(flux%top_dn_dir)) then
        call out_file%put("top_spectral_flux_dn_direct_"//band_name, &
             &  flux%top_dn_dir)
      end if
      if (allocated(flux%roof_in)) then
        call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
             &  flux%roof_in, tmp_spec)
        call out_file%put("roof_spectral_flux_in_"//band_name, tmp_spec)
        call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
             &  flux%roof_net, tmp_spec)
        call out_file%put("roof_spectral_flux_net_"//band_name, tmp_spec)
        call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
             &  flux%wall_in, tmp_spec)
        call out_file%put("wall_spectral_flux_in_"//band_name, tmp_spec)
        call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
             &  flux%wall_net, tmp_spec)
        call out_file%put("wall_spectral_flux_net_"//band_name, tmp_spec)
        if (allocated(flux%roof_in_dir)) then
          call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
               &  flux%roof_in_dir, tmp_spec)
          call out_file%put("roof_spectral_flux_in_direct_"//band_name, tmp_spec)
          call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
               &  flux%wall_in_dir, tmp_spec)
          call out_file%put("wall_spectral_flux_in_direct_"//band_name, tmp_spec)
        end if
      end if
      if (allocated(flux%clear_air_abs)) then
        call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
             &  flux%clear_air_abs, tmp_spec)
        call out_file%put("clear_air_spectral_absorption_"//band_name, tmp_spec)    
      end if
      if (allocated(flux%veg_abs)) then
        call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
             &  flux%veg_abs, tmp_spec)
        call out_file%put("veg_spectral_absorption_"//band_name, tmp_spec)
        call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
             &  flux%veg_air_abs, tmp_spec)
        call out_file%put("veg_air_spectral_absorption_"//band_name, tmp_spec)
      end if
      if (allocated(flux%veg_abs_dir)) then
        call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
             &  flux%veg_abs_dir, tmp_spec)
        call out_file%put("veg_spectral_absorption_direct_"//band_name, tmp_spec)
      end if
      if (allocated(flux%flux_dn_layer_top)) then
        call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
             &  flux%flux_dn_layer_top, tmp_spec)
        call out_file%put("spectral_flux_dn_layer_top_"//band_name, tmp_spec)
        if (allocated(flux%flux_dn_dir_layer_top)) then
          call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
               &  flux%flux_dn_dir_layer_top, tmp_spec)
          call out_file%put("spectral_flux_dn_direct_layer_top_"//band_name, tmp_spec)
        end if
        call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
             &  flux%flux_up_layer_top, tmp_spec)
        call out_file%put("spectral_flux_up_layer_top_"//band_name, tmp_spec)
        call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
             &  flux%flux_dn_layer_base, tmp_spec)
        call out_file%put("spectral_flux_dn_layer_base_"//band_name, tmp_spec)
        if (allocated(flux%flux_dn_dir_layer_base)) then
          call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
               &  flux%flux_dn_dir_layer_base, tmp_spec)
          call out_file%put("spectral_flux_dn_direct_layer_base_"//band_name, tmp_spec)
        end if
        call unpack_variable_spectral(flux%ncol, nmaxlay, nlay, flux%nspec, FillValueFlux, &
             &  flux%flux_up_layer_base, tmp_spec)
        call out_file%put("spectral_flux_up_layer_base_"//band_name, tmp_spec)
      end if
    end if


  end subroutine write_canopy_flux_variables

  subroutine unpack_variable(ncol, nmaxlay, nlay, fill_value, var_in, var_out)

    integer(kind=jpim), intent(in)  :: ncol, nmaxlay, nlay(:)
    real(kind=jprb),    intent(in)  :: fill_value
    real(kind=jprb),    intent(in)  :: var_in(:)
    real(kind=jprb),    intent(out) :: var_out(nmaxlay,ncol)
    
    integer :: jcol, jlay, itotlay

    var_out = fill_value
    
    itotlay = 1
    
    do jcol = 1,ncol
      do jlay = 1,nlay(jcol)
        var_out(jlay,jcol) = var_in(itotlay)
        itotlay = itotlay + 1
      end do
    end do

  end subroutine unpack_variable

  subroutine unpack_variable_broadband(ncol, nmaxlay, nlay, fill_value, var_in, var_out)

    integer(kind=jpim), intent(in)  :: ncol, nmaxlay, nlay(:)
    real(kind=jprb),    intent(in)  :: fill_value
    real(kind=jprb),    intent(in)  :: var_in(:,:)
    real(kind=jprb),    intent(out) :: var_out(nmaxlay,ncol)
    
    integer :: jcol, jlay, itotlay

    var_out = fill_value
    
    itotlay = 1
    
    do jcol = 1,ncol
      do jlay = 1,nlay(jcol)
        var_out(jlay,jcol) = sum(var_in(:,itotlay),1)
        itotlay = itotlay + 1
      end do
    end do

  end subroutine unpack_variable_broadband

  subroutine unpack_variable_spectral(ncol, nmaxlay, nlay, nspec, fill_value, var_in, var_out)

    integer(kind=jpim), intent(in)  :: ncol, nmaxlay, nlay(:), nspec
    real(kind=jprb),    intent(in)  :: fill_value
    real(kind=jprb),    intent(in)  :: var_in(:,:)
    real(kind=jprb),    intent(out) :: var_out(nspec,nmaxlay,ncol)
    
    integer :: jcol, jlay, itotlay

    var_out = fill_value
    
    itotlay = 1
    
    do jcol = 1,ncol
      do jlay = 1,nlay(jcol)
        var_out(:,jlay,jcol) = var_in(:,itotlay)
        itotlay = itotlay + 1
      end do
    end do

  end subroutine unpack_variable_spectral


end module radsurf_save
