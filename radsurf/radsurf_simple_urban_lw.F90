! radsurf_simple_urban_lw.F90 - Longwave solver for unvegetated single-layer urban canopy
!
! (C) Copyright 2021- ECMWF.
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

module radsurf_simple_urban_lw

contains

  ! --------------------------------------------------------
  ! Single-layer solar urban radiative transfer (i.e. the assumption
  ! that all buildings are the same height) using the Harman et
  ! al. (BLM 2004) method of solving a 2x2 matrix problem, but with
  ! the option of two different models for urban geometry: the
  ! original "infinite street" of constant width, and the more recent
  ! "exponential model" for the distribution of wall-to-wall
  ! separation distances
  subroutine simple_urban_lw(config, is_infinite_street, &
       &  nlw, icol, ilay, &
       &  canopy_props, lw_spectral_props, &
       &  top_emissivity, top_emission, &
       &  lw_internal, lw_norm)
    
    use parkind1,                   only : jpim, jprb
    use yomhook,                    only : lhook, dr_hook
    use radiation_io,               only : radiation_abort
    use radsurf_config,             only : config_type
    use radsurf_canopy_properties,  only : canopy_properties_type
    use radsurf_lw_spectral_properties,  only : lw_spectral_properties_type
    use radsurf_canopy_flux,        only : canopy_flux_type
    use radiation_constants,        only : Pi
    use radsurf_norm_perim,         only : calc_norm_perim_urban
    use radsurf_view_factor,        only : calc_view_factors_inf, &
         &                                 calc_view_factors_exp
    use radtool_matrix,             only : solve_vec

    implicit none

    ! Inputs

    ! Algorithm configuration
    type(config_type),             intent(in)  :: config
    ! Do we use the infinite-street or exponential model?
    logical,                       intent(in)  :: is_infinite_street
    ! Number of spectral intervals
    integer(kind=jpim),            intent(in)  :: nlw
    ! Index of current column and layer
    integer(kind=jpim),            intent(in)  :: icol, ilay
    ! Geometric and other spectrally independent properties of the canopy
    type(canopy_properties_type),  intent(in)  :: canopy_props
    ! Spectral properties of the air, vegetation and urban facets
    type(lw_spectral_properties_type),  intent(in)  :: lw_spectral_props

    ! Outputs

    ! Top-of-canopy spectral emissivity and emission (W m-2)
    real(kind=jprb), dimension(nlw),intent(out):: top_emissivity, &
         &                                        top_emission
    ! Flux outputs
    type(canopy_flux_type), intent(inout), optional &
         &  :: lw_internal, & ! LW fluxes from internal emission
         &     lw_norm        ! LW fluxes normalized by top-of-canopy diffuse

    ! Local variables

    ! view_A_B is the fraction of radiaion emanating from facet A that
    ! intercepts facet B
    real(kind=jprb) :: view_ground_sky, view_wall_wall
    real(kind=jprb) :: view_wall_ground, view_ground_wall

    ! Normalized perimeter length between regions (unused), and
    ! between air and walls (m-1)
    real(jprb) :: norm_perim(1), norm_perim_wall(1)

    ! Dummy variables for calc_norm_perim_urban
    real(kind=jprb) :: veg_fraction(1), veg_scale(1), veg_contact_fraction(1)

    ! The "X" in Hogan (BLM 2019a, exponential), metres
    real(kind=jprb) :: building_separation_scale

    ! The street width in the infinite-street assumption, metres
    real(kind=jprb) :: street_width

    ! Fundamentally the Harman et al. (BLM 2004) method solves a 2x2
    ! matrix problem of the form
    ! interaction_matrix*solution_vector=source_vector, which here
    ! includes an additional dimension for the number of spectral
    ! intervals
    real(kind=jprb) :: interaction_matrix(nlw,2,2)
    real(kind=jprb) :: solution_vector(nlw,2), source_vector(nlw,2)

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radsurf_simple_urban_lw:simple_urban_lw',0,hook_handle)

    associate( &
         &  dz                 => canopy_props%dz(ilay), &
         &  building_fraction  => canopy_props%building_fraction(ilay), &
         &  building_scale     => canopy_props%building_scale(ilay), &
         &  ground_emissivity => lw_spectral_props%ground_emissivity(:,icol), &
         &  ground_emission   => lw_spectral_props%ground_emission(:,icol), &
         &  roof_emissivity    => lw_spectral_props%roof_emissivity(:,ilay), &
         &  roof_emission      => lw_spectral_props%roof_emission(:,ilay), &
         &  wall_emissivity    => lw_spectral_props%wall_emissivity(:,ilay), &
         &  wall_emission      => lw_spectral_props%wall_emission(:,ilay) )
      
      ! Compute normalized perimeter length of walls, noting that
      ! calc_norm_perim_urban accepts vectors of inputs but we only
      ! want to compute a single value, and dummy values are entered
      ! for the unused vegetation variables
      veg_fraction = 0.0_jprb
      veg_scale = 1.0_jprb
      veg_contact_fraction = 0.0_jprb
      call calc_norm_perim_urban(config, 1, 1, spread(building_fraction,1,1), &
           &  spread(building_scale,1,1), veg_fraction, veg_scale, &
           &  veg_contact_fraction, norm_perim, norm_perim_wall)

      ! Compute length scales and view factors
      if (is_infinite_street) then
        ! Hogan (BLM 2019b, radiative transfer), Eq. 7
        street_width = 2.0_jprb * (1.0_jprb - building_fraction) / norm_perim_wall(1)
        call calc_view_factors_inf(dz / street_width, &
             &  view_ground_sky, view_wall_wall);
      else
        ! Hogan (BLM 2019b, radiative transfer), Eq. 8
        building_separation_scale = Pi * (1.0_jprb - building_fraction) &
             &                    / norm_perim_wall(1)
        call calc_view_factors_exp(dz / building_separation_scale, &
             &  view_ground_sky, view_wall_wall);
      end if

      ! Compute extra view factors
      view_wall_ground = 0.5_jprb * (1.0_jprb - view_wall_wall)
      view_ground_wall = 1.0_jprb - view_ground_sky

      ! Set to the flux components to zero initially
      call lw_norm%zero(icol, ilay, ilay)
      call lw_internal%zero( icol, ilay, ilay)

      ! First the fluxes due to internal emission

      ! The elements of the interaction matrix are common for direct
      ! and diffuse input fluxes
      interaction_matrix(:,1,1) = 1.0_jprb
      interaction_matrix(:,1,2) = -view_wall_ground*(1.0_jprb-wall_emissivity)
      interaction_matrix(:,2,1) = -view_ground_wall*(1.0_jprb-ground_emissivity)
      interaction_matrix(:,2,2) = 1.0_jprb - view_wall_wall*(1.0_jprb-ground_emissivity)

      ! Incoming radiation at ground and walls due to emission from
      ! the other facet
      source_vector(:,1) = view_wall_ground * wall_emission * norm_perim_wall(1) * dz
      source_vector(:,2) = view_ground_wall * ground_emission * (1.0_jprb-building_fraction) &
           &             + view_wall_wall * wall_emission * norm_perim_wall(1) * dz

      ! Solve 2x2 matrix problem
      solution_vector = solve_vec(nlw,nlw,2,interaction_matrix,source_vector)

      ! Ground fluxes
      lw_internal%ground_dn(:,icol) = solution_vector(:,1)
      lw_internal%ground_net(:,icol) &
           &  = solution_vector(:,1) * ground_emissivity & 
           &  - ground_emission * (1.0_jprb - building_fraction)
      ! vertical flux???

      ! Roof fluxes
      lw_internal%roof_in(:,ilay)  = 0.0_jprb
      lw_internal%roof_net(:,ilay) = - building_fraction * roof_emission

      ! Wall fluxes
      lw_internal%wall_in(:,ilay) = solution_vector(:,2)
      lw_internal%wall_net(:,ilay) &
           &  = solution_vector(:,2) * wall_emissivity &
           &  - wall_emission * norm_perim_wall(1) * dz

      ! Top-of-canopy fluxes
      lw_internal%top_dn(:,icol)     = 0.0_jprb
      lw_internal%top_net(:,icol)    = - building_fraction*roof_emission &
           &  - (lw_internal%ground_dn(:,ilay) - lw_internal%ground_net(:,ilay)) &
           &    * view_ground_sky &
           &  - (lw_internal%wall_in(:,ilay) - lw_internal%wall_net(:,ilay)) &
           &    * view_wall_ground

      ! Flux "profiles"
      if (allocated(lw_internal%flux_dn_layer_top)) then
        lw_internal%flux_dn_layer_top(:,ilay) = 0.0_jprb
        lw_internal%flux_up_layer_top(:,ilay) &
             &  = (lw_internal%ground_dn(:,ilay) - lw_internal%ground_net(:,ilay)) &
             &    * view_ground_sky &
             &  + (lw_internal%wall_in(:,ilay) - lw_internal%wall_net(:,ilay)) &
             &    * view_wall_ground
        lw_internal%flux_dn_layer_base(:,ilay) = lw_internal%ground_dn(:,ilay)
        lw_internal%flux_up_layer_base(:,ilay) &
             &  = lw_internal%ground_dn(:,ilay) - lw_internal%ground_net(:,ilay)
      end if

      ! Second the fluxes normalized by the diffuse downwelling flux
      ! at canopy top.

      ! Incoming radiation at ground and walls due to incoming diffuse
      ! radiation at top-of-canopy
      source_vector(:,1) = view_ground_sky * (1.0_jprb-building_fraction)
      source_vector(:,2) = view_ground_wall* (1.0_jprb-building_fraction)

      ! Solve 2x2 matrix problem
      solution_vector = solve_vec(nlw,nlw,2,interaction_matrix,source_vector)

      ! Ground fluxes
      lw_norm%ground_dn(:,ilay)     = solution_vector(:,1)
      lw_norm%ground_net(:,ilay) &
           &  = lw_norm%ground_dn(:,ilay) * ground_emissivity
      ! vertical flux???

      ! Roof fluxes
      lw_norm%roof_in(:,ilay)  = building_fraction
      lw_norm%roof_net(:,ilay) = building_fraction * roof_emissivity

      ! Wall fluxes
      lw_norm%wall_in(:,ilay) = solution_vector(:,2)
      lw_norm%wall_net(:,ilay) &
           &  =  lw_norm%wall_in(:,ilay) * wall_emissivity

      ! Top-of-canopy fluxes
      lw_norm%top_dn(:,icol)     = 1.0_jprb
      lw_norm%top_net(:,icol)    = 1.0_jprb - building_fraction*(1.0_jprb-roof_emissivity) &
           &  - (lw_norm%ground_dn(:,ilay) - lw_norm%ground_net(:,ilay)) &
           &    * view_ground_sky &
           &  - (lw_norm%wall_in(:,ilay) - lw_norm%wall_net(:,ilay)) &
           &    * view_wall_ground

      ! Flux "profiles"
      if (allocated(lw_norm%flux_dn_layer_top)) then
        lw_norm%flux_dn_layer_top(:,ilay) = 1.0_jprb - building_fraction
        lw_norm%flux_up_layer_top(:,ilay) &
             &  = (lw_norm%ground_dn(:,ilay) - lw_norm%ground_net(:,ilay)) &
             &    * view_ground_sky &
             &  + (lw_norm%wall_in(:,ilay) - lw_norm%wall_net(:,ilay)) &
             &    * view_wall_ground
        lw_norm%flux_dn_layer_base(:,ilay) = lw_norm%ground_dn(:,ilay)
        lw_norm%flux_up_layer_base(:,ilay) &
             &  = lw_norm%ground_dn(:,ilay) - lw_norm%ground_net(:,ilay)
      end if

    end associate

    if (lhook) call dr_hook('radsurf_simple_urban_lw:simple_urban_lw',1,hook_handle)

  end subroutine simple_urban_lw


end module radsurf_simple_urban_lw
