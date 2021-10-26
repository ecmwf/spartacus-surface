! radsurf_simple_urban_sw.F90 - Shortwave solver for unvegetated single-layer urban canopy
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

module radsurf_simple_urban_sw

contains

  ! --------------------------------------------------------
  ! Single-layer solar urban radiative transfer (i.e. the assumption
  ! that all buildings are the same height) using the Harman et
  ! al. (BLM 2004) method of solving a 2x2 matrix problem, but with
  ! the option of two different models for urban geometry: the
  ! original "infinite street" of constant width, and the more recent
  ! "exponential model" for the distribution of wall-to-wall
  ! separation distances
  subroutine simple_urban_sw(config, is_infinite_street, &
       &  nsw, icol, ilay, cos_sza, &
       &  canopy_props, sw_spectral_props, &
       &  ground_albedo_diff, ground_albedo_dir, &
       &  top_albedo_diff, top_albedo_dir, &
       &  sw_norm_dir, sw_norm_diff)
    
    use parkind1,                   only : jpim, jprb
    use yomhook,                    only : lhook, dr_hook
    use radiation_io,               only : radiation_abort
    use radsurf_config,             only : config_type
    use radsurf_canopy_properties,  only : canopy_properties_type
    use radsurf_sw_spectral_properties,  only : sw_spectral_properties_type
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
    integer(kind=jpim),            intent(in)  :: nsw
    ! Index of current column and layer
    integer(kind=jpim),            intent(in)  :: icol, ilay
    ! Cosine of the solar zenith angle
    real(kind=jprb),               intent(in)  :: cos_sza
    ! Geometric and other spectrally independent properties of the canopy
    type(canopy_properties_type),  intent(in)  :: canopy_props
    ! Spectral properties of the air, vegetation and urban facets
    type(sw_spectral_properties_type),  intent(in)  :: sw_spectral_props
    ! Spectral albedo of the ground to diffuse and direct radiation
    real(kind=jprb), dimension(nsw),intent(in) :: ground_albedo_diff, &
         &                                        ground_albedo_dir

    ! Outputs

    ! Top-of-canopy spectral albedo to diffuse and direct radiation
    real(kind=jprb), dimension(nsw),intent(out):: top_albedo_diff, &
         &                                        top_albedo_dir
    ! Flux outputs
    type(canopy_flux_type), intent(inout), optional &
         &  :: sw_norm_dir, &  ! SW fluxes normalized by top-of-canopy direct
         &     sw_norm_diff    ! SW fluxes normalized by top-of-canopy diffuse

    ! Local variables

    ! view_A_B is the fraction of radiaion emanating from facet A that
    ! intercepts facet B, where "dir" is direct radiation emanating
    ! from the sky and all other "A" facets refer to diffuse radiation
    real(kind=jprb) :: view_dir_ground, view_dir_wall
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

    ! Tangent of solar zenith angle
    real(kind=jprb) :: tan_sza

    ! Fundamentally the Harman et al. (BLM 2004) method solves a 2x2
    ! matrix problem of the form
    ! interaction_matrix*solution_vector=source_vector, which here
    ! includes an additional dimension for the number of spectral
    ! intervals
    real(kind=jprb) :: interaction_matrix(nsw,2,2)
    real(kind=jprb) :: solution_vector(nsw,2), source_vector(nsw,2)

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radsurf_simple_urban_sw:simple_urban_sw',0,hook_handle)

    associate( &
         &  dz                 => canopy_props%dz(ilay), &
         &  building_fraction  => canopy_props%building_fraction(ilay), &
         &  building_scale     => canopy_props%building_scale(ilay), &
         &  roof_albedo        => sw_spectral_props%roof_albedo(:,ilay), &
         &  wall_albedo        => sw_spectral_props%wall_albedo(:,ilay), &
         &  wall_specular_frac => sw_spectral_props%wall_specular_frac(:,ilay))
      
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
             &  view_ground_sky, view_wall_wall, &
             &  cos_sza=cos_sza, view_dir_ground=view_dir_ground);
      else
        ! Hogan (BLM 2019b, radiative transfer), Eq. 8
        building_separation_scale = Pi * (1.0_jprb - building_fraction) &
             &                    / norm_perim_wall(1)
        call calc_view_factors_exp(dz / building_separation_scale, &
             &  view_ground_sky, view_wall_wall, &
             &  cos_sza=cos_sza, view_dir_ground=view_dir_ground);
      end if

      ! Compute extra view factors
      view_dir_wall = 1.0_jprb - view_dir_ground
      view_wall_ground = 0.5_jprb * (1.0_jprb - view_wall_wall)
      view_ground_wall = 1.0_jprb - view_ground_sky

      ! Set to the flux components to zero initially
      call sw_norm_diff%zero(icol, ilay, ilay)
      call sw_norm_dir%zero( icol, ilay, ilay)

      ! First the fluxes normalized by the direct downwelling flux at
      ! canopy top.

      ! The elements of the interaction matrix are common for direct
      ! and diffuse input fluxes
      interaction_matrix(:,1,1) = 1.0_jprb
      interaction_matrix(:,1,2) = -view_wall_ground*wall_albedo
      interaction_matrix(:,2,1) = -view_ground_wall*ground_albedo_diff
      interaction_matrix(:,2,2) = 1.0_jprb - view_wall_wall*wall_albedo

      ! Incoming radiation at ground and walls due to incoming direct
      ! radiation at top-of-canopy. Note that because the ground
      ! albedo is separated into diffuse and direct parts, but the
      ! wall albedo is not, we cannot solve for the sum of these terms
      ! at the ground. Rather, solution_vector(:,1) is the diffuse
      ! downwelling flux at the ground, while solution_vector(:,2) is
      ! the total flux (direct plus diffuse) into the walls. Thus, the
      ! source into the walls (source_vector(:,2)) contains the direct
      ! flux into the walls plus the direct flux incident on the
      ! ground and scattered once into the walls.
      source_vector(:,1) = 0.0_jprb
      source_vector(:,2) = (view_dir_wall + ground_albedo_dir*view_dir_ground*view_ground_wall) &
           * (1.0_jprb-building_fraction)

      ! Solve 2x2 matrix problem
      solution_vector = solve_vec(nsw,nsw,2,interaction_matrix,source_vector)

      ! Ground fluxes
      sw_norm_dir%ground_dn_dir(:,icol) = view_dir_ground &
           * (1.0_jprb-building_fraction)
      sw_norm_dir%ground_dn(:,icol) = sw_norm_dir%ground_dn_dir(:,icol) + solution_vector(:,1)
      sw_norm_dir%ground_net(:,icol) &
           &  = sw_norm_dir%ground_dn_dir(:,ilay) * (1.0_jprb-ground_albedo_dir) &
           &       + solution_vector(:,1) * (1.0_jprb-ground_albedo_diff)
      sw_norm_dir%ground_sunlit_frac(icol) = view_dir_ground
      ! vertical flux???

      ! Roof fluxes
      sw_norm_dir%roof_in_dir(:,ilay) = building_fraction
      sw_norm_dir%roof_in(:,ilay)     = building_fraction
      sw_norm_dir%roof_net(:,ilay)    = building_fraction * (1.0_jprb-roof_albedo)
      sw_norm_dir%roof_sunlit_frac(ilay) = 1.0_jprb

      ! Wall fluxes
      sw_norm_dir%wall_in_dir(:,ilay) = view_dir_wall &
           * (1.0_jprb-building_fraction)
      sw_norm_dir%wall_in(:,ilay) = solution_vector(:,2)
      sw_norm_dir%wall_net(:,ilay) &
           &  =  sw_norm_dir%wall_in(:,ilay) * (1.0_jprb-wall_albedo)
      tan_sza = sqrt(1.0_jprb / (cos_sza*cos_sza) - 1.0_jprb)
      sw_norm_dir%wall_sunlit_frac(ilay) = 0.5_jprb * view_dir_wall &
           &  / (max(tan_sza,1.0e-6_jprb) * norm_perim_wall(1)*dz &
           &     / (Pi*(1.0_jprb - building_fraction)))

      ! Top-of-canopy fluxes
      sw_norm_dir%top_dn_dir(:,icol) = 1.0_jprb
      sw_norm_dir%top_dn(:,icol)     = 1.0_jprb
      sw_norm_dir%top_net(:,icol)    = 1.0_jprb - building_fraction*roof_albedo &
           &  - (sw_norm_dir%ground_dn(:,ilay) - sw_norm_dir%ground_net(:,ilay)) &
           &    * view_ground_sky &
           &  - (sw_norm_dir%wall_in(:,ilay) - sw_norm_dir%wall_net(:,ilay)) &
           &    * view_wall_ground

      ! Flux "profiles"
      if (allocated(sw_norm_dir%flux_dn_layer_top)) then
        sw_norm_dir%flux_dn_dir_layer_top(:,ilay) = (1.0_jprb-building_fraction)
        sw_norm_dir%flux_dn_layer_top(:,ilay) = (1.0_jprb-building_fraction)
        sw_norm_dir%flux_up_layer_top(:,ilay) &
             &  = (sw_norm_dir%ground_dn(:,ilay) - sw_norm_dir%ground_net(:,ilay)) &
             &    * view_ground_sky &
             &  + (sw_norm_dir%wall_in(:,ilay) - sw_norm_dir%wall_net(:,ilay)) &
             &    * view_wall_ground
        sw_norm_dir%flux_dn_dir_layer_base(:,ilay) = sw_norm_dir%ground_dn_dir(:,ilay)
        sw_norm_dir%flux_dn_layer_base(:,ilay) = sw_norm_dir%ground_dn(:,ilay)
        sw_norm_dir%flux_up_layer_base(:,ilay) &
             &  = sw_norm_dir%ground_dn(:,ilay) - sw_norm_dir%ground_net(:,ilay)
      end if

      ! Second the fluxes normalized by the diffuse downwelling flux
      ! at canopy top.

      ! Incoming radiation at ground and walls due to incoming diffuse
      ! radiation at top-of-canopy.  This time we do not need to
      ! separate the incoming and scattered radiation at the surface,
      ! so solution_vector(:,1) is the total flux into the ground.
      source_vector(:,1) = view_ground_sky * (1.0_jprb-building_fraction)
      source_vector(:,2) = view_ground_wall* (1.0_jprb-building_fraction)

      ! Solve 2x2 matrix problem
      solution_vector = solve_vec(nsw,nsw,2,interaction_matrix,source_vector)

      ! Ground fluxes
      sw_norm_diff%ground_dn_dir(:,ilay) = 0.0_jprb
      sw_norm_diff%ground_dn(:,ilay)     = solution_vector(:,1)
      sw_norm_diff%ground_net(:,ilay) &
           &  = sw_norm_diff%ground_dn(:,ilay) * (1.0_jprb-ground_albedo_diff)
      ! vertical flux???

      ! Roof fluxes
      sw_norm_diff%roof_in(:,ilay)  = building_fraction
      sw_norm_diff%roof_net(:,ilay) = building_fraction * (1.0_jprb-roof_albedo)

      ! Wall fluxes
      sw_norm_diff%wall_in(:,ilay) = solution_vector(:,2)
      sw_norm_diff%wall_net(:,ilay) &
           &  =  sw_norm_diff%wall_in(:,ilay) * (1.0_jprb-wall_albedo)

      ! Top-of-canopy fluxes
      sw_norm_diff%top_dn_dir(:,icol) = 0.0_jprb
      sw_norm_diff%top_dn(:,icol)     = 1.0_jprb
      sw_norm_diff%top_net(:,icol)    = 1.0_jprb - building_fraction*roof_albedo &
           &  - (sw_norm_diff%ground_dn(:,ilay) - sw_norm_diff%ground_net(:,ilay)) &
           &    * view_ground_sky &
           &  - (sw_norm_diff%wall_in(:,ilay) - sw_norm_diff%wall_net(:,ilay)) &
           &    * view_wall_ground

      ! Flux "profiles"
      if (allocated(sw_norm_diff%flux_dn_layer_top)) then
        sw_norm_diff%flux_dn_layer_top(:,ilay) = (1.0_jprb-building_fraction)
        sw_norm_diff%flux_up_layer_top(:,ilay) &
             &  = (sw_norm_diff%ground_dn(:,ilay) - sw_norm_diff%ground_net(:,ilay)) &
             &    * view_ground_sky &
             &  + (sw_norm_diff%wall_in(:,ilay) - sw_norm_diff%wall_net(:,ilay)) &
             &    * view_wall_ground
        sw_norm_diff%flux_dn_layer_base(:,ilay) = sw_norm_diff%ground_dn(:,ilay)
        sw_norm_diff%flux_up_layer_base(:,ilay) &
             &  = sw_norm_diff%ground_dn(:,ilay) - sw_norm_diff%ground_net(:,ilay)
      end if

    end associate

    if (lhook) call dr_hook('radsurf_simple_urban_sw:simple_urban_sw',1,hook_handle)

  end subroutine simple_urban_sw


end module radsurf_simple_urban_sw
