! radsurf_simple_spectrum.f90 - Simplest calculation of spectral properties
!
! Copyright (C) 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radsurf_simple_spectrum

contains
  
  subroutine calc_simple_spectrum_lw(config, canopy_props, facet_props, &
       &  volume_props, istartcol, iendcol)
    
    use parkind1,                   only : jpim, jprb
    use yomhook,                    only : lhook, dr_hook
    use radiation_io,               only : radiation_abort
    use radsurf_config,             only : config_type
    use radsurf_canopy_properties,  only : canopy_properties_type
    use radsurf_facet_properties,   only : facet_properties_type
    use radsurf_volume_properties,  only : volume_properties_type
    use radiation_constants,        only : StefanBoltzmann
    
    implicit none
    
    type(config_type),             intent(in) :: config
    type(canopy_properties_type),  intent(in) :: canopy_props
    type(facet_properties_type),   intent(inout) :: facet_props
    type(volume_properties_type),  intent(inout) :: volume_props
    integer(kind=jpim),            intent(in) :: istartcol, iendcol
    
    ! Start and end layers for current column
    integer(kind=jpim) :: ilay1, ilay2

    ilay1 = canopy_props%istartlay(istartcol)
    ilay2 = canopy_props%istartlay(iendcol) + canopy_props%nlay(iendcol) - 1

    if (facet_props%nlw > 1) then
      call radiation_abort('Simple longwave spectrum only possible with one input spectral interval')
    end if
    
    facet_props%ground_lw_emission(1,istartcol:iendcol) = StefanBoltzmann &
         &  * facet_props%ground_lw_emissivity(1,istartcol:iendcol) &
         &  * canopy_props%ground_temperature(istartcol:iendcol) ** 4
    if (ilay2 >= ilay1) then
      facet_props%roof_lw_emission(1,ilay1:ilay2) = StefanBoltzmann &
           &  * facet_props%roof_lw_emissivity(1,ilay1:ilay2) &
           &  * canopy_props%roof_temperature(ilay1:ilay2) ** 4
      facet_props%wall_lw_emission(1,ilay1:ilay2) = StefanBoltzmann &
           &  * facet_props%wall_lw_emissivity(1,ilay1:ilay2) &
           &  * canopy_props%wall_temperature(ilay1:ilay2) ** 4
      volume_props%air_lw_planck(1,ilay1:ilay2) = StefanBoltzmann &
           &  * canopy_props%air_temperature(ilay1:ilay2) ** 4
      if (config%do_vegetation) then
        volume_props%veg_lw_planck(1,ilay1:ilay2) = StefanBoltzmann &
             &  * canopy_props%veg_temperature(ilay1:ilay2) ** 4
      end if
    end if

  end subroutine calc_simple_spectrum_lw

end module radsurf_simple_spectrum
