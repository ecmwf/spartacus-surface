SOURCES = \
	radsurf_config.F90 \
	radsurf_canopy_properties.F90 \
	radsurf_boundary_conds_out.F90 \
	radsurf_canopy_flux.F90 \
	radsurf_interface.F90 \
	radsurf_forest_sw.F90 \
	radsurf_forest_lw.F90 \
	radsurf_urban_sw.F90 \
	radsurf_urban_lw.F90 \
	radsurf_overlap.F90 \
	radsurf_save.F90 \
	radsurf_simple_spectrum.F90 \
	radsurf_sw_spectral_properties.F90 \
	radsurf_lw_spectral_properties.F90 \
	radsurf_norm_perim.F90 \
	radsurf_view_factor.F90 \
	radsurf_simple_urban_sw.F90

OBJECTS := $(SOURCES:.F90=.o)
LIBSURF = ../lib/libradsurf.a

all: $(LIBSURF)

$(LIBSURF): $(OBJECTS)
	ar r $(LIBSURF) $(OBJECTS)

%.o: %.F90 ../mod/radtool_*.mod
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(LIBSURF)

radsurf_canopy_properties.o radsurf_norm_perim.o: radsurf_config.o
radsurf_lw_spectral_properties.o: radsurf_canopy_properties.o radsurf_config.o
radsurf_sw_spectral_properties.o: radsurf_canopy_properties.o radsurf_config.o
radsurf_interface.o: radsurf_canopy_properties.o \
	radsurf_sw_spectral_properties.o radsurf_lw_spectral_properties.o \
	radsurf_boundary_conds_out.o radsurf_simple_urban_sw.o \
	radsurf_canopy_flux.o radsurf_config.o radsurf_forest_sw.o \
	radsurf_forest_lw.o radsurf_urban_sw.o radsurf_urban_lw.o 
radsurf_forest_sw.o: radsurf_config.o radsurf_canopy_properties.o \
	radsurf_sw_spectral_properties.o radsurf_boundary_conds_out.o \
	radsurf_canopy_flux.o radsurf_overlap.o radsurf_norm_perim.o
radsurf_forest_lw.o: radsurf_config.o radsurf_canopy_properties.o \
	radsurf_boundary_conds_out.o radsurf_norm_perim.o \
	radsurf_canopy_flux.o radsurf_overlap.o radsurf_lw_spectral_properties.o
radsurf_urban_sw.o: radsurf_config.o radsurf_canopy_properties.o \
	radsurf_boundary_conds_out.o radsurf_norm_perim.o \
	radsurf_canopy_flux.o radsurf_overlap.o radsurf_sw_spectral_properties.o
radsurf_urban_lw.o: radsurf_config.o radsurf_canopy_properties.o \
	radsurf_boundary_conds_out.o radsurf_norm_perim.o \
	radsurf_canopy_flux.o radsurf_overlap.o radsurf_lw_spectral_properties.o
radsurf_save.o: radsurf_config.o radsurf_canopy_properties.o \
	radsurf_canopy_flux.o
radsurf_simple_spectrum.o: radsurf_config.o radsurf_canopy_properties.o \
	radsurf_lw_spectral_properties.o 
radsurf_canopy_flux.o: radsurf_canopy_properties.o
radsurf_simple_urban_sw.o: radsurf_config.o radsurf_canopy_properties.o \
	radsurf_sw_spectral_properties.o radsurf_canopy_flux.o \
	radsurf_norm_perim.o radsurf_view_factor.o
