SOURCES = \
	radsurf_config.F90 \
	radsurf_canopy_properties.F90 \
	radsurf_facet_properties.F90 \
	radsurf_volume_properties.F90 \
	radsurf_boundary_conds_out.F90 \
	radsurf_canopy_flux.F90 \
	radsurf_interface.F90 \
	radsurf_forest_sw.F90 \
	radsurf_overlap.F90

OBJECTS := $(SOURCES:.F90=.o)
LIBSURF = ../lib/libradsurf.a

all: $(LIBSURF)

$(LIBSURF): $(OBJECTS)
	ar r $(LIBSURF) $(OBJECTS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(LIBSURF)

radsurf_canopy_properties.o: radsurf_config.o
radsurf_facet_properties.o: radsurf_canopy_properties.o radsurf_config.o
radsurf_volume_properties.o: radsurf_canopy_properties.o radsurf_config.o
radsurf_interface.o: radsurf_canopy_properties.o radsurf_facet_properties.o \
	radsurf_volume_properties.o radsurf_boundary_conds_out.o \
	radsurf_canopy_flux.o radsurf_config.o radsurf_forest_sw.o
radsurf_forest_sw.o: radsurf_config.o radsurf_canopy_properties.o \
	radsurf_volume_properties.o radsurf_boundary_conds_out.o \
	radsurf_canopy_flux.o radsurf_overlap.o
