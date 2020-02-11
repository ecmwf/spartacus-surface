# Makefile to test SPARTACUS-Surface on the RAMI4PILPS vegetation
# cases

SPSURF = ../../bin/spartacus_surface
CONFIG = config.nam
CHANGE_NL = ./change_namelist.sh
BASENAME = rami4pilps
BASE_PROFILE = $(BASENAME)_base_profile.nc
INPUT = $(BASENAME).nc

# The command "make BAND-SURFACE-COVER[-NVEGREG[-NSTREAM]]" will run a
# SPARTACUS-Surface on a forest scene defined by the target, where
# BAND is the spectral region and may be "vis" or "nir", SURFACE is
# the surface type and may be "med" (bare soil) or "snw" (snow
# covered), COVER is the vegetation cover between 0 and 1, the values
# used in RAMI4PILPS being 0.1, 0.3 and 0.5, NVEGREG is the optional
# number of regions used to describe the vegetated part of the canopy
# and may be 1 or 2, while NSTREAM is the optional number of diffuse
# streams to use in each hemisphere and may be any integer 1 or
# larger.  The default values for NSTREAM and NVEGREG are given in
# config.nam.

# The default scene used by "make test_method"
BASESCEN=vis-snw-0.3
#BASESCEN=nir-med-0.3

test: test_method test_scenes

# Test the SPARTACUS-Surface method for different numbers of streams
# and regions, using the BASESCEN scene
test_method: $(BASESCEN)-1-1_out $(BASESCEN)-1-2_out $(BASESCEN)-1-4_out $(BASESCEN)-1-8_out \
	     $(BASESCEN)-2-1_out $(BASESCEN)-2-2_out $(BASESCEN)-2-4_out $(BASESCEN)-2-8_out 

# Test the various RAMI4PILPS Open Forest Canopy scenarios (with the
# default number of streams and regions)
test_scenes: vis-med-0.1_out vis-med-0.3_out vis-med-0.5_out \
	     vis-snw-0.1_out vis-snw-0.3_out vis-snw-0.5_out \
	     nir-med-0.1_out nir-med-0.3_out nir-med-0.5_out \
	     nir-snw-0.1_out nir-snw-0.3_out nir-snw-0.5_out

# Duplicate the profiles in rami4pilps_base_profile.nc to add a range
# of solar zenith angles, storing the result in rami4pilps.nc
$(INPUT): $(BASE_PROFILE)
	./duplicate_profiles.sh $(BASE_PROFILE) $(INPUT)

# Use the expand_args.sh shell script to convert the make target into
# namelist arguments, which are then inserted into
# config_tmp.nam. SPARTACUS-Surface is then run with this namelist
# file.
%_out: $(INPUT)
	echo Computing $@
	$(CHANGE_NL) $(CONFIG) config_tmp.nam $(shell ./expand_args.sh $@)
	$(SPSURF) config_tmp.nam $(INPUT) $(BASENAME)_$@.nc 

.PHONY: *_out $(BASENAME) test*

clean:
	rm -f $(BASENAME).nc $(BASENAME)_*_out.nc $(BASENAME)_*_out.log config_tmp.nam
