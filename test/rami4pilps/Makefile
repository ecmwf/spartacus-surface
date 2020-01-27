SPSURF = ../../bin/spartacus_surface
CONFIG = config.nam
CHANGE_NL = ./change_namelist.sh
BASENAME = rami4pilps
BASE_PROFILE = $(BASENAME)_base_profile.nc
INPUT = $(BASENAME).nc
BASESCEN=vis-snw-0.3

test: test_method

test_method: $(BASESCEN)-1-1_out $(BASESCEN)-1-2_out \
	     $(BASESCEN)-2-1_out $(BASESCEN)-2-2_out 

test_scenes: vis-med-0.1_out vis-med-0.3_out vis-med-0.5_out \
	     vis-snw-0.1_out vis-snw-0.3_out vis-snw-0.5_out \
	     nir-med-0.1_out nir-med-0.3_out nir-med-0.5_out \
	     nir-snw-0.1_out nir-snw-0.3_out nir-snw-0.5_out

$(INPUT): $(BASE_PROFILE)
	./duplicate_profiles.sh $(BASE_PROFILE) $(INPUT)

%_out: $(INPUT)
	echo Computing $@
	$(CHANGE_NL) $(CONFIG) config_tmp.nam $(shell ./expand_args.sh $@)
	$(SPSURF) config_tmp.nam $(INPUT) $(BASENAME)_$@.nc | tee $(BASENAME)_$@.log

.PHONY: *_out $(BASENAME) test*

clean:
	rm -f $(BASENAME)_*_out.nc $(BASENAME)_*_out.log config_tmp.nam
