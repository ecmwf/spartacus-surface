# Makefile to test SPARTACUS-Surface on the RAMI4PILPS vegetation
# cases

SPSURF = ../../bin/spartacus_surface
CONFIG = config.nam
CONFIG_TMP = config_tmp.nam
CHANGE_NL = ./change_namelist.sh
BASENAME = russell_square
BASE_PROFILE = $(BASENAME).nc
INPUT = $(BASENAME)_sza.nc

test: test_streams

test_streams: $(INPUT)
	$(CHANGE_NL) $(CONFIG) $(CONFIG_TMP) n_stream_sw_urban=1
	$(SPSURF) $(CONFIG_TMP) $(INPUT) $(BASENAME)_1stream_out.nc
	$(CHANGE_NL) $(CONFIG) $(CONFIG_TMP) n_stream_sw_urban=2
	$(SPSURF) $(CONFIG_TMP) $(INPUT) $(BASENAME)_2stream_out.nc
	$(CHANGE_NL) $(CONFIG) $(CONFIG_TMP) n_stream_sw_urban=4
	$(SPSURF) $(CONFIG_TMP) $(INPUT) $(BASENAME)_4stream_out.nc

# Duplicate the profiles in russell_square.nc to add a range of solar
# zenith angles, storing the result in russell_square_sza.nc
$(INPUT): $(BASE_PROFILE)
	./duplicate_profiles.sh $(BASE_PROFILE) $(INPUT)

.PHONY: *_out $(BASENAME) test*

clean:
	rm -f $(INPUT) $(BASENAME)_*_out.nc $(CONFIG_TMP)
