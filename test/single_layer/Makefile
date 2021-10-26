SPSURF = ../../bin/spartacus_surface
CONFIG = config.nam
CONFIG_TMP = config_tmp.nam
CHANGE_NL = ./change_namelist.sh
BASENAME = test_single_layer
INPUT = $(BASENAME).nc

test: test_sp test_exp test_inf

test_sp:
	$(CHANGE_NL) $(CONFIG) $(CONFIG_TMP) isurfacetype=2
	$(SPSURF) $(CONFIG_TMP) $(INPUT) $(BASENAME)_sp_out.nc

test_exp:
	$(CHANGE_NL) $(CONFIG) $(CONFIG_TMP) isurfacetype=4
	$(SPSURF) $(CONFIG_TMP) $(INPUT) $(BASENAME)_exp_out.nc

test_inf:
	$(CHANGE_NL) $(CONFIG) $(CONFIG_TMP) isurfacetype=5
	$(SPSURF) $(CONFIG_TMP) $(INPUT) $(BASENAME)_inf_out.nc

clean:
	rm -f $(BASENAME)_*_out.nc $(CONFIG_TMP)
