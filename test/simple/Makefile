ECRAD = ../../bin/spartacus_surface
CONFIG = config.nam

test: test_surfaces_out.nc test_consistency_out.nc

%_out.nc: %_in.nc
	$(ECRAD) $(CONFIG) $< $@

clean:
	rm -f *_out.nc

.PHONY: *_in.nc
