SPSURF = ../../bin/spartacus_surface
CONFIG = config.nam

# "make" or "make test" runs the two test cases outlined in the README
# file
test: test_surfaces_out.nc test_consistency_out.nc

# This target runs SPARTACUS-Surface on any file named *_in.nc,
# producing an output file *_out.nc
%_out.nc: %_in.nc
	$(SPSURF) $(CONFIG) $< $@

clean:
	rm -f *_out.nc

.PHONY: *_in.nc
