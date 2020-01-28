ECRAD = ../../bin/spartacus_surface
CONFIG = config.nam
INPUTBASE = test_surface
INPUT = $(INPUTBASE).nc
OUTPUT= $(INPUTBASE)_out.nc

test:
	$(ECRAD) $(CONFIG) $(INPUT) $(OUTPUT)

clean:
	rm -f $(OUTPUT)
