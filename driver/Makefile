TEST = test_sw

SOURCES = spartacus_surface_config.F90 spartacus_surface_read_input.F90 spartacus_surface_driver.F90
OBJECTS := $(SOURCES:.F90=.o)

EXECUTABLE = ../bin/spartacus_surface

all: $(EXECUTABLE)

test_code: $(TEST)

$(EXECUTABLE): $(OBJECTS) ../lib/*.a
	$(FC) $(OBJECTS) $(LIBS) -o $(EXECUTABLE)

$(TEST): $(TEST).F90 ../mod/rad*.mod
	$(FC) $(FCFLAGS) $(TEST).F90 $(LIBS) -o $(TEST)

%.o: %.F90 ../mod/rad*.mod
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(TEST) $(EXECUTABLE)

spartacus_surface_driver.o: spartacus_surface_read_input.o
spartacus_surface_driver.o spartacus_surface_read_input.o: spartacus_surface_config.o
