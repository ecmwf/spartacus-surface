TEST = test_sw

SOURCES = spartacus_surface_config.F90 spartacus_surface_read_input.F90 spartacus_surface_driver.F90
OBJECTS := $(SOURCES:.F90=.o)

EXECUTABLE = ../bin/spartacus_surface

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) ../lib/*.a
	$(FC) $(OBJECTS) $(LIBS) -o $(EXECUTABLE)

$(TEST): $(TEST).F90 $(OBJECTS) ../lib/*.a
	$(FC) $(FCFLAGS) $(TEST).F90 $(OBJECTS) $(LIBS) -o $(TEST)

%.o: %.F90 ../mod/*.mod
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(TEST)

spartacus_surface_driver.o: spartacus_surface_read_input.o
spartacus_surface_driver.o spartacus_surface_read_input.o: spartacus_surface_config.o
