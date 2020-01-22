SOURCES = \
	radiation_constants.F90 \
	radtool_legendre_gauss.F90 \
	radtool_matrix.F90 \
	radtool_eigen_decomposition.F90 \
	radtool_calc_matrices_lw_eig.F90 \
	radtool_calc_matrices_sw_eig.F90 \
	radtool_schur.F90

OBJECTS := $(SOURCES:.F90=.o)
LIBRADTOOL = ../lib/libradtool.a

all: $(LIBRADTOOL)

$(LIBRADTOOL): $(OBJECTS)
	ar r $(LIBRADTOOL) $(OBJECTS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(LIBRADTOOL)

radtool_legendre_gauss.o: radiation_constants.o
radtool_calc_matrices_sw_eig.o radtool_calc_matrices_lw_eig.o: \
	radtool_matrix.o radtool_eigen_decomposition.o
radtool_schur.o: radtool_matrix.o
