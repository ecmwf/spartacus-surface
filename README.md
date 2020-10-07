# SPARTACUS-Surface - Canopy radiation scheme

Robin Hogan <r.j.hogan@ecmwf.int>

This document was last updated 7 October 2020


## INTRODUCTION

This package contains the offline version of a scheme for computing
radiative transfer in complex 3D surface canopies, such as forests and
cities. It makes use of the "SPARTACUS" technique, originally
developed for clouds. By treating trees and buildings as randomly
distributed in the horizontal plane, their 3D structure may be
described statistically by a handful of geometric properties, such as
the tree or building cover, and the length of the interface between
tree crowns or building walls and clear air.  See Hogan et al. (2018)
and Hogan (2019) for full details of the algorithm, and their
validation.


## PACKAGE OVERVIEW

The subdirectories are as follows:

  radsurf - source code for radiative transfer in complex surfaces
       such as vegetation and urban areas

  radtool - source code for mathematical routines needed by the
       SPARTACUS algorithm

  utilities - source code for useful utilities, such as reading NetCDF
       files

  driver - the source code for the offline driver program

  mod - where Fortran module files are written

  lib - where the static libraries are written

  bin - where the executable ecrad is written

  test - test cases including Matlab code to plot the outputs

  doc - LaTeX source for SPARTACUS-Surface documentation

## TO COMPILE

1. Ensure you have a Fortran compiler supporting the 2003
standard. Ensure you have the Fortran netCDF library installed
(version 3 or 4) for this compiler.  The command "nc-config --fc", if
available on your system, will tell you the Fortran compiler that your
netCDF library was compiled for.

2. You can compile the code using 

  make PROFILE=<prof>

where `<prof>` is one of gfortran, pgi or intel.  This will read the
compiler-specific configurations from the file
Makefile_include.<prof>.  If you omit the PROFILE= option then
gfortran will be assumed. If you have a compiler other than these
three then create such a file for your compiler following the example
in Makefile_include.gfortran.

If the compile is successful then static libraries should appear in
the lib directory, and then the executable bin/spartacus_surface.

3. To clean-up, type "make clean".  To build an unoptimized version
for debugging, you can do

  make PROFILE=<prof> DEBUG=1

or you can specifically override the variables in Makefile_include.<prof>
using, for example

  make PROFILE=<prof> OPTFLAGS=-O0 DEBUGFLAGS="-g -pg"


## TO TEST

The offline driver is run via

    spartacus_surface <namelist.nam> <input_file.nc> <output_file.nc>

where the radiation scheme is configured using the Fortran namelist
`<namelist.nam>`, and the inputs and outputs are in netCDF format.  

To run the tests, type

   make test

from the top-level directory.  This runs "make" in each of the
subdirectories of the test directory, which in turn runs
spartacus_surface on different input files.

The "test/simple" directory contains a minimal test file of the four
surface types that SPARTACUS-Surface supports: flat, forest,
unvegetated urban and vegetated urban.

The "test/rami4pilps" directory contains a configuration to run the
RAMI4PILPS test cases that were used by Hogan et al. (2018).

The "test/urban" directory contains an urban profile from Fig. 1 of
Hogan (2019).


## LICENSE

(C) Copyright 2019- ECMWF.

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

In applying this licence, ECMWF does not waive the privileges and immunities
granted to it by virtue of its status as an intergovernmental organisation
nor does it submit to any jurisdiction.
Copyright statements are given in the file NOTICE.


## PUBLICATIONS

Hogan, R. J., T. Quaife and R. Braghiere, 2018: Fast matrix treatment
of 3-D radiative transfer in vegetation canopies: SPARTACUS-Vegetation
1.1. Geosci. Model Dev., 11, 339-350.

Hogan, R. J., 2019b: Flexible treatment of radiative transfer in
complex urban canopies for use in weather and climate
models. Boundary-Layer Meteorol., 173, 53-78.


## CONTACT

Any queries or bug fixes, please email Robin Hogan <r.j.hogan@ecmwf.int>



