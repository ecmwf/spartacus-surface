! (C) Copyright 2019- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

program test_lw

  use parkind1, only : jprb
  use radtool_calc_matrices_lw_eig
  use print_matrix_mod, only : print_matrix, print_vector
  use legendre_gauss, only : calc_legendre_gauss
  
  integer,    parameter :: nmat = 1
  integer,    parameter :: norder = 2
  real(jprb), parameter :: dz = 50.0
  
  real(jprb), dimension(nmat,norder,norder) :: gamma1, gamma2, int_flux, &
       &   reflectance, transmittance
  real(jprb), dimension(nmat,norder) :: int_flux_source, source, emiss_rate

  integer, parameter :: nnode = 8
  real(jprb), dimension(nnode) :: xnode(nnode), weight(nnode)

  integer :: jnode
  
  ! gamma1(1,:,:) = reshape([ &
  !    -0.0078528,   4.4823e-05,  0.00011571,  0.00058123, &
  !     6.0502e-05, -0.02199,     0.00043988,  0.0022095, &
  !     7.6933e-05,  0.00021667, -0.05665,     0.0028095, &
  !      .3367e-05,  0.00012213,  0.0003153,  -0.28577], &
  !      [norder, norder])

  ! gamma2(1,:,:) = reshape([ &
  !     1.5916e-05,  4.4823e-05,  0.00011571,  0.00058123, &
  !     6.0502e-05,  0.00017039,  0.00043988,  0.0022095, &
  !     7.6933e-05,  0.00021667,  0.00055934,  0.0028095, &
  !     4.3367e-05,  0.00012213,  0.0003153,   0.0015837], &
  !     [norder, norder])
  
  ! emiss_rate(1,:) = [0.70539, 2.6815, 3.4097, 1.922]

  gamma1(1,:,:) = transpose(reshape([-0.01544, 0.00089297, &
       0.00023927,   -0.091084],[norder,norder]))
  gamma2(1,:,:) = transpose(reshape([0.0001505,  0.00089297, &
       0.00023927,   0.0014196],[norder,norder]))
  emiss_rate(1,:) = [3.4053, 5.4136]

  call calc_matrices_lw_eig(nmat, norder, dz, gamma1, gamma2, emiss_rate, &
       &  reflectance, transmittance, source, int_flux, int_flux_source)

  call print_matrix('R',reflectance(1,:,:))
  call print_matrix('T',transmittance(1,:,:))
  call print_vector('source',source(1,:))
  call print_matrix('int_flux',int_flux(1,:,:))
  call print_vector('int_flux_source',int_flux_source(1,:))

  do jnode = 1,nnode
     write(6,*) 'Legendre-Gauss quadrature, n = ', jnode
     call calc_legendre_gauss(jnode, 0.0_jprb, 1.0_jprb, xnode, weight)
     call print_vector('x     ',xnode(1:jnode))
     call print_vector('weight',weight(1:jnode))
  end do
  
end program test_lw
