! (C) Copyright 2019- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

program test_sw

  use parkind1, only : jprb
  use radtool_calc_matrices_sw_eig
  use print_matrix_mod, only : print_matrix, print_vector
  use radtool_schur, only : schur_invert_sw
  use radtool_matrix, only: invert

  implicit none
    
  integer,    parameter :: nmat  = 1
  integer,    parameter :: ndiff = 2
  integer,    parameter :: ndir  = 2
  real(jprb), parameter :: dz    = 2.0
  real(jprb), parameter :: mu0   = 0.25882 ! 75 degrees
  
  real(jprb), dimension(nmat,ndir,ndir)  :: gamma0
  real(jprb), dimension(nmat,ndiff,ndiff):: gamma1, gamma2
  real(jprb), dimension(nmat,ndiff,ndir) :: gamma3
  real(jprb), dimension(nmat,ndiff*2+ndir,ndiff*2+ndir) :: gamma
  real(jprb), dimension(nmat,ndiff,ndiff) &
       & ::  reflectance, transmittance
  real(jprb), dimension(nmat,ndiff,ndir) :: s_up, s_dn
  real(jprb), dimension(nmat,ndir,ndir)  :: trans_dir
  real(jprb), dimension(nmat,ndir,ndir)  :: gamma0i
  real(jprb), dimension(nmat,ndiff,ndiff):: gamma1i, gamma2i
  real(jprb), dimension(nmat,ndiff,ndir) :: gamma3i
  real(jprb), dimension(nmat,ndiff*2+ndir,ndiff*2+ndir) :: gammai

  gamma0(1,:,:) = transpose(reshape([-0.32359, 0.3849, 0.096224, -1.1987],[ndir,ndir]))
  gamma1(1,:,:) = transpose(reshape([-0.141379, 0.0446576, 0.178631, -0.586107],[ndiff,ndiff]))
  gamma2(1,:,:) = transpose(reshape([0.00880136, 0.0, 0.0, 0.0115973],[ndiff,ndiff]))
  gamma3(1,:,:) = transpose(reshape([0.0189622, 0.0, 0.0, 0.0227582],[ndir,ndiff]))

  gamma = 0.0_jprb
  gamma(1,1:ndiff,1:ndiff) = -gamma1(1,:,:)
  gamma(1,1:ndiff,ndiff+1:2*ndiff) = -gamma2(1,:,:)
  gamma(1,1:ndiff,2*ndiff+1:2*ndiff+ndir) = -gamma3(1,:,:)
  gamma(1,ndiff+1:2*ndiff,1:ndiff) = gamma2(1,:,:)
  gamma(1,ndiff+1:2*ndiff,ndiff+1:2*ndiff) = gamma1(1,:,:)
  gamma(1,ndiff+1:2*ndiff,2*ndiff+1:2*ndiff+ndir) = gamma3(1,:,:)
  gamma(1,2*ndiff+1:2*ndiff+ndir,2*ndiff+1:2*ndiff+ndir) = gamma0(1,:,:)

  call calc_matrices_sw_eig(nmat, ndiff, ndir, dz, mu0, gamma0, gamma1, gamma2, gamma3, &
       &  reflectance, transmittance, s_up, s_dn, trans_dir)

  call schur_invert_sw(nmat,ndir,ndiff,gamma0,gamma1,gamma2,gamma3,gamma0i,gamma1i,gamma2i,gamma3i)

  gammai = invert(1,1,2*ndiff+ndir,gamma)

  call print_matrix('gamma0',gamma0(1,:,:))
  call print_matrix('gamma1',gamma1(1,:,:))
  call print_matrix('gamma2',gamma2(1,:,:))
  call print_matrix('gamma3',gamma3(1,:,:))
  call print_matrix('gamma0i',gamma0i(1,:,:))
  call print_matrix('gamma1i',gamma1i(1,:,:))
  call print_matrix('gamma2i',gamma2i(1,:,:))
  call print_matrix('gamma3i',gamma3i(1,:,:))
  call print_matrix('gamma',gamma(1,:,:))
  call print_matrix('gammai',gammai(1,:,:))

  call print_matrix('R',reflectance(1,:,:))
  call print_matrix('T',transmittance(1,:,:))
  call print_matrix('E',trans_dir(1,:,:))
  write(6,*) s_dn
  call print_matrix('S_dn',s_dn(1,:,:))
  call print_matrix('S_up',s_up(1,:,:))

end program test_sw
