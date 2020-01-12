program test_sw

  use parkind1, only : jprb
  use radtool_calc_matrices_sw_eig
  use print_matrix_mod, only : print_matrix, print_vector

  implicit none
    
  integer,    parameter :: nmat  = 1
  integer,    parameter :: ndiff = 2
  integer,    parameter :: ndir  = 2
  real(jprb), parameter :: dz    = 2.0
  real(jprb), parameter :: mu0   = 0.25882 ! 75 degrees
  
  real(jprb), dimension(nmat,ndir,ndir)  :: gamma0
  real(jprb), dimension(nmat,ndiff,ndiff):: gamma1, gamma2
  real(jprb), dimension(nmat,ndiff,ndir) :: gamma3
  real(jprb), dimension(nmat,ndiff,ndiff) &
       & ::  reflectance, transmittance
  real(jprb), dimension(nmat,ndiff,ndir) :: s_up, s_dn
  real(jprb), dimension(nmat,ndir,ndir)  :: trans_dir

  gamma0(1,:,:) = transpose(reshape([-0.32359, 0.3849, 0.096224, -1.1987],[ndir,ndir]))
  gamma1(1,:,:) = transpose(reshape([-0.141379, 0.0446576, 0.178631, -0.586107],[ndiff,ndiff]))
  gamma2(1,:,:) = transpose(reshape([0.00880136, 0.0, 0.0, 0.0115973],[ndiff,ndiff]))
  gamma3(1,:,:) = transpose(reshape([0.0189622, 0.0, 0.0, 0.0227582],[ndir,ndiff]))

  call calc_matrices_sw_eig(nmat, ndiff, ndir, dz, mu0, gamma0, gamma1, gamma2, gamma3, &
       &  reflectance, transmittance, s_up, s_dn, trans_dir)

  call print_matrix('R',reflectance(1,:,:))
  call print_matrix('T',transmittance(1,:,:))
  call print_matrix('E',trans_dir(1,:,:))
  write(6,*) s_dn
  call print_matrix('S_dn',s_dn(1,:,:))
  call print_matrix('S_up',s_up(1,:,:))

end program test_sw
