! radtool_calc_matrices_sw_eig.F90 - Compute reflectance matrix and friends for shortwave
!
! Copyright (C) 2019 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radtool_calc_matrices_sw_eig

contains

  ! This routine implements computes various properties in a single
  ! layer of an urban or vegetation canopy in the shortwave: the
  ! diffuse reflectance and transmittance matrices, the direct
  ! transmittance "trans_dir" and the rate at which direct radiation
  ! is scattered up out of the top of the layer "s_up" and down out of
  ! the base of the layer "s_dn".  The first dimension of each array
  ! is spectral interval, so multiple spectral intervals are computed
  ! at once.  Equation numbers mentioned in the comments below are
  ! from Hogan (2019, Boundary Layer Meteorology: Flexible treatment
  ! of radiative transfer in complex urban canopies).
  subroutine calc_matrices_sw_eig(nmat, ndiff, ndir, dz, mu0, &
       &  gamma0, gamma1, gamma2, gamma3, &
       &  reflectance, transmittance, s_up, s_dn, trans_dir)

    use parkind1, only : jprb

    use radtool_matrix, only : mat_x_mat, rect_mat_x_mat, mat_x_vec, &
         &                     invert, solve_mat, solve_vec, solve_rect_mat

    use radtool_eigen_decomposition, only : eigen_decomposition_real

    implicit none

    
    ! Inputs

    ! Number of input matrices
    integer, intent(in) :: nmat

    ! Size of diffuse gamma matrices
    integer, intent(in) :: ndiff

    ! Size of direct gamma matrices
    integer, intent(in) :: ndir

    ! Layer thickness (m)
    real(kind=jprb), intent(in) :: dz

    ! Cosine of solar zenith angle
    real(kind=jprb), intent(in) :: mu0

    ! Exchange matrices; the unique parts of the matrix "Gamma" that
    ! is not explicitly constructed:
    !   Gamma = [-gamma1 -gamma2 -gamma3]
    !           [+gamma2 +gamma1 +gamma3]
    !           [                 gamma0]
    real(kind=jprb), intent(in) :: gamma0(nmat,ndir,ndir)
    real(kind=jprb), intent(in) :: gamma1(nmat,ndiff,ndiff)
    real(kind=jprb), intent(in) :: gamma2(nmat,ndiff,ndiff)
    real(kind=jprb), intent(in) :: gamma3(nmat,ndiff,ndir)

    ! Outputs
    
    ! Diffuse reflectance and transmittance matrices, "R" and "T" in
    ! Eq. 49
    real(kind=jprb), intent(out), dimension(nmat,ndiff,ndiff) &
         &  :: reflectance, transmittance
    
    ! Radiation emerging from top and bottom of layer due to
    ! scattering by the direct beam, and the transmittance of the
    ! layer to direct radiation with no scattering on the way
    real(kind=jprb), intent(out), dimension(nmat,ndiff, ndiff) &
         &  :: s_up, s_dn, trans_dir

    
    ! Local variables

    ! Difference between gamma1 and gamma2
    real(kind=jprb) :: gamma_diff(nmat,ndiff,ndiff)

    ! The matrix product (gamma1-gamma2)*(gamma1+gamma2)
    real(kind=jprb) :: gamma_product(nmat,ndiff,ndiff)

    ! Eigenvalues/eigenvectors of gamma_product
    real(kind=jprb) :: eigenval_prod(nmat,ndiff), &
         &             eigenvec_prod(nmat,ndiff,ndiff)

    ! Eigenvalues of gamma0
    real(kind=jprb) :: eigenval_dir(nmat,ndir)
    
    ! The unique submatrices of the eigenvector matrix of the full
    ! Gamma matrix, "G_0", "G_1" etc. in Eq. 45
    real(kind=jprb), dimension(nmat,ndir,ndir)   :: g0
    real(kind=jprb), dimension(nmat,ndiff,ndiff) :: g1, g2
    real(kind=jprb), dimension(nmat,ndiff,ndir)  :: g3, g4
    
    ! g1 and g2 multiplied by the diagonal matrix "D" defined after
    ! Eq. 46
    real(kind=jprb), dimension(nmat,ndiff,ndiff) :: g1_d, g2_d

    ! The square upper and lower parts of the rectangular matrix "C'"
    ! in Eq. 47
    real(kind=jprb), dimension(nmat,ndiff,ndiff) :: cprime_upper, cprime_lower

    ! Similar but for computing s_up and s_dn
    real(kind=jprb), dimension(nmat,ndiff,ndir) :: cdir_upper, cdir_lower

    ! Positive eigenvalues of the full Gamma matrix
    real(kind=jprb) :: lambda(nmat,ndiff)

    ! exp(-lambda*dz)
    real(kind=jprb) :: exp_lambda_dz(nmat,ndiff)

    ! gamma2*g0
    real(kind=jprb), dimension(nmat,ndiff,ndir) :: gamma3_g0

    ! Inverse of g0
    real(kind=jprb), dimension(nmat,ndiff,ndir) :: g0_inv

    ! DEFINE
    real(kind=jprb), dimension(nmat,ndiff,ndir) :: gamma1_d

    real(kind=jprb), dimension(nmat,ndiff,ndir) :: gamma2_inv_gamma1_d

    real(kind=jprb), dimension(nmat,ndiff,ndir) :: g3_d_inv_g0, g4_inv_g0
    
    ! Temporary vector and matrix
    real(kind=jprb), dimension(nmat,ndiff) :: tmp_vec
    real(kind=jprb), dimension(nmat,ndiff,ndiff) :: tmp_mat

    ! Loop index over matrix elements
    integer :: jd, jo

    
    ! Section 1: Compute the eigenvectors and eigenvalues of the
    ! diffuse part of the Gamma matrix using the clever DISORT method
    ! exploiting its regular structure

    ! Compute the eigenvectors and eigenvalues of "gamma_product"
    gamma_diff = gamma1-gamma2;
    gamma_product = mat_x_mat(nmat,nmat,ndiff,gamma_diff,gamma1+gamma2)

    call eigen_decomposition_real(ndiff,nmat,gamma_product, &
         &                        eigenval_prod, eigenvec_prod)

    ! Eigenvalues of the diffuse part of the full Gamma matrix are
    ! (+/-)sqrt(eigenval_prod)
    lambda = sqrt(max(0.0_jprb,eigenval_prod))
    exp_lambda_dz = exp(-lambda*dz)

    ! Compute the submatrices g1 and g2 of the eigenvector matrix of
    ! the full Gamma matrix
    tmp_mat = -solve_mat(nmat,nmat,ndiff,gamma_diff,eigenvec_prod)
    tmp_mat = tmp_mat * spread(lambda,2,ndiff)
    g1 = eigenvec_prod + tmp_mat
    g2 = eigenvec_prod - tmp_mat
    !g1 = -0.5*(eigenvec_prod + tmp_mat)
    !g2 = -0.5*(eigenvec_prod - tmp_mat)

    
    ! Section 2: Compute diffuse reflectance and transmittance
    ! matrices

    ! Various equations require product of g1,g2 and exp(-lambda*dz)
    g1_d = g1 * spread(exp_lambda_dz,2,ndiff)
    g2_d = g2 * spread(exp_lambda_dz,2,ndiff)

    ! Solve Eq. 48 to compute the upper and lower part of "C'", using
    ! the Schur complement to exploit the regular structure of the
    ! matrix on the left hand side of Eq. 48.
    cprime_lower = invert(nmat,nmat,ndiff,(g1 &
         &  -mat_x_mat(nmat,nmat,ndiff,g2_d, &
         &             solve_mat(nmat,nmat,ndiff,g1,g2_d))))
    cprime_upper = -solve_mat(nmat,nmat,ndiff,g1, &
         &   mat_x_mat(nmat,nmat,ndiff,g2_d, cprime_lower))

    ! Apply Eq. 49
    reflectance   = mat_x_mat(nmat,nmat,ndiff,g1_d,cprime_upper) &
         &        + mat_x_mat(nmat,nmat,ndiff,g2,  cprime_lower)
    transmittance = mat_x_mat(nmat,nmat,ndiff,g2  ,cprime_upper) &
         &        + mat_x_mat(nmat,nmat,ndiff,g1_d,cprime_lower)


    ! Section 3: Direct transmittance matrix: the matrix exponential
    ! of gamma0, which we compute using eigen-decomposition
    call eigen_decomposition_real(ndir,nmat,gamma0, &
         &                        eigenval_dir, g0)
    g0_inv = invert(nmat,nmat,ndir,g0)
    trans_dir = mat_x_mat(nmat,nmat,ndir, &
         &  g0*spread(exp(eigenval_dir*dz),2,ndir), g0_inv)
   
    ! Section 4: Mixed direct-diffuse part
    gamma3_g0 = rect_mat_x_mat(nmat,ndiff,ndiff,ndir,gamma3,g0)
    gamma1_d = gamma1
    do jd = 1,ndir
       do jo = 1,ndiff
          gamma1_d(:,jo,jo) = gamma1(:,jo,jo) + eigenval_dir(:,jd)
       end do
       gamma2_inv_gamma1_d = mat_x_mat(nmat,nmat,ndiff,gamma2,&
            &                        invert(nmat,nmat,ndiff,gamma1_d))
       tmp_mat = gamma1 - mat_x_mat(nmat,nmat,ndiff,&
            &                       gamma2_inv_gamma1_d, gamma2)
       do jo = 1,ndiff
          tmp_mat(:,jo,jo) = tmp_mat(:,jo,jo) - eigenval_dir(:,jd)
       end do
       ! Subtract identity matrix
       do jo = 1,ndiff
          gamma2_inv_gamma1_d(:,jo,jo) = gamma2_inv_gamma1_d(:,jo,jo) - 1.0_jprb
       end do
       g4(:,:,jd) = solve_vec(nmat,nmat,ndiff,tmp_mat, &
            &  mat_x_vec(nmat,nmat,ndiff,gamma2_inv_gamma1_d,gamma3_g0(:,:,jd)))
       g3(:,:,jd) = -solve_vec(nmat,nmat,ndiff,gamma1_d, &
            &  gamma3_g0(:,:,jd)+mat_x_vec(nmat,nmat,ndiff,gamma2,g4(:,:,jd)))
    end do

    call direct_diffuse_part(nmat, ndiff, ndir, mu0, &
         &  exp_lambda_dz, exp(eigenval_dir*dz), &
         &  g0, g1, g2, g3, g4, &
         &  s_up, s_dn)

#ifdef FAST_BUT_INCORRECT

    g3_d_inv_g0 = rect_mat_x_mat(nmat,ndiff,ndir,ndir,g3*spread(exp(eigenval_dir*dz),2,ndir), &
         &                  g0_inv)
    g4_inv_g0 = rect_mat_x_mat(nmat,ndiff,ndir,ndir,g4, g0_inv)

    ! Schur complement
    tmp_mat = g1 - mat_x_mat(nmat,nmat,ndiff,g2_d, &
         &    solve_mat(nmat,nmat,ndiff,g1,g2_d))
    cdir_upper = -solve_rect_mat(nmat,ndiff,ndir,g1,g3_d_inv_g0 &
         &  + rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g2_d, &
         &  solve_rect_mat(nmat,ndiff,ndir,tmp_mat, &
         &  rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g2_d, &
         &  rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g1,g3_d_inv_g0) &
         &  - g4_inv_g0))))
    cdir_lower = solve_rect_mat(nmat,ndiff,ndir,tmp_mat, &
         &  rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g2_d, &
         &  solve_rect_mat(nmat,ndiff,ndir,g1,g3_d_inv_g0) &
         &  + g4_inv_g0))
    s_up = mu0 * (rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g1_d,cdir_upper) &
         &      + rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g2,  cdir_lower) &
         &      + rect_mat_x_mat(nmat,ndiff,ndir, ndir,g3,  g0_inv))
    s_dn = mu0 * (rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g2,  cdir_upper) &
         &      + rect_mat_x_mat(nmat,ndiff,ndiff,ndir,g1_d,cdir_lower) &
         &      + rect_mat_x_mat(nmat,ndiff,ndir, ndir, &
         &          g4*spread(exp(eigenval_dir*dz),2,ndir),  g0_inv))

#endif
    
  end subroutine calc_matrices_sw_eig


  ! Brute-force method, without any clever optimizations using the
  ! Schur complement
  subroutine direct_diffuse_part(nmat, ndiff, ndir, mu0, &
       &  exp_lambda_dz, exp_eigenval_dir_dz, &
       &  g0, g1, g2, g3, g4, &
       &  s_up, s_dn)
    
    use parkind1, only : jprb
    use radtool_matrix, only : rect_mat_x_mat, solve_rect_mat

    ! Number of input matrices
    integer, intent(in) :: nmat

    ! Size of diffuse gamma matrices
    integer, intent(in) :: ndiff

    ! Size of direct gamma matrices
    integer, intent(in) :: ndir

    ! Cosine of solar zenith angle
    real(kind=jprb), intent(in) :: mu0

    ! exp(-lambda*dz)
    real(kind=jprb), intent(in) :: exp_lambda_dz(nmat,ndiff)

    ! exp([Eigenvalues of gamma0]*dz)
    real(kind=jprb), intent(in) :: exp_eigenval_dir_dz(nmat,ndir)

    ! The unique submatrices of the eigenvector matrix of the full
    ! Gamma matrix, "G_0", "G_1" etc. in Eq. 45
    real(kind=jprb), intent(in), dimension(nmat,ndir,ndir)   :: g0
    real(kind=jprb), intent(in), dimension(nmat,ndiff,ndiff) :: g1, g2
    real(kind=jprb), intent(in), dimension(nmat,ndiff,ndir)  :: g3, g4

    ! Radiation emerging from top and bottom of layer due to
    ! scattering by the direct beam
    real(kind=jprb), intent(out), dimension(nmat,ndiff, ndiff) &
         &  :: s_up, s_dn

    real(kind=jprb) :: g_d(nmat,2*ndiff+ndir,2*ndiff+ndir)
    real(kind=jprb) :: g_d2(nmat,ndiff,2*ndiff+ndir)
    real(kind=jprb) :: cprime_dir(nmat,2*ndiff+ndir,ndir)
    real(kind=jprb) :: rhs(nmat,2*ndiff+ndir,ndir)

    integer :: jj

    g_d = 0.0_jprb
    ! Fill diffuse part of eigenvector matrix
    g_d(:,1:ndiff,1:ndiff) = g1
    g_d(:,ndiff+1:2*ndiff,1:ndiff) = g2 * spread(exp_lambda_dz,2,ndiff)
    g_d(:,1:ndiff,ndiff+1:2*ndiff) = g_d(:,ndiff+1:2*ndiff,1:ndiff)
    g_d(:,ndiff+1:2*ndiff,ndiff+1:2*ndiff) = g1
    ! Fill direct part
    g_d(:,2*ndiff+1:2+ndiff+ndir,2*ndiff+1:2*ndiff+ndir) = g0
    ! Mixed part

    print *, g_d(:,1:ndiff,2*ndiff+1:2*ndiff+ndir)
!    print *, g3
    print *, spread(exp_eigenval_dir_dz,2,ndiff)

    g_d(:,1:ndiff,2*ndiff+1:2*ndiff+ndir) = g3 * spread(exp_eigenval_dir_dz,2,ndiff)
    g_d(:,ndiff+1:2*ndiff,2*ndiff+1:2*ndiff+ndir) = g4

    rhs = 0.0_jprb
    !rhs = 1.0e-6_jprb
    do jj = 1,ndir
      rhs(:,2*ndiff+jj,jj) = 1.0_jprb
    end do

    cprime_dir = solve_rect_mat(nmat, 2*ndiff+ndir, ndir, g_d, rhs)

    g_d2(:,1:ndiff,1:ndiff) = g1 * spread(exp_lambda_dz,2,ndiff)
    g_d2(:,1:ndiff,ndiff+1:2*ndiff) = g2
    g_d2(:,1:ndiff,2*ndiff+1:2*ndiff+ndir) = g3

    s_up = rect_mat_x_mat(nmat, ndiff, 2*ndiff+ndir, ndir, g_d2, cprime_dir)
    s_up = s_up * mu0

    g_d2(:,1:ndiff,1:ndiff) = g2
    g_d2(:,1:ndiff,ndiff+1:2*ndiff) = g1 * spread(exp_lambda_dz,2,ndiff)
    g_d2(:,1:ndiff,2*ndiff+1:2*ndiff+ndir) = g4 * spread(exp_eigenval_dir_dz,2,ndiff)

    s_dn = rect_mat_x_mat(nmat, ndiff, 2*ndiff+ndir, ndir, g_d2, cprime_dir)
    s_dn = s_dn * mu0

  end subroutine direct_diffuse_part

end module radtool_calc_matrices_sw_eig