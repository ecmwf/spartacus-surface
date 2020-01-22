! radtool_schur.f90 - Matrix manipulation using the Schur complement
!
! Copyright (C) 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radtool_schur


contains

  !---------------------------------------------------------------------
  ! Use the Schur complement method to invert the Gamma matrix
  ! occurring in the shortwave SPARTACUS method.  If 
  !        g = [ -g1  -g2  -g3  ]
  !            [  g2   g1   g3  ]
  !            [            g0  ]
  ! then its inverse is
  !   inv(g) = [ -g1i -g2i -g3i ]
  !            [  g2i  g1i -g3i ]
  !            [            g0i ]
  !
  subroutine schur_invert_sw(nmat, n0, n1, g0, g1, g2, g3, g0i, g1i, g2i, g3i)
    
    use parkind1, only : jprb, jpim
    use radtool_matrix, only : invert, mat_x_mat, solve_mat, rect_mat_x_mat

    implicit none

    integer(kind=jpim), intent(in)  :: nmat, n0, n1
    real(kind=jprb),    intent(in)  :: g0(nmat,n0,n0), g1(nmat,n1,n1)
    real(kind=jprb),    intent(in)  :: g2(nmat,n1,n1), g3(nmat,n1,n0)
    real(kind=jprb),    intent(out) :: g0i(nmat,n0,n0), g1i(nmat,n1,n1)
    real(kind=jprb),    intent(out) :: g2i(nmat,n1,n1), g3i(nmat,n1,n0)

    g0i = invert(nmat,nmat,n0,g0)
    g1i = invert(nmat,nmat,n1,g1 - mat_x_mat(nmat,nmat,n1,g2, &
         &                           solve_mat(nmat,nmat,n1,g1,g2)))
    g2i = mat_x_mat(nmat,nmat,n1,g1i,mat_x_mat(nmat,nmat,n1,g2, &
         &                                  invert(nmat,nmat,n1,g1)))
    g3i = rect_mat_x_mat(nmat,n1,n1,n0,g1i-g2i, &
         &               rect_mat_x_mat(nmat,n1,n0,n0,g3,g0i))

  end subroutine schur_invert_sw


end module radtool_schur
