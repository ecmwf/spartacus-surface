! radtool_legendre_gauss.F90 - Derived type for holding Legendre-Gauss coefficients
!
! Copyright (C) 2019-2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radtool_legendre_gauss

  use parkind1, only : jpim, jprb

  implicit none

  !---------------------------------------------------------------------
  ! Derived type storing the points and weights for Legendre-Gauss
  ! quadrature with a particular resolution
  type legendre_gauss_type

    real(kind=jprb), allocatable, dimension(:) &
         &  :: mu, &       ! Cosine of zenith angle of each stream
         &     sin_ang, &  ! Sine of zenith angle of each stream
         &     tan_ang, &  ! Tangent of zenith angle of each stream
         &     weight, &   ! Weight of each stream in each hemisphere
         &     hweight     ! Weight for isotropic emission by a horizontal surface

    ! Number of streams in one hemisphere
    integer(kind=jpim) :: nstream = 0

  contains
    
    procedure :: initialize => initialize_legendre_gauss
    procedure :: deallocate => deallocate_legendre_gauss

  end type legendre_gauss_type

contains

  subroutine initialize_legendre_gauss(this, nstream)

    class(legendre_gauss_type), intent(inout) :: this
    integer(kind=jpim),         intent(in)    :: nstream

    call this%deallocate()

    this%nstream = nstream

    ! Calculate N-stream coefficients
    allocate(this%mu(nstream))
    allocate(this%sin_ang(nstream))
    allocate(this%tan_ang(nstream))
    allocate(this%weight(nstream))
    allocate(this%hweight(nstream))

    call calc_legendre_gauss(nstream, 0.0_jprb, 1.0_jprb, this%mu, this%weight)

    this%sin_ang = sqrt(1.0_jprb - this%mu*this%mu)
    this%tan_ang = this%sin_ang / this%mu
    this%hweight = this%weight * this%mu
    this%hweight = this%hweight / sum(this%hweight)
    
  end subroutine initialize_legendre_gauss

  subroutine deallocate_legendre_gauss(this)

    class(legendre_gauss_type), intent(inout) :: this
    
    if (allocated(this%mu))      deallocate(this%mu)
    if (allocated(this%sin_ang)) deallocate(this%sin_ang)
    if (allocated(this%tan_ang)) deallocate(this%tan_ang)
    if (allocated(this%weight))  deallocate(this%weight)
    if (allocated(this%hweight)) deallocate(this%hweight)

    this%nstream = 0

  end subroutine deallocate_legendre_gauss


  ! Compute Legendre-Gauss nodes and weights
  subroutine calc_legendre_gauss(nnode, x1, x2, xnode, weight)

    use parkind1, only : jprb

    use radiation_constants, only : Pi
    
    implicit none

    ! Number of nodes required...
    integer,         intent(in) :: nnode

    ! ...between x1 and x2
    real(kind=jprb), intent(in)  :: x1, x2

    ! Location of nodes and their weights
    real(kind=jprb), intent(out) :: xnode(nnode)
    real(kind=jprb), intent(out) :: weight(nnode)

    real(kind=jprb) :: ynode(nnode), ynode0(nnode)

    ! Legendre-Gauss Vandermonde matrix
    real(kind=jprb) :: lgvm(nnode,nnode+1), lgvm_deriv(nnode)
    
    integer jn

    do jn = 1,nnode
       ynode(jn) = cos((2*(jn-1)+1)*Pi/(2*nnode)) &
            &    + (0.27/nnode)*sin(Pi*(-1.0_jprb+2.0_jprb*jn)/(nnode+1))
       ynode0(jn) = 2.0_jprb
    end do

    do while (maxval(abs(ynode-ynode0)) > epsilon(1.0_jprb))
       lgvm(:,1) = 1.0_jprb
       lgvm(:,2) = ynode
       do jn = 2,nnode
          lgvm(:,jn+1) = ((2*jn-1)*ynode * lgvm(:,jn)-(jn-1)*lgvm(:,jn-1)) &
               &       / jn
       end do
       lgvm_deriv = (nnode+1) * (lgvm(:,nnode) - ynode*lgvm(:,nnode+1)) &
            &     / (1.0_jprb - ynode*ynode)

       ynode0 = ynode
       ynode  = ynode0 - lgvm(:,nnode+1) / lgvm_deriv
       
    end do

    ! Linear map from [-1,1] to [x1,x2]
    xnode  = 0.5_jprb * (x1*(1.0_jprb-ynode) + x2*(1.0_jprb-ynode))
    weight = ((nnode+1)*(nnode+1)/real(nnode*nnode,jprb)) * (x2-x1) &
         & / ((1.0_jprb-ynode*ynode) * lgvm_deriv*lgvm_deriv)

  end subroutine calc_legendre_gauss

end module radtool_legendre_gauss
