! (C) Copyright 2020- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

module print_matrix_mod
contains
  subroutine print_vector(name, vec)
    use parkind1, only : jprb
    character(*), intent(in) :: name
    real(jprb),   intent(in) :: vec(:)
    integer    :: i
    write(6,'(a,a)') name, '='
    do i = 1,size(vec,1)
       write(6,'(f16.8,$)') vec(i)
       write(6,'(x)')
    end do
  end subroutine print_vector

  subroutine print_matrix(name, mat)
    use parkind1, only : jprb
    character(*), intent(in) :: name
    real(jprb),   intent(in) :: mat(:,:)
    integer    :: i, j
    write(6,'(a,a)') name, '='
    do i = 1,size(mat,1)
       do j = 1,size(mat,2)
          write(6,'(f16.8,$)') mat(i,j)
       end do
       write(6,'(x)')
    end do
  end subroutine print_matrix
  
  subroutine print_array3(name, mat)
    use parkind1, only : jprb
    character(*), intent(in) :: name
    real(jprb),   intent(in) :: mat(:,:,:)
    integer    :: i, j, k
    write(6,'(a,a)') name, '='
    do k = 1,size(mat,3)
      do i = 1,size(mat,1)
        do j = 1,size(mat,2)
             write(6,'(f16.8,$)') mat(i,j,k)
          end do
          write(6,'(x)')
       end do
       write(6,'(x)')
    end do
  end subroutine print_array3
  

end module print_matrix_mod
