! radiation_io.F90 - Provides logging and abort functionality
!
! Copyright (C) 2015, 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
!  This file provides file units used for logging (nulout and nulerr)
!  and for reading data files (nulrad), as well as an abort routine
!  that should do clean-up appropriate for the environment in which
!  the radiation code is embedded.
!
!  This version is for the offline package, and should be rewritten as
!  appropriate if the radiation code is to be embedded in a larger
!  model.

module radiation_io

  use parkind1, only : jpim

  implicit none
  
  ! In the IFS, nulout is equivalent to standard output but only
  ! output from the primary node will be logged, while nulerr is
  ! equivalent to standard error, and text sent to this unit from any
  ! node will be logged. Normally, nulerr should only be used before
  ! calling radiation_abort.
  integer(kind=jpim), parameter :: nulout = 6
  integer(kind=jpim), parameter :: nulerr = 0

  ! This unit may be used for reading radiation configuration files,
  ! but should be closed as soon as the file is read
  integer(kind=jpim), parameter :: nulrad = 25

contains

  ! Abort the program with optional error message. Normally you would
  ! log details of the error to nulerr before calling this subroutine.
  subroutine radiation_abort(text)
    character(len=*), intent(in), optional :: text
    if (present(text)) then
      write(nulerr, '(a)') text
    else
      write(nulerr, '(a)') 'Error in radiation calculation'
    end if
    call abort
  end subroutine radiation_abort

end module radiation_io
