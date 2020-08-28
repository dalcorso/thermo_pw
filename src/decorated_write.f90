!
! Copyright (C) 2020 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE decorated_write(message)
!--------------------------------------------------------------------------

USE io_global, ONLY : stdout

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: message

WRITE(stdout,'(/,5x,40("%"))') 
WRITE(stdout,'(a)') TRIM(message)
WRITE(stdout,'(5x,40("%"),/)') 

RETURN
END SUBROUTINE decorated_write

!--------------------------------------------------------------------------
SUBROUTINE decorated1_write(message)
!--------------------------------------------------------------------------

USE io_global, ONLY : stdout

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: message

WRITE(stdout,'(/,2x,76("+"))') 
WRITE(stdout,'(a)') TRIM(message)
WRITE(stdout,'(2x,76("+"),/)') 

RETURN
END SUBROUTINE decorated1_write
