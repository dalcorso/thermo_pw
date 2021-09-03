!
! Copyright (C) 2021 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
PROGRAM template
!--------------------------------------------------------------------
!
!
USE kinds, ONLY : DP
USE mp_global,        ONLY : mp_startup, mp_global_end
USE environment,      ONLY : environment_start, environment_end

USE io_global,        ONLY : stdout

IMPLICIT NONE
CHARACTER(LEN=9) :: code='TEMPLATE'


CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

WRITE(stdout,'("Hello world")')

CALL environment_end( code )
CALL mp_global_end ()

END PROGRAM template

