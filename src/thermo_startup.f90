!
! Copyright (C) 2020 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------
SUBROUTINE thermo_startup(code)
!-------------------------------------------------
!
! Initialize MPI, clocks, pools, bands, images, etc. print initial messages
!
USE mp_global,        ONLY : mp_startup
USE environment,      ONLY : environment_start

IMPLICIT NONE
CHARACTER(LEN=*) :: code
  !
  CALL mp_startup( start_images=.TRUE. )
  !
  CALL environment_start ( code )
  !
  CALL start_clock( 'PWSCF' )
  !
  RETURN
END SUBROUTINE thermo_startup
!
!-------------------------------------------------
SUBROUTINE thermo_end(code)
!-------------------------------------------------
!
!  This subroutine deallocates all the variables of thermo_pw, close
!  the files and closes MPI
!
USE mp_global,        ONLY : mp_global_end
USE environment,      ONLY : environment_end

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: code
   !
   CALL deallocate_thermo()
   !
   CALL laxlib_end()
   !
   CALL environment_end( code )
   !
   CALL mp_global_end ()
   !
   RETURN
END SUBROUTINE thermo_end
