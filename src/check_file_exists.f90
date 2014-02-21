!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
LOGICAL FUNCTION check_file_exists(filename)
!
!   This function checks if a file exists and return true or false.
!   It has to be called by all processors when they are syncronized.
!   Only the meta_ionode checks the file existence and send the information
!   to all processors
!
USE io_global, ONLY : meta_ionode, meta_ionode_id
USE mp_world,  ONLY : world_comm
USE mp,        ONLY : mp_bcast
IMPLICIT NONE

CHARACTER(LEN=256), INTENT(IN) :: filename
LOGICAL :: exst

IF (meta_ionode) INQUIRE(FILE=TRIM(filename),EXIST=exst)
CALL mp_bcast(exst, meta_ionode_id, world_comm)

check_file_exists=exst
write(6,*) 'file ', TRIM(filename), ' exists', exst

RETURN
END FUNCTION check_file_exists

LOGICAL FUNCTION check_dyn_file_exists(filename)
!
!   This function checks if all the dynamical matrices file at a given geometry
!   exist. Only in this case the calculations are not redone.
!   It has to be called by all processors when they are syncronized.
!   Only the meta_ionode checks the file existence and send the information
!   to all processors
!
USE io_global, ONLY : meta_ionode, meta_ionode_id
USE mp_world,  ONLY : world_comm
USE mp,        ONLY : mp_bcast
USE disp,      ONLY : nqs
USE control_ph, ONLY : ldisp, xmldyn
IMPLICIT NONE

CHARACTER (LEN=256), INTENT(IN) :: filename
CHARACTER (LEN=256) :: fildyn
CHARACTER (LEN=6) :: int_to_char
LOGICAL :: exst, exst_all
INTEGER :: iq

exst_all=.TRUE.
fildyn = TRIM( filename ) 
DO iq=1,nqs
   IF (exst_all) THEN
      IF (ldisp) fildyn = TRIM( filename ) // TRIM( int_to_char( iq ) )
      IF (xmldyn) fildyn = TRIM(fildyn) // '.xml'
      IF (meta_ionode) INQUIRE(FILE=TRIM(fildyn),EXIST=exst)
      CALL mp_bcast(exst, meta_ionode_id, world_comm)
      exst_all=exst_all.AND.exst
   ENDIF
END DO

check_dyn_file_exists=exst_all
write(6,*) 'file ', TRIM(filename), ' exists', exst_all

RETURN
END FUNCTION check_dyn_file_exists
