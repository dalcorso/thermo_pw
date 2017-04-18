!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
LOGICAL FUNCTION check_file_exists(filename)
!
!   This function checks if a file exists and returns true or false.
!   It has to be called by all processors when they are synchronized.
!   Only meta_ionode checks that the file exists and sends the information
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

RETURN
END FUNCTION check_file_exists

LOGICAL FUNCTION check_dyn_file_exists(filename)
!
!   This function checks if all the dynamical matrix files at a given geometry
!   exist. In this case the phonon calculation at that geometry is not redone.
!   It has to be called by all processors when they are synchronized.
!   Only the meta_ionode checks that the file exists and sends the information
!   to all processors.
!   If after_disp is true but the dynamical matrix files are not found, 
!   this routine stops with an error.
!
USE io_global, ONLY : meta_ionode, meta_ionode_id
USE mp_world,  ONLY : world_comm
USE mp,        ONLY : mp_bcast
USE control_thermo, ONLY : after_disp
USE control_ph, ONLY : ldisp, xmldyn
IMPLICIT NONE

CHARACTER (LEN=256), INTENT(IN) :: filename
CHARACTER (LEN=256) :: fildyn
CHARACTER (LEN=6) :: int_to_char
LOGICAL :: exst, exst_all
INTEGER :: iq, nq1, nq2, nq3, nqs, iudyn
INTEGER :: find_free_unit

exst_all=.TRUE.
nqs=0
fildyn = TRIM( filename ) // '0'
IF (meta_ionode) INQUIRE(FILE=TRIM(fildyn),EXIST=exst)
CALL mp_bcast(exst, meta_ionode_id, world_comm)
exst_all=exst_all.AND.exst
IF (exst_all) THEN
   IF (meta_ionode) THEN
      iudyn=find_free_unit()
      OPEN (UNIT=iudyn, FILE=TRIM(fildyn), STATUS='old',FORM='formatted')
      READ(iudyn,*) nq1, nq2, nq3
      READ(iudyn,*) nqs
      CLOSE(iudyn)
   END IF
   CALL mp_bcast(nq1, meta_ionode_id, world_comm)
   CALL mp_bcast(nq2, meta_ionode_id, world_comm)
   CALL mp_bcast(nq3, meta_ionode_id, world_comm)
   CALL mp_bcast(nqs, meta_ionode_id, world_comm)
END IF

DO iq=1,nqs
   IF (exst_all) THEN
      IF (ldisp) fildyn = TRIM( filename ) // TRIM( int_to_char( iq ) )
      IF (xmldyn) fildyn = TRIM(fildyn) // '.xml'
      IF (meta_ionode) INQUIRE(FILE=TRIM(fildyn),EXIST=exst)
      CALL mp_bcast(exst, meta_ionode_id, world_comm)
      exst_all=exst_all.AND.exst
   ENDIF
END DO

IF (.NOT. exst_all .AND. after_disp) &
      CALL errore('check_dyn_file_exists','after_disp but dynamical matrix &
                      &files not found', 1)

check_dyn_file_exists=exst_all

RETURN
END FUNCTION check_dyn_file_exists
