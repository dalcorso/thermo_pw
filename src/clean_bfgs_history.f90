!
! Copyright (C) 2013-2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE clean_bfgs_history()
!-----------------------------------------------------------------------

USE io_files,    ONLY : tmp_dir, prefix, seqopn
USE io_global,   ONLY : ionode

IMPLICIT NONE
INTEGER :: iunupdate
INTEGER :: find_free_unit
LOGICAL :: exst
CHARACTER(LEN=256) :: filename

IF (ionode) THEN
   !
   !  clean the bfgs history
   !
   iunupdate=find_free_unit()
   CALL seqopn( iunupdate, 'update', 'FORMATTED', exst )
   CLOSE(iunupdate, STATUS='DELETE')
   filename = TRIM( tmp_dir ) // TRIM( prefix ) // '.bfgs'
   OPEN( iunupdate, FILE=TRIM(filename), FORM='FORMATTED')
   CLOSE(iunupdate, STATUS='DELETE')
ENDIF
RETURN
END SUBROUTINE clean_bfgs_history
