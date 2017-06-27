!
! Copyright (C) 2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE check_all_geometry_done(all_geometry_done)

USE thermo_mod, ONLY : tot_ngeo
USE output, ONLY : fildyn

IMPLICIT NONE

LOGICAL, INTENT(OUT) :: all_geometry_done

INTEGER :: igeom
LOGICAL  :: check_dyn_file_exists
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=256) :: auxdyn

all_geometry_done=.TRUE.
DO igeom=1,tot_ngeo
   auxdyn=TRIM(fildyn)//'.g'//TRIM(int_to_char(igeom))//'.'
   WRITE(6,*) 'auxdyn ', TRIM(auxdyn)
   IF (all_geometry_done) all_geometry_done=all_geometry_done.AND. &
                                         check_dyn_file_exists(auxdyn)
   WRITE(6,*) 'igeom all_geometry_done', igeom, all_geometry_done
ENDDO

RETURN
END SUBROUTINE check_all_geometry_done

