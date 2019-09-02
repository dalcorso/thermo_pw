!
! Copyright (C) 2019 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE check_phgeo_on_file()
!--------------------------------------------------------------------------
!
!   This routine analyses the dynamical_matrices directory and sets
!   phgeo_on_file(igeom) to .TRUE. if all the dynamical matrices for
!   this geometry are on file. 
!

USE thermo_mod,       ONLY : start_geometry, last_geometry, phgeo_on_file, &
                             tot_ngeo, no_ph
USE control_thermo,   ONLY : after_disp
USE output,           ONLY : fildyn

IMPLICIT NONE

INTEGER  :: igeom
CHARACTER(LEN=6) :: int_to_char
CHARACTER (LEN=256) :: auxdyn=' '
LOGICAL :: check_dyn_file_exists

ALLOCATE(phgeo_on_file(tot_ngeo))
phgeo_on_file=.FALSE.
!
!  loop on all the geometries calculated in this run
!
DO igeom=start_geometry,last_geometry
   IF (no_ph(igeom)) CYCLE
   auxdyn=TRIM(fildyn)//'.g'//TRIM(int_to_char(igeom))//'.'
   IF (auxdyn(1:18)/='dynamical_matrices') &
          auxdyn='dynamical_matrices/'//TRIM(auxdyn)
   IF (check_dyn_file_exists(auxdyn)) phgeo_on_file(igeom)=.TRUE.
   IF (after_disp.AND..NOT.phgeo_on_file(igeom)) &
      CALL errore('check_phgeo_on_file','after_disp but dynamical matrix &
                      &files not found', 1)
ENDDO

RETURN
END SUBROUTINE check_phgeo_on_file
