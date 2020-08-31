!
! Copyright (C) 2019 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE check_dynmat_all_geo_on_file()
!--------------------------------------------------------------------------
!
!   This routine analyses the dynamical_matrices directory and sets
!   dynmat_on_file(igeom) to .TRUE. if all the dynamical matrices for
!   this geometry are on file. 
!
!   NB: this routine works only if after_disp=.TRUE. and/or fildyn has been
!       given in input to thermo_pw. 
!

USE thermo_mod,       ONLY : start_geometry, last_geometry, dynmat_on_file, &
                             tot_ngeo, no_ph
USE control_thermo,   ONLY : after_disp

IMPLICIT NONE

INTEGER  :: igeom
!
!  Allocate the dynmat_on_file variable. This routine can be called only
!  once
!
ALLOCATE(dynmat_on_file(tot_ngeo))
dynmat_on_file=.FALSE.
!
!  loop on all the geometries calculated in this run
!
DO igeom=start_geometry,last_geometry
   IF (no_ph(igeom)) CYCLE
   CALL set_fildyn_name(igeom)
   CALL check_dynmat_on_file_1g(igeom)
ENDDO

RETURN
END SUBROUTINE check_dynmat_all_geo_on_file
!
 
!--------------------------------------------------------------
SUBROUTINE check_dynmat_on_file_1g(igeom)
!--------------------------------------------------------------

USE thermo_mod,     ONLY : dynmat_on_file
USE control_thermo, ONLY : after_disp
USE output,         ONLY : fildyn

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom
CHARACTER(LEN=6)    :: int_to_char
LOGICAL             :: check_dyn_file_exists

IF (check_dyn_file_exists(fildyn)) dynmat_on_file(igeom)=.TRUE.
IF (after_disp.AND..NOT.dynmat_on_file(igeom)) &
   CALL errore('check_dynmat_on_file','after_disp but dynamical matrix &
                      &files not found', 1)
RETURN
END SUBROUTINE check_dynmat_on_file_1g
