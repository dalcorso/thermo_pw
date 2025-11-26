!
! Copyright (C) 2020 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE manage_syn_2(nwork)
!-------------------------------------------------------------------------
!
!   This subroutine manages the synchronous calculations made in
!   part2 after the run of many pw.x calculations.
!
USE kinds, ONLY : DP
USE control_elastic_constants, ONLY : ngeom
USE control_thermo,   ONLY : lelastic_const, lpiezoelectric_tensor, &
                             lpolarization
USE control_eldos,    ONLY : lel_free_energy

IMPLICIT NONE
INTEGER :: nwork
!
! here we return synchronized and calculate the elastic constants 
! from energy or stress 
!
IF (lel_free_energy) CALL manage_el_free_energy(nwork)
!
IF (lelastic_const) CALL manage_elastic_cons(nwork, 1)
!
IF (lpiezoelectric_tensor) CALL manage_piezo_tensor(nwork, ngeom)
!
IF (lpolarization) CALL manage_polarization(nwork)

RETURN

END SUBROUTINE manage_syn_2

