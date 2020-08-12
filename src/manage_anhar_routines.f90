!
! Copyright (C) 2020 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE manage_anhar_routines()
!
!   This subroutine manages the calculation of the anharmonic properties
!   and selects the routines that make the isotropic (where volume is a 
!   basic variable) or anisotropic (where strain is a basic variable)
!   thermodynamic.
!
USE control_mur, ONLY : lmurn
USE mp_world,    ONLY : world_comm
USE mp,          ONLY : mp_barrier

IMPLICIT NONE
!
!  We resinchronize the processors here, to be sure that all phonon
!  calculations finish.
!
   CALL mp_barrier(world_comm)
   IF (lmurn) THEN
      CALL manage_anhar()
   ELSE
      CALL manage_anhar_anis()
   ENDIF
   RETURN
END SUBROUTINE manage_anhar_routines
