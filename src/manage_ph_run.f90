!
! Copyright (C) 2020 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------
SUBROUTINE manage_ph_run
!--------------------------------------------------
!
!   This is a simple driver that decides if to compute the phonon of all
!   geometries together of one by one
!
USE control_thermo, ONLY : all_geometries_together
IMPLICIT NONE

   IF (all_geometries_together) THEN
      CALL manage_all_geometries_ph()
   ELSE
      CALL manage_ph()
   ENDIF
   RETURN
END SUBROUTINE manage_ph_run

