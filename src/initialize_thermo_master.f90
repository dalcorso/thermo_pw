!
! Copyright (C) 2013 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE initialize_thermo_master(nwork, part)
  !-----------------------------------------------------------------------
  !
  !  This routine decides how much work there is to do and calculates
  !  the priority of the different jobs
  !
  USE kinds,      ONLY : DP
  USE thermo_mod, ONLY : what, ngeo
  USE thermo_priority, ONLY : npriority, priority, max_priority
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nwork, part
  INTEGER :: igeom
  !
  IF (nwork == 0) RETURN
  !
  IF ( part==1 ) THEN
     SELECT CASE (TRIM(what))
        CASE ('scf', 'scf_ke', 'scf_nk', 'scf_bands', 'scf_ph', 'scf_disp', &
              'mur_lc', 'mur_lc_bands', 'mur_lc_ph', 'mur_lc_disp', &
              'mur_lc_t')
           max_priority=1
        CASE DEFAULT
           CALL errore('initialize_thermo_master','unknown what',1)
     END SELECT
  ELSE IF ( part==2 ) THEN
     SELECT CASE (TRIM(what))
        CASE ('scf_ph', 'scf_bands', 'scf_disp', 'mur_lc', 'mur_lc_bands', &
              'mur_lc_ph', 'mur_lc_t')
           max_priority=1
     END SELECT
  ELSE
     CALL errore('initialize_thermo_master','unknown part',1)
  END IF

  ALLOCATE(npriority(nwork))
  ALLOCATE(priority(nwork,max_priority))

  npriority=0
  priority=0
  IF ( part == 1) THEN
     SELECT CASE (TRIM(what))
        CASE ('scf', 'scf_bands','scf_ph', 'scf_disp', 'mur_lc', &
                     'mur_lc_bands', 'mur_lc_ph', 'mur_lc_disp', 'mur_lc_t')
     END SELECT
  ELSE IF ( part == 2) THEN
     SELECT CASE (TRIM(what))
        CASE ('scf', 'scf_ph', 'scf_disp', 'mur_lc', 'mur_lc_ph', &
              'mur_lc_disp', 'mur_lc_t')
        CASE DEFAULT
           CALL errore('initialize_thermo_master','unknown what',1)
     END SELECT
  ELSE
     CALL errore('initialize_thermo_master','unknown part',1)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE initialize_thermo_master
