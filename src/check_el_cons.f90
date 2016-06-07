! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE check_el_cons()
  !-----------------------------------------------------------------------
  !
  !  This routine tries to read the elastic constants from file.
  !  If it finds them, it prints a few information on output and
  !  set the flag el_cons_available.
  !
  USE kinds,         ONLY : DP
  USE mp_images,     ONLY : my_image_id, root_image
  USE io_global,     ONLY : stdout
  USE cell_base,     ONLY : ibrav
  USE elastic_constants, ONLY : read_elastic, el_con, el_compliances, &
                           print_macro_elasticity, print_elastic_constants, &
                           print_elastic_compliances
  USE control_elastic_constants, ONLY : el_cons_available, frozen_ions
  USE control_macro_elasticity, ONLY : macro_el
  USE data_files,    ONLY : fl_el_cons
  !
  IMPLICIT NONE
  !
  LOGICAL  :: exst
  !
  IF ( my_image_id /= root_image ) RETURN
  !
  !  First check if the elastic constants are on file
  !
  CALL read_elastic(fl_el_cons, exst)
  !
  !  If the elastic constants are not available, the routine exits.
  !
  IF (.NOT. exst) THEN
     WRITE(stdout,'(/,5x,"The elastic constants are needed to compute ")')
     WRITE(stdout,'(5x,"thermal expansions from Gruneisen parameters")')
     RETURN
  ELSE
     WRITE(stdout,'(/,5x,"Elastic constants found on file ",a)') &
                                                            TRIM(fl_el_cons)
     el_cons_available=.TRUE.
     CALL print_elastic_constants(el_con, frozen_ions)
     CALL print_elastic_compliances(el_compliances, frozen_ions)
     CALL print_macro_elasticity(ibrav,el_con,el_compliances,macro_el)
  ENDIF

  RETURN
END SUBROUTINE check_el_cons
