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
  !  If it finds one file, it prints a few information on output and
  !  set the flag el_cons_available. Only elastic constants assumed
  !  independent from temperature will be available.
  !  If it finds one file for each geometry it will set the flag
  !  el_cons_t_available. The full temperature dependence of the 
  !  elastic constants will be available.
  !
  !
  USE kinds,         ONLY : DP
  USE thermo_mod,    ONLY : ibrav_geo, celldm_geo
  USE io_global,     ONLY : stdout
  USE cell_base,     ONLY : ibrav
  USE elastic_constants, ONLY : read_elastic, el_con, el_compliances, &
                           print_macro_elasticity, print_elastic_constants, &
                           print_elastic_compliances
  USE control_elastic_constants, ONLY : el_cons_available, frozen_ions, &
                                        el_cons_t_available, el_con_geo
  USE control_macro_elasticity, ONLY : macro_el
  USE control_grun,  ONLY : lb0_t
  USE thermo_mod,    ONLY : tot_ngeo, no_ph
  USE data_files,    ONLY : fl_el_cons
  !
  IMPLICIT NONE
  CHARACTER(LEN=256) :: filelastic
  CHARACTER(LEN=6) :: int_to_char
  INTEGER :: igeo, central_geo
  !
  LOGICAL  :: exst
  !
  ALLOCATE( el_con_geo(6,6,tot_ngeo) )
  IF (.NOT.lb0_t) CALL find_central_geo(tot_ngeo,no_ph,central_geo)

  DO igeo = 1, tot_ngeo
     filelastic='elastic_constants/'//TRIM(fl_el_cons)//'.g'//&
                                              TRIM(int_to_char(igeo))
     CALL read_elastic(filelastic, exst)
     IF (.NOT.exst.AND.igeo==1) THEN
        WRITE(stdout,'(/,5x,"No elastic constants file found.")')
        WRITE(stdout,'(5x,"Some plots might be missing")')
        RETURN
     ELSEIF(.NOT.exst.AND.igeo>1) THEN
         WRITE(stdout,'(/,5x,"File ",i5," with elastic constants missing")')&
                                       igeo
         WRITE(stdout,'(5x,"The dependence on temperature will be neglected")')
         RETURN
     END IF
     el_con_geo(:,:,igeo)=el_con(:,:)
     WRITE(stdout,'(/,5x,"Geometry number",i5)') igeo
     CALL print_elastic_constants(el_con, frozen_ions)
     CALL print_elastic_compliances(el_compliances, frozen_ions)
     CALL print_macro_elasticity(ibrav,el_con,el_compliances,macro_el,.TRUE.)
     IF (igeo==1) el_cons_available=.TRUE.
     IF (.NOT.lb0_t.AND.igeo==central_geo) EXIT
  ENDDO
  el_cons_t_available=.TRUE.

  RETURN
END SUBROUTINE check_el_cons
