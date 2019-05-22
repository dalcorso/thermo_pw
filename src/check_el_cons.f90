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
  !  This routine tries to read the elastic constants from file for all
  !  the geometries considered in the calculation
  !  If it finds one geometry it sets el_cons_available to true.
  !  If it finds all geometries it sets also el_cons_t_available to true.
  !  In the variables el_cons and el_compliances and macro_el it sets
  !  the elastic properties of the central geometry if available or of
  !  the last geometry read if not available.
  !
  USE kinds,             ONLY : DP
  USE io_global,         ONLY : stdout
  USE cell_base,         ONLY : ibrav
  USE elastic_constants, ONLY : read_elastic, el_con, el_compliances, &
                                compute_elastic_compliances,  &
                                print_macro_elasticity,  &
                                print_elastic_constants, &
                                print_elastic_compliances
  USE control_elastic_constants, ONLY : el_cons_available, frozen_ions, &
                                        el_cons_t_available, el_con_geo, &
                                        el_cons_qha_available
  USE control_macro_elasticity, ONLY : macro_el
  USE control_grun,     ONLY : lb0_t
  USE thermo_mod,       ONLY : tot_ngeo, no_ph
  USE data_files,       ONLY : fl_el_cons
  !
  IMPLICIT NONE
  CHARACTER(LEN=256) :: filelastic
  CHARACTER(LEN=6) :: int_to_char
  INTEGER :: igeo, central_geo
  !
  LOGICAL  :: exst, found(tot_ngeo)
  !
  CALL read_elastic_qha()
  IF (el_cons_qha_available) RETURN

  ALLOCATE( el_con_geo(6,6,tot_ngeo) )
  el_con_geo=0.0_DP
  found=.FALSE.
  el_cons_available=.FALSE.
  el_cons_t_available=.FALSE.
  DO igeo = 1, tot_ngeo
     filelastic='elastic_constants/'//TRIM(fl_el_cons)//'.g'//&
                                              TRIM(int_to_char(igeo))
     CALL read_elastic(filelastic, exst)
     IF (.NOT.exst) CYCLE
     found(igeo)=.TRUE.
     el_con_geo(:,:,igeo)=el_con(:,:)
     WRITE(stdout,'(/,5x,"Geometry number",i5," elastic constants found")') igeo
     el_cons_available=.TRUE.
  ENDDO
  IF (.NOT.el_cons_available) RETURN
!
!  If there are the elastic constants of the central geometry we
!  set the el_cons and el_compliances of that geometry. Otherwise
!  they are those of the last geometry read.
!
  el_cons_t_available=.TRUE.
  CALL find_central_geo(tot_ngeo,no_ph,central_geo) 
  IF (found(central_geo)) THEN
     el_con(:,:)=el_con_geo(:,:,central_geo)
     WRITE(stdout,'(/,5x,"Central geometry is number",i5)') central_geo
  ENDIF
  CALL compute_elastic_compliances(el_con,el_compliances)
  CALL print_elastic_constants(el_con, frozen_ions)
  CALL print_elastic_compliances(el_compliances, frozen_ions)
  CALL print_macro_elasticity(ibrav,el_con,el_compliances,macro_el,.TRUE.)
  IF (.NOT.lb0_t) RETURN
!
!  If the code arrives here to have the temperature dependent elastic
!  constants we need the elastic constants at all geometries
! 
  DO igeo=1,tot_ngeo
     el_cons_t_available=el_cons_t_available.AND.found(igeo)
  ENDDO     

  RETURN
END SUBROUTINE check_el_cons
