! Copyright (C) 2019 C. Malica
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE check_el_cons_qha()
  !-----------------------------------------------------------------------
  !
  !  This routine tries to read the temperature dependent elastic constants 
  !  from file obtained with the option 'elastic_constants_qha' for all the 
  !  geometries considered in the calculation. 
  !  If it finds one geometry it sets el_cons_qha_available to true.
  !  If it finds all geometries it sets also el_cons_qha_t_available to
  !  true.
  !  In the variables el_cons and el_compliances and macro_el it sets
  !  the elastic properties of the central geometry if available or of
  !  the last geometry read if not available.
  !
  USE kinds,             ONLY : DP
  USE io_global,         ONLY : stdout
  USE cell_base,         ONLY : ibrav
  USE temperature,       ONLY : ntemp, temp
  USE elastic_constants, ONLY : read_elastic, el_con, el_compliances, &
                                compute_elastic_compliances,  &
                                print_macro_elasticity,  &
                                print_elastic_constants, &
                                print_elastic_compliances

  USE anharmonic,         ONLY : el_cons_t, el_comp_t, b0_t, el_con_geo_t
  USE ph_freq_anharmonic, ONLY : el_consf_t, el_compf_t, b0f_t, el_con_geo_t_ph

  USE control_elastic_constants, ONLY : frozen_ions, &
                                        el_cons_qha_available,     &
                                        el_cons_qha_available_ph,  &
                                        el_cons_qha_geo_available, &
                                        el_cons_qha_geo_available_ph
 
  USE control_thermo, ONLY : ltherm_dos, ltherm_freq 
  USE data_files,     ONLY : flanhar 
  USE control_macro_elasticity, ONLY : macro_el
  USE control_grun,     ONLY : lb0_t
  USE thermo_mod,       ONLY : tot_ngeo, no_ph
  USE data_files,       ONLY : fl_el_cons
  !
  IMPLICIT NONE
  CHARACTER(LEN=256) :: filelastic, filelastic_ph 
  CHARACTER(LEN=6) :: int_to_char
  INTEGER :: igeo, central_geo, i, j
  !
  LOGICAL  :: check_file_exists
  LOGICAL  :: exst, found_dos(tot_ngeo), found_ph(tot_ngeo) 
  !
  el_cons_qha_available=.FALSE.
  el_cons_qha_available_ph=.FALSE.
  found_ph=.FALSE.  
  found_dos=.FALSE.  

  IF (ltherm_dos) THEN
     el_con_geo_t=0.0_DP
     DO igeo = 1, tot_ngeo
        filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons'//&
                                                 TRIM(int_to_char(igeo))
        exst=check_file_exists(filelastic)
        IF (.NOT.exst) CYCLE

        !The loop on temperatures is inside read_el_cons_from_file
        CALL read_el_cons_from_file(temp, ntemp, ibrav, &
                  el_con_geo_t(:,:,:,igeo), b0_t(:), filelastic)

        found_dos(igeo)=.TRUE.
        el_cons_t(:,:,:)=el_con_geo_t(:,:,:,igeo)

        WRITE(stdout,'(/,5x,"Geometry number",i5," of elastic constants &
                   computed via ltherm_dos found")') igeo
        
        el_cons_qha_available=.TRUE.
     END DO
  END IF

  IF (ltherm_freq) THEN
     el_con_geo_t_ph=0.0_DP
     DO igeo=1, tot_ngeo
        filelastic_ph='anhar_files/'//TRIM(flanhar)//'.el_cons_ph'//&
                                                 TRIM(int_to_char(igeo))
        exst=check_file_exists(filelastic_ph)
        IF (.NOT.exst) CYCLE

        !The loop on temperatures is inside read_el_cons_from_file
        CALL read_el_cons_from_file(temp, ntemp, ibrav,& 
                el_con_geo_t_ph(:,:,:,igeo), b0f_t(:), filelastic_ph)
      
        found_ph(igeo)=.TRUE.
        el_consf_t(:,:,:)=el_con_geo_t_ph(:,:,:,igeo)

        WRITE(stdout,'(/,5x,"Geometry number",i5," of elastic constants &
                   computed via ltherm_freq found")') igeo
    
        el_cons_qha_available_ph=.TRUE.
     END DO 
  END IF
 
  IF ((.NOT.el_cons_qha_available).AND.(.NOT.el_cons_qha_available_ph)) RETURN
!
!
!  If there are the elastic constants of the central geometry we
!  set the el_cons and el_compliances of that geometry. Otherwise
!  they are those of the last geometry read.
!
!
  CALL find_central_geo(tot_ngeo,no_ph,central_geo)
  WRITE(stdout,'(/,5x,"Central geometry is number",i5)') central_geo

  IF (ltherm_dos) THEN
     IF (found_dos(central_geo)) el_cons_t(:,:,:)=&
                                    el_con_geo_t(:,:,:,central_geo)
     CALL compute_elastic_compliances_t(el_cons_t,el_comp_t, b0_t)
  ENDIF

  IF (ltherm_freq) THEN
     IF (found_ph(central_geo)) el_consf_t(:,:,:)= &
                                    el_con_geo_t_ph(:,:,:,central_geo)
     CALL compute_elastic_compliances_t(el_consf_t,el_compf_t, b0f_t)
  ENDIF
!
!  If the code arrives here to have the temperature dependent elastic
!  constants we need the elastic constants at all geometries
! 
  el_cons_qha_geo_available=.TRUE.
  el_cons_qha_geo_available_ph=.TRUE.
  DO igeo=1,tot_ngeo
     el_cons_qha_geo_available=el_cons_qha_geo_available.AND.found_dos(igeo)
     el_cons_qha_geo_available_ph=el_cons_qha_geo_available_ph.AND.&
                                                              found_ph(igeo)
  ENDDO
!
  RETURN
END SUBROUTINE check_el_cons_qha

SUBROUTINE compute_elastic_compliances_t(el_cons_t, el_comp_t, b0_t)

 USE kinds,              ONLY : DP
 USE cell_base,          ONLY : ibrav
 USE temperature,        ONLY : ntemp
 USE elastic_constants,  ONLY : compute_elastic_compliances, &
                                print_macro_elasticity
 USE mp_world,           ONLY : world_comm
 USE mp,                 ONLY : mp_sum

 IMPLICIT NONE

 REAL(DP) :: el_cons_t(6,6,ntemp), el_comp_t(6,6,ntemp), b0_t(ntemp)
 REAL(DP) :: macro_el(8)
 INTEGER  :: itemp, startt, lastt
 INTEGER  :: i, j

 CALL divide(world_comm, ntemp, startt, lastt)
 el_comp_t=0.0_DP
 b0_t=0.0_DP
 DO itemp=startt,lastt
    IF (itemp==1.OR.itemp==ntemp) CYCLE
    CALL compute_elastic_compliances(el_cons_t(:,:,itemp),el_comp_t(:,:,itemp))
    CALL print_macro_elasticity(ibrav,el_cons_t(:,:,itemp), &
                          el_comp_t(:,:,itemp),macro_el,.FALSE.)
    b0_t(itemp)=macro_el(5)
 ENDDO

 CALL mp_sum(el_comp_t, world_comm)
 CALL mp_sum(b0_t, world_comm)

 RETURN
END SUBROUTINE compute_elastic_compliances_t
