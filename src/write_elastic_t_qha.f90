!
! Copyright (C) 2019 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_elastic_t_qha()

USE kinds,             ONLY : DP
USE thermo_mod,        ONLY : energy_geo, tot_ngeo
USE control_elastic_constants, ONLY : ngeo_strain, frozen_ions, elcpvar
USE control_elastic_constants_qha, ONLY : ngeom, work_base, omega0_qha
USE initial_conf,      ONLY : ibrav_save
USE thermo_sym,        ONLY : laue
USE elastic_constants, ONLY : print_elastic_constants, epsilon_geo,        &
                              el_con, el_compliances,                      &
                              compute_elastic_compliances,                 &
                              print_macro_elasticity,                      &
                              compute_elastic_constants_ene
USE equilibrium_conf,  ONLY : omega0
USE control_macro_elasticity, ONLY : macro_el
USE temperature,       ONLY : ntemp, temp
USE data_files,        ONLY : flanhar
USE io_global,         ONLY : stdout
USE mp_world,          ONLY : world_comm
USE mp,                ONLY : mp_sum
USE thermodynamics,        ONLY : ph_free_ener
USE ph_freq_thermodynamics,  ONLY : phf_free_ener
USE anharmonic,        ONLY : el_cons_t, el_comp_t, b0_t, lelastic
USE ph_freq_anharmonic,ONLY : el_consf_t, el_compf_t, b0f_t, lelasticf
USE control_thermo,    ONLY : ltherm_dos, ltherm_freq

IMPLICIT NONE

REAL(DP), ALLOCATABLE :: free_energy_geo(:), epsilon_geo_loc(:,:,:)
INTEGER :: itemp, startt, lastt, igeom, base_ind
CHARACTER(LEN=256)  :: filelastic
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: all_geometry_done

CALL check_all_geometry_done(all_geometry_done)
IF (.NOT.all_geometry_done) RETURN
!
!  the elastic constants are calculated here
!  for each temperature
!
CALL divide(world_comm, ntemp, startt, lastt)
ALLOCATE(free_energy_geo(work_base))
ALLOCATE(epsilon_geo_loc(3,3,work_base))

DO igeom=1, ngeom 

   base_ind=(igeom-1)*work_base
   el_cons_t=0.0_DP
   el_comp_t=0.0_DP
   b0_t=0.0_DP
   el_consf_t=0.0_DP
   el_compf_t=0.0_DP
   b0f_t=0.0_DP
   DO itemp = startt, lastt
      WRITE(stdout,'(5x, 70("-"))') 
      WRITE(stdout,*) 'Computing elastic constants at temperature', itemp, &
                                                                temp(itemp)
      IF (ltherm_dos) THEN
         free_energy_geo(:)=energy_geo(base_ind+1:base_ind+work_base) &
                          +ph_free_ener(itemp,base_ind+1:base_ind+work_base)
         epsilon_geo_loc(:,:,:)=epsilon_geo(:,:,base_ind+1:base_ind+work_base)
         CALL compute_elastic_constants_ene(free_energy_geo, epsilon_geo_loc, &
                            tot_ngeo, ngeo_strain, ibrav_save, laue,          &
                            omega0_qha(igeom), elcpvar)
         CALL compute_elastic_compliances(el_con,el_compliances)
         CALL print_macro_elasticity(ibrav_save,el_con, el_compliances,&
                                                   macro_el,.FALSE.)
         el_cons_t(:,:,itemp) = el_con(:,:)
         el_comp_t(:,:,itemp) = el_compliances(:,:)
         b0_t(itemp)=macro_el(5)
         CALL print_elastic_constants(el_con, frozen_ions)
      ENDIF

      IF (ltherm_freq) THEN
         free_energy_geo(:)=energy_geo(base_ind+1:base_ind+work_base) &
                          +ph_free_ener(itemp,base_ind+1:base_ind+work_base)
         epsilon_geo_loc(:,:,:)=epsilon_geo(:,:,base_ind+1:base_ind+work_base)
         CALL compute_elastic_constants_ene(free_energy_geo, epsilon_geo_loc, &
                         tot_ngeo, ngeo_strain, ibrav_save, laue,       &
                         omega0_qha(igeom), elcpvar)
         CALL compute_elastic_compliances(el_con,el_compliances)
         CALL print_macro_elasticity(ibrav_save,el_con, el_compliances, &
                                                macro_el,.FALSE.)
         el_consf_t(:,:,itemp) = el_con(:,:)
         el_compf_t(:,:,itemp) = el_compliances(:,:)
         b0f_t(itemp)=macro_el(5)
         CALL print_elastic_constants(el_con, frozen_ions)
      ENDIF 
   ENDDO
!
!  Now collect the results on all processors and write them on file.
!
   IF (ltherm_dos) THEN
      CALL mp_sum(el_cons_t, world_comm)
      CALL mp_sum(el_comp_t, world_comm)
      CALL mp_sum(b0_t, world_comm)
      lelastic=.TRUE.
      filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons.g'//&
                                                  TRIM(int_to_char(igeom))
      CALL write_el_cons_on_file(temp, ntemp, ibrav_save, el_cons_t, b0_t, &
                                                       filelastic, 0)
      filelastic='anhar_files/'//TRIM(flanhar)//'.el_comp.g'//&
                                                  TRIM(int_to_char(igeom))
      CALL write_el_cons_on_file(temp, ntemp, ibrav_save, el_comp_t, b0_t, &
                                                       filelastic, 1)
   ENDIF

   IF (ltherm_freq) THEN
      CALL mp_sum(el_consf_t, world_comm)
      CALL mp_sum(el_compf_t, world_comm)
      CALL mp_sum(b0f_t, world_comm)
      lelasticf=.TRUE.
      filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons.g'//&
                                           TRIM(int_to_char(igeom))//'_ph'
      CALL write_el_cons_on_file(temp, ntemp, ibrav_save, el_consf_t, b0f_t, &
                                                          filelastic, 0)

      filelastic='anhar_files/'//TRIM(flanhar)//'.el_comp.g'//&
                                           TRIM(int_to_char(igeom))//'_ph'
      CALL write_el_cons_on_file(temp, ntemp, ibrav_save, el_compf_t, b0f_t, &
                                                           filelastic,1)
   ENDIF
ENDDO

DEALLOCATE(free_energy_geo)
DEALLOCATE(epsilon_geo_loc)

RETURN
END SUBROUTINE write_elastic_t_qha
