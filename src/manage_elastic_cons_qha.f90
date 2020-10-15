!
! Copyright (C) 2019 Cristiano Malica and Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------
SUBROUTINE manage_elastic_cons_qha()
!----------------------------------------------------------
!
!  This routine is similar to manage_elastic_cons, but it
!  works for the finite temperature case.
!
USE kinds,             ONLY : DP
USE thermo_mod,        ONLY : energy_geo, tot_ngeo 
USE control_elastic_constants, ONLY : ngeo_strain, elcpvar, ngeom, &
                              work_base, el_con_omega_geo,         &
                              start_geometry_qha, last_geometry_qha
USE initial_conf,      ONLY : ibrav_save
USE thermo_sym,        ONLY : laue
USE elastic_constants, ONLY : epsilon_geo, el_con, el_compliances,         &
                              compute_elastic_constants_ene,               &
                              write_el_cons_on_file  
USE thermodynamics,    ONLY : ph_free_ener
USE ph_freq_thermodynamics, ONLY : phf_free_ener
USE anharmonic,        ONLY : el_cons_t, el_comp_t, b0_t, lelastic
USE ph_freq_anharmonic,ONLY : el_consf_t, el_compf_t, b0f_t, lelasticf
USE control_thermo,    ONLY : ltherm_dos, ltherm_freq
USE temperature,       ONLY : ntemp, temp
USE data_files,        ONLY : flanhar
USE io_global,         ONLY : stdout
USE mp_world,          ONLY : world_comm
USE mp,                ONLY : mp_sum

IMPLICIT NONE

REAL(DP), ALLOCATABLE :: free_energy_geo(:), epsilon_geo_loc(:,:,:)
INTEGER :: itemp, startt, lastt, igeom, base_ind
CHARACTER(LEN=256)  :: filelastic, filename
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: all_geometry_done, exst, check_file_exists

CALL check_all_geometries_done(all_geometry_done)
IF (.NOT.all_geometry_done) RETURN

filename='anhar_files/'//TRIM(flanhar)//'.celldm'
exst=check_file_exists(filename)

IF (exst) THEN
   CALL manage_elastic_cons_qha_2()
ELSE
!
!  the elastic constants are calculated here
!  for each temperature
!
CALL divide(world_comm, ntemp, startt, lastt)
ALLOCATE(free_energy_geo(work_base))
ALLOCATE(epsilon_geo_loc(3,3,work_base))

DO igeom=start_geometry_qha, last_geometry_qha
   base_ind=(igeom-1)*work_base
   el_cons_t=0.0_DP
   el_consf_t=0.0_DP
   DO itemp = startt, lastt
      WRITE(stdout,'(5x, 70("-"))') 
      WRITE(stdout,'(5x,"Computing elastic constants at temperature",&
                                           &i5, f15.5)') itemp, temp(itemp)
      IF (ltherm_dos) THEN
         WRITE(stdout,'(5x,"From vibrational density of states")')
         free_energy_geo(:)=energy_geo(base_ind+1:base_ind+work_base) &
                          +ph_free_ener(itemp,base_ind+1:base_ind+work_base)
         epsilon_geo_loc(:,:,:)=epsilon_geo(:,:,base_ind+1:base_ind+work_base)
         CALL compute_elastic_constants_ene(free_energy_geo, epsilon_geo_loc, &
                            work_base, ngeo_strain, ibrav_save, laue,         &
                            el_con_omega_geo(igeom), elcpvar)
         el_cons_t(:,:,itemp) = el_con(:,:)
      ENDIF
      IF (ltherm_freq) THEN
         WRITE(stdout, '(5x, "Using Brillouin zone integrals")')
         free_energy_geo(:)=energy_geo(base_ind+1:base_ind+work_base) &
                          +phf_free_ener(itemp,base_ind+1:base_ind+work_base)
         epsilon_geo_loc(:,:,:)=epsilon_geo(:,:,base_ind+1:base_ind+work_base)
         CALL compute_elastic_constants_ene(free_energy_geo, epsilon_geo_loc, &
                         work_base, ngeo_strain, ibrav_save, laue,       &
                         el_con_omega_geo(igeom), elcpvar)
         el_consf_t(:,:,itemp) = el_con(:,:)
      ENDIF 
   ENDDO
!
!  Now collect the results on all processors and write them on file.
!
   IF (ltherm_dos) THEN
      CALL mp_sum(el_cons_t, world_comm)
      CALL compute_el_comp_t(el_cons_t,el_comp_t,b0_t)
      CALL add_geometry_number('anhar_files/', TRIM(flanhar)//'.el_cons', &
                                              filelastic, igeom)
      CALL write_el_cons_on_file(temp, ntemp, ibrav_save, laue, el_cons_t, &
                                              b0_t, filelastic, 0)
      CALL add_geometry_number('anhar_files/', TRIM(flanhar)//'.el_comp', &
                                              filelastic, igeom)
      CALL write_el_cons_on_file(temp, ntemp, ibrav_save, laue, el_comp_t, &
                                              b0_t, filelastic, 1)
      lelastic=.TRUE.
   ENDIF

   IF (ltherm_freq) THEN
      CALL mp_sum(el_consf_t, world_comm)
      CALL compute_el_comp_t(el_consf_t,el_compf_t,b0f_t)
      CALL add_geometry_number('anhar_files/', TRIM(flanhar)//'.el_cons', &
                                filelastic, igeom)
      filelastic=TRIM(filelastic)//'_ph'
      CALL write_el_cons_on_file(temp, ntemp, ibrav_save, laue, el_consf_t, &
                                                    b0f_t, filelastic, 0)

      CALL add_geometry_number('anhar_files/', TRIM(flanhar)//'.el_comp', &
                                filelastic, igeom)
      filelastic=TRIM(filelastic)//'_ph'
      CALL write_el_cons_on_file(temp, ntemp, ibrav_save, laue, el_compf_t, &
                                                    b0f_t, filelastic, 1)
      lelasticf=.TRUE.
   ENDIF
ENDDO

DEALLOCATE(free_energy_geo)
DEALLOCATE(epsilon_geo_loc)

CALL plot_elastic_t(0,.FALSE.)
CALL plot_elastic_t(1,.FALSE.)

ENDIF

RETURN
END SUBROUTINE manage_elastic_cons_qha
