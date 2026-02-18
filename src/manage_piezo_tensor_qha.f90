
! Copyright (C) 2025 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!----------------------------------------------------------
SUBROUTINE manage_piezo_tensor_qha()
!----------------------------------------------------------
!
!  This routine is similar to manage_piezo_tensor, but it
!  works for the finite temperature case. Internal degrees of 
!  freedom are calculated as a function of temperature minimizing
!  the free energy if the phonons as a function of the internal
!  degrees of freedom are available.
!
USE kinds,             ONLY : DP
USE constants,         ONLY : electron_si, bohr_radius_si
USE thermo_mod,        ONLY : energy_geo, tot_ngeo_eos
USE control_elastic_constants, ONLY : ngeo_strain, ngeom,                 &
                              work_base, el_con_omega_geo,                &
                              start_geometry_qha, last_geometry_qha,      &
                              lelastic, lelasticf, all_geometry_done_geo, &
                              epsil_geo, min_y_t, el_con_at_geo,          &
                              el_con_celldm_geo, el_con_omega_geo,        &
                              min_yf_t, nstep_ec, nint_var_ec, stypec
USE initial_conf,      ONLY : ibrav_save
USE thermo_sym,        ONLY : code_group_save
USE control_atomic_pos, ONLY : max_nint_var
USE elastic_constants, ONLY : epsilon_geo, el_con, el_compliances,         &
                              compute_elastic_constants_ene,               &
                              write_el_cons_on_file  
USE piezoelectric_tensor, ONLY : polar_strain, tot_b_phase, e_piezo_tensor, &
                                 compute_proper_piezo_tensor,               &
                                 write_piezo_tensor_on_file,                &
                                 print_piezo_tensor
USE control_piezoelectric_tensor,  ONLY : lpiezo, lpiezof
USE thermodynamics,       ONLY : e_piezo_tensor_eos_t
USE ph_freq_thermodynamics, ONLY : e_piezo_tensorf_eos_t
USE thermodynamics,    ONLY : ph_free_ener
USE ph_freq_thermodynamics, ONLY : phf_free_ener
USE control_thermo,    ONLY : ltherm_dos, ltherm_freq
USE temperature,       ONLY : ntemp, temp
USE data_files,        ONLY : flanhar
USE io_global,         ONLY : stdout
USE mp_world,          ONLY : world_comm
USE mp,                ONLY : mp_sum

IMPLICIT NONE

REAL(DP), ALLOCATABLE :: free_energy_geo(:), epsilon_geo_loc(:,:,:),        &
                         free_energy_geo_eff(:), epsilon_geo_eff(:,:,:),    &
                         polar_strain_eff(:,:), tot_b_phase_eff(:,:),       &
                         epsil_geo_eff(:), min_y_aux(:,:)
INTEGER :: itemp, startt, lastt, igeom, base_ind, work_base_eff, ipol, jpol, &
           istep, igeo
CHARACTER(LEN=256)  :: filename
CHARACTER(LEN=80)   :: label
CHARACTER(LEN=6) :: int_to_char
REAL(DP) :: fact
LOGICAL :: all_geometry_done

CALL check_all_geometries_done(all_geometry_done)
!IF (.NOT.all_geometry_done) RETURN
!
!  the elastic constants are calculated here
DO igeom=start_geometry_qha, last_geometry_qha
   IF (all_geometry_done_geo(igeom)) GOTO 100
ENDDO
RETURN
100 CONTINUE
!
!  the elastic constants are calculated here
!  for each temperature
!
CALL divide(world_comm, ntemp, startt, lastt)
ALLOCATE(free_energy_geo(work_base))
ALLOCATE(epsilon_geo_loc(3,3,work_base))
ALLOCATE(free_energy_geo_eff(work_base))
ALLOCATE(polar_strain_eff(3,work_base))
ALLOCATE(tot_b_phase_eff(3,work_base))
ALLOCATE(epsilon_geo_eff(3,3,work_base))
ALLOCATE(epsil_geo_eff(work_base))
ALLOCATE(min_y_aux(max_nint_var,ntemp))

min_y_t=0.0_DP
min_yf_t=0.0_DP
e_piezo_tensor_eos_t=0.0_DP
e_piezo_tensorf_eos_t=0.0_DP

DO itemp = startt, lastt
   WRITE(stdout,'(5x, 70("-"))') 
   WRITE(stdout,'(5x,"Computing piezoelectric tensor at temperature",&
                                           &i5, f15.5)') itemp, temp(itemp)
   DO igeom=start_geometry_qha, last_geometry_qha
      base_ind=(igeom-1)*work_base
      IF (.NOT.all_geometry_done_geo(igeom)) CYCLE

      IF (ltherm_dos) THEN
         free_energy_geo(:)=energy_geo(base_ind+1:base_ind+work_base)    &
                          +ph_free_ener(itemp,base_ind+1:base_ind+work_base)
         epsilon_geo_loc(:,:,:)=epsilon_geo(:,:,base_ind+1:base_ind+work_base)
         CALL redefine_energies_qua_t(free_energy_geo, epsilon_geo_loc,  &
                       epsil_geo(base_ind+1), work_base, free_energy_geo_eff,&
                       epsilon_geo_eff, min_y_t(1,1,1,igeom,itemp),      &
                       work_base_eff, igeom)

         CALL redefine_polar_qha_t(polar_strain(:,base_ind+1:            &
                       base_ind+work_base), tot_b_phase(:,base_ind+1:    &
                       base_ind+work_base), work_base, polar_strain_eff, &
                       tot_b_phase_eff, epsil_geo, epsil_geo_eff,        &
                       min_y_t(1,1,1,igeom,itemp))

         CALL compute_proper_piezo_tensor(tot_b_phase_eff,               &
              epsilon_geo_eff, work_base_eff, ngeo_strain,               &
              ibrav_save, code_group_save, el_con_at_geo(:,:,igeom) )
         e_piezo_tensor=e_piezo_tensor * el_con_celldm_geo(1,igeom) / &
                                           el_con_omega_geo(igeom)
         label="Proper total piezoelectric tensor gamma_ij [ C/m^2 ]"
         fact= electron_si / (bohr_radius_si)**2
         CALL print_piezo_tensor(e_piezo_tensor, fact, label, .FALSE.)

         e_piezo_tensor_eos_t(:,:,itemp,igeom) = e_piezo_tensor(:,:)
      ENDIF

      IF (ltherm_freq) THEN
         free_energy_geo(:)=energy_geo(base_ind+1:base_ind+work_base)    &
                          +phf_free_ener(itemp,base_ind+1:base_ind+work_base)
         epsilon_geo_loc(:,:,:)=epsilon_geo(:,:,base_ind+1:base_ind+work_base)
         CALL redefine_energies_qua_t(free_energy_geo, epsilon_geo_loc,  &
                       epsil_geo(base_ind+1), work_base, free_energy_geo_eff,&
                       epsilon_geo_eff, min_y_t(1,1,1,igeom,itemp),      &
                       work_base_eff, igeom)

         CALL redefine_polar_qha_t(polar_strain(:,base_ind+1:            &
                       base_ind+work_base), tot_b_phase(:,base_ind+1:    &
                       base_ind+work_base), work_base, polar_strain_eff, &
                       tot_b_phase_eff, epsil_geo, epsil_geo_eff,        &
                       min_y_t(1,1,1,igeom,itemp))

         CALL compute_proper_piezo_tensor(tot_b_phase_eff,               &
              epsilon_geo_eff, work_base_eff, ngeo_strain,               &
              ibrav_save, code_group_save, el_con_at_geo(:,:,igeom) )


         e_piezo_tensor=e_piezo_tensor * el_con_celldm_geo(1,igeom) / &
                                           el_con_omega_geo(igeom)
         label="Proper total piezoelectric tensor gamma_ij [ C/m^2 ]"
         fact= electron_si / (bohr_radius_si)**2
         CALL print_piezo_tensor(e_piezo_tensor, fact, label, .FALSE.)

         e_piezo_tensorf_eos_t(:,:,itemp,igeom) = e_piezo_tensor(:,:)

      ENDIF 
   ENDDO
ENDDO
!
!  Now collect the results on all processors and write them on file.
!
DO igeom=start_geometry_qha, last_geometry_qha
   IF (ltherm_dos) THEN
      CALL mp_sum(e_piezo_tensor_eos_t(:,:,:,igeom), world_comm)
      CALL add_geometry_number('anhar_files/', TRIM(flanhar)//&
                             '.e_piezo', filename, igeom)
      lpiezo=.TRUE.
      CALL write_piezo_tensor_on_file(temp, ntemp, ibrav_save,          &
                       code_group_save, e_piezo_tensor_eos_t(:,:,:,igeom),  &
                                              filename, 0, 1)
      DO istep=1, nstep_ec
         DO igeo=1, ngeo_strain
            IF (stypec(istep)) THEN
               min_y_aux(:,:)=min_y_t(:,igeo,istep,igeom,:)
               CALL mp_sum(min_y_aux, world_comm)
               CALL add_geometry_number('anhar_files/', TRIM(flanhar)//      &
                             '.uint.ffem', filename, igeom)
               filename=TRIM(filename)//'.'//TRIM(int_to_char(igeo))//'.'//&
                   TRIM(int_to_char(istep))
               CALL write_only_uint(min_y_aux, max_nint_var, nint_var_ec(istep),&
                                                  temp, ntemp, filename, 0 )
            ENDIF
         ENDDO
      ENDDO
   ENDIF

   IF (ltherm_freq) THEN
      CALL mp_sum(e_piezo_tensorf_eos_t(:,:,:,igeom), world_comm)
      CALL add_geometry_number('anhar_files/', TRIM(flanhar)//&
                              '.e_piezo', filename, igeom)
      filename=TRIM(filename)//'_ph'
      lpiezof=.TRUE.
      CALL write_piezo_tensor_on_file(temp, ntemp, ibrav_save,               &
                       code_group_save, e_piezo_tensorf_eos_t(:,:,:,igeom),  &
                                              filename, 0, 1)
      DO istep=1, nstep_ec
         DO igeo=1, ngeo_strain
            IF (stypec(istep)) THEN
               min_y_aux(:,:)=min_yf_t(:,igeo,istep,igeom,:)
               CALL mp_sum(min_y_aux, world_comm)
               CALL add_geometry_number('anhar_files/', TRIM(flanhar)//      &
                             '.uint.ffem', filename, igeom)
               filename=TRIM(filename)//'.'//TRIM(int_to_char(igeo))//'.'//&
                   TRIM(int_to_char(istep))//'_ph'
               CALL write_only_uint(min_y_aux, max_nint_var, nint_var_ec(istep),&
                                     temp, ntemp, filename, 0 )
            ENDIF
         ENDDO
      ENDDO
   ENDIF
ENDDO

DEALLOCATE(min_y_aux)
DEALLOCATE(free_energy_geo)
DEALLOCATE(epsilon_geo_loc)
DEALLOCATE(free_energy_geo_eff)
DEALLOCATE(epsilon_geo_eff)
DEALLOCATE(epsil_geo_eff)
DEALLOCATE(polar_strain_eff)
DEALLOCATE(tot_b_phase_eff)

CALL plot_piezo_t(0)

RETURN
END SUBROUTINE manage_piezo_tensor_qha

