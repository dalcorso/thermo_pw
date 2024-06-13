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
                              start_geometry_qha, last_geometry_qha, &
                              lelastic, lelasticf, all_geometry_done_geo, &
                              epsil_geo, min_y_t, stype
USE initial_conf,      ONLY : ibrav_save
USE thermo_sym,        ONLY : laue
USE elastic_constants, ONLY : epsilon_geo, el_con, el_compliances,         &
                              compute_elastic_constants_ene,               &
                              write_el_cons_on_file  
USE thermodynamics,    ONLY : ph_free_ener
USE control_thermo,  ONLY : lstress
USE ph_freq_thermodynamics, ONLY : phf_free_ener
USE el_thermodynamics, ONLY : el_ener, el_free_ener, el_entr, el_ce
USE anharmonic,        ONLY : el_cons_t, el_comp_t, b0_t
USE ph_freq_anharmonic,ONLY : el_consf_t, el_compf_t, b0f_t
USE control_thermo,    ONLY : ltherm_dos, ltherm_freq
USE control_eldos,     ONLY : lel_free_energy
USE temperature,       ONLY : ntemp, temp
USE data_files,        ONLY : flanhar, fleltherm
USE io_global,         ONLY : stdout
USE mp_world,          ONLY : world_comm
USE mp,                ONLY : mp_sum

IMPLICIT NONE

REAL(DP), ALLOCATABLE :: free_energy_geo(:), epsilon_geo_loc(:,:,:), &
                         free_energy_geo_eff(:), epsilon_geo_eff(:,:,:)
INTEGER :: itemp, startt, lastt, igeom, base_ind, iwork, work_base_eff
CHARACTER(LEN=256)  :: filelastic, filename, filedata
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: all_geometry_done, all_el_free, exst, ldummy, check_file_exists, &
           all_found, run

CALL check_all_geometries_done(all_geometry_done)
!IF (.NOT.all_geometry_done) RETURN
IF (lel_free_energy) THEN
   CALL check_all_el_free_ener_done(all_el_free)
   IF (.NOT.all_el_free) CALL errore('manage_anhar',&
                        'missing electron thermodynamics',1)
   DO igeom=1, tot_ngeo
      CALL set_el_files_names(igeom)
      filedata="therm_files/"//TRIM(fleltherm)
      CALL read_thermo(ntemp, temp, el_ener(:,igeom),           &
                       el_free_ener(:,igeom), el_entr(:,igeom), &
                       el_ce(:,igeom), ldummy, filedata)
   ENDDO
   CALL restore_el_file_names()
ENDIF


filename='anhar_files/'//TRIM(flanhar)//'.celldm'
exst=check_file_exists(filename)

IF (exst.AND..NOT.ANY(stype)) THEN
   CALL manage_elastic_cons_qha_2()
ELSE
!
!  the elastic constants are calculated here
!  for each temperature
!
CALL divide(world_comm, ntemp, startt, lastt)
ALLOCATE(free_energy_geo(work_base))
ALLOCATE(epsilon_geo_loc(3,3,work_base))
ALLOCATE(free_energy_geo_eff(work_base))
ALLOCATE(epsilon_geo_eff(3,3,work_base))
ALLOCATE(lstress(work_base*ngeom))
lstress=.FALSE.
min_y_t=0.0_DP

DO igeom=1, ngeom
   base_ind=(igeom-1)*work_base
   el_cons_t=0.0_DP
   el_consf_t=0.0_DP
   IF (.NOT.all_geometry_done_geo(igeom)) CYCLE
!
!  If for some geometry the energy is not available, try to read it
!  from the restart file. If not available for some geometry the elastic
!  constants are not computed for this igeom.
!
   all_found=.TRUE.
   DO iwork=1,work_base
      run=.FALSE.
      IF (energy_geo(base_ind+iwork)==0.0_DP) &
           CALL check_existence(base_ind+iwork,1,run)
      all_found=all_found.AND..NOT.run
   ENDDO
   IF (.NOT.all_found) CYCLE
   DO itemp = startt, lastt
      WRITE(stdout,'(5x, 70("-"))') 
      WRITE(stdout,'(5x,"Computing elastic constants at temperature",&
                                           &i5, f15.5)') itemp, temp(itemp)
      IF (ltherm_dos) THEN
         WRITE(stdout,'(5x,"From vibrational density of states")')
         free_energy_geo(:)=energy_geo(base_ind+1:base_ind+work_base) &
                          +ph_free_ener(itemp,base_ind+1:base_ind+work_base)
         IF (lel_free_energy) free_energy_geo(:)=free_energy_geo(:) + &
                           el_free_ener(itemp,base_ind+1:base_ind+work_base)
         epsilon_geo_loc(:,:,:)=epsilon_geo(:,:,base_ind+1:base_ind+work_base)
         CALL redefine_energies_qua_t(free_energy_geo, epsilon_geo_loc, &
                       epsil_geo(base_ind+1), work_base, free_energy_geo_eff, &
                       epsilon_geo_eff, work_base_eff, igeom, itemp)
         CALL compute_elastic_constants_ene(free_energy_geo_eff,  &
                            epsilon_geo_eff, work_base_eff, ngeo_strain, &
                            ibrav_save, laue, el_con_omega_geo(igeom), elcpvar)
         el_cons_t(:,:,itemp) = el_con(:,:)
      ENDIF
      IF (ltherm_freq) THEN
         WRITE(stdout, '(5x, "Using Brillouin zone integrals")')
         free_energy_geo(:)=energy_geo(base_ind+1:base_ind+work_base) &
                          +phf_free_ener(itemp,base_ind+1:base_ind+work_base)
         IF (lel_free_energy) free_energy_geo(:)=free_energy_geo(:) + &
                           el_free_ener(itemp,base_ind+1:base_ind+work_base)
         epsilon_geo_loc(:,:,:)=epsilon_geo(:,:,base_ind+1:base_ind+work_base)
         CALL redefine_energies_qua_t(free_energy_geo, epsilon_geo_loc, &
                       epsil_geo(base_ind+1), work_base, free_energy_geo_eff, &
                       epsilon_geo_eff, work_base_eff, igeom, itemp)
         CALL compute_elastic_constants_ene(free_energy_geo_eff, &
                         epsilon_geo_eff, work_base_eff, ngeo_strain, & 
                         ibrav_save, laue, el_con_omega_geo(igeom), elcpvar)
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
IF (ANY(stype)) THEN
   CALL mp_sum(min_y_t, world_comm)
   CALL write_min_y()
ENDIF

DEALLOCATE(free_energy_geo)
DEALLOCATE(epsilon_geo_loc)
DEALLOCATE(free_energy_geo_eff)
DEALLOCATE(epsilon_geo_eff)
DEALLOCATE(lstress)

CALL plot_elastic_t(0,.FALSE.)
CALL plot_elastic_t(1,.FALSE.)

ENDIF

RETURN
END SUBROUTINE manage_elastic_cons_qha
