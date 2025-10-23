!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------
SUBROUTINE manage_elastic_cons(nwork,ngeom)
!----------------------------------------------------------
!
!   This routine coodinates the work to compute the T=0 K
!   elastic constants after the scf calculations of the total
!   energy or stress of all the strained geometries have been
!   done. ngeom is the number of unperturbed geometries for which
!   we calculate the elastic constants. 
!
!   It coordinates also the calculation of material quantities
!   such as the sound velocity which depend on the elastic
!   constant.
!

USE kinds,             ONLY : DP
USE thermo_mod,        ONLY : energy_geo, density, start_geometry,         &
                              last_geometry
USE control_elastic_constants, ONLY : ngeo_strain, frozen_ions,            &
                              elastic_algorithm, rot_mat, elcpvar,         &
                              el_con_omega_geo, elalgen, start_geometry_qha, &
                              last_geometry_qha, stype, min_y, nstep_ec, &
                              lcm_ec, atom_dir, move_at, epsil_y, epsil_geo, &
                              old_ec
USE initial_conf,      ONLY : ibrav_save
USE thermo_sym,        ONLY : laue
USE elastic_constants, ONLY : print_elastic_constants,                     &
                              compute_elastic_constants, epsilon_geo,      &
                              sigma_geo, el_con, el_compliances,           &
                              compute_elastic_compliances,                 &
                              print_elastic_compliances,                   &
                              write_elastic, print_macro_elasticity,       &
                              compute_elastic_constants_ene,               &
                              print_sound_velocities
USE rotate,            ONLY : rotate_tensors2
USE control_macro_elasticity, ONLY : macro_el, vp, vb, vg, approx_debye_t
USE debye_module,      ONLY : compute_debye_temperature,                   &
                              compute_debye_temperature_poisson
USE control_debye,     ONLY : debye_t
USE data_files,        ONLY : fl_el_cons
USE ions_base,         ONLY : nat
USE cell_base,         ONLY : omega
USE io_global,         ONLY : stdout
USE mp,                ONLY : mp_sum
USE mp_images,         ONLY : my_image_id, root_image, nproc_image
USE mp_world,          ONLY : world_comm

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork, ngeom

REAL(DP) :: poisson, bulkm
REAL(DP), ALLOCATABLE :: sigma_geo_aux(:,:,:), epsilon_geo_eff(:,:,:), &
                                               energy_geo_eff(:)
INTEGER :: iwork, work_base, igeom, base_ind, nwork_eff, istep, igeo
CHARACTER(LEN=6)    :: int_to_char
CHARACTER(LEN=256)  :: filelastic
!
!  First collect the total energies
!
CALL mp_sum(energy_geo, world_comm)
energy_geo=energy_geo / nproc_image
!
!  the elastic constants are calculated here if we have the energies
!  of all geometries
!
lreturn=.FALSE.
DO iwork=1,nwork
   lreturn=lreturn.OR.(ABS(energy_geo(iwork))<1.D-10)
ENDDO
IF (lreturn) RETURN

ALLOCATE(energy_geo_eff(nwork))
ALLOCATE(epsilon_geo_eff(3,3,nwork))

CALL redefine_energies_qua(energy_geo, epsilon_geo, epsil_geo, nwork,  &
                       energy_geo_eff, epsilon_geo_eff, nwork_eff)
!
!  Then collect the stress if it has been calculated
!
IF (.NOT.elalgen) THEN
   CALL mp_sum(sigma_geo, world_comm)
   sigma_geo=sigma_geo / nproc_image
ENDIF
!
!  work_base is the nwork to compute the elastic constants of one geometry
!
work_base=nwork_eff/ngeom
ALLOCATE(sigma_geo_aux(3,3,work_base))
DO igeom=start_geometry_qha, last_geometry_qha
   base_ind=(igeom-1)*work_base
   IF (elastic_algorithm=='standard') THEN
      sigma_geo_aux(1:3,1:3,1:work_base)=sigma_geo(1:3,1:3, &
                                         base_ind+1:base_ind+work_base)
   ELSEIF (elastic_algorithm=='advanced') THEN
      DO iwork=1,work_base
         CALL rotate_tensors2(rot_mat(1,1,base_ind+iwork), 1, &
                 sigma_geo(1,1,base_ind+iwork), sigma_geo_aux(1,1,iwork),-1)
      ENDDO
   ENDIF

   IF (.NOT.elalgen) THEN
      CALL compute_elastic_constants(sigma_geo_aux, &
                       epsilon_geo(1,1,base_ind+1), &
                       work_base, ngeo_strain, ibrav_save, laue, elcpvar)
   ELSE
      CALL compute_elastic_constants_ene(energy_geo_eff(base_ind+1), &
           epsilon_geo_eff(1,1,base_ind+1), work_base, ngeo_strain,  &
                    ibrav_save, laue, el_con_omega_geo(igeom), elcpvar, old_ec)
   END IF
   CALL print_elastic_constants(el_con, frozen_ions)
!
!  now compute the elastic compliances and prints them
!
   CALL compute_elastic_compliances(el_con,el_compliances)
   CALL print_elastic_compliances(el_compliances, frozen_ions)
   CALL print_macro_elasticity(ibrav_save,el_con,el_compliances,&
                                                       macro_el,.TRUE.)
!
!  here compute the sound velocities, using the density of the solid and
!  the elastic constants
!
   CALL print_sound_velocities( ibrav_save, el_con, el_compliances, &
                                       density, vp, vb, vg, .TRUE.)
!
!  here we compute the Debye temperature approximatively from the
!  poisson ratio and the bulk modulus
!
   poisson=(macro_el(4)+macro_el(8) ) * 0.5_DP
   bulkm=(macro_el(1)+macro_el(5) ) * 0.5_DP
   CALL compute_debye_temperature_poisson(poisson, bulkm, &
                            density, nat, omega, approx_debye_t)
!
!  compute the Debye temperature and the thermodynamic quantities
!  within the Debye model
!
   CALL compute_debye_temperature(el_con, density, nat, omega, debye_t)
   CALL write_thermo_debye(igeom)
   CALL plot_thermo_debye(igeom)
!
!  save elastic constants and compliances on file
!
   filelastic='elastic_constants/'//TRIM(fl_el_cons)//'.g'//&
                                                     TRIM(int_to_char(igeom))
   IF (my_image_id==root_image) CALL write_elastic(filelastic)

   IF (ANY(stype)) THEN
      WRITE(stdout,'(/,5x,"Calculated internal relaxations (a.u.) &
                                           &for geometry",i5)') igeom 
      DO istep=1, nstep_ec
         IF (stype(istep)) THEN
            WRITE(stdout,'(5x,"Strain step number ",i5)') istep
            WRITE(stdout,'(5x,"Atom that moves",i5)') move_at(istep)
            WRITE(stdout,'(5x,"Movement direction",3f12.5)') atom_dir(:,istep)
            IF (lcm_ec) WRITE(stdout,'(5x,"Conserving the center of mass")')
            DO igeo=1,ngeo_strain
               WRITE(stdout,'(2f20.8)') epsil_y(igeo,istep,igeom), &
                                          min_y(1,igeo,istep,igeom)  
            ENDDO
         ENDIF
      ENDDO
   ENDIF
ENDDO

DEALLOCATE(sigma_geo_aux)
DEALLOCATE(energy_geo_eff)
DEALLOCATE(epsilon_geo_eff)

RETURN
END SUBROUTINE manage_elastic_cons
