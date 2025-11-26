!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE manage_piezo_tensor(nwork,ngeom)
!----------------------------------------------------------------------------
!
USE initial_conf,         ONLY : ibrav_save
USE thermo_mod,           ONLY : energy_geo
USE thermo_sym,           ONLY : code_group_save
USE control_elastic_constants, ONLY : ngeo_strain, frozen_ions, epsil_geo, &
                                 el_con_ibrav_geo, el_con_celldm_geo,      &
                                 el_con_at_geo, el_con_celldm_geo,         &
                                 el_con_omega_geo, epsil_geo,              &
                                 start_geometry_qha, last_geometry_qha

USE control_piezoelectric_tensor, ONLY : g_piezo_tensor_geo, &
                                 eg_piezo_tensor_geo, e_piezo_tensor_geo, &
                                 d_piezo_tensor_geo, polar0_geo
USE piezoelectric_tensor, ONLY : compute_improper_piezo_tensor, &
                                 compute_d_piezo_tensor,        &
                                 polar_strain, print_d_piezo_tensor,         &
                                 print_g_piezo_tensor, print_e_piezo_tensor, &
                                 e_piezo_tensor, tot_b_phase,     &
                                 eg_piezo_tensor, g_piezo_tensor, &
                                 d_piezo_tensor,                  &
                                 compute_proper_piezo_tensor,     &
                                 compute_polarization_equil,      &
                                 proper_improper_piezo,           &
                                 print_eg_piezo_tensor, clean_piezo_tensor, &  
                                 write_piezo_tensor
                                    
USE elastic_constants,    ONLY : epsilon_geo, el_con, el_compliances, &
                                 read_elastic
USE data_files,           ONLY : fl_el_cons, fl_piezo

USE mp_world,             ONLY : world_comm
USE mp_images,            ONLY : my_image_id, root_image, nproc_image
USE mp,                   ONLY : mp_bcast, mp_sum
USE io_global,            ONLY : stdout, meta_ionode_id 

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork, ngeom
LOGICAL :: exst

INTEGER :: iwork, igeom, base_ind, work_base
CHARACTER(LEN=256) :: filepiezo, filelastic
CHARACTER(LEN=6)   :: int_to_char
LOGICAL :: lreturn
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
work_base = nwork / ngeom
DO igeom=start_geometry_qha,last_geometry_qha
   base_ind=(igeom-1)*work_base
   DO iwork=1,work_base
      lreturn=lreturn.OR.(ABS(energy_geo(base_ind+iwork))<1.D-10)
   ENDDO
ENDDO
IF (lreturn) RETURN

!
!  First collect the polarization among all images
!
CALL mp_sum(polar_strain, world_comm)
polar_strain=polar_strain / nproc_image
CALL mp_sum(tot_b_phase, world_comm)
tot_b_phase=tot_b_phase / nproc_image

DO igeom=start_geometry_qha, last_geometry_qha

   WRITE(stdout,'(2x,76("*"),/)')
   WRITE(stdout,'(5x,"Computing the piezoelectric tensor for the equilibrium &
                                         &geometry=",i4)') igeom
   WRITE(stdout,'(5x,i3,6f10.5)') el_con_ibrav_geo(igeom), &
                                   el_con_celldm_geo(:,igeom)
   WRITE(stdout,'(2x,76("*"),/)')
   base_ind= (igeom-1)*work_base
   !
   !  First compute the polarization of the unstrained state
   !
   WRITE(stdout,'(/,2x,76("-"))')
   WRITE(stdout,'(5x,"Polarization of the equilibrium geometry")')
   WRITE(stdout,'(2x,76("-"),/)')
   CALL compute_polarization_equil(polar_strain(:,base_ind+1), &
          epsil_geo(base_ind+1), polar0_geo(:,igeom), work_base, ngeo_strain)
   WRITE(stdout,'(/,2x,76("-"))')
   WRITE(stdout,'(5x,"Computing the improper piezoelectric tensor")')
   WRITE(stdout,'(2x,76("-"),/)')
!
!  the piezoelectric tensor is calculated here. First the improper one
!
   CALL compute_improper_piezo_tensor(polar_strain(:,base_ind+1), &
                epsilon_geo(:,:,base_ind+1), work_base, ngeo_strain, &
                ibrav_save, code_group_save)
   g_piezo_tensor_geo(:,:,igeom)=g_piezo_tensor(:,:)
   CALL print_g_piezo_tensor(frozen_ions)
   CALL proper_improper_piezo(polar0_geo(:,igeom), g_piezo_tensor, &
                                                      eg_piezo_tensor, -1)
   CALL clean_piezo_tensor(eg_piezo_tensor, ibrav_save, code_group_save)
   eg_piezo_tensor_geo(:,:,igeom)=eg_piezo_tensor(:,:)
   CALL print_eg_piezo_tensor(frozen_ions)

   WRITE(stdout,'(/,2x,76("-"))')
   WRITE(stdout,'(5x,"Computing the proper piezoelectric tensor")')
   WRITE(stdout,'(2x,76("-"),/)')
!
!  and then the proper one (this can be compared with experiments)
!
   CALL compute_proper_piezo_tensor(tot_b_phase(:,base_ind+1), &
                      epsilon_geo(:,:,base_ind+1), work_base,  &
           ngeo_strain, ibrav_save, code_group_save,           &
           el_con_at_geo(:,:,igeom) )

   e_piezo_tensor=e_piezo_tensor * el_con_celldm_geo(1,igeom) / &
                                           el_con_omega_geo(igeom)
   e_piezo_tensor_geo(:,:,igeom)=e_piezo_tensor(:,:)
   CALL print_e_piezo_tensor(frozen_ions)
!
!  If a file with the elastic constants exists read them and compute the
!  strain piezoelectric tensor.
!
   filelastic='elastic_constants/'//TRIM(fl_el_cons)//'.g'//&
                                                     TRIM(int_to_char(igeom))
   IF (my_image_id==root_image) CALL read_elastic(filelastic, exst)
   CALL mp_bcast(exst, meta_ionode_id, world_comm)
   IF (exst) THEN
      CALL mp_bcast(el_con, meta_ionode_id, world_comm)
      CALL mp_bcast(el_compliances, meta_ionode_id, world_comm)
      CALL compute_d_piezo_tensor(el_compliances)
      d_piezo_tensor_geo(:,:,igeom)=d_piezo_tensor(:,:)
      CALL print_d_piezo_tensor(frozen_ions)
   ENDIF
!
!   Now write all piezo tensors on file
!
   filepiezo='elastic_constants/'//TRIM(fl_piezo)//'.g'//&
                                                     TRIM(int_to_char(igeom))
   IF (my_image_id==root_image) CALL write_piezo_tensor(filepiezo,&
                                              polar0_geo(:,igeom))
ENDDO
RETURN
END SUBROUTINE manage_piezo_tensor
