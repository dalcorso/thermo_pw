!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE manage_piezo_tensor(nwork)
!----------------------------------------------------------------------------
!
USE initial_conf,         ONLY : ibrav_save
USE thermo_sym,           ONLY : code_group_save
USE control_elastic_constants, ONLY : ngeo_strain, frozen_ions, epsil_geo
USE piezoelectric_tensor, ONLY : compute_improper_piezo_tensor, &
                                 compute_d_piezo_tensor,        &
                                 polar_geo, print_d_piezo_tensor,            &
                                 print_g_piezo_tensor, print_e_piezo_tensor, &
                                 e_piezo_tensor, tot_b_phase,     &
                                 eg_piezo_tensor, g_piezo_tensor, &
                                 compute_proper_piezo_tensor,     &
                                 compute_polarization_equil,      &
                                 proper_improper_piezo,           &
                                 print_eg_piezo_tensor, clean_piezo_tensor, &  
                                 write_piezo_tensor
                                    
USE elastic_constants,    ONLY : epsilon_geo, el_con, el_compliances, &
                                 read_elastic
USE equilibrium_conf,     ONLY : at0, celldm0, omega0, polar0
USE data_files,           ONLY : fl_el_cons, fl_piezo

USE mp_world,             ONLY : world_comm
USE mp_images,            ONLY : my_image_id, root_image, nproc_image
USE mp,                   ONLY : mp_bcast, mp_sum
USE io_global,            ONLY : stdout, meta_ionode_id 

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork
LOGICAL :: exst

INTEGER :: igeom
CHARACTER(LEN=256) :: filepiezo, filelastic
CHARACTER(LEN=6)   :: int_to_char
!
!  First collect the polarization among all images
!
igeom=1
CALL mp_sum(polar_geo, world_comm)
polar_geo=polar_geo / nproc_image
CALL mp_sum(tot_b_phase, world_comm)
tot_b_phase=tot_b_phase / nproc_image

WRITE(stdout,'(/,2x,76("-"))')
WRITE(stdout,'(5x,"Computing the polarization of the equilibrium structure")')
WRITE(stdout,'(2x,76("-"),/)')
!
!  First compute the polarization of the unstrained state
!
CALL compute_polarization_equil(polar_geo, epsil_geo, polar0, nwork, ngeo_strain)
WRITE(stdout,'(/,2x,76("-"))')
WRITE(stdout,'(5x,"Computing the improper piezoelectric tensor")')
WRITE(stdout,'(2x,76("-"),/)')
!
!  the piezoelectric tensor is calculated here. First the improper one
!
CALL compute_improper_piezo_tensor(polar_geo, epsilon_geo, nwork, &
                               ngeo_strain, ibrav_save, code_group_save)
CALL print_g_piezo_tensor(frozen_ions)
CALL proper_improper_piezo(polar0, g_piezo_tensor, eg_piezo_tensor, -1)
CALL clean_piezo_tensor(eg_piezo_tensor, ibrav_save, code_group_save)
CALL print_eg_piezo_tensor(frozen_ions)

WRITE(stdout,'(/,2x,76("-"))')
WRITE(stdout,'(5x,"Computing the proper piezoelectric tensor")')
WRITE(stdout,'(2x,76("-"),/)')
!
!  and then the proper one (this can be compared with experiments)
!
CALL compute_proper_piezo_tensor(tot_b_phase, epsilon_geo, nwork, &
           ngeo_strain, ibrav_save, code_group_save, at0 )
e_piezo_tensor=e_piezo_tensor * celldm0(1) / omega0
CALL print_e_piezo_tensor(frozen_ions)
filepiezo='elastic_constants/'//TRIM(fl_piezo)//'.g'//&
                                                     TRIM(int_to_char(igeom))
IF (my_image_id==root_image) CALL write_piezo_tensor(filepiezo,polar0)

filelastic='elastic_constants/'//TRIM(fl_el_cons)//'.g'//&
                                                     TRIM(int_to_char(igeom))
IF (my_image_id==root_image) CALL read_elastic(filelastic, exst)
CALL mp_bcast(exst, meta_ionode_id, world_comm)
IF (exst) THEN
   CALL mp_bcast(el_con, meta_ionode_id, world_comm)
   CALL mp_bcast(el_compliances, meta_ionode_id, world_comm)
   CALL compute_d_piezo_tensor(el_compliances)
   CALL print_d_piezo_tensor(frozen_ions)
ENDIF

RETURN
END SUBROUTINE manage_piezo_tensor
