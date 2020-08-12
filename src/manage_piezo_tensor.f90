!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE manage_piezo_tensor(nwork)

USE initial_conf,         ONLY : ibrav_save
USE thermo_sym,           ONLY : code_group_save
USE control_elastic_constants, ONLY : ngeo_strain, frozen_ions
USE piezoelectric_tensor, ONLY : compute_piezo_tensor, compute_d_piezo_tensor,&
                                 polar_geo, print_d_piezo_tensor,             &
                                 print_g_piezo_tensor
USE elastic_constants,    ONLY : epsilon_geo, el_con, el_compliances, &
                                 read_elastic
USE data_files,           ONLY : fl_el_cons

USE mp_world,             ONLY : world_comm
USE mp_images,            ONLY : my_image_id, root_image, nproc_image
USE mp,                   ONLY : mp_bcast, mp_sum
USE io_global,            ONLY : meta_ionode_id 

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork
LOGICAL :: exst

!
!  First collect the polarization among all images
!
CALL mp_sum(polar_geo, world_comm)
polar_geo=polar_geo / nproc_image
!
!  the piezoelectric tensor is calculated here
!
CALL compute_piezo_tensor(polar_geo, epsilon_geo, nwork, &
                               ngeo_strain, ibrav_save, code_group_save)
CALL print_g_piezo_tensor(frozen_ions)

IF (my_image_id==root_image) CALL read_elastic(fl_el_cons, exst)
CALL mp_bcast(exst, meta_ionode_id, world_comm)
IF (exst) THEN
   CALL mp_bcast(el_con, meta_ionode_id, world_comm)
   CALL mp_bcast(el_compliances, meta_ionode_id, world_comm)
   CALL compute_d_piezo_tensor(el_compliances)
   CALL print_d_piezo_tensor(frozen_ions)
ENDIF

RETURN
END SUBROUTINE manage_piezo_tensor
