!
! Copyright (C) 2013-2017 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE collect_everything(auxdyn, igeom)
!
!  This routine calls ph.x another time in order to collect the results 
!  calculated by all images. The work is divided among images. Each image 
!  collects a set of q points. Which one is controlled by the routine 
!  check_stc when all_geometries_together is .TRUE.. Otherwise
!  the q points are distributed sequentially to each image.
!
USE mp_images,   ONLY : my_image_id, nimage
USE mp_world,    ONLY : world_comm
USE mp,          ONLY : mp_barrier
USE io_global,   ONLY : stdout
USE control_ph,  ONLY : recover
USE control_qe,  ONLY : tcollect_all
USE control_thermo, ONLY : all_geometries_together
USE ions_base,   ONLY : nat
USE disp,        ONLY : comp_iq, nqs
USE grid_irr_iq, ONLY : comp_irr_iq, irr_iq

IMPLICIT NONE
CHARACTER (LEN=256), INTENT(IN) :: auxdyn
INTEGER :: irr, iq, igeom
LOGICAL :: save_recover, std
!
WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(5x,"Collecting the results and writing fildyn files")') 
WRITE(stdout,'(2x,76("+"),/)')
!
!  All quantities should be already in phsave
!  Each image receives a subset of q points and collects the results
!  for those q points.
!
save_recover=recover
recover=.TRUE.
tcollect_all=.TRUE.

IF (all_geometries_together) THEN
   CALL check_stc(nqs, igeom, nat, nimage, my_image_id, comp_iq, &
                                                        comp_irr_iq, std)
ELSE
   comp_irr_iq=.TRUE.
   comp_iq=.FALSE.
   std=.FALSE.
   DO iq=1,nqs
      IF (MOD(iq-1,nimage)==my_image_id) THEN
         comp_iq(iq)=.TRUE.
         DO irr=0, irr_iq(iq)
            comp_irr_iq(:,iq)=.TRUE.
         ENDDO
         std=.TRUE.
      ENDIF
   ENDDO
ENDIF

IF (std) CALL do_phonon_tpw(auxdyn)

tcollect_all=.FALSE.
recover=save_recover

RETURN
END SUBROUTINE collect_everything

SUBROUTINE check_stc(nqs, igeom, nat, nimage, my_image_id, comp_iq, &
                                                     comp_irr_iq, std)
!
!  This subroutine checks if this image has something to collect
!  in the current geometry and sets the appropriate comp_iq and comp_irr_iq
!
USE thermo_mod,   ONLY : tot_ngeo, start_geometry, last_geometry
USE initial_conf, ONLY : collect_info_save
IMPLICIT NONE
INTEGER, INTENT(IN) :: nqs, nat, igeom, nimage, my_image_id
LOGICAL, INTENT(INOUT) :: comp_iq(nqs), comp_irr_iq(0:3*nat,nqs)
LOGICAL, INTENT(OUT) :: std

INTEGER :: irr, iq, jgeom, task
INTEGER :: start_task(nimage), last_task(nimage)

CALL compute_total_task_and_divide(start_task, last_task, nimage)

comp_irr_iq=.TRUE.
comp_iq=.FALSE.
std=.FALSE.
task=0
DO jgeom=start_geometry, last_geometry
   DO iq=1,collect_info_save(jgeom)%nqs
      task = task+1
      IF (task >=start_task(my_image_id+1).AND.task<=last_task(my_image_id+1)&
                                                     .AND.igeom==jgeom) THEN
         comp_iq(iq)=.TRUE.
         DO irr=0, collect_info_save(igeom)%irr_iq(iq)
            comp_irr_iq(irr,iq)=.TRUE.
         ENDDO
         std=.TRUE.
      ENDIF
   ENDDO
ENDDO

RETURN
END SUBROUTINE check_stc 

SUBROUTINE check_stc_g(igeom, nimage, my_image_id, std)
!
!  This subroutine checks if this image has something to collect
!  in the current geometry, but does not set the comp_irr_iq and comp_iq 
!  flags
!
USE thermo_mod,   ONLY : tot_ngeo, start_geometry, last_geometry
USE initial_conf, ONLY : collect_info_save
IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom, nimage, my_image_id
LOGICAL, INTENT(OUT) :: std

INTEGER :: iq, jgeom, task
INTEGER :: start_task(nimage), last_task(nimage)

CALL compute_total_task_and_divide(start_task, last_task, nimage)

std=.FALSE.
task=0
DO jgeom=start_geometry, last_geometry
   DO iq=1,collect_info_save(jgeom)%nqs
      task = task+1
      IF (task >=start_task(my_image_id+1).AND.task<=last_task(my_image_id+1)&
                                             .AND.igeom==jgeom) std=.TRUE.
   ENDDO
ENDDO

RETURN
END SUBROUTINE check_stc_g

SUBROUTINE compute_total_task_and_divide(start_task, last_task, nimage)
!
!  This routine distributes all the collection tasks to the available images. 
!  When there are more tasks than images, it gives the q points belonging
!  to the same geometry to the same image as much as possible
!  in order to minimize the number of geometries that are initialized 
!  again by each image. 
!
USE thermo_mod,   ONLY : tot_ngeo, start_geometry, last_geometry
USE initial_conf, ONLY : collect_info_save
IMPLICIT NONE

INTEGER, INTENT(IN) :: nimage
INTEGER, INTENT(INOUT) :: start_task(nimage), last_task(nimage)

INTEGER :: iq, jgeom, image, total_task, task_per_image, resto, task

total_task=0
DO jgeom=start_geometry, last_geometry
   DO iq=1,collect_info_save(jgeom)%nqs
      total_task = total_task+1
   ENDDO
ENDDO

task_per_image=total_task/nimage
resto=total_task - task_per_image * nimage

start_task(1)=1
last_task(1)=task_per_image
IF (resto >=1) last_task(1)=task_per_image+1
DO image=2, nimage
   start_task(image)=last_task(image-1)+1
   last_task(image)=start_task(image)+task_per_image-1
   IF (resto >=1) last_task(image)=last_task(image)+1
ENDDO

RETURN
END SUBROUTINE compute_total_task_and_divide
