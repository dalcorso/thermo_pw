!
! Copyright (C) 2013-2017 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE manage_collection(auxdyn,igeom)

USE control_ph, ONLY : trans
USE control_lr, ONLY : lgamma
USE freq_ph,    ONLY : fpol
USE mp_asyn,    ONLY : with_asyn_images

IMPLICIT NONE
CHARACTER (LEN=256), INTENT(IN) :: auxdyn
INTEGER :: igeom

IF (trans) THEN
   IF (with_asyn_images) CALL collect_everything(auxdyn, igeom)
ELSEIF (fpol) THEN
   IF (lgamma) THEN
      IF (with_asyn_images) CALL collect_all_epsilon()
      CALL plot_epsilon_omega_opt()
   ELSE
      IF (with_asyn_images) CALL collect_all_chi()
      CALL plot_epsilon_omega_q()
   ENDIF
ENDIF

RETURN
END SUBROUTINE manage_collection

SUBROUTINE collect_everything(auxdyn, igeom)
!
!  This routine calls ph.x another time in order to collect the results 
!  calculated by all images. The work is divided among images. Each image 
!  collects a set of q points. Which one is controlled by the routine 
!  check_stc when all_geometries_together is .TRUE.. Otherwise
!  the q points are distributed sequentially to each image.
!
USE io_global,   ONLY : stdout
USE control_ph,  ONLY : recover
USE control_qe,  ONLY : tcollect_all
USE ions_base,   ONLY : nat
USE disp,        ONLY : comp_iq, nqs
USE grid_irr_iq, ONLY : comp_irr_iq, irr_iq

USE distribute_collection, ONLY : me_igeom_iq

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
!
!   set the work that has to be done by the present image
!
comp_irr_iq=.FALSE.
comp_iq=.FALSE.
std=.FALSE.
DO iq=1,nqs
   IF (me_igeom_iq(igeom,iq)) THEN
      comp_iq(iq)=.TRUE.
      DO irr=0, irr_iq(iq)
         comp_irr_iq(irr,iq)=.TRUE.
      ENDDO
      std=.TRUE.
   ENDIF
ENDDO
!
!  and do the collection work
!
IF (std) CALL do_phonon_tpw(auxdyn)

tcollect_all=.FALSE.
recover=save_recover

RETURN
END SUBROUTINE collect_everything

SUBROUTINE compute_total_task_and_divide(start_task, last_task, nimage)
!
!  This routine distributes all the collection tasks to the available images. 
!  When there are more tasks than images, it gives the q points belonging
!  to the same geometry to the same image as much as possible
!  in order to minimize the number of geometries that are initialized 
!  again by each image. 
!
USE thermo_mod,   ONLY : no_ph, start_geometry, last_geometry, phgeo_on_file
USE initial_conf, ONLY : collect_info_save
IMPLICIT NONE

INTEGER, INTENT(IN) :: nimage
INTEGER, INTENT(INOUT) :: start_task(nimage), last_task(nimage)

INTEGER :: iq, igeom, image, total_task, task_per_image, resto

total_task=0
DO igeom=start_geometry, last_geometry
   IF (no_ph(igeom).OR.phgeo_on_file(igeom)) CYCLE
   DO iq=1,collect_info_save(igeom)%nqs
      total_task = total_task+1
   ENDDO
ENDDO

task_per_image=total_task/nimage
resto=total_task - task_per_image * nimage

start_task(1)=1
last_task(1)=task_per_image
IF (resto > 0) last_task(1)=task_per_image+1
DO image=2, nimage
   start_task(image)=last_task(image-1)+1
   last_task(image)=start_task(image)+task_per_image-1
   IF (resto >= image) last_task(image)=last_task(image)+1
ENDDO

RETURN
END SUBROUTINE compute_total_task_and_divide

SUBROUTINE divide_all_collection_work()
!
!   This routine, that is called only when all_geometry_together=.TRUE.
!   sets two arrays (different for each image). 
!   me_igeom : for each geometry is .true. if the present image must
!              do some collection work in this geometry
!   me_igeom_iq : is .true. if the present image must collect this
!                 q point of this geometry
!   
!
USE thermo_mod,     ONLY : no_ph, phgeo_on_file, start_geometry, last_geometry
USE distribute_collection, ONLY : me_igeom, me_igeom_iq
USE initial_conf,   ONLY : collect_info_save
USE mp_images,      ONLY : nimage, my_image_id

IMPLICIT NONE

INTEGER :: nqsx, start_task(nimage), last_task(nimage), task, igeom, iq

nqsx=MAXVAL(collect_info_save(:)%nqs)
ALLOCATE(me_igeom(start_geometry:last_geometry))
ALLOCATE(me_igeom_iq(start_geometry:last_geometry,nqsx))

me_igeom=.FALSE.
me_igeom_iq=.FALSE.

CALL compute_total_task_and_divide(start_task, last_task, nimage)

task=0
DO igeom = start_geometry, last_geometry
   IF (no_ph(igeom).OR.phgeo_on_file(igeom)) CYCLE
   DO iq = 1, collect_info_save(igeom)%nqs
      task = task+1
      IF (task >=start_task(my_image_id+1).AND. &
                           task<=last_task(my_image_id+1)) THEN
         me_igeom(igeom)=.TRUE.
         me_igeom_iq(igeom,iq)=.TRUE.                                    
      ENDIF
   ENDDO
ENDDO

RETURN
END SUBROUTINE divide_all_collection_work

SUBROUTINE divide_collection_work(igeom)
!
!  This routine, similar to the previous one, works for the
!  case all_geometry_together=.FALSE.. It sets only
!  me_igeom_iq
!  
!
USE distribute_collection, ONLY : me_igeom_iq
USE mp_images,    ONLY : nimage, my_image_id

USE disp,         ONLY : nqs
IMPLICIT NONE

INTEGER :: igeom
INTEGER :: iq

ALLOCATE(me_igeom_iq(igeom:igeom,nqs))
DO iq=1,nqs
   me_igeom_iq(igeom,iq)=(MOD(iq-1,nimage)==my_image_id)
ENDDO

RETURN
END SUBROUTINE divide_collection_work

SUBROUTINE clean_collection_work()

USE distribute_collection, ONLY : me_igeom, me_igeom_iq
IMPLICIT NONE

IF (ALLOCATED(me_igeom)) DEALLOCATE(me_igeom)
IF (ALLOCATED(me_igeom_iq)) DEALLOCATE(me_igeom_iq)

RETURN
END SUBROUTINE clean_collection_work
