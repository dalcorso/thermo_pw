!
! Copyright (C) 2020 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------
SUBROUTINE collect_the_info(info)
!-------------------------------------------------------------
!
!  This routine must be called after the call to check_initial_status
!  to save on the info structure the status information on the phonon
!  grid calculation. It acts as an interface between the phonon
!  and the infor structure
!
USE collect_info,     ONLY : collect_info_type, init_collect_info, &
                             save_collect_info, comm_collect_info
USE mp_images,        ONLY : nimage, my_image_id, inter_image_comm
USE ions_base,        ONLY : nat
USE disp,             ONLY : nqs, comp_iq, done_iq
USE grid_irr_iq,      ONLY : comp_irr_iq, done_irr_iq, irr_iq
USE control_qe,       ONLY : use_ph_images

IMPLICIT NONE
TYPE(collect_info_type), INTENT(INOUT) :: info
INTEGER :: pos, nima

IF (use_ph_images) THEN
   pos=my_image_id+1
   nima=nimage
ELSE
   pos=1
   nima=1
ENDIF
CALL init_collect_info(info, nqs, nat, nima, irr_iq)
CALL save_collect_info(info, nqs, nat, pos, comp_irr_iq, done_irr_iq, &
                                                         comp_iq, done_iq)
!
!  If the phonon images are used, all images must have the same information 
!  on what must be calculated by each image, so the collect_info structure 
!  is shared between all images for each geometry. 
!  Otherwise here there is a syncronization of all processors.
!
IF (use_ph_images) CALL comm_collect_info(info, inter_image_comm)

RETURN
END SUBROUTINE collect_the_info
