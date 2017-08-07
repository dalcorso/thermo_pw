!
! Copyright (C) 2013 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE collect_everything(auxdyn)
!
!  This routine calls phonon a last time in order to collect the results
!  of all the images. The work is divided among images. Each image collects
!  a set of q points if nimage<nqs, or the first nqs images collect something
!  if nimage>nqs.
!
USE mp_images,   ONLY : my_image_id, nimage
USE mp_world,    ONLY : world_comm
USE mp,          ONLY : mp_barrier
USE io_global,   ONLY : stdout
USE control_ph,  ONLY : recover, tmp_dir_ph
USE control_qe,  ONLY : tcollect_all
USE disp,        ONLY : comp_iq, nqs
USE grid_irr_iq, ONLY : comp_irr_iq

IMPLICIT NONE
CHARACTER (LEN=256), INTENT(IN) :: auxdyn
INTEGER :: iq
!
! We cannot start the collection before all images finished
!
CALL mp_barrier(world_comm)
!
WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(5x,"Collecting the results and writing fildyn files")') 
WRITE(stdout,'(2x,76("+"),/)')
!
!  All quantities should be already in phsave
!  Each image receives a subset of q points and collects the results
!
comp_irr_iq=.TRUE.
comp_iq=.FALSE.
recover=.TRUE.
tcollect_all=.TRUE.

DO iq=1,nqs
   IF (MOD(iq,nimage)==my_image_id) comp_iq(iq)=.TRUE.
ENDDO

CALL do_phonon_tpw(auxdyn)

tcollect_all=.FALSE.
!
! resyncronize all images
!
CALL mp_barrier(world_comm)

RETURN
END SUBROUTINE collect_everything
