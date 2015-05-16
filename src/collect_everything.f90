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
!  of all the images. Only the root image calls the phonon.
!
USE mp_images,   ONLY : my_image_id, root_image
USE io_global,   ONLY : stdout
USE control_ph,  ONLY : recover
USE control_qe,  ONLY : tcollect_all
USE disp,        ONLY : comp_iq
USE grid_irr_iq, ONLY : comp_irr_iq

IMPLICIT NONE
CHARACTER (LEN=256), INTENT(IN) :: auxdyn
!
!  Only the first image collects the results if there are several images
!
IF (  my_image_id /= root_image ) RETURN
!
WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(5x,"Collecting the results and writing fildyn files")') 
WRITE(stdout,'(2x,76("+"),/)')
!
!  Try to compute everything. All quantities should be already in phsave
!
comp_irr_iq=.TRUE.
comp_iq=.TRUE.
recover=.TRUE.
tcollect_all=.TRUE.

CALL do_phonon_tpw(auxdyn)

tcollect_all=.FALSE.

RETURN
END SUBROUTINE collect_everything
