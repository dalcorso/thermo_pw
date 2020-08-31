!
! Copyright (C) 2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE prepare_do_phonon(igeom)
!
!   This routine initializes the phonon calculation by calling 
!   fast_phq_readin that reads the xml file produced by pw.x and 
!   allocates the variables usually allocated by phq_readin.
!   Then it calls check_initial_status to initialize the phonon.
!
USE io_global,        ONLY : stdout

USE control_phrun,    ONLY : auxdyn

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom
CHARACTER(LEN=80) :: message
LOGICAL :: ldcs

WRITE(message,'(5x,"Computing geometry ", i5)') igeom
CALL decorated_write(message)
CALL set_outdir_name(igeom)
!
! ... reads the xml file produced by pw.x for this geometry
!
CALL fast_phq_readin(.TRUE., igeom)
!
CALL set_fildyn_name(igeom)
!
!  no need to recheck the phsave directory here. Call check_initial_status_tpw
!  with iflag=1 to skip the recheck.
!
CALL check_initial_status_tpw(auxdyn, 1)
!
RETURN
END SUBROUTINE prepare_do_phonon

SUBROUTINE close_ph_geometry(close_ph)
!
!  This routine frees the variables allocated by phq_readin, closes
!  all the files, and reset the control variables of the phonon to their
!  input values. Note that it frees also the pw.x variables allocated
!  by read_file.
!
USE ph_restart,       ONLY : destroy_status_run
USE save_ph,          ONLY : clean_input_variables
USE control_ph,       ONLY : done_epsil, done_zeu, done_zue, done_start_zstar, &
                             rec_code, rec_code_read, epsil, zeu, zue, &
                             start_q, last_q
USE initial_conf,     ONLY : epsil_save, zeu_save, zue_save, start_q_save, &
                             last_q_save
USE ramanm,           ONLY : done_lraman, done_elop

IMPLICIT NONE

LOGICAL, INTENT(IN) :: close_ph

CALL clean_pw(.TRUE.)
IF (close_ph) CALL close_phq(.FALSE.)
CALL clean_input_variables()
CALL destroy_status_run()
CALL deallocate_ph_tpw()
CALL deallocate_part()

done_epsil=.FALSE.
done_zue=.FALSE.
done_zeu=.FALSE.
done_start_zstar=.FALSE.
done_lraman=.FALSE.
done_elop=.FALSE.
epsil=epsil_save
zeu=zeu_save
zue=zue_save
start_q=start_q_save
last_q=last_q_save
rec_code_read=-1000
rec_code=-1000

RETURN
END SUBROUTINE close_ph_geometry
