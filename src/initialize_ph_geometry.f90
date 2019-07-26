!
! Copyright (C) 2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE initialize_ph_geometry(igeom, auxdyn)

USE input_parameters, ONLY : outdir
USE control_thermo,   ONLY : outdir_thermo
USE io_global,        ONLY : stdout

USE control_ph,       ONLY : rec_code_read, low_directory_check
USE grid_irr_iq,      ONLY : done_irr_iq
USE disp,             ONLY : done_iq

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom
CHARACTER (LEN=256), INTENT(INOUT) :: auxdyn
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: ldcs

WRITE(stdout,'(/,5x,40("%"))')
WRITE(stdout,'(5x,"Computing geometry ", i5)') igeom
WRITE(stdout,'(5x,40("%"),/)')
outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(igeom))//'/'
!
! ... reads the xml file produced by pw.x for this geometry
!
CALL initialize_geometry_and_ph(.TRUE., igeom, auxdyn)
!
!  no need to recheck the phsave directory here
!
ldcs=low_directory_check
low_directory_check=.TRUE.
CALL check_initial_status_tpw(auxdyn)
low_directory_check=ldcs
!
!  The recover is managed after this routine
!
done_irr_iq=.FALSE.
done_iq=.FALSE.
!
!  When there are many geometries the possible recover file is not used.
!  Actually the rec_code_read read from file is the one written in the
!  previous task of this image. We do all the tasks of a given image if
!  the dynamical matrix is not already on file. In this case thermo_pw
!  do not call the phonon and rec_code_read is not relevant.
!
rec_code_read=-40
!
RETURN
END SUBROUTINE initialize_ph_geometry

SUBROUTINE close_ph_geometry(close_ph)
!
!  This routine frees the variables allocated by the phq_readin, closes
!  all the files, and reset the control variables of the phonon to their
!  default values.
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
