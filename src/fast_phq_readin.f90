!
! Copyright (C) 2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE fast_phq_readin(recover_, igeom)
!--------------------------------------------------------------------------
!
!  This routine substitutes phq_readin. It reads the output of pw.x and
!  initializes the same variables initialized by phq_readin, but does not 
!  read the ph_control input file.
!  It must be called by all processors of an image.
!
USE input_parameters, ONLY : outdir
USE control_ph,       ONLY : tmp_dir_ph, tmp_dir_phq, recover
USE control_lr,       ONLY : rec_code_read
USE ph_restart,       ONLY : ph_readfile
USE save_ph,          ONLY : save_ph_input_variables, tmp_dir_save
USE io_files,         ONLY : tmp_dir, check_tempdir
USE mp_images,        ONLY : my_image_id
USE mp_pools,         ONLY : kunit
USE ions_base,        ONLY : nat

IMPLICIT NONE
LOGICAL, INTENT(IN) :: recover_
INTEGER, INTENT(IN) :: igeom

LOGICAL :: exst, parallelfs
INTEGER :: ierr, ierr1
CHARACTER(LEN=256), EXTERNAL :: trimcheck
CHARACTER(LEN=6) :: int_to_char
!
!  Set the tmp_dir_ph for the phonon depending on the image number
!  Set also tmp_dir where pw.x is supposed to have written the xml file
!
kunit=1
tmp_dir=trimcheck(outdir)
tmp_dir_save=tmp_dir
tmp_dir_ph= TRIM (tmp_dir) // '_ph' // TRIM(int_to_char(my_image_id)) //'/'
CALL check_tempdir ( tmp_dir_ph, exst, parallelfs )
tmp_dir_phq=tmp_dir_ph
!
!   Read the xml file produced by pw.x
!
CALL read_file()
!
!   And allocate the variables usually allocated by phq_readin
!
CALL allocate_part(nat)
CALL allocate_ph_tpw()
CALL save_ph_input_variables()
tmp_dir=tmp_dir_save

rec_code_read=-1000
IF (recover_) THEN
   recover=.TRUE.
   CALL ph_readfile('init', 0, 0, ierr)
   CALL ph_readfile('status_ph', 0, 0, ierr1)
ENDIF

RETURN
END SUBROUTINE fast_phq_readin
