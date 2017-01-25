!
! Copyright (C) 2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE deallocate_ph_tpw

USE optical, ONLY : fru
USE images_omega, ONLY : comp_f
IMPLICIT NONE
    IF (ALLOCATED(fru)) DEALLOCATE(fru)
    IF (ALLOCATED(comp_f)) DEALLOCATE(comp_f)
RETURN
END SUBROUTINE deallocate_ph_tpw
