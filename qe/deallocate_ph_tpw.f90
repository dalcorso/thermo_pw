!
! Copyright (C) 2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE allocate_ph_tpw
!------------------------------------------------------------------------
!
!  This allocation is done by phq_readin_tpw. To avoid crash when only
!  the xml file is read we allocate these variables here
!
USE kinds, ONLY : DP
USE optical, ONLY : fru
USE images_omega, ONLY : comp_f
USE freq_ph, ONLY : fiu
IMPLICIT NONE
    IF (.NOT. ALLOCATED(fru)) ALLOCATE(fru(1))
    IF (.NOT. ALLOCATED(fiu)) ALLOCATE(fiu(1))
    IF (.NOT. ALLOCATED(comp_f)) ALLOCATE(comp_f(1))
    fiu=0.0_DP
    fru=0.0_DP
    comp_f=.TRUE.
RETURN
END SUBROUTINE allocate_ph_tpw

!------------------------------------------------------------------------
SUBROUTINE deallocate_ph_tpw
!------------------------------------------------------------------------

USE optical, ONLY : fru
USE images_omega, ONLY : comp_f
USE many_k_ph_mod, ONLY : becp1k_d, alphak_d, becptk_d, alphatk_d
IMPLICIT NONE

IF (ALLOCATED(fru)) DEALLOCATE(fru)
IF (ALLOCATED(comp_f)) DEALLOCATE(comp_f)
IF (ALLOCATED(becp1k_d)) DEALLOCATE(becp1k_d)
IF (ALLOCATED(alphak_d)) DEALLOCATE(alphak_d)
IF (ALLOCATED(becptk_d)) DEALLOCATE(becptk_d)
IF (ALLOCATED(alphatk_d)) DEALLOCATE(alphatk_d)

RETURN
END SUBROUTINE deallocate_ph_tpw
