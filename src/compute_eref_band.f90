!
! Copyright (C) 2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
SUBROUTINE compute_eref_band(e, nbnd, eref, print_eref)
!--------------------------------------------------------------------
USE kinds,            ONLY : DP
USE constants,        ONLY : rytoev
USE ener,             ONLY : ef
USE klist,            ONLY : degauss, nelec
USE noncollin_module, ONLY : noncolin
IMPLICIT NONE
INTEGER, INTENT(IN) :: nbnd
REAL(DP), INTENT(IN)  :: e(nbnd)
REAL(DP), INTENT(OUT) :: eref
LOGICAL, INTENT(OUT) :: print_eref

INTEGER :: nband_occ, ibnd

print_eref=.FALSE.
eref=-1d20
IF (degauss > 0.0_DP) THEN
   eref=ef * rytoev
   print_eref=.TRUE.
ELSE
   nband_occ=NINT(nelec/2.0_DP)
   IF (noncolin) nband_occ=nband_occ*2        
   DO ibnd=1, nband_occ
      IF (e(ibnd) > eref) eref=e(ibnd)
   END DO
END IF
RETURN
END SUBROUTINE compute_eref_band

