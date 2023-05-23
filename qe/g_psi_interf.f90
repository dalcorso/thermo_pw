!
! Copyright (C) 2023- Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
INTERFACE g_psii_interf
!-----------------------------------------------------------------------
    ATTRIBUTES(GLOBAL) SUBROUTINE g_psii( lda, outk, npw, nveck, nb1k, &
                                stx, npol, psi, e, nvecx, nset )
  !-----------------------------------------------------------------------
  !
   USE cudafor
   USE util_param, ONLY : DP
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN), VALUE :: lda
   INTEGER, INTENT(IN), VALUE :: nset, nvecx
   LOGICAL, INTENT(IN), DEVICE :: outk(nset)
   INTEGER, INTENT(IN), DEVICE :: npw(nset)
   INTEGER, INTENT(IN), DEVICE :: nveck(nset)
   INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
   INTEGER, INTENT(IN), DEVICE :: stx(nset)
   INTEGER, INTENT(IN), VALUE :: npol
   COMPLEX(DP), INTENT(INOUT), DEVICE :: psi(lda, npol, nvecx * nset)
   REAL(DP), INTENT(IN), DEVICE :: e(nvecx,nset)

   END SUBROUTINE g_psii

END INTERFACE g_psii_interf
#endif
