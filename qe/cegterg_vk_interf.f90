!
! Copyright (C) 2023- Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
INTERFACE cegterg_interf
!-----------------------------------------------------------------------
  ATTRIBUTES(GLOBAL) SUBROUTINE cegterg_init (outk_d, nbasek_d, kdimk_d, &
                     stx_d, hpsi, psi, hc, sc, kdmx, nvecx, nk)
!-----------------------------------------------------------------------
    USE cudafor
    USE util_param, ONLY : DP
    IMPLICIT NONE

    INTEGER, VALUE :: kdmx, nvecx, nk
    LOGICAL, DEVICE :: outk_d(nk)
    INTEGER, DEVICE :: nbasek_d(nk)
    INTEGER, DEVICE :: kdimk_d(nk)
    INTEGER, DEVICE :: stx_d(nk)
    COMPLEX(DP), DEVICE :: hpsi(kdmx,nvecx*nk)
    COMPLEX(DP), DEVICE :: psi(kdmx,nvecx*nk) 
    COMPLEX(DP), DEVICE :: hc(nvecx,nvecx,nk)
    COMPLEX(DP), DEVICE :: sc(nvecx,nvecx,nk)

   END SUBROUTINE cegterg_init

  !-----------------------------------------------------------------------
   ATTRIBUTES(GLOBAL) SUBROUTINE cegterg_init_us (outk_d, nbasek_d, kdimk_d, &
               stx_d, hpsi, spsi, psi, hc, sc, kdmx, nvecx, nk)
  !-----------------------------------------------------------------------
   USE cudafor
   USE util_param, ONLY : DP
   IMPLICIT NONE

   INTEGER, VALUE :: kdmx, nvecx, nk
   LOGICAL, DEVICE :: outk_d(nk)
   INTEGER, DEVICE :: nbasek_d(nk)
   INTEGER, DEVICE :: kdimk_d(nk)
   INTEGER, DEVICE :: stx_d(nk)
   COMPLEX(DP), DEVICE :: hpsi(kdmx,nvecx*nk)
   COMPLEX(DP), DEVICE :: psi(kdmx,nvecx*nk)
   COMPLEX(DP), DEVICE :: spsi(kdmx,nvecx*nk)
   COMPLEX(DP), DEVICE :: hc(nvecx,nvecx,nk)
   COMPLEX(DP), DEVICE :: sc(nvecx,nvecx,nk)
   END SUBROUTINE cegterg_init_us
  !-----------------------------------------------------------------------
   ATTRIBUTES(GLOBAL) SUBROUTINE cegterg_overlap(outk_d, nbasek_d, kdimk_d, &
         notcnvk_d, nb1k_d, stx_d, hpsi, psi, hc, sc, kdmx, nvecx, nk)
  !-----------------------------------------------------------------------
   USE cudafor
   USE util_param, ONLY : DP
   IMPLICIT NONE

   INTEGER, VALUE :: kdmx, nvecx, nk
   LOGICAL, DEVICE :: outk_d(nk)
   INTEGER, DEVICE :: nbasek_d(nk)
   INTEGER, DEVICE :: kdimk_d(nk)
   INTEGER, DEVICE :: notcnvk_d(nk)
   INTEGER, DEVICE :: nb1k_d(nk)
   INTEGER, DEVICE :: stx_d(nk)
   COMPLEX(DP), DEVICE :: hpsi(kdmx,nvecx*nk)
   COMPLEX(DP), DEVICE :: psi(kdmx,nvecx*nk)
   COMPLEX(DP), DEVICE :: hc(nvecx,nvecx,nk)
   COMPLEX(DP), DEVICE :: sc(nvecx,nvecx,nk)
   END SUBROUTINE cegterg_overlap
  !
  !-----------------------------------------------------------------------
   ATTRIBUTES(GLOBAL) SUBROUTINE cegterg_overlap_us(outk_d, nbasek_d, &
                  kdimk_d, notcnvk_d, nb1k_d, stx_d, hpsi, spsi, psi, &
                               hc, sc, kdmx, nvecx, nk)
  !-----------------------------------------------------------------------
   USE cudafor
   USE util_param, ONLY : DP
   IMPLICIT NONE

   INTEGER, VALUE :: kdmx, nvecx, nk
   LOGICAL, DEVICE :: outk_d(nk)
   INTEGER, DEVICE :: nbasek_d(nk)
   INTEGER, DEVICE :: kdimk_d(nk)
   INTEGER, DEVICE :: notcnvk_d(nk)
   INTEGER, DEVICE :: nb1k_d(nk)
   INTEGER, DEVICE :: stx_d(nk)
   COMPLEX(DP), DEVICE :: hpsi(kdmx,nvecx*nk)
   COMPLEX(DP), DEVICE :: psi(kdmx,nvecx*nk)
   COMPLEX(DP), DEVICE :: spsi(kdmx,nvecx*nk)
   COMPLEX(DP), DEVICE :: hc(nvecx,nvecx,nk)
   COMPLEX(DP), DEVICE :: sc(nvecx,nvecx,nk)
   END SUBROUTINE cegterg_overlap_us

  !-----------------------------------------------------------------------
   ATTRIBUTES(GLOBAL) SUBROUTINE cegterg_upd0(outk_d, nbasek_d, st_d, aux_d, &
                       conv_d, nb1k_d, vc, ew, e, nvecx, nvec, kter, nk)
  !-----------------------------------------------------------------------
    USE util_param, ONLY : DP
    IMPLICIT NONE
    INTEGER, VALUE :: nk, nvecx, nvec, kter
    LOGICAL, DEVICE :: outk_d(nk)
    INTEGER, DEVICE :: nbasek_d(nk)
    INTEGER, DEVICE :: st_d(nk)
    INTEGER, DEVICE :: aux_d(nk)
    LOGICAL, DEVICE :: conv_d(nvec*nk)
    INTEGER, DEVICE :: nb1k_d(nk)
    COMPLEX(DP), DEVICE :: vc(nvecx,nvecx,nk)
    REAL(DP), DEVICE :: ew(nvecx,nk), e(nvec)
    END SUBROUTINE cegterg_upd0

!-------------------------------------------------------------------------
  ATTRIBUTES(GLOBAL) SUBROUTINE compute_dot_ew (outk_d, npw, nbasek_d, &
             notcnvk_d, stx_d, psi, ew, npol, npwx, nvecx, nk)
!-------------------------------------------------------------------------
USE cudafor
USE util_param, ONLY : DP
IMPLICIT NONE

INTEGER, VALUE :: npwx, nvecx, npol, nk
LOGICAL, DEVICE :: outk_d(nk)
INTEGER, DEVICE :: nbasek_d(nk)
INTEGER, DEVICE :: npw(nk)
INTEGER, DEVICE :: notcnvk_d(nk)
INTEGER, DEVICE :: stx_d(nk)
COMPLEX(DP), DEVICE :: psi(npwx*npol,nvecx*nk)
REAL(DP), DEVICE :: ew(nvecx,nk)

END SUBROUTINE compute_dot_ew

!-------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE cegterg_upd3(outk_d, nbasek_d, kdimk_d, &
        notcnvk_d, nb1k_d, stx_d, hpsi, psi, vc, kdmx, nvecx, nk)
!-------------------------------------------------------------------------
USE cudafor
USE util_param, ONLY : DP
IMPLICIT NONE

INTEGER, VALUE :: kdmx, nvecx, nk
LOGICAL, DEVICE :: outk_d(nk)
INTEGER, DEVICE :: nbasek_d(nk)
INTEGER, DEVICE :: kdimk_d(nk)
INTEGER, DEVICE :: notcnvk_d(nk)
INTEGER, DEVICE :: nb1k_d(nk)
INTEGER, DEVICE :: stx_d(nk)
COMPLEX(DP), DEVICE :: hpsi(kdmx,nvecx*nk) 
COMPLEX(DP), DEVICE :: psi(kdmx,nvecx*nk)
COMPLEX(DP), DEVICE :: vc(nvecx,nvecx,nk)
END SUBROUTINE cegterg_upd3

!-------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE cegterg_herm(outk_d, nbasek_d, nb1k_d, &
                              hc, sc, nvecx, nk)
!-------------------------------------------------------------------------
USE cudafor
USE util_param, ONLY : DP
IMPLICIT NONE

INTEGER, VALUE :: nvecx, nk
LOGICAL, DEVICE :: outk_d(nk)
INTEGER, DEVICE :: nbasek_d(nk)
INTEGER, DEVICE :: nb1k_d(nk)
COMPLEX(DP), DEVICE :: hc(nvecx,nvecx,nk)
COMPLEX(DP), DEVICE :: sc(nvecx,nvecx,nk)
END SUBROUTINE cegterg_herm

  !-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE copy_psi_gpu( lda, outk, enter, kdim, stx, st, &
                              npol, psi_d, evc, nvec, nvecx, nk )
  !-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  INTEGER, INTENT(IN), VALUE :: nk, nvec, nvecx
  LOGICAL, INTENT(IN), DEVICE :: outk(nk), enter(nk)
  INTEGER, INTENT(IN), DEVICE :: kdim(nk)
  INTEGER, INTENT(IN), DEVICE :: stx(nk), st(nk)
  INTEGER, INTENT(IN), VALUE :: npol
  COMPLEX(DP), DEVICE :: psi_d(lda*npol, nvecx * nk)
  COMPLEX(DP), DEVICE :: evc(lda*npol, nvec * nk)
END SUBROUTINE copy_psi_gpu

END INTERFACE cegterg_interf
#endif
