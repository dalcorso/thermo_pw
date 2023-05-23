!
! Copyright (C) 2023- Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
INTERFACE h_psi_gpu_interf
!
!-----------------------------------------------------------------------
    ATTRIBUTES(GLOBAL) SUBROUTINE h_psi_ke( lda, outk, npw, nveck, nb1k, &
                             stx, ikblk, npol, psi_d, hpsi_d, nvecx, nset )
  !-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP

  INTEGER, INTENT(IN), VALUE :: lda
  INTEGER, INTENT(IN), VALUE :: nset, nvecx
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: npw(nset)
  INTEGER, INTENT(IN), DEVICE :: ikblk(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  INTEGER, INTENT(IN), VALUE :: npol
  COMPLEX(DP), DEVICE :: psi_d(lda*npol, nvecx * nset)
  COMPLEX(DP), DEVICE :: hpsi_d(lda*npol, nvecx * nset)

END SUBROUTINE h_psi_ke
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE h_psi_calbec( lda, outk, npw, nveck, &
                             nb1k, stx, ikblk, npol, psi_d, nvecx, nset )
  !-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  INTEGER, INTENT(IN), VALUE :: nset, nvecx
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: ikblk(nset)
  INTEGER, INTENT(IN), DEVICE :: npw(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  INTEGER, INTENT(IN), VALUE :: npol
  COMPLEX(DP), DEVICE :: psi_d(lda*npol, nvecx * nset)

  END SUBROUTINE h_psi_calbec
  !
    !-----------------------------------------------------------------------
  ATTRIBUTES(GLOBAL) SUBROUTINE copy_s_gpu( lda, outk, kdimk, nveck, &
                        nb1k, stx, npol, psi_d, spsi_d, nvecx, nset )
  !-----------------------------------------------------------------------
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  INTEGER, INTENT(IN), VALUE :: nset, nvecx
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: kdimk(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  INTEGER, INTENT(IN), VALUE :: npol
  COMPLEX(DP), DEVICE :: psi_d(lda*npol, nvecx * nset)
  COMPLEX(DP), DEVICE :: spsi_d(lda*npol, nvecx * nset)
  END SUBROUTINE copy_s_gpu

!-----------------------------------------------------------------------
  ATTRIBUTES(GLOBAL) SUBROUTINE add_vnlpsi_us_gpu( lda, outk, npw, nveck, &
                     nb1k, stx, ikblk, npol, hpsi_d, spsi_d, nvecx, nset )
  !-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  INTEGER, INTENT(IN), VALUE :: nset, nvecx
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: ikblk(nset)
  INTEGER, INTENT(IN), DEVICE :: npw(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  INTEGER, INTENT(IN), VALUE :: npol
  COMPLEX(DP), DEVICE :: hpsi_d(lda*npol, nvecx * nset)
  COMPLEX(DP), DEVICE :: spsi_d(lda*npol, nvecx * nset)
 
  END  SUBROUTINE add_vnlpsi_us_gpu
!-----------------------------------------------------------------------
  ATTRIBUTES(GLOBAL) SUBROUTINE add_vnlpsi_gpu( lda, outk, npw, nveck, &
                        nb1k, stx, ikblk, npol, hpsi_d, nvecx, nset )
  !-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  INTEGER, INTENT(IN), VALUE :: nset, nvecx
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: ikblk(nset)
  INTEGER, INTENT(IN), DEVICE :: npw(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  INTEGER, INTENT(IN), VALUE :: npol
  COMPLEX(DP), DEVICE :: hpsi_d(lda*npol, nvecx * nset)
  END  SUBROUTINE add_vnlpsi_gpu

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE vlocpsi_gpu_vp(outk, nveck, st, ikt, npol, &
                              psicr, nvec, nnr, nset)
!-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP

  INTEGER, INTENT(IN), VALUE  :: nset, nnr, nvec
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  INTEGER, INTENT(IN), DEVICE :: ikt(nset)
  INTEGER, INTENT(IN), VALUE  :: npol
  REAL(DP), INTENT(INOUT), DEVICE :: psicr(2, nnr, npol, nvec * nset)

END SUBROUTINE vlocpsi_gpu_vp

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE vlocpsi_gpu_setpsi_zero( outk, nveck, st, &
                   npol, psicr, nvec, nnr, nset)
!-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE  :: nset, nnr, nvec
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  INTEGER, INTENT(IN), VALUE  :: npol
  REAL(DP), INTENT(INOUT), DEVICE :: psicr(2, nnr, npol, nvec * nset)
END SUBROUTINE vlocpsi_gpu_setpsi_zero
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE put_psi_on_grid( lda, outk, npw, nveck,  &
                nb1k, st, stx, ikt, npol, psi_d, psicr, nvec, nvecx, nnr, &
                nset)
!-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN), VALUE :: lda
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nvecx, nvec
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: ikt(nset)
  INTEGER, INTENT(IN), DEVICE :: npw(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  INTEGER, INTENT(IN), VALUE :: npol
  COMPLEX(DP), INTENT(IN), DEVICE :: psi_d(lda*npol, nvecx * nset)
  REAL(DP), INTENT(INOUT), DEVICE :: psicr(2, nnr, npol, nvec * nset)

END SUBROUTINE put_psi_on_grid
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE add_grid_to_hpsi(lda, outk, npw, nveck,   &
                nb1k, st, stx, ikt, npol, hpsi_d, psicr, nvec, nvecx, nnr, &
                nset)
!-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN), VALUE :: lda
  INTEGER, INTENT(IN), VALUE :: nset, nvecx, nnr, nvec
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: ikt(nset)
  INTEGER, INTENT(IN), DEVICE :: npw(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  INTEGER, INTENT(IN), VALUE :: npol
  COMPLEX(DP), INTENT(INOUT), DEVICE :: hpsi_d(lda*npol, nvecx * nset)
  REAL(DP), INTENT(IN), DEVICE :: psicr(2, nnr, npol, nvec * nset)

END SUBROUTINE add_grid_to_hpsi

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE fft1inv_dev(outk, nveck, st, npol, psicr, &
                   nvec, nr1, nr2, nr3, nr1x, nr2x, nnr, nset, adim)
!-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nr1, nr2, nr3, nr1x, nr2x, nvec, &
                                adim
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  INTEGER, INTENT(IN), VALUE :: npol
  REAL(DP), INTENT(INOUT), DEVICE :: psicr(2, nnr, npol, nvec * nset)
END SUBROUTINE fft1inv_dev

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE fft2inv_dev(outk, nveck, st, npol, psicr, &
                   nvec, nr1, nr2, nr3, nr1x, nr2x, nnr, nset, adim)
!-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nr1, nr2, nr3, nr1x, nr2x, nvec, &
                                adim
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  INTEGER, INTENT(IN), VALUE :: npol
  REAL(DP), INTENT(INOUT), DEVICE :: psicr(2, nnr, npol, nvec * nset)
END SUBROUTINE fft2inv_dev
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE fft3inv_dev(outk, nveck, st, npol, psicr, &
                   nvec, nr1, nr2, nr3, nr1x, nr2x, nnr, nset, adim)
!-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nr1, nr2, nr3, nr1x, nr2x, nvec, &
                                adim
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  INTEGER, INTENT(IN), VALUE :: npol
  REAL(DP), INTENT(INOUT), DEVICE :: psicr(2, nnr, npol, nvec * nset)
END SUBROUTINE fft3inv_dev
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE fft1fwd_dev(outk, nveck, st, npol, psicr, &
                   nvec, nr1, nr2, nr3, nr1x, nr2x, nnr, nset, adim)
!-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN), VALUE  :: nset, nnr, nr1, nr2, nr3, nr1x, nr2x, nvec, &
                                 adim
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  INTEGER, INTENT(IN), VALUE  :: npol
  REAL(DP), INTENT(INOUT), DEVICE :: psicr(2, nnr, npol, nvec * nset)

END SUBROUTINE fft1fwd_dev
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE fft2fwd_dev(outk, nveck, st, npol, psicr, &
                   nvec, nr1, nr2, nr3, nr1x, nr2x, nnr, nset, adim)
  !-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nr1, nr2, nr3, nr1x, nr2x, nvec, &
                                adim
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  INTEGER, INTENT(IN), VALUE :: npol
  REAL(DP), INTENT(INOUT), DEVICE :: psicr(2, nnr, npol, nvec * nset)
END SUBROUTINE fft2fwd_dev
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE fft3fwd_dev(outk, nveck, st, npol, psicr, &
                   nvec, nr1, nr2, nr3, nr1x, nr2x, nnr, nset, adim)
!-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nr1, nr2, nr3, nr1x, nr2x, nvec, &
                                adim
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  INTEGER, INTENT(IN), VALUE  :: npol
  REAL(DP), INTENT(INOUT), DEVICE :: psicr(2, nnr, npol, nvec * nset)
END SUBROUTINE fft3fwd_dev

!-------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE compute_ps_gpu(outk, nveck, ikt, nset)
  !-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nset
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  INTEGER, INTENT(IN), DEVICE :: ikt(nset)

  END SUBROUTINE compute_ps_gpu
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE compute_ps_nc_gpu(outk, nveck, nset)
  !-----------------------------------------------------------------------
  !
  ! This routine computes becp for all k points.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nset
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)

END SUBROUTINE compute_ps_nc_gpu

END INTERFACE h_psi_gpu_interf
#endif
