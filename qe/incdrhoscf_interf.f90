!
! Copyright (C) 2023 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
#if defined(__CUDA)
INTERFACE incdrhoscf_interf
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE setdpsi_zero(nbndk, st, npol, dpsicr, &
                                           nbnd, nnrs, nk, npe)
!-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, npe, nnrs, nbnd
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe)
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe)
  INTEGER, INTENT(IN), VALUE :: npol
  REAL(DP), DEVICE, INTENT(INOUT) :: dpsicr(2, nnrs, npol, nbnd * nk * npe)

  END SUBROUTINE setdpsi_zero

! -----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE set_dpsi_grid( lda, nbndk, st, npol, dpsik, &
                                        dpsicr, nbnd, nnrs, nk, npe, ikb)
!-----------------------------------------------------------------------
!
  USE cudafor
  USE util_param,  ONLY : DP

  IMPLICIT NONE

  INTEGER, INTENT(IN), VALUE :: lda
  INTEGER, INTENT(IN), VALUE :: nk, npe, nnrs, nbnd, ikb
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe)
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe)
  INTEGER, INTENT(IN), VALUE :: npol
  COMPLEX(DP), DEVICE, INTENT(IN) :: dpsik(lda*npol, nbnd * nk * npe)
  REAL(DP), DEVICE, INTENT(INOUT) :: dpsicr(2, nnrs, npol, nbnd * nk * npe)
  !
  END SUBROUTINE set_dpsi_grid

! -----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE incdrho_dev( nbndk, st, npol, drhoscf, &
     psicr, dpsicr, nnrs, nnr, nk, npe, nbnd, ikb, nspin_mag, omega)
!-----------------------------------------------------------------------
!
  USE cudafor
  USE kinds, ONLY : DP
  IMPLICIT NONE

  REAL(DP), INTENT(IN), VALUE :: omega
  INTEGER, INTENT(IN), VALUE :: nk, npe, nnrs, nnr, nbnd, ikb, nspin_mag
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe)
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe)
  INTEGER, INTENT(IN), VALUE :: npol
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: drhoscf(nnr, nspin_mag, npe)
  REAL(DP), DEVICE, INTENT(IN) :: psicr(2, nnrs, npol, nbnd * nk * npe)
  REAL(DP), DEVICE, INTENT(IN) :: dpsicr(2, nnrs, npol, nbnd * nk * npe)

  END SUBROUTINE incdrho_dev
 !
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE incdrho_calbec(ndmx, st, nbndk, npwk, &
                           dpsi, current_ikb_ph, npol, nk, npe, nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, npol, nbnd, current_ikb_ph

 INTEGER, DEVICE :: nbndk(nk*npe)
 INTEGER, DEVICE :: st(nk*npe)
 INTEGER, DEVICE :: npwk(nk*npe)

 COMPLEX(DP), DEVICE :: dpsi(ndmx*npol,nbnd*nk*npe)

 END SUBROUTINE incdrho_calbec
 !
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE incdrho_addusdbec(st, nbndk, current_ikb_ph, &
                               npol, nk, npe)
 !--------------------------------------------------------------------------
 !
 USE cudafor
 USE util_param,    ONLY : DP

 IMPLICIT NONE
 INTEGER, VALUE :: nk, npe, npol, current_ikb_ph
 INTEGER, DEVICE :: nbndk(nk*npe)
 INTEGER, DEVICE :: st(nk*npe)

END SUBROUTINE incdrho_addusdbec

END INTERFACE incdrhoscf_interf
#endif
