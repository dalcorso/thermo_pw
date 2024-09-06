!
! Copyright (C) 2024 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
INTERFACE ke_g2kin_interf
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE ke_g2kin( ikb, nk, g2kink_d, npwx)
!-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE klist,          ONLY : ngk => ngk_d, igk_k => igk_k_d
  USE many_k_mod,     ONLY : qcutz => qcutz_d, ecfixed => ecfixed_d, &
                             q2sigma => q2sigma_d, tpiba => tpiba_d, &
                             xk => xk_d, g => g_d
  USE many_k_ph_mod,  ONLY : ikks=>ikks_d, ikqs=>ikqs_d, ikmks=>ikmks_d,    &
                             ikmkmqs=>ikmkmqs_d, startkb_ph => startkb_ph_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, npwx, ikb
  REAL(DP), DEVICE, INTENT(INOUT) :: g2kink_d(npwx, nk)

END SUBROUTINE ke_g2kin
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE ke_eprec( ikb, nk, g2kink_d, evqk_d, &
                                        eprec_d, nbnd, npol, npwx, nsolv)
!-----------------------------------------------------------------------
  !
  !  This routine computes the array eprec for all k points
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE klist,          ONLY : ngk => ngk_d
  USE many_k_mod,     ONLY : noncolin_d
  USE many_k_ph_mod,  ONLY : ikks=>ikks_d, ikqs=>ikqs_d,    &
                             startkb_ph => startkb_ph_d, &
                             nbnd_occ => nbnd_occ_d
  !
IMPLICIT NONE
!
INTEGER, INTENT(IN), VALUE :: nk, npwx, npol, ikb, nbnd, nsolv
REAL(DP), DEVICE, INTENT(IN) :: g2kink_d(npwx, nk)
COMPLEX(DP), DEVICE, INTENT(IN) :: evqk_d(npwx*npol, nbnd*nk)
REAL(DP), DEVICE, INTENT(INOUT) :: eprec_d(nbnd*nk)

END SUBROUTINE ke_eprec

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE ke_hprec( st, ikb, nk, g2kink_d, h_diag_ph_d, &
                                        eprec_d, nbnd, npe, nsolv, npol, npwx)
!-----------------------------------------------------------------------
  USE cudafor
  USE util_param,     ONLY : DP
  USE klist,          ONLY : ngk => ngk_d 
  USE many_k_mod,     ONLY : noncolin_d
  USE many_k_ph_mod,  ONLY : ikks=>ikks_d, ikqs=>ikqs_d,   &
                             startkb_ph => startkb_ph_d,   &
                             nbnd_occ => nbnd_occ_d 
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, npwx, npol, ikb, nbnd, npe, nsolv
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  REAL(DP), DEVICE, INTENT(IN) :: g2kink_d(npwx, nk)
  REAL(DP), DEVICE, INTENT(INOUT) :: h_diag_ph_d(npwx*npol, &
                                                     nbnd*nk*npe*nsolv)
  REAL(DP), DEVICE, INTENT(IN) :: eprec_d(nbnd*nk)
END SUBROUTINE ke_hprec
END INTERFACE ke_g2kin_interf
#endif
