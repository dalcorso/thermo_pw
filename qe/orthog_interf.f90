!
! Copyright (C) 2024 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
INTERFACE orthog_interf
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE orthog_ps( st, nbndk, ikb, nk, npe, nsolv, &
                    dvpsik_d, evqk_d, ortho_ps, npol, npwx, nbnd, nksbx_ph )
!-----------------------------------------------------------------------
  !
  !  This routine computes the scalar product of the vector evqk_d and
  !  dvpsik_d.
  !  Each thread computes one k point, one perturbation, and one band
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE klist,          ONLY : ngk => ngk_d
  USE many_k_mod,     ONLY : lgauss => lgauss_d, ltetra => ltetra_d
  USE many_k_ph_mod,  ONLY : ikks=>ikks_d, iksq=>ikqs_d, ikmks=>ikmks_d, &
                             ikmkmqs=>ikmkmqs_d, startkb_ph => startkb_ph_d, &
                             nbnd_occ => nbnd_occ_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, npe, nsolv, npwx, nbnd, ikb, nksbx_ph
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  INTEGER, INTENT(IN), VALUE :: npol
  COMPLEX(DP), DEVICE, INTENT(IN) :: evqk_d(npwx*npol, nbnd*nk*nsolv)
  COMPLEX(DP), DEVICE, INTENT(IN) :: dvpsik_d(npwx*npol, nbnd*nk*npe*nsolv)
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: ortho_ps(nk*nbnd*nsolv, &
                                                 nbnd*nk*npe*nsolv)
END SUBROUTINE orthog_ps

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE orthog_last( st, nbndk, ikb, nk, npe, nsolv, &
                      dvpsik_d, sevqk_d, ortho_ps, npol, npwx, nbnd, nksbx_ph )
!-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE klist,          ONLY : ngk => ngk_d
  USE many_k_mod,     ONLY : lgauss => lgauss_d, ltetra => ltetra_d,    &
                             degauss=> degauss_d, ngauss => ngauss_d
  USE many_k_ph_mod,  ONLY : ikks=>ikks_d, ikqs=>ikqs_d, ikmks=>ikmks_d,    &
                             ikmkmqs=>ikmkmqs_d, startkb_ph => startkb_ph_d,&
                             nbnd_occ => nbnd_occ_d, ef => ef_d, et => et_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, npe, nsolv, npwx, nbnd, ikb, nksbx_ph
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  INTEGER, INTENT(IN), VALUE :: npol
  COMPLEX(DP), DEVICE, INTENT(IN) :: sevqk_d(npwx*npol, nbnd*nk*nsolv)
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: dvpsik_d(npwx*npol, nbnd*nk*npe*nsolv)
  COMPLEX(DP), DEVICE, INTENT(IN) :: ortho_ps(nksbx_ph*nbnd*nsolv, &
                                                 nbnd*nksbx_ph*npe*nsolv)

END SUBROUTINE orthog_last

END INTERFACE orthog_interf
#endif
