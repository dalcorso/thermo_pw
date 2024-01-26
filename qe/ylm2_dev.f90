!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! Copyright (C) 2023 Andrea Dal Corso for the device routine
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
  !
  !----------------------------------------------------------------------
  ATTRIBUTES(GLOBAL) SUBROUTINE ylm2_dev(ylm_d, lmaxkb, npwx, nk, ikt)
  !----------------------------------------------------------------------
  !
    USE kinds,        ONLY : DP
    USE many_k_mod,   ONLY : xk=>xk_d, g=>g_d
    USE klist,        ONLY : igk=>igk_k_d, ngk=>ngk_d

    !
    IMPLICIT NONE

#include<ylmr2_interf.f90>
    !
    INTEGER, INTENT(IN), VALUE :: lmaxkb, npwx, nk
    INTEGER, INTENT(IN), DEVICE :: ikt(nk)
    !! number of PWs 
    REAL(DP), INTENT(OUT), DEVICE :: ylm_d((lmaxkb + 1) **2, npwx, nk)
    !! beta functions 
    !
    !     Local variables
    !
    INTEGER :: ig, iv_d, ik, ik1, np
    REAL(DP) :: q1, q2, q3, rv_d

    REAL(DP) :: gk_d(3), qg_d
    !
    IF (lmaxkb<0) RETURN
    !
    ik1=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
    IF (ik1>nk) RETURN
    ik=ikt(ik1)
    np=ngk(ik)
    ig=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
    IF (ig>np) RETURN

    q1 = xk(1,ik)
    q2 = xk(2,ik)
    q3 = xk(3,ik)

    iv_d = igk(ig,ik)
    gk_d (1) = q1 + g(1, iv_d )
    gk_d (2) = q2 + g(2, iv_d )
    gk_d (3) = q3 + g(3, iv_d )
    qg_d = gk_d(1)*gk_d(1) + &
           gk_d(2)*gk_d(2) + &
           gk_d(3)*gk_d(3)
    !
    CALL ylmr2_dev ((lmaxkb+1)**2, gk_d, qg_d, ylm_d(1,ig,ik1))
    !
    END SUBROUTINE ylm2_dev
#endif
  !
  !
