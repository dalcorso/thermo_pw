!
! Copyright (C) 2023 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
  !
  !----------------------------------------------------------------------
#if defined(__CUDA)
INTERFACE initus2
  ATTRIBUTES(GLOBAL) SUBROUTINE init_us_2_kernel(vkbk_d, ylm_d, vkb1_d, &
                      nhm, lmaxkb, nkb, npwx, nk, ikt )
    !
    USE cudafor
    USE kinds,        ONLY : DP
    !
    IMPLICIT NONE

    INTEGER, INTENT(IN), VALUE :: npwx, lmaxkb, nhm, nkb, nk
    COMPLEX(DP), INTENT(OUT), DEVICE :: vkbk_d(npwx,nkb*nk)
    INTEGER, INTENT(IN), DEVICE :: ikt(nk)
    REAL(DP), INTENT(IN), DEVICE :: ylm_d((lmaxkb + 1) **2, npwx, nk)
    REAL(DP), INTENT(INOUT), DEVICE :: vkb1_d(nhm, npwx, nk)
  END SUBROUTINE init_us_2_kernel
END INTERFACE initus2
#endif
