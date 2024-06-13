!
! Copyright (C) 2023 Andrea Dal Corso  
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
INTERFACE ch_psi_all_interf
 !
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE ch_psi_computeps(ndmx, outk, kdimk, st, &
         nbndk, evqk, spsi, ps, current_ikb_ph, npol, nk,        &
         npe, nsolv, nbnd, my_nbnd, alpha_pv)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, my_nbnd, current_ikb_ph
 REAL(DP), VALUE :: alpha_pv

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: nbndk(nk*npe*nsolv)
 INTEGER, DEVICE :: kdimk(nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: evqk(ndmx*npol,my_nbnd*nk*nsolv)
 COMPLEX(DP), DEVICE :: spsi(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: ps(nbnd,my_nbnd*nk*npe*nsolv)

 END SUBROUTINE ch_psi_computeps
 !
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE ch_psi_ah(ndmx, outk, st, &
         nbndk, ah, hpsi, spsi, eu, current_ikb_ph, npol, nk, &
         npe, nsolv, nbnd, my_nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, my_nbnd, current_ikb_ph

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: nbndk(nk*npe*nsolv)

 REAL(DP), DEVICE :: eu(my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: ah(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: hpsi(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: spsi(ndmx*npol,my_nbnd*nk*npe*nsolv)

 END SUBROUTINE ch_psi_ah
 !
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE ch_psi_sqdpsi(lda, outk, npw, st,  &
         nbndk, hpsi, ah, current_ikb_ph, npol, nk, npe, nsolv, nkb, &
         nbnd, my_nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 IMPLICIT NONE

 INTEGER, VALUE :: lda, nk, npe, nsolv, npol, nbnd, my_nbnd, nkb, &
                   current_ikb_ph

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: npw(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: nbndk(nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: ah(lda*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: hpsi(lda*npol,my_nbnd*nk*npe*nsolv)

 END SUBROUTINE ch_psi_sqdpsi
 !
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE ch_psi_sdpsi(lda, outk, st, nbndk, hpsi, &
         ah, current_ikb_ph, npol, nk, npe, nsolv, nkb, nbnd, my_nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 IMPLICIT NONE

 INTEGER, VALUE :: lda, nk, npe, nsolv, npol, nbnd, my_nbnd, nkb, &
                   current_ikb_ph

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: nbndk(nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: ah(lda*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: hpsi(lda*npol,my_nbnd*nk*npe*nsolv)

 END SUBROUTINE ch_psi_sdpsi
 !
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE ch_psi_lo2(ndmx, outk, st, &
         nbndk, evqk, hpsi, ps, current_ikb_ph, npol, nk,      &
         npe, nsolv, nbnd, my_nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, my_nbnd, current_ikb_ph

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: nbndk(nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: evqk(ndmx*npol,my_nbnd*nk*nsolv)
 COMPLEX(DP), DEVICE :: hpsi(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: ps(nbnd,my_nbnd*nk*npe*nsolv)

 END SUBROUTINE ch_psi_lo2
 !
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE ch_psi_calbec(ndmx, outk, st,  &
         nbndk, npw, hpsi, current_ikb_ph, npol, nk,   &
         npe, nsolv, nkb, nbnd, my_nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, my_nbnd, nkb, &
                   current_ikb_ph

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: nbndk(nk*npe*nsolv)
 INTEGER, DEVICE :: npw(nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: hpsi(ndmx*npol,my_nbnd*nk*npe*nsolv)
 END SUBROUTINE ch_psi_calbec
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE ch_psi_qqps( outk, nbndk, nkb, nk, npe, &
                                            nsolv, npol )
  !-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, nkb, npe, nsolv, npol
  LOGICAL, INTENT(IN), DEVICE :: outk(nk*npe*nsolv)
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  !
  END SUBROUTINE ch_psi_qqps
END INTERFACE ch_psi_all_interf
#endif
