! Copyright (C) 2023 Andrea Dal Corso  
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
INTERFACE
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE cgsolve_all_loop2(ndmx, outk, st, conv, &
         lbndk, h, g, h_diag, rho, current_ikb_ph, npol, nk, &
         npe, nsolv, nbnd, my_nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, my_nbnd, current_ikb_ph

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: lbndk(my_nbnd*nk*npe*nsolv)
 INTEGER, DEVICE :: conv(nbnd*nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: h(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: g(ndmx*npol,my_nbnd*nk*npe*nsolv)
 REAL(DP), DEVICE :: h_diag(ndmx*npol,nbnd*nk*npe*nsolv)
 REAL(DP), DEVICE :: rho(my_nbnd*nk*npe*nsolv)

 END SUBROUTINE cgsolve_all_loop2
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE cgsolve_all_loop4(ndmx, outk, st, conv, &
         lbndk, h, hold, dcgammak, eu, e, current_ikb_ph, npol, nk,    &
         npe, nsolv, nbnd, my_nbnd, iter, nks)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, my_nbnd, &
                   current_ikb_ph, iter, nks

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: lbndk(my_nbnd*nk*npe*nsolv)
 INTEGER, DEVICE :: conv(nbnd*nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: h(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: hold(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: dcgammak(my_nbnd*nk*npe*nsolv)
 REAL(DP), DEVICE :: eu(my_nbnd*nk*npe*nsolv)
 REAL(DP), DEVICE :: e(nbnd,nks)

 END SUBROUTINE cgsolve_all_loop4
 !
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE cgsolve_all_loop5(ndmx, outk, st, conv, &
         lbndk, g, h, t, a, c, current_ikb_ph, npol, nk, &
         npe, nsolv, nbnd, my_nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, my_nbnd, current_ikb_ph

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: lbndk(my_nbnd*nk*npe*nsolv)
 INTEGER, DEVICE :: conv(nbnd*nk*npe*nsolv)

 REAL(DP), DEVICE :: a(my_nbnd*nk*npe*nsolv)
 REAL(DP), DEVICE :: c(my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: g(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: h(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: t(ndmx*npol,my_nbnd*nk*npe*nsolv)

 END SUBROUTINE cgsolve_all_loop5
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE cgsolve_all_loop6(ndmx, outk, st, conv,  &
       lbndk, dpsi, g, h, hold, t, dclambdak, current_ikb_ph, npol, nk, &
       npe, nsolv, nbnd, my_nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, my_nbnd, current_ikb_ph

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: lbndk(my_nbnd*nk*npe*nsolv)
 INTEGER, DEVICE :: conv(nbnd*nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: dpsi(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: g(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: h(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: hold(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: t(ndmx*npol,my_nbnd*nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: dclambdak(my_nbnd*nk*npe*nsolv)

 END SUBROUTINE cgsolve_all_loop6

END INTERFACE
#endif
