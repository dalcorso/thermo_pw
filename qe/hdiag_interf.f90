! Copyright (C) 2023 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
  !
INTERFACE hdiag_interf
  !----------------------------------------------------------------------
#if defined(__CUDA)
  ATTRIBUTES(GLOBAL) SUBROUTINE hdiag_kernel(g2kink, h_diagk, s_diagk, npwx, &
                                        v_of_0, npol, nkb, nk, current_ikb )
    !----------------------------------------------------------------------
    !
    ! Calculate the kinetic energy for several k points on the GPU
    !
    USE cudafor
    USE kinds,        ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER, VALUE, INTENT(IN) :: npwx, nkb, nk, npol, current_ikb
    REAL(DP), VALUE, INTENT(IN) :: v_of_0
    !! number of PWs 
    REAL(DP), DEVICE, INTENT(OUT) :: g2kink(npwx,nk)
    REAL(DP), DEVICE, INTENT(OUT) :: h_diagk(npwx,npol,nk)
    REAL(DP), DEVICE, INTENT(OUT) :: s_diagk(npwx,npol,nk)
    !! beta functions 
    !
    END SUBROUTINE hdiag_kernel
END INTERFACE hdiag_interf
#endif
  !
