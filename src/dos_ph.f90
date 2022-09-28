!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
SUBROUTINE dos_ph (freq0, nmodes, nq, wq, degauss, ngauss, freq, dosg)
  !--------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER :: nq, nmodes, ngauss

  real(DP) :: wq (nq), freq0 (nmodes, nq), degauss, freq, dosg 
  real(DP) :: w0gauss
  REAL(DP), PARAMETER :: thr_ph=1.D-3
  INTEGER :: n, iq
  EXTERNAL w0gauss
  !
  dosg = 0.0d0
  DO iq = 1, nq
     DO n = 1, nmodes
        IF (freq0(n, iq)>thr_ph) &
           dosg = dosg + wq (iq) * w0gauss ( (freq-freq0 (n, iq) ) &
                / degauss, ngauss)
     ENDDO
  ENDDO
  !
  dosg = dosg / degauss
  !
  RETURN
END SUBROUTINE dos_ph
