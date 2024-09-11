!
! Copyright (C) 2024 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
!
INTERFACE wgauss_interf
!-----------------------------------------------------------------------
ATTRIBUTES(DEVICE) SUBROUTINE wgauss_dev (x, n, wgauss)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : pi
  implicit none
  real(DP), DEVICE :: wgauss
  real(DP), VALUE :: x
  integer, VALUE :: n
END SUBROUTINE wgauss_dev
!
!-----------------------------------------------------------------------
ATTRIBUTES(DEVICE) SUBROUTINE w0gauss_dev (x, n, w0gauss)
  !-----------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE constants, ONLY : sqrtpm1
  implicit none
  real(DP), DEVICE :: w0gauss
  real(DP), VALUE :: x
  integer, VALUE :: n
END SUBROUTINE w0gauss_dev

END INTERFACE wgauss_interf

#endif
