!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
FUNCTION compute_omega_geo(ibrav, celldm_)
  !-----------------------------------------------------------------------
  !
  !  This function receives ibrav and celldm and computes the volume
  !  of the unit cell. 
  !
  USE kinds,            ONLY : DP
  !
  IMPLICIT NONE
  REAL(DP) :: compute_omega_geo
  INTEGER, INTENT(IN) :: ibrav
  REAL(DP), INTENT(IN) :: celldm_(6)
  REAL(DP) :: at(3,3), omega

  CALL latgen( ibrav, celldm_, at(1,1), at(1,2), at(1,3), omega )
  compute_omega_geo = omega

  RETURN
END FUNCTION compute_omega_geo

!-----------------------------------------------------------------------
SUBROUTINE compute_celldm_geo(omega, celldm, celldm1, omega1)
  !-----------------------------------------------------------------------
  !
  !  This routine receives the volume, the celldm1 corresponding 
  !  to omega1 and computes the lattice constant corresponding to omega.
  !  It assumes that celldm(2) ... celldm(6) remain constant.
  !
  USE kinds,            ONLY : DP
  !
  IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: celldm(6)
  REAL(DP), INTENT(IN) :: omega, celldm1(6), omega1
  !
  celldm=celldm1
  celldm(1) = celldm1(1) * ( omega / omega1 ) ** ( 1.0_DP / 3.0_DP )
  !
  RETURN
END SUBROUTINE compute_celldm_geo

