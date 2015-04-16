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
  !  This routine receives the celldm and computes the
  !  volume, using ibrav read from input. The celldm read from input is
  !  kept.
  !
  USE kinds,            ONLY : DP
  USE input_parameters, ONLY : celldm, a, b, c, cosab, cosac, cosbc, &
                               trd_ht, rd_ht, cell_units
  USE cell_base,        ONLY : at, cell_base_init, omega, alat
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ibrav
  REAL(DP), INTENT(IN) :: celldm_(6)
  REAL(DP) :: compute_omega_geo, celldm_save(6)
  !
  celldm_save=celldm
  celldm=celldm_
  CALL cell_base_init ( ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                        trd_ht, rd_ht, cell_units )
  celldm=celldm_save
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

