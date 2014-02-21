!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
FUNCTION compute_omega_geo(current_alat)
  !-----------------------------------------------------------------------
  !
  !  This routine receives the lattice constant and computes the
  !  volume, using ibrav and celldm read from input.
  !
  USE kinds,            ONLY : DP
  USE input_parameters, ONLY : ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                               trd_ht, rd_ht, cell_units
  USE cell_base,        ONLY : cell_base_init, omega, alat
  !
  IMPLICIT NONE
  REAL(DP) :: compute_omega_geo, alat_save
  REAL(DP), INTENT(IN) :: current_alat
  !
  celldm(1)=current_alat
  alat_save=alat
  CALL cell_base_init ( ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                        trd_ht, rd_ht, cell_units )
  alat=alat_save
  compute_omega_geo = omega

  RETURN
END FUNCTION compute_omega_geo

!-----------------------------------------------------------------------
FUNCTION compute_alat_geo(omega, alat1, omega1)
  !-----------------------------------------------------------------------
  !
  !  This routine receives the volume, the lattice constant alat1 corresponding 
  !  to omega1 and computes the lattice constant corresponding to omega.
  !
  USE kinds,            ONLY : DP
  !
  IMPLICIT NONE
  REAL(DP) :: compute_alat_geo
  REAL(DP), INTENT(IN) :: omega, alat1, omega1
  !
  compute_alat_geo = alat1 * ( omega / omega1 ) ** ( 1.0_DP / 3.0_DP )
  !
  RETURN
END FUNCTION compute_alat_geo
