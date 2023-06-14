INTERFACE compute_deff_interf
!-------------------------------------------------------------------------
ATTRIBUTES(DEVICE) SUBROUTINE compute_deff_dev(nhm, nat, current_spin, okvan, deff, et )
  !-----------------------------------------------------------------------
  !
  USE cudafor
  USE kinds,       ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER, VALUE :: current_spin, nhm, nat
  LOGICAL, VALUE :: okvan
  !
  REAL(DP), INTENT(IN), VALUE :: et
  !
  REAL(DP), INTENT(OUT), DEVICE :: deff(nhm,nhm,nat)

END SUBROUTINE compute_deff_dev

!-------------------------------------------------------------------------
ATTRIBUTES(DEVICE) SUBROUTINE compute_deff_nc_dev( nhm, nat, nspin, isolv, &
                                      npol, okvan, nsolv, deff_nc, et )
  !-----------------------------------------------------------------------
  !
  USE cudafor
  USE kinds,       ONLY: DP
  IMPLICIT NONE
  !
  INTEGER, VALUE :: nspin, nhm, nat, npol, isolv, nsolv
  LOGICAL, VALUE :: okvan
  !
  REAL(DP), INTENT(IN), VALUE :: et
  COMPLEX(DP), INTENT(OUT), DEVICE :: deff_nc(nhm,nhm,nat,nspin)

END SUBROUTINE compute_deff_nc_dev

END INTERFACE compute_deff_interf

