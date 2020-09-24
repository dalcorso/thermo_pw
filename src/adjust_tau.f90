!
! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE adjust_tau(tau_save, tau, at)
!-----------------------------------------------------------------------
!
!  This routine receives as input the crystal coordinates of the atoms in the
!  original cell and transforms them in cartesian coordinates with the new
!  at. This corresponds a uniform strain of the atomic coordinates  that 
!  should be already in units of the new celldm
!
USE kinds,          ONLY : DP
USE ions_base,      ONLY : nat
IMPLICIT NONE
REAL(DP), INTENT(IN)    :: tau_save(3,nat), at(3,3)
REAL(DP), INTENT(INOUT) :: tau(3,nat)

tau=tau_save
CALL cryst_to_cart( nat, tau, at, 1 )

RETURN
END SUBROUTINE adjust_tau

