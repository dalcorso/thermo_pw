!
! Copyright (C) 2023 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE print_thermo_clocks()
!--------------------------------------------------------------------------
!
!   Print the clocks of thermo_pw
!
USE io_global, ONLY : stdout
IMPLICIT NONE
  !
  WRITE(stdout,'(/,5x,"Times for scf pw.x")')
  CALL print_clock('tpw_scf_pw')
  WRITE(stdout,'(/,5x,"Times for non scf pw.x")')
  CALL print_clock('tpw_nscf_pw')
  WRITE(stdout,'(/,5x,"Times for ph.x ")')
  CALL print_clock('tpw_ph')
  CALL print_clock('tpw_nscf_ph')
  CALL print_clock('tpw_init_ph')
  CALL print_clock('solve_linter')
  CALL print_clock('drhodv')
  !
  RETURN
END SUBROUTINE print_thermo_clocks
