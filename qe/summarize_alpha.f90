!
! Copyright (C) 2019 Andrea Urru
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE summarize_alpha()
!
!  This routine converts the magnetoelectric tensor from atomic units
!  to the SI units (ps/m) and prints the results.
!
  USE kinds,                 ONLY : DP
  USE magnetic_charges,      ONLY : alpha_me
  USE io_global,             ONLY : stdout
  USE constants,             ONLY : h_planck_si, electronmass_si, & 
                                    c_si, bohr_radius_si  

  IMPLICIT NONE

  INTEGER:: jcart
  REAL(DP) :: conv_factor
  REAL(DP) :: alpha(3,3)

  conv_factor = h_planck_si/(electronmass_si*c_si**2*bohr_radius_si)*1.0D+12
  alpha = REAL(alpha_me)*conv_factor

  WRITE(stdout, '(/,5x,"Magnetoelectric tensor in S.I. units (ps/m)", /)') 
  
  WRITE(stdout, '(6x,"E_x (",3f15.5," )")') (alpha (1, jcart), jcart = 1,3)
  WRITE(stdout, '(6x,"E_y (",3f15.5," )")') (alpha (2, jcart), jcart = 1,3)
  WRITE(stdout, '(6x,"E_z (",3f15.5," )")') (alpha (3, jcart), jcart = 1,3)

  RETURN
END SUBROUTINE summarize_alpha
