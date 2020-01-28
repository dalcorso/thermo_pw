!
! Copyright (C) 2019 Andrea Urru
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE sym_and_write_magnetic_charges()
  !
  !  This routine brings the dynamical magnetic charges to the cartesian
  !  basis from the basis of the random modes and writes them on 
  !  output in mu_B/A.
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, atm, ityp
  USE io_global,        ONLY : stdout
  USE modes,            ONLY : u
  USE magnetic_charges, ONLY : mag_charge_mode, mag_charge
  USE constants,        ONLY : bohr_radius_angs
  
  IMPLICIT NONE

  INTEGER :: icart, jcart, na, nu, mu
  ! counter on cartesian coordinates
  ! counter on atoms and modes
  ! counter on modes
  REAL(DP) :: mag_charge_loc(3,nat,3)

  DO jcart = 1, 3
     DO mu = 1, 3 * nat
        na = (mu - 1) / 3 + 1
        icart = mu - 3 * (na - 1)
        DO nu = 1, 3 * nat
           mag_charge (icart, na, jcart) = mag_charge (icart, na, jcart) + &
                u (mu, nu) * mag_charge_mode (nu, jcart)
        ENDDO
     ENDDO
  ENDDO

  mag_charge_loc = REAL(mag_charge) / bohr_radius_angs

  WRITE( stdout, '(/,10x,"Dynamical magnetic charges (dM/du) &
                                     &in cartesian axis (mu_B/A)",/)')
  DO na = 1, nat
     WRITE(stdout, '(/,10x,"atom # ",i6,a6)') na, atm(ityp(na))
     WRITE(stdout, '(6x,"u_x (",3f15.5," )")') (mag_charge_loc(1, na, jcart),& 
              jcart = 1,3)
     WRITE(stdout, '(6x,"u_y (",3f15.5," )")') (mag_charge_loc(2, na, jcart),& 
              jcart = 1,3)
     WRITE(stdout, '(6x,"u_z (",3f15.5," )")') (mag_charge_loc(3, na, jcart),& 
              jcart = 1,3)
  ENDDO

  RETURN
END SUBROUTINE sym_and_write_magnetic_charges
