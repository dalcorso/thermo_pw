!
! Copyright (C) 2019 Andrea Urru
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE rotate_and_write_charges(u,zeta)
  !
  ! This routine rotates the effective charge tensor (given as input)
  ! onto the eigenvectors of the dynamical matrix obtained from the
  ! displacements.
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE ions_base,        ONLY : nat, ityp, amass

  IMPLICIT NONE 

  INTEGER:: nu, ia, ita, icart, jcart
  REAL(DP):: norm
  COMPLEX(DP):: zeta_mode(3)
  COMPLEX(DP), INTENT(IN):: u(3*nat, 3*nat), zeta(3, nat, 3)
  COMPLEX(DP) :: eigen(3*nat, 3*nat)

  DO nu = 1, 3 * nat
     norm = 0.0_DP
     DO ia = 1, nat
        ita = ityp (ia)
        DO icart = 1, 3
           eigen(3 * (ia - 1) + icart, nu) = u(3 * (ia - 1) + icart, nu) &
                                              * SQRT(amass(ita))
           norm = norm + eigen(3 * (ia - 1) + icart, nu) * &
                   CONJG(eigen(3 * (ia - 1) + icart, nu))
        ENDDO
     ENDDO
     eigen(:, nu) = eigen(:, nu) / SQRT(norm)
  ENDDO

  DO nu = 1, 3 * nat
     zeta_mode = (0.0_DP, 0.0_DP)
     DO jcart = 1, 3
        DO icart = 1, 3
           DO ia=1 , nat
              zeta_mode(jcart) = zeta_mode(jcart) + zeta(icart, ia, jcart) &
                               * u(3 * (ia - 1) + icart, nu)
           ENDDO
        ENDDO
     ENDDO
     WRITE(stdout, '(6x,"mode #", i6, " charge (",3f15.5, ")")') nu, &
                                    (REAL(zeta_mode (icart)), icart = 1, 3)
  ENDDO

  RETURN
END SUBROUTINE rotate_and_write_charges
