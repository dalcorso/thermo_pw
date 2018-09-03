!
! Copyright (C) 2018  Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

SUBROUTINE tra_write_zstar(zstareu0_wrk, zstareu, flag)

USE kinds,     ONLY : DP
USE ions_base, ONLY : nat, atm, ityp
USE cell_base, ONLY : bg 
USE modes,     ONLY : u
USE io_global, ONLY : stdout

IMPLICIT NONE

REAL(DP) :: zstareu(3,3,nat)
COMPLEX(DP) :: zstareu0_wrk (3,3*nat)
LOGICAL :: flag

INTEGER :: jpol, icart, mu, nu, na

zstareu (:,:,:) = 0.0_DP
DO jpol = 1, 3
   DO mu = 1, 3 * nat
      na = (mu - 1) / 3 + 1
      icart = mu - 3 * (na - 1)
      DO nu = 1, 3 * nat
         zstareu (jpol, icart, na) = zstareu (jpol, icart, na) + &
              CONJG(u (mu, nu) ) * ( zstareu0_wrk (1,nu) * bg(jpol,1) + &
                                     zstareu0_wrk (2,nu) * bg(jpol,2) + &
                                     zstareu0_wrk (3,nu) * bg(jpol,3) )
      ENDDO
   ENDDO
ENDDO

IF (flag) THEN
   DO na = 1, nat
      WRITE( stdout, '(10x," atom ",i6, a6)') na, atm(ityp(na))
      WRITE( stdout, '(6x,"Ex  (",3f10.5," )")')  (zstareu(1,jpol,na), &
             jpol = 1, 3)
      WRITE( stdout, '(6x,"Ey  (",3f10.5," )")')  (zstareu(2,jpol,na), &
             jpol = 1, 3)
      WRITE( stdout, '(6x,"Ez  (",3f10.5," )")')  (zstareu(3,jpol,na), &
             jpol = 1, 3)
   ENDDO
ENDIF

RETURN
END SUBROUTINE tra_write_zstar
