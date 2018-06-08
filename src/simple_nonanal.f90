!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE simple_nonanal(nat, epsil, q, zeu, omega, dyn )
  !-----------------------------------------------------------------------
  !     add the nonanalytic term with macroscopic electric fields
  !
  USE kinds, ONLY: dp
  USE constants, ONLY: fpi, e2
  USE io_global, ONLY : stdout
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: nat
 !
 COMPLEX(DP), INTENT(INOUT) :: dyn(3,3,nat,nat) ! dynamical matrix
 REAL(DP), INTENT(IN) :: q(3),  &! polarization vector
      &       epsil(3,3),     &! dielectric constant tensor
      &       zeu(3,3,nat),   &! effective charges tensor
      &       omega            ! unit cell volume
 !
 ! local variables
 !
 REAL(DP)::  zag(3),zbg(3),  &! eff. charges  times g-vector
      &      qeq              !  <q| epsil | q>
 INTEGER :: na,nb,           &! counters on atoms
      &  i,j                  ! counters on cartesian coordinates
 !
 qeq = (q(1)*(epsil(1,1)*q(1)+epsil(1,2)*q(2)+epsil(1,3)*q(3))+    &
        q(2)*(epsil(2,1)*q(1)+epsil(2,2)*q(2)+epsil(2,3)*q(3))+    &
        q(3)*(epsil(3,1)*q(1)+epsil(3,2)*q(2)+epsil(3,3)*q(3)))
 !
 IF (qeq < 1.d-8) THEN
    WRITE(stdout,'(5x,"A direction for q was not specified:", &
      &          "TO-LO splitting will be absent")')
    RETURN
 END IF
 !
 DO na = 1,nat
    DO nb = 1,nat
       !
       DO i=1,3
          !
          zag(i) = q(1)*zeu(1,i,na) +  q(2)*zeu(2,i,na) + &
                   q(3)*zeu(3,i,na)
          zbg(i) = q(1)*zeu(1,i,nb) +  q(2)*zeu(2,i,nb) + &
                   q(3)*zeu(3,i,nb)
       ENDDO
       !
       DO i = 1,3
          DO j = 1,3
             dyn(i,j,na,nb) = dyn(i,j,na,nb)+ fpi*e2*zag(i)*zbg(j)/qeq/omega
!             print*, zag(i),zbg(j),qeq, fpi*e2*zag(i)*zbg(j)/qeq/omega
          ENDDO
       ENDDO
    ENDDO
 ENDDO
 !
 RETURN
END SUBROUTINE simple_nonanal

