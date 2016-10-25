!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE rotate
!
!  This module contains all the routines that have to do with rotations.
!  Presently it has the follofing:
!  rotate_vect: takes a set of vectors, a rotation matrix and gives the 
!               rotated vectors.  
!  rotate_tensors2: takes a set of rank two tensors, a rotation matrix
!               and gives the rotated rank two tensors.
!
!  find_su2_euler: takes three euler angles psi, theta, phi and gives the
!               Cayley Klein parameters of the su2 rotation
!
!  find_rotation: given three 3D vectors at(3,3) obtained by rotating three
!                 vectors at1(3,3) and the vectors at1(3,3), the routine 
!                 finds the rotation matrix that brings the at1 into the at.
!
!  is_rotation: a function that receives a 3x3 matrix and returns .true.
!               if it is a rotation matrix
!
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  SAVE


  PUBLIC rotate_vect, rotate_tensors2, euler_to_su2, find_rotation, &
         is_rotation

CONTAINS
   SUBROUTINE rotate_vect(rot, n, a, ra, flag)
!
!  if flag is 1 rotate the n vectors according to rot
!  if flag is -1 rotate the n vectors according to rot^-1
!

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   REAL(DP), INTENT(IN) :: rot(3,3), a(3,n)
   REAL(DP), INTENT(OUT) :: ra(3,n)
   INTEGER, INTENT(IN) :: flag

   INTEGER :: ivect, irot, jrot
  
   ra=0.0_DP
   DO ivect=1,n
      IF (flag==1) THEN
         DO irot=1,3
            DO jrot=1,3
               ra(irot,ivect) = ra(irot,ivect) + rot(irot,jrot) * a(jrot,ivect)
            END DO
         END DO
      ELSE
         DO irot=1,3
            DO jrot=1,3
               ra(irot,ivect) = ra(irot,ivect) + rot(jrot,irot) * a(jrot,ivect)
            END DO
         END DO
      ENDIF
   END DO

   RETURN
   END SUBROUTINE rotate_vect

   SUBROUTINE rotate_tensors2(rot, n, a, ra, flag)
!
!  this routine apply a rotation to n tensors of rank 2
!  if flag is 1 rotate the n tensors according to rot
!  if flag is -1 rotate the n tensors according to rot^-1
!

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   REAL(DP), INTENT(IN) :: rot(3,3), a(3,3,n)
   REAL(DP), INTENT(OUT) :: ra(3,3,n)
   INTEGER, INTENT(IN) :: flag

   INTEGER :: ivect, irot, jrot, krot, lrot
  
   ra=0.0_DP
   DO ivect=1,n
      IF (flag==1) THEN
         DO irot=1,3
            DO jrot=1,3
               DO krot=1,3
                  DO lrot=1,3
                     ra(irot,jrot,ivect)=ra(irot,jrot,ivect)+rot(irot,krot) &
                                         *rot(jrot,lrot)*a(krot,lrot,ivect)
                  END DO
               END DO
            END DO
         END DO
      ELSE
         DO irot=1,3
            DO jrot=1,3
               DO krot=1,3
                  DO lrot=1,3
                     ra(irot,jrot,ivect) = ra(irot,jrot,ivect)+rot(krot,irot) &
                                         *rot(lrot,jrot)*a(krot,lrot,ivect)
                  END DO
               END DO
            END DO
         END DO
      ENDIF
   END DO

   RETURN
   END SUBROUTINE rotate_tensors2


SUBROUTINE euler_to_su2(psi,theta,phi,a,b)

USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP), INTENT(IN) :: psi, theta, phi
COMPLEX(DP), INTENT(OUT) :: a, b

REAL(DP) :: arg

arg = -(phi + psi)*0.5_DP
a=COS(theta/2.0_DP) * CMPLX( COS(arg), SIN(arg), KIND=DP )

arg = (phi - psi)*0.5_DP
b=SIN(theta/2.0_DP) * CMPLX( COS(arg), SIN(arg), KIND=DP ) * (0.0_DP, -1.0_DP)

RETURN
END SUBROUTINE euler_to_su2

SUBROUTINE find_rotation(at,at1,sr)
!
!   This routine receives six vectors at and at1. The three vectors
!   at are obtained by applying a matrix sr to the vectors at1. 
!   The routines gives as output the matrix sr. You can
!   check if it is a rotation with the is_rotation function.
!   at, at1, and sr are all in Cartesian coordinates.
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP), INTENT(IN) :: at(3,3), at1(3,3)
REAL(DP), INTENT(OUT) :: sr(3,3)

REAL(DP) :: bg(3,3)
INTEGER :: ipol, jpol

CALL recips(at1(1,1), at1(1,2), at1(1,3), bg(1,1), bg(1,2), bg(1,3)) 

DO ipol=1,3
   DO jpol=1,3
      sr(ipol,jpol)=at(ipol,1) * bg(jpol,1) +  &
                    at(ipol,2) * bg(jpol,2) +  &
                    at(ipol,3) * bg(jpol,3) 
   END DO
END DO

RETURN
END SUBROUTINE find_rotation

LOGICAL FUNCTION is_rotation(rmat)
!
!  This function receives a 3x3 real matrix  becomes .TRUE. if
!  it is a rotation.
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP), INTENT(IN) :: rmat(3,3)

REAL(DP) :: unity(3,3)
REAL(DP), PARAMETER :: eps1=1.D-8
INTEGER :: ipol, jpol

unity = matmul(rmat, TRANSPOSE(rmat))

is_rotation=.TRUE.
DO ipol=1,3
   is_rotation = is_rotation .AND. ABS(unity(ipol,ipol)-1.0_DP) < eps1
   DO jpol=ipol+1,3
      is_rotation = is_rotation .AND. ( ABS(unity(ipol,jpol)) < eps1 ) .AND. &
                                      ( ABS(unity(jpol,ipol)) < eps1 )
   END DO 
END DO


RETURN
END FUNCTION is_rotation

END MODULE rotate
