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
!  Presently it has the following:
!
!  rotate_vect: takes a set of vectors, a rotation matrix and gives the 
!               rotated vectors.  
!  rotate_tensors2: takes a set of rank two tensors, a rotation matrix
!               and gives the rotated rank two tensors.
!  rotate_tensors3: takes a set of rank three tensors, a rotation matrix
!               and gives the rotated rank three tensor.
!  rotate_tensors4: takes a set of rank four tensors, a rotation matrix
!               and gives the rotated rank four tensors.
!
!  rotate_rot : receive a rotation matrix that describe a rotation of
!               the reference system and a rotation matrix in the old
!               reference system. Gives the rotation in the rotated system
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
!  set_rot_xyz: set a rotation matrix of an angle phi about x, y, or z
!               It sets cos(phi), -sin(phi)
!                       sin(phi),  cos(phi)
!               on the axis that are not fix.
!
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC rotate_vect, rotate_tensors2, rotate_tensors3, rotate_tensors4,  &
         euler_to_su2, find_rotation, is_rotation, rotate_rot, set_rot_xyz

CONTAINS
!--------------------------------------------------------------------
   SUBROUTINE rotate_vect(rot, n, a, ra, flag)
!--------------------------------------------------------------------
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


!--------------------------------------------------------------------
   SUBROUTINE rotate_tensors2(rot, n, a, ra, flag)
!--------------------------------------------------------------------
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

!--------------------------------------------------------------------
   SUBROUTINE rotate_tensors3(rot, n, a, ra, flag)
!--------------------------------------------------------------------
!
!  this routine applies a rotation to n tensors of rank 4
!  if flag is 1 rotate the n tensors according to rot
!  if flag is -1 rotate the n tensors according to rot^-1
!

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   REAL(DP), INTENT(IN) :: rot(3,3), a(3,3,3,n)
   REAL(DP), INTENT(OUT) :: ra(3,3,3,n)
   INTEGER, INTENT(IN) :: flag

   INTEGER :: ivect, irot1, irot2, irot3, jrot1, jrot2, jrot3
  
   ra=0.0_DP
   DO ivect=1,n
      IF (flag==1) THEN
         DO irot1=1,3
          DO irot2=1,3
           DO irot3=1,3
            DO jrot1=1,3
             DO jrot2=1,3
              DO jrot3=1,3
                ra(irot1,irot2,irot3,ivect)=            &
                   ra(irot1,irot2,irot3,ivect)+         &
                       rot(irot1,jrot1)*rot(irot2,jrot2)*     &
                       rot(irot3,jrot3)*a(jrot1,jrot2,jrot3,ivect)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ELSE
         DO irot1=1,3
          DO irot2=1,3
           DO irot3=1,3
            DO jrot1=1,3
             DO jrot2=1,3
              DO jrot3=1,3
                 ra(irot1,irot2,irot3,ivect)=            &
                    ra(irot1,irot2,irot3,ivect)+         &
                        rot(jrot1,irot1)*rot(jrot2,irot2)*     &
                        rot(jrot3,irot3)*    &
                                a(jrot1,jrot2,jrot3,ivect)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
   END DO

   RETURN
   END SUBROUTINE rotate_tensors3

!--------------------------------------------------------------------
   SUBROUTINE rotate_tensors4(rot, n, a, ra, flag)
!--------------------------------------------------------------------
!
!  this routine applies a rotation to n tensors of rank 4
!  if flag is 1 rotate the n tensors according to rot
!  if flag is -1 rotate the n tensors according to rot^-1
!

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   REAL(DP), INTENT(IN) :: rot(3,3), a(3,3,3,3,n)
   REAL(DP), INTENT(OUT) :: ra(3,3,3,3,n)
   INTEGER, INTENT(IN) :: flag

   INTEGER :: ivect, irot1, irot2, irot3, irot4,   & 
                     jrot1, jrot2, jrot3, jrot4
  
   ra=0.0_DP
   DO ivect=1,n
      IF (flag==1) THEN
         DO irot1=1,3
          DO irot2=1,3
           DO irot3=1,3
            DO irot4=1,3
             DO jrot1=1,3
              DO jrot2=1,3
               DO jrot3=1,3
                DO jrot4=1,3
                  ra(irot1,irot2,irot3,irot4,ivect)=            &
                     ra(irot1,irot2,irot3,irot4,ivect)+         &
                         rot(irot1,jrot1)*rot(irot2,jrot2)*     &
                         rot(irot3,jrot3)*rot(irot4,jrot4)*     &
                                 a(jrot1,jrot2,jrot3,jrot4,ivect)
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ELSE
         DO irot1=1,3
          DO irot2=1,3
           DO irot3=1,3
            DO irot4=1,3
             DO jrot1=1,3
              DO jrot2=1,3
               DO jrot3=1,3
                DO jrot4=1,3
                  ra(irot1,irot2,irot3,irot4,ivect)=            &
                     ra(irot1,irot2,irot3,irot4,ivect)+         &
                         rot(jrot1,irot1)*rot(jrot2,irot2)*     &
                         rot(jrot3,irot3)*rot(jrot4,irot4)*     &
                                 a(jrot1,jrot2,jrot3,jrot4,ivect)
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
   END DO

   RETURN
   END SUBROUTINE rotate_tensors4

!--------------------------------------------------------------------
   SUBROUTINE rotate_rot(rot, rot_in, rot_out)
!--------------------------------------------------------------------

   USE kinds, ONLY : DP
   IMPLICIT NONE
   REAL(DP) :: rot(3,3), rot_in(3,3), rot_out(3,3)

   rot_out=MATMUL(TRANSPOSE(rot), MATMUL(rot_in, rot)) 

   RETURN
   END SUBROUTINE rotate_rot


!--------------------------------------------------------------------
SUBROUTINE euler_to_su2(psi,theta,phi,a,b)
!--------------------------------------------------------------------

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

!--------------------------------------------------------------------
SUBROUTINE find_rotation(at,at1,sr)
!--------------------------------------------------------------------
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

!--------------------------------------------------------------------
LOGICAL FUNCTION is_rotation(rmat)
!--------------------------------------------------------------------
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

!--------------------------------------------------------------------
SUBROUTINE set_rot_xyz(rot, axis, angle)
!--------------------------------------------------------------------
!
!  Set a rotation matrix about the axis x, y, or z according to axis (1,2,3).
!  Angle is the rotation angle and it is supposed to be in radiants
!
IMPLICIT NONE
REAL(DP), INTENT(INOUT) :: rot(3,3)
REAL(DP), INTENT(IN) :: angle
INTEGER, INTENT(IN) :: axis

rot=0.0_DP

IF (axis==1) THEN
   rot(1,1)=1.0_DP
   rot(2,2)=COS(angle)
   rot(2,3)=-SIN(angle)
   rot(3,2)=SIN(angle)
   rot(3,3)=COS(angle)
ELSEIF (axis==2) THEN
   rot(2,2)=1.0_DP
   rot(1,1)=COS(angle)
   rot(1,3)=-SIN(angle)
   rot(3,1)=SIN(angle)
   rot(3,3)=COS(angle)
ELSEIF (axis==3) THEN
   rot(3,3)=1.0_DP
   rot(1,1)=COS(angle)
   rot(1,2)=-SIN(angle)
   rot(2,1)=SIN(angle)
   rot(2,2)=COS(angle)
ELSE
   CALL errore('set_rot_xyz','axis not programmed',1)
ENDIF

RETURN
END SUBROUTINE set_rot_xyz

END MODULE rotate
