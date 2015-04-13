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
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  SAVE


  PUBLIC rotate_vect, rotate_tensors2

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

END MODULE rotate
