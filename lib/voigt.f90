!
! Copyright (C) 2017 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE voigt
!
!   this module contains the support routines for the transformation
!   of tensors in voigt notation
!   factors 1/2 must be dealt with by the calling routine depending
!   on the nature of the tensor
!


  USE kinds,     ONLY : DP
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC            &
       voigt_index, & ! given two indeces of a symmetric tensor gives the voigt
                      ! index 
       voigt_extract_indices, & ! given the voigt index gives the two
                      ! indices of a symmetric tensor. Only one couple 
                      ! is given.
       to_voigt2,   & ! transform a tensor of rank 2 to the voigt form
                      ! and viceversa
       to_voigt3,   & ! transform a tensor of rank 3 to the voigt form
                      ! and viceversa
       to_voigt4      ! transform a tensor of rank 4 to the voigt form 
                      ! and viceversa

CONTAINS

!-----------------------------------------------------------------------
SUBROUTINE voigt_index(m, n, mn)
!-----------------------------------------------------------------------
!
!  This routine receives two indeces 1<= m, n <=3 and
!  gives the voigt index 1<=mn<=6 corresponding to these two indices.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: m, n
INTEGER, INTENT(OUT) :: mn
INTEGER :: voigt(3,3)
DATA voigt / 1, 6, 5, 6, 2, 4, 5, 4, 3 / 

IF (m<1.OR.m>3.OR.n<1.OR.n>3) &
   CALL errore('voigt_index','m or n out or range',1)
mn=voigt(m,n) 

RETURN
END SUBROUTINE voigt_index

!-----------------------------------------------------------------------
SUBROUTINE voigt_extract_indices(m, n, mn)
!-----------------------------------------------------------------------
!
!  This routine receives the voigt index 1<=mn<=6 and gives
!  the two indices m, n, that correspond to it. It provide the couple
!  with n>=m.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: mn
INTEGER, INTENT(OUT) :: m, n
INTEGER :: mind(6), nind(6)
DATA mind  / 1, 2, 3, 2, 1, 1 /
DATA nind  / 1, 2, 3, 3, 3, 2 /

IF (mn<1.OR.mn>6) CALL errore('oigt_extract_indeces','mn out of range',1)
m=mind(mn)
n=nind(mn)

RETURN
END SUBROUTINE voigt_extract_indices

!-----------------------------------------------------------------------
SUBROUTINE to_voigt2 (av, a, flag)
!-----------------------------------------------------------------------
!
!  This routine transforms a rank 2 3x3 tensor in a 6 component Voigt array
!  (flag=.true.) or viceversa (flag=.false.)
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP), INTENT(INOUT) :: a(3,3)
REAL(DP), INTENT(INOUT) :: av(6)
LOGICAL, INTENT(IN) :: flag

INTEGER :: ij, i, j

IF (flag) THEN
   av=0.0_DP
   DO ij=1,6
      CALL voigt_extract_indices(i,j,ij)
      av(ij) = a(i,j) 
   ENDDO
ELSE
   a=0.0_DP
   DO i=1,3
      DO j=1,3
         CALL voigt_index(i,j,ij)
         a(i,j) = av(ij)
      ENDDO
   ENDDO
ENDIF

RETURN
END SUBROUTINE to_voigt2

!-----------------------------------------------------------------------
SUBROUTINE to_voigt3(av, a, flag)
!-----------------------------------------------------------------------
!
!  This routine transforms a rank 3 3x3x3 tensor in a 3x6 component Voigt array
!  (flag=.true.) or viceversa (flag=.false.)
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP), INTENT(INOUT) :: a(3,3,3)
REAL(DP), INTENT(INOUT) :: av(3,6)
LOGICAL, INTENT(IN) :: flag

INTEGER :: ij, i, j

IF (flag) THEN
   av=0.0_DP
   DO ij=1,6
      CALL voigt_extract_indices(i,j,ij)
      av(:,ij) = a(:,i,j) 
   ENDDO
ELSE
   a=0.0_DP
   DO i=1,3
      DO j=1,3
         CALL voigt_index(i,j,ij)
         a(:,i,j) = av(:,ij)
      ENDDO
   ENDDO
ENDIF

RETURN
END SUBROUTINE to_voigt3

!-----------------------------------------------------------------------
SUBROUTINE to_voigt4(av, a, flag)
!-----------------------------------------------------------------------
!
!  This routine transforms a rank 4 3x3x3x3 tensor in a 6x6 component Voigt 
!  array (flag=.true.) or viceversa (flag=.false.)
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP), INTENT(INOUT) :: a(3,3,3,3)
REAL(DP), INTENT(INOUT) :: av(6,6)
LOGICAL, INTENT(IN) :: flag

INTEGER :: ij, mn, i, j, m, n

IF (flag) THEN
   av=0.0_DP
   DO ij=1,6
      CALL voigt_extract_indices(i,j,ij)
      DO mn=1,6
         CALL voigt_extract_indices(m,n,mn)
         av(ij,mn) = a(i,j,m,n) 
      ENDDO
   ENDDO
ELSE
   a=0.0_DP
   DO i=1,3
      DO j=1,3
         CALL voigt_index(i,j,ij)
         DO m=1,3
            DO n=1,3
               CALL voigt_index(m,n,mn)
               a(i,j,m,n) = av(ij,mn)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDIF

RETURN
END SUBROUTINE to_voigt4

END MODULE voigt

