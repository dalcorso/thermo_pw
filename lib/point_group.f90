!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE point_group
!
!  This module contains variables and routines to deal with the point group
!  symmetry. It complements the routines that are in find_group.f90,
!  divide_class.f90 and divide_class_so.f90 in the PW/src directory of the
!  QE package.
!  The conventions, such as the code group, the symmetry operation types,
!  the irreducible representations etc. are the same.
!  Presently it has routines to perform the following task:
!  Given two point group, the second a subgroup of the first,
!  an a list of representations of the first point group, it
!  transform it in a list of representations of the second point group.
!  Double groups are supported.
!
!
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  SAVE


  PUBLIC convert_rap, find_aux_ind_two_groups, has_sigma_h, is_right_oriented

CONTAINS
  SUBROUTINE convert_rap(n, list_in, list_out, group_in, group_out, aux_ind, &
                         lspinorb)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: group_in, group_out, n, aux_ind
  INTEGER, INTENT(IN) :: list_in(n)
  INTEGER, INTENT(OUT) :: list_out(n)
  LOGICAL, INTENT(IN) :: lspinorb

  INTEGER :: i, j, ndeg
  LOGICAL :: done(n)
  INTEGER :: rap_list(4)

  done=.FALSE.
  DO i=1,n
     IF (done(i)) CYCLE
     IF (list_in(i)<=0) THEN
         list_out(i)=list_in(i)
         CYCLE
     ENDIF
     IF (lspinorb) THEN
        CALL convert_one_rap_so(list_in(i), ndeg, rap_list, group_in, group_out )
     ELSE
        CALL convert_one_rap(list_in(i), ndeg, rap_list, group_in, group_out, &
                                               aux_ind)
     ENDIF
     list_out(i) = rap_list(1)
     DO j=2, ndeg
        IF (list_in(i+j-1) /= list_in(i)) &
           CALL errore('conver_rap','Problem with degeneracy',1)
        list_out(i+j-1) = rap_list(j)
        done(i+j-1) = .TRUE.
     END DO
  END DO

  RETURN
  END SUBROUTINE convert_rap
 
  SUBROUTINE convert_one_rap(rap, ndeg, rap_list, group_in, group_out, aux_ind)
!
!  This routine sets the subduction table for the group subgroup relationship.
!  This subduction table is organized like this. It is set for each group_in,
!  for all the possibile subgroups. The first index is group_out, the
!  second index is aux_ind (there might be several possibilities for the
!  same subgroup, depending on which operations are selected), the
!  third index is rap, and the fourth index contains in the first position
!  the degeneracy of the rappresentation (1, 2, or 3) and in the three 
!  following positions the indices of the representation.
!  The first part of the routine set the information for all
!  the representations of group_in, and for all the possible aux_ind,
!  and the final instructions copy in ndeg and rap_list only the
!  information for the required representation.
!  The representations numbers are syncronized with those defined in
!  the routine set_irr (in PW/src/divide_class.f90).
!
!  This routine can be used in the scalar relativistic case. For the
!  fully relativistic case see the similar routine convert_one_rap_so
!
!
  USE io_global, ONLY : stdout

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: group_in, group_out, aux_ind, rap
  INTEGER, INTENT(OUT) :: ndeg, rap_list(3)
  INTEGER :: sub_table(32, 6, 12, 4)
  INTEGER :: ideg

  sub_table=0
!
!  C_1  has only representation 1
!
  sub_table(1,:,:,1)=1
  sub_table(1,:,:,2)=1

  SELECT CASE (group_in) 
!
!    
!
     CASE(1,2,3,4,5)

     CASE(6)
!
!  C_4
!
        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2

     CASE(7)
!
!  C_6
!
        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=2
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2
        sub_table(4,1,5,1)=1
        sub_table(4,1,5,2)=1
        sub_table(4,1,6,1)=1
        sub_table(4,1,6,2)=1

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=1
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=2
        sub_table(5,1,4,1)=1
        sub_table(5,1,4,2)=3
        sub_table(5,1,5,1)=1
        sub_table(5,1,5,2)=2
        sub_table(5,1,6,1)=1
        sub_table(5,1,6,2)=3

     CASE(8)
!
! D_2
!
        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=1
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=2
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=2

        sub_table(3,2,1,1)=1
        sub_table(3,2,1,2)=1
        sub_table(3,2,2,1)=1
        sub_table(3,2,2,2)=2
        sub_table(3,2,3,1)=1
        sub_table(3,2,3,2)=1
        sub_table(3,2,4,1)=1
        sub_table(3,2,4,2)=2

        sub_table(3,3,1,1)=1
        sub_table(3,3,1,2)=1
        sub_table(3,3,2,1)=1
        sub_table(3,3,2,2)=2
        sub_table(3,3,3,1)=1
        sub_table(3,3,3,2)=2
        sub_table(3,3,4,1)=1
        sub_table(3,3,4,2)=1

     CASE(9)
!
! D_3
!
        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=2
        sub_table(4,1,3,1)=2
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=1
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=2
        sub_table(5,1,3,3)=3

     CASE(10)
!
!  D_4
!
        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=1
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=1
        sub_table(4,1,5,1)=2
        sub_table(4,1,5,2)=2
        sub_table(4,1,5,3)=2

        sub_table(4,2,1,1)=1
        sub_table(4,2,1,2)=1
        sub_table(4,2,2,1)=1
        sub_table(4,2,2,2)=2
        sub_table(4,2,3,1)=1
        sub_table(4,2,3,2)=1
        sub_table(4,2,4,1)=1
        sub_table(4,2,4,2)=2
        sub_table(4,2,5,1)=2
        sub_table(4,2,5,2)=1
        sub_table(4,2,5,3)=2

        sub_table(4,3,1,1)=1
        sub_table(4,3,1,2)=1
        sub_table(4,3,2,1)=1
        sub_table(4,3,2,2)=2
        sub_table(4,3,3,1)=1
        sub_table(4,3,3,2)=2
        sub_table(4,3,4,1)=1
        sub_table(4,3,4,2)=1
        sub_table(4,3,5,1)=2
        sub_table(4,3,5,2)=1
        sub_table(4,3,5,3)=2

        sub_table(6,1,1,1)=1
        sub_table(6,1,1,2)=1
        sub_table(6,1,2,1)=1
        sub_table(6,1,2,2)=1
        sub_table(6,1,3,1)=1
        sub_table(6,1,3,2)=2
        sub_table(6,1,4,1)=1
        sub_table(6,1,4,2)=2
        sub_table(6,1,5,1)=2
        sub_table(6,1,5,2)=3
        sub_table(6,1,5,3)=4

        sub_table(8,1,1,1)=1
        sub_table(8,1,1,2)=1
        sub_table(8,1,2,1)=1
        sub_table(8,1,2,2)=2
        sub_table(8,1,3,1)=1
        sub_table(8,1,3,2)=1
        sub_table(8,1,4,1)=1
        sub_table(8,1,4,2)=2
        sub_table(8,1,5,1)=2
        sub_table(8,1,5,2)=3
        sub_table(8,1,5,3)=4

        sub_table(8,2,1,1)=1
        sub_table(8,2,1,2)=1
        sub_table(8,2,2,1)=1
        sub_table(8,2,2,2)=2
        sub_table(8,2,3,1)=1
        sub_table(8,2,3,2)=2
        sub_table(8,2,4,1)=1
        sub_table(8,2,4,2)=1
        sub_table(8,2,5,1)=2
        sub_table(8,2,5,2)=3
        sub_table(8,2,5,3)=4

     CASE(11)
!
!  D_6
!
        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2
        sub_table(4,1,5,1)=2
        sub_table(4,1,5,2)=1
        sub_table(4,1,5,3)=1
        sub_table(4,1,6,1)=2
        sub_table(4,1,6,2)=2
        sub_table(4,1,6,3)=2

        sub_table(4,2,1,1)=1
        sub_table(4,2,1,2)=1
        sub_table(4,2,2,1)=1
        sub_table(4,2,2,2)=2
        sub_table(4,2,3,1)=1
        sub_table(4,2,3,2)=1
        sub_table(4,2,4,1)=1
        sub_table(4,2,4,2)=2
        sub_table(4,2,5,1)=2
        sub_table(4,2,5,2)=1
        sub_table(4,2,5,3)=2
        sub_table(4,2,6,1)=2
        sub_table(4,2,6,2)=1
        sub_table(4,2,6,3)=2

        sub_table(4,3,1,1)=1
        sub_table(4,3,1,2)=1
        sub_table(4,3,2,1)=1
        sub_table(4,3,2,2)=2
        sub_table(4,3,3,1)=1
        sub_table(4,3,3,2)=2
        sub_table(4,3,4,1)=1
        sub_table(4,3,4,2)=1
        sub_table(4,3,5,1)=2
        sub_table(4,3,5,2)=1
        sub_table(4,3,5,3)=2
        sub_table(4,3,6,1)=2
        sub_table(4,3,6,2)=1
        sub_table(4,3,6,3)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=1
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=1
        sub_table(5,1,4,1)=1
        sub_table(5,1,4,2)=1
        sub_table(5,1,5,1)=2
        sub_table(5,1,5,2)=2
        sub_table(5,1,5,3)=3
        sub_table(5,1,6,1)=2
        sub_table(5,1,6,2)=2
        sub_table(5,1,6,3)=3

        sub_table(7,1,1,1)=1
        sub_table(7,1,1,2)=1
        sub_table(7,1,2,1)=1
        sub_table(7,1,2,2)=1
        sub_table(7,1,3,1)=1
        sub_table(7,1,3,2)=2
        sub_table(7,1,4,1)=1
        sub_table(7,1,4,2)=2
        sub_table(7,1,5,1)=2
        sub_table(7,1,5,2)=3
        sub_table(7,1,5,3)=4
        sub_table(7,1,6,1)=2
        sub_table(7,1,6,2)=5
        sub_table(7,1,6,3)=6

        sub_table(9,1,1,1)=1
        sub_table(9,1,1,2)=1
        sub_table(9,1,2,1)=1
        sub_table(9,1,2,2)=2
        sub_table(9,1,3,1)=1
        sub_table(9,1,3,2)=1
        sub_table(9,1,4,1)=1
        sub_table(9,1,4,2)=2
        sub_table(9,1,5,1)=2
        sub_table(9,1,5,2)=3
        sub_table(9,1,5,3)=3
        sub_table(9,1,6,1)=2
        sub_table(9,1,6,2)=3
        sub_table(9,1,6,3)=3

        sub_table(9,2,1,1)=1
        sub_table(9,2,1,2)=1
        sub_table(9,2,2,1)=1
        sub_table(9,2,2,2)=2
        sub_table(9,2,3,1)=1
        sub_table(9,2,3,2)=2
        sub_table(9,2,4,1)=1
        sub_table(9,2,4,2)=1
        sub_table(9,2,5,1)=2
        sub_table(9,2,5,2)=3
        sub_table(9,2,5,3)=3
        sub_table(9,2,6,1)=2
        sub_table(9,2,6,2)=3
        sub_table(9,2,6,3)=3

     CASE(12)
!
!   C_2v
!
        !
        !  C_s
        !
        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=1
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=2
        
        sub_table(3,2,1,1)=1
        sub_table(3,2,1,2)=1
        sub_table(3,2,2,1)=1
        sub_table(3,2,2,2)=2
        sub_table(3,2,3,1)=1
        sub_table(3,2,3,2)=2
        sub_table(3,2,4,1)=1
        sub_table(3,2,4,2)=1
        !
        !   C_2
        !
        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=2
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2

     CASE(13)
!
!   C_3v
!
        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=2
        sub_table(3,1,3,2)=1
        sub_table(3,1,3,3)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=1
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=2
        sub_table(5,1,3,3)=2


     CASE(14)
!
!   C_4v
!
        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=1
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=2
        sub_table(3,1,5,1)=2
        sub_table(3,1,5,2)=1
        sub_table(3,1,5,3)=2

        sub_table(3,2,1,1)=1
        sub_table(3,2,1,2)=1
        sub_table(3,2,2,1)=1
        sub_table(3,2,2,2)=2
        sub_table(3,2,3,1)=1
        sub_table(3,2,3,2)=2
        sub_table(3,2,4,1)=1
        sub_table(3,2,4,2)=1
        sub_table(3,2,5,1)=2
        sub_table(3,2,5,2)=1
        sub_table(3,2,5,3)=2

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=1
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=1
        sub_table(4,1,5,1)=2
        sub_table(4,1,5,2)=2
        sub_table(4,1,5,3)=2

        sub_table(6,1,1,1)=1
        sub_table(6,1,1,2)=1
        sub_table(6,1,2,1)=1
        sub_table(6,1,2,2)=1
        sub_table(6,1,3,1)=1
        sub_table(6,1,3,2)=2
        sub_table(6,1,4,1)=1
        sub_table(6,1,4,2)=2
        sub_table(6,1,5,1)=2
        sub_table(6,1,5,2)=3
        sub_table(6,1,5,3)=4

        sub_table(12,1,1,1)=1
        sub_table(12,1,1,2)=1
        sub_table(12,1,2,1)=1
        sub_table(12,1,2,2)=2
        sub_table(12,1,3,1)=1
        sub_table(12,1,3,2)=1
        sub_table(12,1,4,1)=1
        sub_table(12,1,4,2)=2
        sub_table(12,1,5,1)=2
        sub_table(12,1,5,2)=3
        sub_table(12,1,5,3)=4

        sub_table(12,2,1,1)=1
        sub_table(12,2,1,2)=1
        sub_table(12,2,2,1)=1
        sub_table(12,2,2,2)=2
        sub_table(12,2,3,1)=1
        sub_table(12,2,3,2)=2
        sub_table(12,2,4,1)=1
        sub_table(12,2,4,2)=1
        sub_table(12,2,5,1)=2
        sub_table(12,2,5,2)=3
        sub_table(12,2,5,3)=4

     CASE(15)
!
!   C_6v
!
        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=1
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=2
        sub_table(3,1,5,1)=2
        sub_table(3,1,5,2)=1
        sub_table(3,1,5,3)=2
        sub_table(3,1,6,1)=2
        sub_table(3,1,6,2)=1
        sub_table(3,1,6,3)=2

        sub_table(3,2,1,1)=1
        sub_table(3,2,1,2)=1
        sub_table(3,2,2,1)=1
        sub_table(3,2,2,2)=2
        sub_table(3,2,3,1)=1
        sub_table(3,2,3,2)=2
        sub_table(3,2,4,1)=1
        sub_table(3,2,4,2)=1
        sub_table(3,2,5,1)=2
        sub_table(3,2,5,2)=1
        sub_table(3,2,5,3)=2
        sub_table(3,2,6,1)=2
        sub_table(3,2,6,2)=1
        sub_table(3,2,6,3)=2

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2
        sub_table(4,1,5,1)=2
        sub_table(4,1,5,2)=1
        sub_table(4,1,5,3)=1
        sub_table(4,1,6,1)=2
        sub_table(4,1,6,2)=2
        sub_table(4,1,6,3)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=1
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=1
        sub_table(5,1,4,1)=1
        sub_table(5,1,4,2)=1
        sub_table(5,1,5,1)=2
        sub_table(5,1,5,2)=2
        sub_table(5,1,5,3)=3
        sub_table(5,1,6,1)=2
        sub_table(5,1,6,2)=2
        sub_table(5,1,6,3)=3

        sub_table(7,1,1,1)=1
        sub_table(7,1,1,2)=1
        sub_table(7,1,2,1)=1
        sub_table(7,1,2,2)=1
        sub_table(7,1,3,1)=1
        sub_table(7,1,3,2)=2
        sub_table(7,1,4,1)=1
        sub_table(7,1,4,2)=2
        sub_table(7,1,5,1)=2
        sub_table(7,1,5,2)=3
        sub_table(7,1,5,3)=4
        sub_table(7,1,6,1)=2
        sub_table(7,1,6,2)=5
        sub_table(7,1,6,3)=6

        sub_table(12,1,1,1)=1
        sub_table(12,1,1,2)=1
        sub_table(12,1,2,1)=1
        sub_table(12,1,2,2)=2
        sub_table(12,1,3,1)=1
        sub_table(12,1,3,2)=3
        sub_table(12,1,4,1)=1
        sub_table(12,1,4,2)=4
        sub_table(12,1,5,1)=2
        sub_table(12,1,5,2)=3
        sub_table(12,1,5,3)=4
        sub_table(12,1,6,1)=2
        sub_table(12,1,6,2)=1
        sub_table(12,1,6,3)=2

        sub_table(13,1,1,1)=1
        sub_table(13,1,1,2)=1
        sub_table(13,1,2,1)=1
        sub_table(13,1,2,2)=2
        sub_table(13,1,3,1)=1
        sub_table(13,1,3,2)=1
        sub_table(13,1,4,1)=1
        sub_table(13,1,4,2)=2
        sub_table(13,1,5,1)=2
        sub_table(13,1,5,2)=3
        sub_table(13,1,5,3)=3
        sub_table(13,1,6,1)=2
        sub_table(13,1,6,2)=3
        sub_table(13,1,6,3)=3

        sub_table(13,2,1,1)=1
        sub_table(13,2,1,2)=1
        sub_table(13,2,2,1)=1
        sub_table(13,2,2,2)=2
        sub_table(13,2,3,1)=1
        sub_table(13,2,3,2)=2
        sub_table(13,2,4,1)=1
        sub_table(13,2,4,2)=1
        sub_table(13,2,5,1)=2
        sub_table(13,2,5,2)=3
        sub_table(13,2,5,3)=3
        sub_table(13,2,6,1)=2
        sub_table(13,2,6,2)=3
        sub_table(13,2,6,3)=3

     CASE(16)
!
! C_2h
!
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=1
        sub_table(2,1,3,2)=2
        sub_table(2,1,4,1)=1
        sub_table(2,1,4,2)=2

        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=2
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=1

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=2
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=1
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2

     CASE(17)
!
! C_3h
!
        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=1
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=1
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=2
        sub_table(3,1,5,1)=1
        sub_table(3,1,5,2)=2
        sub_table(3,1,6,1)=1
        sub_table(3,1,6,2)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=2
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=3
        sub_table(5,1,4,1)=1
        sub_table(5,1,4,2)=1
        sub_table(5,1,5,1)=1
        sub_table(5,1,5,2)=2
        sub_table(5,1,6,1)=1
        sub_table(5,1,6,2)=3

     CASE(18)
!
!   C_4h
!
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=1
        sub_table(2,1,3,2)=1
        sub_table(2,1,4,1)=1
        sub_table(2,1,4,2)=1
        sub_table(2,1,5,1)=1
        sub_table(2,1,5,2)=2
        sub_table(2,1,6,1)=1
        sub_table(2,1,6,2)=2
        sub_table(2,1,7,1)=1
        sub_table(2,1,7,2)=2
        sub_table(2,1,8,1)=1
        sub_table(2,1,8,2)=2

        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=1
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=2
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=2
        sub_table(3,1,5,1)=1
        sub_table(3,1,5,2)=2
        sub_table(3,1,6,1)=1
        sub_table(3,1,6,2)=2
        sub_table(3,1,7,1)=1
        sub_table(3,1,7,2)=1
        sub_table(3,1,8,1)=1
        sub_table(3,1,8,2)=1

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2
        sub_table(4,1,5,1)=1
        sub_table(4,1,5,2)=1
        sub_table(4,1,6,1)=1
        sub_table(4,1,6,2)=1
        sub_table(4,1,7,1)=1
        sub_table(4,1,7,2)=2
        sub_table(4,1,8,1)=1
        sub_table(4,1,8,2)=2

        sub_table(6,1,1,1)=1
        sub_table(6,1,1,2)=1
        sub_table(6,1,2,1)=1
        sub_table(6,1,2,2)=2
        sub_table(6,1,3,1)=1
        sub_table(6,1,3,2)=3
        sub_table(6,1,4,1)=1
        sub_table(6,1,4,2)=4
        sub_table(6,1,5,1)=1
        sub_table(6,1,5,2)=1
        sub_table(6,1,6,1)=1
        sub_table(6,1,6,2)=2
        sub_table(6,1,7,1)=1
        sub_table(6,1,7,2)=3
        sub_table(6,1,8,1)=1
        sub_table(6,1,8,2)=4

        sub_table(16,1,1,1)=1
        sub_table(16,1,1,2)=1
        sub_table(16,1,2,1)=1
        sub_table(16,1,2,2)=1
        sub_table(16,1,3,1)=1
        sub_table(16,1,3,2)=2
        sub_table(16,1,4,1)=1
        sub_table(16,1,4,2)=2
        sub_table(16,1,5,1)=1
        sub_table(16,1,5,2)=3
        sub_table(16,1,6,1)=1
        sub_table(16,1,6,2)=3
        sub_table(16,1,7,1)=1
        sub_table(16,1,7,2)=4
        sub_table(16,1,8,1)=1
        sub_table(16,1,8,2)=4

        sub_table(26,1,1,1)=1
        sub_table(26,1,1,2)=1
        sub_table(26,1,2,1)=1
        sub_table(26,1,2,2)=2
        sub_table(26,1,3,1)=1
        sub_table(26,1,3,2)=3
        sub_table(26,1,4,1)=1
        sub_table(26,1,4,2)=4
        sub_table(26,1,5,1)=1
        sub_table(26,1,5,2)=2
        sub_table(26,1,6,1)=1
        sub_table(26,1,6,2)=1
        sub_table(26,1,7,1)=1
        sub_table(26,1,7,2)=4
        sub_table(26,1,8,1)=1
        sub_table(26,1,8,2)=3

     CASE(19)
!
!  C_6h
!
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=1
        sub_table(2,1,3,2)=1
        sub_table(2,1,4,1)=1
        sub_table(2,1,4,2)=1
        sub_table(2,1,5,1)=1
        sub_table(2,1,5,2)=1
        sub_table(2,1,6,1)=1
        sub_table(2,1,6,2)=1
        sub_table(2,1,7,1)=1
        sub_table(2,1,7,2)=2
        sub_table(2,1,8,1)=1
        sub_table(2,1,8,2)=2
        sub_table(2,1,9,1)=1
        sub_table(2,1,9,2)=2
        sub_table(2,1,10,1)=1
        sub_table(2,1,10,2)=2
        sub_table(2,1,11,1)=1
        sub_table(2,1,11,2)=2
        sub_table(2,1,12,1)=1
        sub_table(2,1,12,2)=2

        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=2
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=2
        sub_table(3,1,5,1)=1
        sub_table(3,1,5,2)=1
        sub_table(3,1,6,1)=1
        sub_table(3,1,6,2)=1
        sub_table(3,1,7,1)=1
        sub_table(3,1,7,2)=2
        sub_table(3,1,8,1)=1
        sub_table(3,1,8,2)=1
        sub_table(3,1,9,1)=1
        sub_table(3,1,9,2)=1
        sub_table(3,1,10,1)=1
        sub_table(3,1,10,2)=1
        sub_table(3,1,11,1)=1
        sub_table(3,1,11,2)=2
        sub_table(3,1,12,1)=1
        sub_table(3,1,12,2)=2

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=2
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2
        sub_table(4,1,5,1)=1
        sub_table(4,1,5,2)=1
        sub_table(4,1,6,1)=1
        sub_table(4,1,6,2)=1
        sub_table(4,1,7,1)=1
        sub_table(4,1,7,2)=1
        sub_table(4,1,8,1)=1
        sub_table(4,1,8,2)=2
        sub_table(4,1,9,1)=1
        sub_table(4,1,9,2)=2
        sub_table(4,1,10,1)=1
        sub_table(4,1,10,2)=2
        sub_table(4,1,11,1)=1
        sub_table(4,1,11,2)=1
        sub_table(4,1,12,1)=1
        sub_table(4,1,12,2)=1

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=1
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=2
        sub_table(5,1,4,1)=1
        sub_table(5,1,4,2)=3
        sub_table(5,1,5,1)=1
        sub_table(5,1,5,2)=2
        sub_table(5,1,6,1)=1
        sub_table(5,1,6,2)=3
        sub_table(5,1,7,1)=1
        sub_table(5,1,7,2)=1
        sub_table(5,1,8,1)=1
        sub_table(5,1,8,2)=1
        sub_table(5,1,9,1)=1
        sub_table(5,1,9,2)=2
        sub_table(5,1,10,1)=1
        sub_table(5,1,10,2)=3
        sub_table(5,1,11,1)=1
        sub_table(5,1,11,2)=2
        sub_table(5,1,12,1)=1
        sub_table(5,1,12,2)=3

        sub_table(7,1,1,1)=1
        sub_table(7,1,1,2)=1
        sub_table(7,1,2,1)=1
        sub_table(7,1,2,2)=2
        sub_table(7,1,3,1)=1
        sub_table(7,1,3,2)=3
        sub_table(7,1,4,1)=1
        sub_table(7,1,4,2)=4
        sub_table(7,1,5,1)=1
        sub_table(7,1,5,2)=5
        sub_table(7,1,6,1)=1
        sub_table(7,1,6,2)=6
        sub_table(7,1,7,1)=1
        sub_table(7,1,7,2)=1
        sub_table(7,1,8,1)=1
        sub_table(7,1,8,2)=2
        sub_table(7,1,9,1)=1
        sub_table(7,1,9,2)=3
        sub_table(7,1,10,1)=1
        sub_table(7,1,10,2)=4
        sub_table(7,1,11,1)=1
        sub_table(7,1,11,2)=5
        sub_table(7,1,12,1)=1
        sub_table(7,1,12,2)=6

        sub_table(16,1,1,1)=1
        sub_table(16,1,1,2)=1
        sub_table(16,1,2,1)=1
        sub_table(16,1,2,2)=2
        sub_table(16,1,3,1)=1
        sub_table(16,1,3,2)=2
        sub_table(16,1,4,1)=1
        sub_table(16,1,4,2)=2
        sub_table(16,1,5,1)=1
        sub_table(16,1,5,2)=1
        sub_table(16,1,6,1)=1
        sub_table(16,1,6,2)=1
        sub_table(16,1,7,1)=1
        sub_table(16,1,7,2)=3
        sub_table(16,1,8,1)=1
        sub_table(16,1,8,2)=4
        sub_table(16,1,9,1)=1
        sub_table(16,1,9,2)=4
        sub_table(16,1,10,1)=1
        sub_table(16,1,10,2)=4
        sub_table(16,1,11,1)=1
        sub_table(16,1,11,2)=3
        sub_table(16,1,12,1)=1
        sub_table(16,1,12,2)=3

        sub_table(17,1,1,1)=1
        sub_table(17,1,1,2)=1
        sub_table(17,1,2,1)=1
        sub_table(17,1,2,2)=4
        sub_table(17,1,3,1)=1
        sub_table(17,1,3,2)=5
        sub_table(17,1,4,1)=1
        sub_table(17,1,4,2)=6
        sub_table(17,1,5,1)=1
        sub_table(17,1,5,2)=2
        sub_table(17,1,6,1)=1
        sub_table(17,1,6,2)=3
        sub_table(17,1,7,1)=1
        sub_table(17,1,7,2)=4
        sub_table(17,1,8,1)=1
        sub_table(17,1,8,2)=1
        sub_table(17,1,9,1)=1
        sub_table(17,1,9,2)=2
        sub_table(17,1,10,1)=1
        sub_table(17,1,10,2)=3
        sub_table(17,1,11,1)=1
        sub_table(17,1,11,2)=5
        sub_table(17,1,12,1)=1
        sub_table(17,1,12,2)=6

        sub_table(27,1,1,1)=1
        sub_table(27,1,1,2)=1
        sub_table(27,1,2,1)=1
        sub_table(27,1,2,2)=1
        sub_table(27,1,3,1)=1
        sub_table(27,1,3,2)=2
        sub_table(27,1,4,1)=1
        sub_table(27,1,4,2)=3
        sub_table(27,1,5,1)=1
        sub_table(27,1,5,2)=2
        sub_table(27,1,6,1)=1
        sub_table(27,1,6,2)=3
        sub_table(27,1,7,1)=1
        sub_table(27,1,7,2)=4
        sub_table(27,1,8,1)=1
        sub_table(27,1,8,2)=4
        sub_table(27,1,9,1)=1
        sub_table(27,1,9,2)=5
        sub_table(27,1,10,1)=1
        sub_table(27,1,10,2)=6
        sub_table(27,1,11,1)=1
        sub_table(27,1,11,2)=5
        sub_table(27,1,12,1)=1
        sub_table(27,1,12,2)=6
 
     CASE(20)
!
!  D_2h
!
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=1
        sub_table(2,1,3,2)=1
        sub_table(2,1,4,1)=1
        sub_table(2,1,4,2)=1
        sub_table(2,1,5,1)=1
        sub_table(2,1,5,2)=2
        sub_table(2,1,6,1)=1
        sub_table(2,1,6,2)=2
        sub_table(2,1,7,1)=1
        sub_table(2,1,7,2)=2
        sub_table(2,1,8,1)=1
        sub_table(2,1,8,2)=2

        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=2
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=1
        sub_table(3,1,5,1)=1
        sub_table(3,1,5,2)=2
        sub_table(3,1,6,1)=1
        sub_table(3,1,6,2)=1
        sub_table(3,1,7,1)=1
        sub_table(3,1,7,2)=1
        sub_table(3,1,8,1)=1
        sub_table(3,1,8,2)=2

        sub_table(3,2,1,1)=1
        sub_table(3,2,1,2)=1
        sub_table(3,2,2,1)=1
        sub_table(3,2,2,2)=2
        sub_table(3,2,3,1)=1
        sub_table(3,2,3,2)=1
        sub_table(3,2,4,1)=1
        sub_table(3,2,4,2)=2
        sub_table(3,2,5,1)=1
        sub_table(3,2,5,2)=2
        sub_table(3,2,6,1)=1
        sub_table(3,2,6,2)=1
        sub_table(3,2,7,1)=1
        sub_table(3,2,7,2)=2
        sub_table(3,2,8,1)=1
        sub_table(3,2,8,2)=1

        sub_table(3,3,1,1)=1
        sub_table(3,3,1,2)=1
        sub_table(3,3,2,1)=1
        sub_table(3,3,2,2)=1
        sub_table(3,3,3,1)=1
        sub_table(3,3,3,2)=2
        sub_table(3,3,4,1)=1
        sub_table(3,3,4,2)=2
        sub_table(3,3,5,1)=1
        sub_table(3,3,5,2)=2
        sub_table(3,3,6,1)=1
        sub_table(3,3,6,2)=2
        sub_table(3,3,7,1)=1
        sub_table(3,3,7,2)=1
        sub_table(3,3,8,1)=1
        sub_table(3,3,8,2)=1

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=2
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=1
        sub_table(4,1,5,1)=1
        sub_table(4,1,5,2)=1
        sub_table(4,1,6,1)=1
        sub_table(4,1,6,2)=2
        sub_table(4,1,7,1)=1
        sub_table(4,1,7,2)=2
        sub_table(4,1,8,1)=1
        sub_table(4,1,8,2)=1

        sub_table(4,2,1,1)=1
        sub_table(4,2,1,2)=1
        sub_table(4,2,2,1)=1
        sub_table(4,2,2,2)=2
        sub_table(4,2,3,1)=1
        sub_table(4,2,3,2)=1
        sub_table(4,2,4,1)=1
        sub_table(4,2,4,2)=2
        sub_table(4,2,5,1)=1
        sub_table(4,2,5,2)=1
        sub_table(4,2,6,1)=1
        sub_table(4,2,6,2)=2
        sub_table(4,2,7,1)=1
        sub_table(4,2,7,2)=1
        sub_table(4,2,8,1)=1
        sub_table(4,2,8,2)=2

        sub_table(4,3,1,1)=1
        sub_table(4,3,1,2)=1
        sub_table(4,3,2,1)=1
        sub_table(4,3,2,2)=1
        sub_table(4,3,3,1)=1
        sub_table(4,3,3,2)=2
        sub_table(4,3,4,1)=1
        sub_table(4,3,4,2)=2
        sub_table(4,3,5,1)=1
        sub_table(4,3,5,2)=1
        sub_table(4,3,6,1)=1
        sub_table(4,3,6,2)=1
        sub_table(4,3,7,1)=1
        sub_table(4,3,7,2)=2
        sub_table(4,3,8,1)=1
        sub_table(4,3,8,2)=2

        sub_table(8,1,1,1)=1
        sub_table(8,1,1,2)=1
        sub_table(8,1,2,1)=1
        sub_table(8,1,2,2)=2
        sub_table(8,1,3,1)=1
        sub_table(8,1,3,2)=3
        sub_table(8,1,4,1)=1
        sub_table(8,1,4,2)=4
        sub_table(8,1,5,1)=1
        sub_table(8,1,5,2)=1
        sub_table(8,1,6,1)=1
        sub_table(8,1,6,2)=2
        sub_table(8,1,7,1)=1
        sub_table(8,1,7,2)=3
        sub_table(8,1,8,1)=1
        sub_table(8,1,8,2)=4

        sub_table(12,1,1,1)=1
        sub_table(12,1,1,2)=1
        sub_table(12,1,2,1)=1
        sub_table(12,1,2,2)=2
        sub_table(12,1,3,1)=1
        sub_table(12,1,3,2)=3
        sub_table(12,1,4,1)=1
        sub_table(12,1,4,2)=4
        sub_table(12,1,5,1)=1
        sub_table(12,1,5,2)=2
        sub_table(12,1,6,1)=1
        sub_table(12,1,6,2)=1
        sub_table(12,1,7,1)=1
        sub_table(12,1,7,2)=4
        sub_table(12,1,8,1)=1
        sub_table(12,1,8,2)=3

        sub_table(12,2,1,1)=1
        sub_table(12,2,1,2)=1
        sub_table(12,2,2,1)=1
        sub_table(12,2,2,2)=3
        sub_table(12,2,3,1)=1
        sub_table(12,2,3,2)=2
        sub_table(12,2,4,1)=1
        sub_table(12,2,4,2)=4
        sub_table(12,2,5,1)=1
        sub_table(12,2,5,2)=2
        sub_table(12,2,6,1)=1
        sub_table(12,2,6,2)=4
        sub_table(12,2,7,1)=1
        sub_table(12,2,7,2)=1
        sub_table(12,2,8,1)=1
        sub_table(12,2,8,2)=3

        sub_table(12,3,1,1)=1
        sub_table(12,3,1,2)=1
        sub_table(12,3,2,1)=1
        sub_table(12,3,2,2)=3
        sub_table(12,3,3,1)=1
        sub_table(12,3,3,2)=4
        sub_table(12,3,4,1)=1
        sub_table(12,3,4,2)=2
        sub_table(12,3,5,1)=1
        sub_table(12,3,5,2)=2
        sub_table(12,3,6,1)=1
        sub_table(12,3,6,2)=4
        sub_table(12,3,7,1)=1
        sub_table(12,3,7,2)=3
        sub_table(12,3,8,1)=1
        sub_table(12,3,8,2)=1

        sub_table(16,1,1,1)=1
        sub_table(16,1,1,2)=1
        sub_table(16,1,2,1)=1
        sub_table(16,1,2,2)=2
        sub_table(16,1,3,1)=1
        sub_table(16,1,3,2)=2
        sub_table(16,1,4,1)=1
        sub_table(16,1,4,2)=1
        sub_table(16,1,5,1)=1
        sub_table(16,1,5,2)=3
        sub_table(16,1,6,1)=1
        sub_table(16,1,6,2)=4
        sub_table(16,1,7,1)=1
        sub_table(16,1,7,2)=4
        sub_table(16,1,8,1)=1
        sub_table(16,1,8,2)=3

        sub_table(16,2,1,1)=1
        sub_table(16,2,1,2)=1
        sub_table(16,2,2,1)=1
        sub_table(16,2,2,2)=2
        sub_table(16,2,3,1)=1
        sub_table(16,2,3,2)=1
        sub_table(16,2,4,1)=1
        sub_table(16,2,4,2)=2
        sub_table(16,2,5,1)=1
        sub_table(16,2,5,2)=3
        sub_table(16,2,6,1)=1
        sub_table(16,2,6,2)=4
        sub_table(16,2,7,1)=1
        sub_table(16,2,7,2)=3
        sub_table(16,2,8,1)=1
        sub_table(16,2,8,2)=4

        sub_table(16,3,1,1)=1
        sub_table(16,3,1,2)=1
        sub_table(16,3,2,1)=1
        sub_table(16,3,2,2)=1
        sub_table(16,3,3,1)=1
        sub_table(16,3,3,2)=2
        sub_table(16,3,4,1)=1
        sub_table(16,3,4,2)=2
        sub_table(16,3,5,1)=1
        sub_table(16,3,5,2)=3
        sub_table(16,3,6,1)=1
        sub_table(16,3,6,2)=3
        sub_table(16,3,7,1)=1
        sub_table(16,3,7,2)=4
        sub_table(16,3,8,1)=1
        sub_table(16,3,8,2)=4

     CASE(21)
!
!  D_3h
!
        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=1
        sub_table(3,1,3,1)=2
        sub_table(3,1,3,2)=1
        sub_table(3,1,3,3)=1
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=2
        sub_table(3,1,5,1)=1
        sub_table(3,1,5,2)=2
        sub_table(3,1,6,1)=2
        sub_table(3,1,6,2)=2
        sub_table(3,1,6,3)=2

        sub_table(3,2,1,1)=1
        sub_table(3,2,1,2)=1
        sub_table(3,2,2,1)=1
        sub_table(3,2,2,2)=2
        sub_table(3,2,3,1)=2
        sub_table(3,2,3,2)=1
        sub_table(3,2,3,3)=2
        sub_table(3,2,4,1)=1
        sub_table(3,2,4,2)=2
        sub_table(3,2,5,1)=1
        sub_table(3,2,5,2)=1
        sub_table(3,2,6,1)=2
        sub_table(3,2,6,2)=1
        sub_table(3,2,6,3)=2

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=2
        sub_table(4,1,3,1)=2
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=1
        sub_table(4,1,5,1)=1
        sub_table(4,1,5,2)=2
        sub_table(4,1,6,1)=2
        sub_table(4,1,6,2)=1
        sub_table(4,1,6,3)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=1
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=2
        sub_table(5,1,3,3)=3
        sub_table(5,1,4,1)=1
        sub_table(5,1,4,2)=1
        sub_table(5,1,5,1)=1
        sub_table(5,1,5,2)=1
        sub_table(5,1,6,1)=2
        sub_table(5,1,6,2)=2
        sub_table(5,1,6,3)=3

        sub_table(9,1,1,1)=1
        sub_table(9,1,1,2)=1
        sub_table(9,1,2,1)=1
        sub_table(9,1,2,2)=2
        sub_table(9,1,3,1)=2
        sub_table(9,1,3,2)=3
        sub_table(9,1,3,3)=3
        sub_table(9,1,4,1)=1
        sub_table(9,1,4,2)=1
        sub_table(9,1,5,1)=1
        sub_table(9,1,5,2)=2
        sub_table(9,1,6,1)=2
        sub_table(9,1,6,2)=3
        sub_table(9,1,6,3)=3

        sub_table(12,1,1,1)=1
        sub_table(12,1,1,2)=1
        sub_table(12,1,2,1)=1
        sub_table(12,1,2,2)=3
        sub_table(12,1,3,1)=2
        sub_table(12,1,3,2)=1
        sub_table(12,1,3,3)=3
        sub_table(12,1,4,1)=1
        sub_table(12,1,4,2)=2
        sub_table(12,1,5,1)=1
        sub_table(12,1,5,2)=4
        sub_table(12,1,6,1)=2
        sub_table(12,1,6,2)=2
        sub_table(12,1,6,3)=4

        sub_table(13,1,1,1)=1
        sub_table(13,1,1,2)=1
        sub_table(13,1,2,1)=1
        sub_table(13,1,2,2)=2
        sub_table(13,1,3,1)=2
        sub_table(13,1,3,2)=3
        sub_table(13,1,3,3)=3
        sub_table(13,1,4,1)=1
        sub_table(13,1,4,2)=2
        sub_table(13,1,5,1)=1
        sub_table(13,1,5,2)=1
        sub_table(13,1,6,1)=2
        sub_table(13,1,6,2)=3
        sub_table(13,1,6,3)=3

        sub_table(17,1,1,1)=1
        sub_table(17,1,1,2)=1
        sub_table(17,1,2,1)=1
        sub_table(17,1,2,2)=1
        sub_table(17,1,3,1)=2
        sub_table(17,1,3,2)=2
        sub_table(17,1,3,3)=3
        sub_table(17,1,4,1)=1
        sub_table(17,1,4,2)=4
        sub_table(17,1,5,1)=1
        sub_table(17,1,5,2)=4
        sub_table(17,1,6,1)=2
        sub_table(17,1,6,2)=5
        sub_table(17,1,6,3)=6
 
     CASE(22)
!
!  D_4h
!
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=1
        sub_table(2,1,3,2)=1
        sub_table(2,1,4,1)=1
        sub_table(2,1,4,2)=1
        sub_table(2,1,5,1)=2
        sub_table(2,1,5,2)=1
        sub_table(2,1,5,3)=1
        sub_table(2,1,6,1)=1
        sub_table(2,1,6,2)=2
        sub_table(2,1,7,1)=1
        sub_table(2,1,7,2)=2
        sub_table(2,1,8,1)=1
        sub_table(2,1,8,2)=2
        sub_table(2,1,9,1)=1
        sub_table(2,1,9,2)=2
        sub_table(2,1,10,1)=2
        sub_table(2,1,10,2)=2
        sub_table(2,1,10,3)=2

        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=1
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=1
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=1
        sub_table(3,1,5,1)=2
        sub_table(3,1,5,2)=2
        sub_table(3,1,5,3)=2
        sub_table(3,1,6,1)=1
        sub_table(3,1,6,2)=2
        sub_table(3,1,7,1)=1
        sub_table(3,1,7,2)=2
        sub_table(3,1,8,1)=1
        sub_table(3,1,8,2)=2
        sub_table(3,1,9,1)=1
        sub_table(3,1,9,2)=2
        sub_table(3,1,10,1)=2
        sub_table(3,1,10,2)=1
        sub_table(3,1,10,3)=1

        sub_table(3,2,1,1)=1
        sub_table(3,2,1,2)=1
        sub_table(3,2,2,1)=1
        sub_table(3,2,2,2)=2
        sub_table(3,2,3,1)=1
        sub_table(3,2,3,2)=1
        sub_table(3,2,4,1)=1
        sub_table(3,2,4,2)=2
        sub_table(3,2,5,1)=2
        sub_table(3,2,5,2)=1
        sub_table(3,2,5,3)=2
        sub_table(3,2,6,1)=1
        sub_table(3,2,6,2)=2
        sub_table(3,2,7,1)=1
        sub_table(3,2,7,2)=1
        sub_table(3,2,8,1)=1
        sub_table(3,2,8,2)=2
        sub_table(3,2,9,1)=1
        sub_table(3,2,9,2)=1
        sub_table(3,2,10,1)=2
        sub_table(3,2,10,2)=1
        sub_table(3,2,10,3)=2

        sub_table(3,3,1,1)=1
        sub_table(3,3,1,2)=1
        sub_table(3,3,2,1)=1
        sub_table(3,3,2,2)=2
        sub_table(3,3,3,1)=1
        sub_table(3,3,3,2)=2
        sub_table(3,3,4,1)=1
        sub_table(3,3,4,2)=1
        sub_table(3,3,5,1)=2
        sub_table(3,3,5,2)=1
        sub_table(3,3,5,3)=2
        sub_table(3,3,6,1)=1
        sub_table(3,3,6,2)=2
        sub_table(3,3,7,1)=1
        sub_table(3,3,7,2)=1
        sub_table(3,3,8,1)=1
        sub_table(3,3,8,2)=1
        sub_table(3,3,9,1)=1
        sub_table(3,3,9,2)=2
        sub_table(3,3,10,1)=2
        sub_table(3,3,10,2)=1
        sub_table(3,3,10,3)=2

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=1
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=1
        sub_table(4,1,5,1)=2
        sub_table(4,1,5,2)=2
        sub_table(4,1,5,3)=2
        sub_table(4,1,6,1)=1
        sub_table(4,1,6,2)=1
        sub_table(4,1,7,1)=1
        sub_table(4,1,7,2)=1
        sub_table(4,1,8,1)=1
        sub_table(4,1,8,2)=1
        sub_table(4,1,9,1)=1
        sub_table(4,1,9,2)=1
        sub_table(4,1,10,1)=2
        sub_table(4,1,10,2)=2
        sub_table(4,1,10,3)=2

        sub_table(4,2,1,1)=1
        sub_table(4,2,1,2)=1
        sub_table(4,2,2,1)=1
        sub_table(4,2,2,2)=2
        sub_table(4,2,3,1)=1
        sub_table(4,2,3,2)=1
        sub_table(4,2,4,1)=1
        sub_table(4,2,4,2)=2
        sub_table(4,2,5,1)=2
        sub_table(4,2,5,2)=1
        sub_table(4,2,5,3)=2
        sub_table(4,2,6,1)=1
        sub_table(4,2,6,2)=1
        sub_table(4,2,7,1)=1
        sub_table(4,2,7,2)=2
        sub_table(4,2,8,1)=1
        sub_table(4,2,8,2)=1
        sub_table(4,2,9,1)=1
        sub_table(4,2,9,2)=2
        sub_table(4,2,10,1)=2
        sub_table(4,2,10,2)=1
        sub_table(4,2,10,3)=2

        sub_table(4,3,1,1)=1
        sub_table(4,3,1,2)=1
        sub_table(4,3,2,1)=1
        sub_table(4,3,2,2)=2
        sub_table(4,3,3,1)=1
        sub_table(4,3,3,2)=2
        sub_table(4,3,4,1)=1
        sub_table(4,3,4,2)=1
        sub_table(4,3,5,1)=2
        sub_table(4,3,5,2)=1
        sub_table(4,3,5,3)=2
        sub_table(4,3,6,1)=1
        sub_table(4,3,6,2)=1
        sub_table(4,3,7,1)=1
        sub_table(4,3,7,2)=2
        sub_table(4,3,8,1)=1
        sub_table(4,3,8,2)=2
        sub_table(4,3,9,1)=1
        sub_table(4,3,9,2)=1
        sub_table(4,3,10,1)=2
        sub_table(4,3,10,2)=1
        sub_table(4,3,10,3)=2

        sub_table(6,1,1,1)=1
        sub_table(6,1,1,2)=1
        sub_table(6,1,2,1)=1
        sub_table(6,1,2,2)=1
        sub_table(6,1,3,1)=1
        sub_table(6,1,3,2)=2
        sub_table(6,1,4,1)=1
        sub_table(6,1,4,2)=2
        sub_table(6,1,5,1)=2
        sub_table(6,1,5,2)=3
        sub_table(6,1,5,3)=4
        sub_table(6,1,6,1)=1
        sub_table(6,1,6,2)=1
        sub_table(6,1,7,1)=1
        sub_table(6,1,7,2)=1
        sub_table(6,1,8,1)=1
        sub_table(6,1,8,2)=2
        sub_table(6,1,9,1)=1
        sub_table(6,1,9,2)=2
        sub_table(6,1,10,1)=2
        sub_table(6,1,10,2)=3
        sub_table(6,1,10,3)=4

        sub_table(8,1,1,1)=1
        sub_table(8,1,1,2)=1
        sub_table(8,1,2,1)=1
        sub_table(8,1,2,2)=2
        sub_table(8,1,3,1)=1
        sub_table(8,1,3,2)=1
        sub_table(8,1,4,1)=1
        sub_table(8,1,4,2)=2
        sub_table(8,1,5,1)=2
        sub_table(8,1,5,2)=3
        sub_table(8,1,5,3)=4
        sub_table(8,1,6,1)=1
        sub_table(8,1,6,2)=1
        sub_table(8,1,7,1)=1
        sub_table(8,1,7,2)=2
        sub_table(8,1,8,1)=1
        sub_table(8,1,8,2)=1
        sub_table(8,1,9,1)=1
        sub_table(8,1,9,2)=2
        sub_table(8,1,10,1)=2
        sub_table(8,1,10,2)=3
        sub_table(8,1,10,3)=4

        sub_table(8,2,1,1)=1
        sub_table(8,2,1,2)=1
        sub_table(8,2,2,1)=1
        sub_table(8,2,2,2)=2
        sub_table(8,2,3,1)=1
        sub_table(8,2,3,2)=2
        sub_table(8,2,4,1)=1
        sub_table(8,2,4,2)=1
        sub_table(8,2,5,1)=2
        sub_table(8,2,5,2)=3
        sub_table(8,2,5,3)=4
        sub_table(8,2,6,1)=1
        sub_table(8,2,6,2)=1
        sub_table(8,2,7,1)=1
        sub_table(8,2,7,2)=2
        sub_table(8,2,8,1)=1
        sub_table(8,2,8,2)=2
        sub_table(8,2,9,1)=1
        sub_table(8,2,9,2)=1
        sub_table(8,2,10,1)=2
        sub_table(8,2,10,2)=3
        sub_table(8,2,10,3)=4

        sub_table(10,1,1,1)=1
        sub_table(10,1,1,2)=1
        sub_table(10,1,2,1)=1
        sub_table(10,1,2,2)=2
        sub_table(10,1,3,1)=1
        sub_table(10,1,3,2)=3
        sub_table(10,1,4,1)=1
        sub_table(10,1,4,2)=4
        sub_table(10,1,5,1)=2
        sub_table(10,1,5,2)=5
        sub_table(10,1,5,3)=5
        sub_table(10,1,6,1)=1
        sub_table(10,1,6,2)=1
        sub_table(10,1,7,1)=1
        sub_table(10,1,7,2)=2
        sub_table(10,1,8,1)=1
        sub_table(10,1,8,2)=3
        sub_table(10,1,9,1)=1
        sub_table(10,1,9,2)=4
        sub_table(10,1,10,1)=2
        sub_table(10,1,10,2)=5
        sub_table(10,1,10,3)=5
        !
        !  C_2v four types
        !
        sub_table(12,1,1,1)=1
        sub_table(12,1,1,2)=1
        sub_table(12,1,2,1)=1
        sub_table(12,1,2,2)=2
        sub_table(12,1,3,1)=1
        sub_table(12,1,3,2)=1
        sub_table(12,1,4,1)=1
        sub_table(12,1,4,2)=2
        sub_table(12,1,5,1)=2
        sub_table(12,1,5,2)=3
        sub_table(12,1,5,3)=4
        sub_table(12,1,6,1)=1
        sub_table(12,1,6,2)=2
        sub_table(12,1,7,1)=1
        sub_table(12,1,7,2)=1
        sub_table(12,1,8,1)=1
        sub_table(12,1,8,2)=2
        sub_table(12,1,9,1)=1
        sub_table(12,1,9,2)=1
        sub_table(12,1,10,1)=2
        sub_table(12,1,10,2)=3
        sub_table(12,1,10,3)=4

        sub_table(12,2,1,1)=1
        sub_table(12,2,1,2)=1
        sub_table(12,2,2,1)=1
        sub_table(12,2,2,2)=2
        sub_table(12,2,3,1)=1
        sub_table(12,2,3,2)=2
        sub_table(12,2,4,1)=1
        sub_table(12,2,4,2)=1
        sub_table(12,2,5,1)=2
        sub_table(12,2,5,2)=3
        sub_table(12,2,5,3)=4
        sub_table(12,2,6,1)=1
        sub_table(12,2,6,2)=2
        sub_table(12,2,7,1)=1
        sub_table(12,2,7,2)=1
        sub_table(12,2,8,1)=1
        sub_table(12,2,8,2)=1
        sub_table(12,2,9,1)=1
        sub_table(12,2,9,2)=2
        sub_table(12,2,10,1)=2
        sub_table(12,2,10,2)=3
        sub_table(12,2,10,3)=4

        sub_table(12,3,1,1)=1
        sub_table(12,3,1,2)=1
        sub_table(12,3,2,1)=1
        sub_table(12,3,2,2)=3
        sub_table(12,3,3,1)=1
        sub_table(12,3,3,2)=1
        sub_table(12,3,4,1)=1
        sub_table(12,3,4,2)=3
        sub_table(12,3,5,1)=2
        sub_table(12,3,5,2)=2
        sub_table(12,3,5,3)=4
        sub_table(12,3,6,1)=1
        sub_table(12,3,6,2)=2
        sub_table(12,3,7,1)=1
        sub_table(12,3,7,2)=4
        sub_table(12,3,8,1)=1
        sub_table(12,3,8,2)=2
        sub_table(12,3,9,1)=1
        sub_table(12,3,9,2)=4
        sub_table(12,3,10,1)=2
        sub_table(12,3,10,2)=1
        sub_table(12,3,10,3)=3

        sub_table(12,4,1,1)=1
        sub_table(12,4,1,2)=1
        sub_table(12,4,2,1)=1
        sub_table(12,4,2,2)=3
        sub_table(12,4,3,1)=1
        sub_table(12,4,3,2)=3
        sub_table(12,4,4,1)=1
        sub_table(12,4,4,2)=1
        sub_table(12,4,5,1)=2
        sub_table(12,4,5,2)=2
        sub_table(12,4,5,3)=4
        sub_table(12,4,6,1)=1
        sub_table(12,4,6,2)=2
        sub_table(12,4,7,1)=1
        sub_table(12,4,7,2)=4
        sub_table(12,4,8,1)=1
        sub_table(12,4,8,2)=4
        sub_table(12,4,9,1)=1
        sub_table(12,4,9,2)=2
        sub_table(12,4,10,1)=2
        sub_table(12,4,10,2)=1
        sub_table(12,4,10,3)=3

        sub_table(12,5,1,1)=1
        sub_table(12,5,1,2)=1
        sub_table(12,5,2,1)=1
        sub_table(12,5,2,2)=4
        sub_table(12,5,3,1)=1
        sub_table(12,5,3,2)=1
        sub_table(12,5,4,1)=1
        sub_table(12,5,4,2)=4
        sub_table(12,5,5,1)=2
        sub_table(12,5,5,2)=2
        sub_table(12,5,5,3)=3
        sub_table(12,5,6,1)=1
        sub_table(12,5,6,2)=2
        sub_table(12,5,7,1)=1
        sub_table(12,5,7,2)=3
        sub_table(12,5,8,1)=1
        sub_table(12,5,8,2)=2
        sub_table(12,5,9,1)=1
        sub_table(12,5,9,2)=3
        sub_table(12,5,10,1)=2
        sub_table(12,5,10,2)=1
        sub_table(12,5,10,3)=4

        sub_table(12,6,1,1)=1
        sub_table(12,6,1,2)=1
        sub_table(12,6,2,1)=1
        sub_table(12,6,2,2)=4
        sub_table(12,6,3,1)=1
        sub_table(12,6,3,2)=4
        sub_table(12,6,4,1)=1
        sub_table(12,6,4,2)=1
        sub_table(12,6,5,1)=2
        sub_table(12,6,5,2)=2
        sub_table(12,6,5,3)=3
        sub_table(12,6,6,1)=1
        sub_table(12,6,6,2)=2
        sub_table(12,6,7,1)=1
        sub_table(12,6,7,2)=3
        sub_table(12,6,8,1)=1
        sub_table(12,6,8,2)=2
        sub_table(12,6,9,1)=1
        sub_table(12,6,9,2)=3
        sub_table(12,6,10,1)=2
        sub_table(12,6,10,2)=1
        sub_table(12,6,10,3)=4

        sub_table(14,1,1,1)=1
        sub_table(14,1,1,2)=1
        sub_table(14,1,2,1)=1
        sub_table(14,1,2,2)=2
        sub_table(14,1,3,1)=1
        sub_table(14,1,3,2)=3
        sub_table(14,1,4,1)=1
        sub_table(14,1,4,2)=4
        sub_table(14,1,5,1)=2
        sub_table(14,1,5,2)=5
        sub_table(14,1,5,3)=5
        sub_table(14,1,6,1)=1
        sub_table(14,1,6,2)=2
        sub_table(14,1,7,1)=1
        sub_table(14,1,7,2)=1
        sub_table(14,1,8,1)=1
        sub_table(14,1,8,2)=4
        sub_table(14,1,9,1)=1
        sub_table(14,1,9,2)=3
        sub_table(14,1,10,1)=2
        sub_table(14,1,10,2)=5
        sub_table(14,1,10,3)=5

        sub_table(16,1,1,1)=1
        sub_table(16,1,1,2)=1
        sub_table(16,1,2,1)=1
        sub_table(16,1,2,2)=1
        sub_table(16,1,3,1)=1
        sub_table(16,1,3,2)=1
        sub_table(16,1,4,1)=1
        sub_table(16,1,4,2)=1
        sub_table(16,1,5,1)=2
        sub_table(16,1,5,2)=2
        sub_table(16,1,5,3)=2
        sub_table(16,1,6,1)=1
        sub_table(16,1,6,2)=3
        sub_table(16,1,7,1)=1
        sub_table(16,1,7,2)=3
        sub_table(16,1,8,1)=1
        sub_table(16,1,8,2)=3
        sub_table(16,1,9,1)=1
        sub_table(16,1,9,2)=3
        sub_table(16,1,10,1)=2
        sub_table(16,1,10,2)=4
        sub_table(16,1,10,3)=4

        sub_table(16,2,1,1)=1
        sub_table(16,2,1,2)=1
        sub_table(16,2,2,1)=1
        sub_table(16,2,2,2)=2
        sub_table(16,2,3,1)=1
        sub_table(16,2,3,2)=1
        sub_table(16,2,4,1)=1
        sub_table(16,2,4,2)=2
        sub_table(16,2,5,1)=2
        sub_table(16,2,5,2)=1
        sub_table(16,2,5,3)=2
        sub_table(16,2,6,1)=1
        sub_table(16,2,6,2)=3
        sub_table(16,2,7,1)=1
        sub_table(16,2,7,2)=4
        sub_table(16,2,8,1)=1
        sub_table(16,2,8,2)=3
        sub_table(16,2,9,1)=1
        sub_table(16,2,9,2)=4
        sub_table(16,2,10,1)=2
        sub_table(16,2,10,2)=3
        sub_table(16,2,10,3)=4

        sub_table(16,3,1,1)=1
        sub_table(16,3,1,2)=1
        sub_table(16,3,2,1)=1
        sub_table(16,3,2,2)=2
        sub_table(16,3,3,1)=1
        sub_table(16,3,3,2)=2
        sub_table(16,3,4,1)=1
        sub_table(16,3,4,2)=1
        sub_table(16,3,5,1)=2
        sub_table(16,3,5,2)=1
        sub_table(16,3,5,3)=2
        sub_table(16,3,6,1)=1
        sub_table(16,3,6,2)=3
        sub_table(16,3,7,1)=1
        sub_table(16,3,7,2)=4
        sub_table(16,3,8,1)=1
        sub_table(16,3,8,2)=4
        sub_table(16,3,9,1)=1
        sub_table(16,3,9,2)=3
        sub_table(16,3,10,1)=2
        sub_table(16,3,10,2)=3
        sub_table(16,3,10,3)=4

        sub_table(18,1,1,1)=1
        sub_table(18,1,1,2)=1
        sub_table(18,1,2,1)=1
        sub_table(18,1,2,2)=1
        sub_table(18,1,3,1)=1
        sub_table(18,1,3,2)=2
        sub_table(18,1,4,1)=1
        sub_table(18,1,4,2)=2
        sub_table(18,1,5,1)=2
        sub_table(18,1,5,2)=3
        sub_table(18,1,5,3)=4
        sub_table(18,1,6,1)=1
        sub_table(18,1,6,2)=5
        sub_table(18,1,7,1)=1
        sub_table(18,1,7,2)=5
        sub_table(18,1,8,1)=1
        sub_table(18,1,8,2)=6
        sub_table(18,1,9,1)=1
        sub_table(18,1,9,2)=6
        sub_table(18,1,10,1)=2
        sub_table(18,1,10,2)=7
        sub_table(18,1,10,3)=8

        sub_table(20,1,1,1)=1
        sub_table(20,1,1,2)=1
        sub_table(20,1,2,1)=1
        sub_table(20,1,2,2)=2
        sub_table(20,1,3,1)=1
        sub_table(20,1,3,2)=1
        sub_table(20,1,4,1)=1
        sub_table(20,1,4,2)=2
        sub_table(20,1,5,1)=2
        sub_table(20,1,5,2)=3
        sub_table(20,1,5,3)=4
        sub_table(20,1,6,1)=1
        sub_table(20,1,6,2)=5
        sub_table(20,1,7,1)=1
        sub_table(20,1,7,2)=6
        sub_table(20,1,8,1)=1
        sub_table(20,1,8,2)=5
        sub_table(20,1,9,1)=1
        sub_table(20,1,9,2)=6
        sub_table(20,1,10,1)=2
        sub_table(20,1,10,2)=7
        sub_table(20,1,10,3)=8

        sub_table(20,2,1,1)=1
        sub_table(20,2,1,2)=1
        sub_table(20,2,2,1)=1
        sub_table(20,2,2,2)=2
        sub_table(20,2,3,1)=1
        sub_table(20,2,3,2)=2
        sub_table(20,2,4,1)=1
        sub_table(20,2,4,2)=1
        sub_table(20,2,5,1)=2
        sub_table(20,2,5,2)=3
        sub_table(20,2,5,3)=4
        sub_table(20,2,6,1)=1
        sub_table(20,2,6,2)=5
        sub_table(20,2,7,1)=1
        sub_table(20,2,7,2)=6
        sub_table(20,2,8,1)=1
        sub_table(20,2,8,2)=6
        sub_table(20,2,9,1)=1
        sub_table(20,2,9,2)=5
        sub_table(20,2,10,1)=2
        sub_table(20,2,10,2)=7
        sub_table(20,2,10,3)=8

        sub_table(24,1,1,1)=1
        sub_table(24,1,1,2)=1
        sub_table(24,1,2,1)=1
        sub_table(24,1,2,2)=2
        sub_table(24,1,3,1)=1
        sub_table(24,1,3,2)=3
        sub_table(24,1,4,1)=1
        sub_table(24,1,4,2)=4
        sub_table(24,1,5,1)=2
        sub_table(24,1,5,2)=5
        sub_table(24,1,5,3)=5
        sub_table(24,1,6,1)=1
        sub_table(24,1,6,2)=3
        sub_table(24,1,7,1)=1
        sub_table(24,1,7,2)=4
        sub_table(24,1,8,1)=1
        sub_table(24,1,8,2)=1
        sub_table(24,1,9,1)=1
        sub_table(24,1,9,2)=2
        sub_table(24,1,10,1)=2
        sub_table(24,1,10,2)=5
        sub_table(24,1,10,3)=5

        sub_table(24,2,1,1)=1
        sub_table(24,2,1,2)=1
        sub_table(24,2,2,1)=1
        sub_table(24,2,2,2)=2
        sub_table(24,2,3,1)=1
        sub_table(24,2,3,2)=4
        sub_table(24,2,4,1)=1
        sub_table(24,2,4,2)=3
        sub_table(24,2,5,1)=2
        sub_table(24,2,5,2)=5
        sub_table(24,2,5,3)=5
        sub_table(24,2,6,1)=1
        sub_table(24,2,6,2)=3
        sub_table(24,2,7,1)=1
        sub_table(24,2,7,2)=4
        sub_table(24,2,8,1)=1
        sub_table(24,2,8,2)=2
        sub_table(24,2,9,1)=1
        sub_table(24,2,9,2)=1
        sub_table(24,2,10,1)=2
        sub_table(24,2,10,2)=5
        sub_table(24,2,10,3)=5

        sub_table(26,1,1,1)=1
        sub_table(26,1,1,2)=1
        sub_table(26,1,2,1)=1
        sub_table(26,1,2,2)=1
        sub_table(26,1,3,1)=1
        sub_table(26,1,3,2)=2
        sub_table(26,1,4,1)=1
        sub_table(26,1,4,2)=2
        sub_table(26,1,5,1)=2
        sub_table(26,1,5,2)=3
        sub_table(26,1,5,3)=4
        sub_table(26,1,6,1)=1
        sub_table(26,1,6,2)=2
        sub_table(26,1,7,1)=1
        sub_table(26,1,7,2)=2
        sub_table(26,1,8,1)=1
        sub_table(26,1,8,2)=1
        sub_table(26,1,9,1)=1
        sub_table(26,1,9,2)=1
        sub_table(26,1,10,1)=2
        sub_table(26,1,10,2)=3
        sub_table(26,1,10,3)=4

     CASE(23)
!
! D_6h
!
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=1
        sub_table(2,1,3,2)=1
        sub_table(2,1,4,1)=1
        sub_table(2,1,4,2)=1
        sub_table(2,1,5,1)=2
        sub_table(2,1,5,2)=1
        sub_table(2,1,5,3)=1
        sub_table(2,1,6,1)=2
        sub_table(2,1,6,2)=1
        sub_table(2,1,6,3)=1
        sub_table(2,1,7,1)=1
        sub_table(2,1,7,2)=2
        sub_table(2,1,8,1)=1
        sub_table(2,1,8,2)=2
        sub_table(2,1,9,1)=1
        sub_table(2,1,9,2)=2
        sub_table(2,1,10,1)=1
        sub_table(2,1,10,2)=2
        sub_table(2,1,11,1)=2
        sub_table(2,1,11,2)=2
        sub_table(2,1,11,3)=2
        sub_table(2,1,12,1)=2
        sub_table(2,1,12,2)=2
        sub_table(2,1,12,3)=2

        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=1
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=2
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=2
        sub_table(3,1,5,1)=2
        sub_table(3,1,5,2)=2
        sub_table(3,1,5,3)=2
        sub_table(3,1,6,1)=2
        sub_table(3,1,6,2)=1
        sub_table(3,1,6,3)=1
        sub_table(3,1,7,1)=1
        sub_table(3,1,7,2)=2
        sub_table(3,1,8,1)=1
        sub_table(3,1,8,2)=2
        sub_table(3,1,9,1)=1
        sub_table(3,1,9,2)=1
        sub_table(3,1,10,1)=1
        sub_table(3,1,10,2)=1
        sub_table(3,1,11,1)=2
        sub_table(3,1,11,2)=1
        sub_table(3,1,11,3)=1
        sub_table(3,1,12,1)=2
        sub_table(3,1,12,2)=2
        sub_table(3,1,12,3)=2

        sub_table(3,2,1,1)=1
        sub_table(3,2,1,2)=1
        sub_table(3,2,2,1)=1
        sub_table(3,2,2,2)=2
        sub_table(3,2,3,1)=1
        sub_table(3,2,3,2)=1
        sub_table(3,2,4,1)=1
        sub_table(3,2,4,2)=2
        sub_table(3,2,5,1)=2
        sub_table(3,2,5,2)=1
        sub_table(3,2,5,3)=2
        sub_table(3,2,6,1)=2
        sub_table(3,2,6,2)=1
        sub_table(3,2,6,3)=2
        sub_table(3,2,7,1)=1
        sub_table(3,2,7,2)=2
        sub_table(3,2,8,1)=1
        sub_table(3,2,8,2)=1
        sub_table(3,2,9,1)=1
        sub_table(3,2,9,2)=2
        sub_table(3,2,10,1)=1
        sub_table(3,2,10,2)=1
        sub_table(3,2,11,1)=2
        sub_table(3,2,11,2)=1
        sub_table(3,2,11,3)=2
        sub_table(3,2,12,1)=2
        sub_table(3,2,12,2)=1
        sub_table(3,2,12,3)=2

        sub_table(3,3,1,1)=1
        sub_table(3,3,1,2)=1
        sub_table(3,3,2,1)=1
        sub_table(3,3,2,2)=2
        sub_table(3,3,3,1)=1
        sub_table(3,3,3,2)=2
        sub_table(3,3,4,1)=1
        sub_table(3,3,4,2)=1
        sub_table(3,3,5,1)=2
        sub_table(3,3,5,2)=1
        sub_table(3,3,5,3)=2
        sub_table(3,3,6,1)=2
        sub_table(3,3,6,2)=1
        sub_table(3,3,6,3)=2
        sub_table(3,3,7,1)=1
        sub_table(3,3,7,2)=2
        sub_table(3,3,8,1)=1
        sub_table(3,3,8,2)=1
        sub_table(3,3,9,1)=1
        sub_table(3,3,9,2)=1
        sub_table(3,3,10,1)=1
        sub_table(3,3,10,2)=2
        sub_table(3,3,11,1)=2
        sub_table(3,3,11,2)=1
        sub_table(3,3,11,3)=2
        sub_table(3,3,12,1)=2
        sub_table(3,3,12,2)=1
        sub_table(3,3,12,3)=2

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2
        sub_table(4,1,5,1)=2
        sub_table(4,1,5,2)=2
        sub_table(4,1,5,3)=2
        sub_table(4,1,6,1)=2
        sub_table(4,1,6,2)=1
        sub_table(4,1,6,3)=1
        sub_table(4,1,7,1)=1
        sub_table(4,1,7,2)=1
        sub_table(4,1,8,1)=1
        sub_table(4,1,8,2)=1
        sub_table(4,1,9,1)=1
        sub_table(4,1,9,2)=2
        sub_table(4,1,10,1)=1
        sub_table(4,1,10,2)=2
        sub_table(4,1,11,1)=2
        sub_table(4,1,11,2)=2
        sub_table(4,1,11,3)=2
        sub_table(4,1,12,1)=2
        sub_table(4,1,12,2)=1
        sub_table(4,1,12,3)=1

        sub_table(4,2,1,1)=1
        sub_table(4,2,1,2)=1
        sub_table(4,2,2,1)=1
        sub_table(4,2,2,2)=2
        sub_table(4,2,3,1)=1
        sub_table(4,2,3,2)=1
        sub_table(4,2,4,1)=1
        sub_table(4,2,4,2)=2
        sub_table(4,2,5,1)=2
        sub_table(4,2,5,2)=1
        sub_table(4,2,5,3)=2
        sub_table(4,2,6,1)=2
        sub_table(4,2,6,2)=1
        sub_table(4,2,6,3)=2
        sub_table(4,2,7,1)=1
        sub_table(4,2,7,2)=1
        sub_table(4,2,8,1)=1
        sub_table(4,2,8,2)=2
        sub_table(4,2,9,1)=1
        sub_table(4,2,9,2)=1
        sub_table(4,2,10,1)=1
        sub_table(4,2,10,2)=2
        sub_table(4,2,11,1)=2
        sub_table(4,2,11,2)=1
        sub_table(4,2,11,3)=2
        sub_table(4,2,12,1)=2
        sub_table(4,2,12,2)=1
        sub_table(4,2,12,3)=2

        sub_table(4,3,1,1)=1
        sub_table(4,3,1,2)=1
        sub_table(4,3,2,1)=1
        sub_table(4,3,2,2)=2
        sub_table(4,3,3,1)=1
        sub_table(4,3,3,2)=2
        sub_table(4,3,4,1)=1
        sub_table(4,3,4,2)=1
        sub_table(4,3,5,1)=2
        sub_table(4,3,5,2)=1
        sub_table(4,3,5,3)=2
        sub_table(4,3,6,1)=2
        sub_table(4,3,6,2)=1
        sub_table(4,3,6,3)=2
        sub_table(4,3,7,1)=1
        sub_table(4,3,7,2)=1
        sub_table(4,3,8,1)=1
        sub_table(4,3,8,2)=2
        sub_table(4,3,9,1)=1
        sub_table(4,3,9,2)=2
        sub_table(4,3,10,1)=1
        sub_table(4,3,10,2)=1
        sub_table(4,3,11,1)=2
        sub_table(4,3,11,2)=1
        sub_table(4,3,11,3)=2
        sub_table(4,3,12,1)=2
        sub_table(4,3,12,2)=1
        sub_table(4,3,12,3)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=1
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=1
        sub_table(5,1,4,1)=1
        sub_table(5,1,4,2)=1
        sub_table(5,1,5,1)=2
        sub_table(5,1,5,2)=2
        sub_table(5,1,5,3)=3
        sub_table(5,1,6,1)=2
        sub_table(5,1,6,2)=2
        sub_table(5,1,6,3)=3
        sub_table(5,1,7,1)=1
        sub_table(5,1,7,2)=1
        sub_table(5,1,8,1)=1
        sub_table(5,1,8,2)=1
        sub_table(5,1,9,1)=1
        sub_table(5,1,9,2)=1
        sub_table(5,1,10,1)=1
        sub_table(5,1,10,2)=1
        sub_table(5,1,11,1)=2
        sub_table(5,1,11,2)=2
        sub_table(5,1,11,3)=3
        sub_table(5,1,12,1)=2
        sub_table(5,1,12,2)=2
        sub_table(5,1,12,3)=3

        sub_table(7,1,1,1)=1
        sub_table(7,1,1,2)=1
        sub_table(7,1,2,1)=1
        sub_table(7,1,2,2)=1
        sub_table(7,1,3,1)=1
        sub_table(7,1,3,2)=2
        sub_table(7,1,4,1)=1
        sub_table(7,1,4,2)=2
        sub_table(7,1,5,1)=2
        sub_table(7,1,5,2)=3
        sub_table(7,1,5,3)=4
        sub_table(7,1,6,1)=2
        sub_table(7,1,6,2)=5
        sub_table(7,1,6,3)=6
        sub_table(7,1,7,1)=1
        sub_table(7,1,7,2)=1
        sub_table(7,1,8,1)=1
        sub_table(7,1,8,2)=1
        sub_table(7,1,9,1)=1
        sub_table(7,1,9,2)=2
        sub_table(7,1,10,1)=1
        sub_table(7,1,10,2)=2
        sub_table(7,1,11,1)=2
        sub_table(7,1,11,2)=3
        sub_table(7,1,11,3)=4
        sub_table(7,1,12,1)=2
        sub_table(7,1,12,2)=5
        sub_table(7,1,12,3)=6

        sub_table(8,1,1,1)=1
        sub_table(8,1,1,2)=1
        sub_table(8,1,2,1)=1
        sub_table(8,1,2,2)=2
        sub_table(8,1,3,1)=1
        sub_table(8,1,3,2)=3
        sub_table(8,1,4,1)=1
        sub_table(8,1,4,2)=4
        sub_table(8,1,5,1)=2
        sub_table(8,1,5,2)=3
        sub_table(8,1,5,3)=4
        sub_table(8,1,6,1)=2
        sub_table(8,1,6,2)=1
        sub_table(8,1,6,3)=2
        sub_table(8,1,7,1)=1
        sub_table(8,1,7,2)=1
        sub_table(8,1,8,1)=1
        sub_table(8,1,8,2)=2
        sub_table(8,1,9,1)=1
        sub_table(8,1,9,2)=3
        sub_table(8,1,10,1)=1
        sub_table(8,1,10,2)=4
        sub_table(8,1,11,1)=2
        sub_table(8,1,11,2)=3
        sub_table(8,1,11,3)=4
        sub_table(8,1,12,1)=2
        sub_table(8,1,12,2)=1
        sub_table(8,1,12,3)=2

        sub_table(9,1,1,1)=1
        sub_table(9,1,1,2)=1
        sub_table(9,1,2,1)=1
        sub_table(9,1,2,2)=2
        sub_table(9,1,3,1)=1
        sub_table(9,1,3,2)=1
        sub_table(9,1,4,1)=1
        sub_table(9,1,4,2)=2
        sub_table(9,1,5,1)=2
        sub_table(9,1,5,2)=3
        sub_table(9,1,5,3)=3
        sub_table(9,1,6,1)=2
        sub_table(9,1,6,2)=3
        sub_table(9,1,6,3)=3
        sub_table(9,1,7,1)=1
        sub_table(9,1,7,2)=1
        sub_table(9,1,8,1)=1
        sub_table(9,1,8,2)=2
        sub_table(9,1,9,1)=1
        sub_table(9,1,9,2)=1
        sub_table(9,1,10,1)=1
        sub_table(9,1,10,2)=2
        sub_table(9,1,11,1)=2
        sub_table(9,1,11,2)=3
        sub_table(9,1,11,3)=3
        sub_table(9,1,12,1)=2
        sub_table(9,1,12,2)=3
        sub_table(9,1,12,3)=3

        sub_table(9,2,1,1)=1
        sub_table(9,2,1,2)=1
        sub_table(9,2,2,1)=1
        sub_table(9,2,2,2)=2
        sub_table(9,2,3,1)=1
        sub_table(9,2,3,2)=2
        sub_table(9,2,4,1)=1
        sub_table(9,2,4,2)=1
        sub_table(9,2,5,1)=2
        sub_table(9,2,5,2)=3
        sub_table(9,2,5,3)=3
        sub_table(9,2,6,1)=2
        sub_table(9,2,6,2)=3
        sub_table(9,2,6,3)=3
        sub_table(9,2,7,1)=1
        sub_table(9,2,7,2)=1
        sub_table(9,2,8,1)=1
        sub_table(9,2,8,2)=2
        sub_table(9,2,9,1)=1
        sub_table(9,2,9,2)=2
        sub_table(9,2,10,1)=1
        sub_table(9,2,10,2)=1
        sub_table(9,2,11,1)=2
        sub_table(9,2,11,2)=3
        sub_table(9,2,11,3)=3
        sub_table(9,2,12,1)=2
        sub_table(9,2,12,2)=3
        sub_table(9,2,12,3)=3

        sub_table(11,1,1,1)=1
        sub_table(11,1,1,2)=1
        sub_table(11,1,2,1)=1
        sub_table(11,1,2,2)=2
        sub_table(11,1,3,1)=1
        sub_table(11,1,3,2)=3
        sub_table(11,1,4,1)=1
        sub_table(11,1,4,2)=4
        sub_table(11,1,5,1)=2
        sub_table(11,1,5,2)=5
        sub_table(11,1,5,3)=5
        sub_table(11,1,6,1)=2
        sub_table(11,1,6,2)=6
        sub_table(11,1,6,3)=6
        sub_table(11,1,7,1)=1
        sub_table(11,1,7,2)=1
        sub_table(11,1,8,1)=1
        sub_table(11,1,8,2)=2
        sub_table(11,1,9,1)=1
        sub_table(11,1,9,2)=3
        sub_table(11,1,10,1)=1
        sub_table(11,1,10,2)=4
        sub_table(11,1,11,1)=2
        sub_table(11,1,11,2)=5
        sub_table(11,1,11,3)=5
        sub_table(11,1,12,1)=2
        sub_table(11,1,12,2)=6
        sub_table(11,1,12,3)=6

        sub_table(12,1,1,1)=1
        sub_table(12,1,1,2)=1
        sub_table(12,1,2,1)=1
        sub_table(12,1,2,2)=2
        sub_table(12,1,3,1)=1
        sub_table(12,1,3,2)=3
        sub_table(12,1,4,1)=1
        sub_table(12,1,4,2)=4
        sub_table(12,1,5,1)=2
        sub_table(12,1,5,2)=3
        sub_table(12,1,5,3)=4
        sub_table(12,1,6,1)=2
        sub_table(12,1,6,2)=1
        sub_table(12,1,6,3)=2
        sub_table(12,1,7,1)=1
        sub_table(12,1,7,2)=2
        sub_table(12,1,8,1)=1
        sub_table(12,1,8,2)=1
        sub_table(12,1,9,1)=1
        sub_table(12,1,9,2)=4
        sub_table(12,1,10,1)=1
        sub_table(12,1,10,2)=3
        sub_table(12,1,11,1)=2
        sub_table(11,1,11,2)=3
        sub_table(12,1,11,3)=4
        sub_table(12,1,12,1)=2
        sub_table(12,1,12,2)=1
        sub_table(12,1,12,3)=2

        sub_table(12,2,1,1)=1
        sub_table(12,2,1,2)=1
        sub_table(12,2,2,1)=1
        sub_table(12,2,2,2)=3
        sub_table(12,2,3,1)=1
        sub_table(12,2,3,2)=2
        sub_table(12,2,4,1)=1
        sub_table(12,2,4,2)=4
        sub_table(12,2,5,1)=2
        sub_table(12,2,5,2)=2
        sub_table(12,2,5,3)=4
        sub_table(12,2,6,1)=2
        sub_table(12,2,6,2)=1
        sub_table(12,2,6,3)=3
        sub_table(12,2,7,1)=1
        sub_table(12,2,7,2)=2
        sub_table(12,2,8,1)=1
        sub_table(12,2,8,2)=4
        sub_table(12,2,9,1)=1
        sub_table(12,2,9,2)=1
        sub_table(12,2,10,1)=1
        sub_table(12,2,10,2)=3
        sub_table(12,2,11,1)=2
        sub_table(11,2,11,2)=1
        sub_table(12,2,11,3)=3
        sub_table(12,2,12,1)=2
        sub_table(12,2,12,2)=2
        sub_table(12,2,12,3)=4

        sub_table(12,3,1,1)=1
        sub_table(12,3,1,2)=1
        sub_table(12,3,2,1)=1
        sub_table(12,3,2,2)=3
        sub_table(12,3,3,1)=1
        sub_table(12,3,3,2)=4
        sub_table(12,3,4,1)=1
        sub_table(12,3,4,2)=2
        sub_table(12,3,5,1)=2
        sub_table(12,3,5,2)=2
        sub_table(12,3,5,3)=4
        sub_table(12,3,6,1)=2
        sub_table(12,3,6,2)=1
        sub_table(12,3,6,3)=3
        sub_table(12,3,7,1)=1
        sub_table(12,3,7,2)=2
        sub_table(12,3,8,1)=1
        sub_table(12,3,8,2)=4
        sub_table(12,3,9,1)=1
        sub_table(12,3,9,2)=3
        sub_table(12,3,10,1)=1
        sub_table(12,3,10,2)=1
        sub_table(12,3,11,1)=2
        sub_table(11,3,11,2)=1
        sub_table(12,3,11,3)=3
        sub_table(12,3,12,1)=2
        sub_table(12,3,12,2)=2
        sub_table(12,3,12,3)=4

        sub_table(13,1,1,1)=1
        sub_table(13,1,1,2)=1
        sub_table(13,1,2,1)=1
        sub_table(13,1,2,2)=2
        sub_table(13,1,3,1)=1
        sub_table(13,1,3,2)=1
        sub_table(13,1,4,1)=1
        sub_table(13,1,4,2)=2
        sub_table(13,1,5,1)=2
        sub_table(13,1,5,2)=3
        sub_table(13,1,5,3)=3
        sub_table(13,1,6,1)=2
        sub_table(13,1,6,2)=3
        sub_table(13,1,6,3)=3
        sub_table(13,1,7,1)=1
        sub_table(13,1,7,2)=2
        sub_table(13,1,8,1)=1
        sub_table(13,1,8,2)=1
        sub_table(13,1,9,1)=1
        sub_table(13,1,9,2)=2
        sub_table(13,1,10,1)=1
        sub_table(13,1,10,2)=1
        sub_table(13,1,11,1)=2
        sub_table(13,1,11,2)=3
        sub_table(13,1,11,3)=3
        sub_table(13,1,12,1)=2
        sub_table(13,1,12,2)=3
        sub_table(13,1,12,3)=3

        sub_table(13,2,1,1)=1
        sub_table(13,2,1,2)=1
        sub_table(13,2,2,1)=1
        sub_table(13,2,2,2)=2
        sub_table(13,2,3,1)=1
        sub_table(13,2,3,2)=2
        sub_table(13,2,4,1)=1
        sub_table(13,2,4,2)=1
        sub_table(13,2,5,1)=2
        sub_table(13,2,5,2)=3
        sub_table(13,2,5,3)=3
        sub_table(13,2,6,1)=2
        sub_table(13,2,6,2)=3
        sub_table(13,2,6,3)=3
        sub_table(13,2,7,1)=1
        sub_table(13,2,7,2)=2
        sub_table(13,2,8,1)=1
        sub_table(13,2,8,2)=1
        sub_table(13,2,9,1)=1
        sub_table(13,2,9,2)=1
        sub_table(13,2,10,1)=1
        sub_table(13,2,10,2)=2
        sub_table(13,2,11,1)=2
        sub_table(13,2,11,2)=3
        sub_table(13,2,11,3)=3
        sub_table(13,2,12,1)=2
        sub_table(13,2,12,2)=3
        sub_table(13,2,12,3)=3

        sub_table(15,1,1,1)=1
        sub_table(15,1,1,2)=1
        sub_table(15,1,2,1)=1
        sub_table(15,1,2,2)=2
        sub_table(15,1,3,1)=1
        sub_table(15,1,3,2)=3
        sub_table(15,1,4,1)=1
        sub_table(15,1,4,2)=4
        sub_table(15,1,5,1)=2
        sub_table(15,1,5,2)=5
        sub_table(15,1,5,3)=5
        sub_table(15,1,6,1)=2
        sub_table(15,1,6,2)=6
        sub_table(15,1,6,3)=6
        sub_table(15,1,7,1)=1
        sub_table(15,1,7,2)=2
        sub_table(15,1,8,1)=1
        sub_table(15,1,8,2)=1
        sub_table(15,1,9,1)=1
        sub_table(15,1,9,2)=4
        sub_table(15,1,10,1)=1
        sub_table(15,1,10,2)=3
        sub_table(15,1,11,1)=2
        sub_table(15,1,11,2)=5
        sub_table(15,1,11,3)=5
        sub_table(15,1,12,1)=2
        sub_table(15,1,12,2)=6
        sub_table(15,1,12,3)=6

        sub_table(16,1,1,1)=1
        sub_table(16,1,1,2)=1
        sub_table(16,1,2,1)=1
        sub_table(16,1,2,2)=1
        sub_table(16,1,3,1)=1
        sub_table(16,1,3,2)=2
        sub_table(16,1,4,1)=1
        sub_table(16,1,4,2)=2
        sub_table(16,1,5,1)=2
        sub_table(16,1,5,2)=2
        sub_table(16,1,5,3)=2
        sub_table(16,1,6,1)=2
        sub_table(16,1,6,2)=1
        sub_table(16,1,6,3)=1
        sub_table(16,1,7,1)=1
        sub_table(16,1,7,2)=3
        sub_table(16,1,8,1)=1
        sub_table(16,1,8,2)=3
        sub_table(16,1,9,1)=1
        sub_table(16,1,9,2)=4
        sub_table(16,1,10,1)=1
        sub_table(16,1,10,2)=4
        sub_table(16,1,11,1)=2
        sub_table(16,1,11,2)=4
        sub_table(16,1,11,3)=4
        sub_table(16,1,12,1)=2
        sub_table(16,1,12,2)=3
        sub_table(16,1,12,3)=3

        sub_table(16,2,1,1)=1
        sub_table(16,2,1,2)=1
        sub_table(16,2,2,1)=1
        sub_table(16,2,2,2)=2
        sub_table(16,2,3,1)=1
        sub_table(16,2,3,2)=1
        sub_table(16,2,4,1)=1
        sub_table(16,2,4,2)=2
        sub_table(16,2,5,1)=2
        sub_table(16,2,5,2)=1
        sub_table(16,2,5,3)=2
        sub_table(16,2,6,1)=2
        sub_table(16,2,6,2)=1
        sub_table(16,2,6,3)=2
        sub_table(16,2,7,1)=1
        sub_table(16,2,7,2)=3
        sub_table(16,2,8,1)=1
        sub_table(16,2,8,2)=4
        sub_table(16,2,9,1)=1
        sub_table(16,2,9,2)=3
        sub_table(16,2,10,1)=1
        sub_table(16,2,10,2)=4
        sub_table(16,2,11,1)=2
        sub_table(16,2,11,2)=3
        sub_table(16,2,11,3)=4
        sub_table(16,2,12,1)=2
        sub_table(16,2,12,2)=3
        sub_table(16,2,12,3)=4

        sub_table(16,3,1,1)=1
        sub_table(16,3,1,2)=1
        sub_table(16,3,2,1)=1
        sub_table(16,3,2,2)=2
        sub_table(16,3,3,1)=1
        sub_table(16,3,3,2)=2
        sub_table(16,3,4,1)=1
        sub_table(16,3,4,2)=1
        sub_table(16,3,5,1)=2
        sub_table(16,3,5,2)=1
        sub_table(16,3,5,3)=2
        sub_table(16,3,6,1)=2
        sub_table(16,3,6,2)=1
        sub_table(16,3,6,3)=2
        sub_table(16,3,7,1)=1
        sub_table(16,3,7,2)=3
        sub_table(16,3,8,1)=1
        sub_table(16,3,8,2)=4
        sub_table(16,3,9,1)=1
        sub_table(16,3,9,2)=3
        sub_table(16,3,10,1)=1
        sub_table(16,3,10,2)=4
        sub_table(16,3,11,1)=2
        sub_table(16,3,11,2)=3
        sub_table(16,3,11,3)=4
        sub_table(16,3,12,1)=2
        sub_table(16,3,12,2)=3
        sub_table(16,3,12,3)=4
 
        sub_table(17,1,1,1)=1
        sub_table(17,1,1,2)=1
        sub_table(17,1,2,1)=1
        sub_table(17,1,2,2)=1
        sub_table(17,1,3,1)=1
        sub_table(17,1,3,2)=4
        sub_table(17,1,4,1)=1
        sub_table(17,1,4,2)=4
        sub_table(17,1,5,1)=2
        sub_table(17,1,5,2)=5
        sub_table(17,1,5,3)=6
        sub_table(17,1,6,1)=2
        sub_table(17,1,6,2)=2
        sub_table(17,1,6,3)=3
        sub_table(17,1,7,1)=1
        sub_table(17,1,7,2)=4
        sub_table(17,1,8,1)=1
        sub_table(17,1,8,2)=4
        sub_table(17,1,9,1)=1
        sub_table(17,1,9,2)=1
        sub_table(17,1,10,1)=1
        sub_table(17,1,10,2)=1
        sub_table(17,1,11,1)=2
        sub_table(17,1,11,2)=2
        sub_table(17,1,11,3)=3
        sub_table(17,1,12,1)=2
        sub_table(17,1,12,2)=5
        sub_table(17,1,12,3)=6
 
        sub_table(19,1,1,1)=1
        sub_table(19,1,1,2)=1
        sub_table(19,1,2,1)=1
        sub_table(19,1,2,2)=1
        sub_table(19,1,3,1)=1
        sub_table(19,1,3,2)=2
        sub_table(19,1,4,1)=1
        sub_table(19,1,4,2)=2
        sub_table(19,1,5,1)=2
        sub_table(19,1,5,2)=3
        sub_table(19,1,5,3)=4
        sub_table(19,1,6,1)=2
        sub_table(19,1,6,2)=5
        sub_table(19,1,6,3)=6
        sub_table(19,1,7,1)=1
        sub_table(19,1,7,2)=7
        sub_table(19,1,8,1)=1
        sub_table(19,1,8,2)=7
        sub_table(19,1,9,1)=1
        sub_table(19,1,9,2)=8
        sub_table(19,1,10,1)=1
        sub_table(19,1,10,2)=8
        sub_table(19,1,11,1)=2
        sub_table(19,1,11,2)=9
        sub_table(19,1,11,3)=10
        sub_table(19,1,12,1)=2
        sub_table(19,1,12,2)=11
        sub_table(19,1,12,3)=12

        sub_table(20,1,1,1)=1
        sub_table(20,1,1,2)=1
        sub_table(20,1,2,1)=1
        sub_table(20,1,2,2)=2
        sub_table(20,1,3,1)=1
        sub_table(20,1,3,2)=3
        sub_table(20,1,4,1)=1
        sub_table(20,1,4,2)=4
        sub_table(20,1,5,1)=2
        sub_table(20,1,5,2)=3
        sub_table(20,1,5,3)=4
        sub_table(20,1,6,1)=2
        sub_table(20,1,6,2)=1
        sub_table(20,1,6,3)=2
        sub_table(20,1,7,1)=1
        sub_table(20,1,7,2)=5
        sub_table(20,1,8,1)=1
        sub_table(20,1,8,2)=6
        sub_table(20,1,9,1)=1
        sub_table(20,1,9,2)=7
        sub_table(20,1,10,1)=1
        sub_table(20,1,10,2)=8
        sub_table(20,1,11,1)=2
        sub_table(20,1,11,2)=7
        sub_table(20,1,11,3)=8
        sub_table(20,1,12,1)=2
        sub_table(20,1,12,2)=5
        sub_table(20,1,12,3)=6

        sub_table(21,1,1,1)=1
        sub_table(21,1,1,2)=1
        sub_table(21,1,2,1)=1
        sub_table(21,1,2,2)=2
        sub_table(21,1,3,1)=1
        sub_table(21,1,3,2)=4
        sub_table(21,1,4,1)=1
        sub_table(21,1,4,2)=5
        sub_table(21,1,5,1)=2
        sub_table(21,1,5,2)=6
        sub_table(21,1,5,3)=6
        sub_table(21,1,6,1)=2
        sub_table(21,1,6,2)=3
        sub_table(21,1,6,3)=3
        sub_table(21,1,7,1)=1
        sub_table(21,1,7,2)=4
        sub_table(21,1,8,1)=1
        sub_table(21,1,8,2)=5
        sub_table(21,1,9,1)=1
        sub_table(21,1,9,2)=1
        sub_table(21,1,10,1)=1
        sub_table(21,1,10,2)=2
        sub_table(21,1,11,1)=2
        sub_table(21,1,11,2)=3
        sub_table(21,1,11,3)=3
        sub_table(21,1,12,1)=2
        sub_table(21,1,12,2)=6
        sub_table(21,1,12,3)=6

        sub_table(21,2,1,1)=1
        sub_table(21,2,1,2)=1
        sub_table(21,2,2,1)=1
        sub_table(21,2,2,2)=2
        sub_table(21,2,3,1)=1
        sub_table(21,2,3,2)=5
        sub_table(21,2,4,1)=1
        sub_table(21,2,4,2)=4
        sub_table(21,2,5,1)=2
        sub_table(21,2,5,2)=6
        sub_table(21,2,5,3)=6
        sub_table(21,2,6,1)=2
        sub_table(21,2,6,2)=3
        sub_table(21,2,6,3)=3
        sub_table(21,2,7,1)=1
        sub_table(21,2,7,2)=4
        sub_table(21,2,8,1)=1
        sub_table(21,2,8,2)=5
        sub_table(21,2,9,1)=1
        sub_table(21,2,9,2)=2
        sub_table(21,2,10,1)=1
        sub_table(21,2,10,2)=1
        sub_table(21,2,11,1)=2
        sub_table(21,2,11,2)=3
        sub_table(21,2,11,3)=3
        sub_table(21,2,12,1)=2
        sub_table(21,2,12,2)=6
        sub_table(21,2,12,3)=6

        sub_table(25,1,1,1)=1
        sub_table(25,1,1,2)=1
        sub_table(25,1,2,1)=1
        sub_table(25,1,2,2)=2
        sub_table(25,1,3,1)=1
        sub_table(25,1,3,2)=1
        sub_table(25,1,4,1)=1
        sub_table(25,1,4,2)=2
        sub_table(25,1,5,1)=2
        sub_table(25,1,5,2)=3
        sub_table(25,1,5,3)=3
        sub_table(25,1,6,1)=2
        sub_table(25,1,6,2)=3
        sub_table(25,1,6,3)=3
        sub_table(25,1,7,1)=1
        sub_table(25,1,7,2)=4
        sub_table(25,1,8,1)=1
        sub_table(25,1,8,2)=5
        sub_table(25,1,9,1)=1
        sub_table(25,1,9,2)=4
        sub_table(25,1,10,1)=1
        sub_table(25,1,10,2)=5
        sub_table(25,1,11,1)=2
        sub_table(25,1,11,2)=6
        sub_table(25,1,11,3)=6
        sub_table(25,1,12,1)=2
        sub_table(25,1,12,2)=6
        sub_table(25,1,12,3)=6

        sub_table(25,2,1,1)=1
        sub_table(25,2,1,2)=1
        sub_table(25,2,2,1)=1
        sub_table(25,2,2,2)=2
        sub_table(25,2,3,1)=1
        sub_table(25,2,3,2)=2
        sub_table(25,2,4,1)=1
        sub_table(25,2,4,2)=1
        sub_table(25,2,5,1)=2
        sub_table(25,2,5,2)=3
        sub_table(25,2,5,3)=3
        sub_table(25,2,6,1)=2
        sub_table(25,2,6,2)=3
        sub_table(25,2,6,3)=3
        sub_table(25,2,7,1)=1
        sub_table(25,2,7,2)=4
        sub_table(25,2,8,1)=1
        sub_table(25,2,8,2)=5
        sub_table(25,2,9,1)=1
        sub_table(25,2,9,2)=5
        sub_table(25,2,10,1)=1
        sub_table(25,2,10,2)=4
        sub_table(25,2,11,1)=2
        sub_table(25,2,11,2)=6
        sub_table(25,2,11,3)=6
        sub_table(25,2,12,1)=2
        sub_table(25,2,12,2)=6
        sub_table(25,2,12,3)=6

     CASE(24)
!
! D_2d
!
        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=2
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=1
        sub_table(3,1,5,1)=2
        sub_table(3,1,5,2)=1
        sub_table(3,1,5,3)=2

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=1
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=1
        sub_table(4,1,5,1)=2
        sub_table(4,1,5,2)=2
        sub_table(4,1,5,3)=2

        sub_table(4,2,1,1)=1
        sub_table(4,2,1,2)=1
        sub_table(4,2,2,1)=1
        sub_table(4,2,2,2)=2
        sub_table(4,2,3,1)=1
        sub_table(4,2,3,2)=1
        sub_table(4,2,4,1)=1
        sub_table(4,2,4,2)=2
        sub_table(4,2,5,1)=2
        sub_table(4,2,5,2)=1
        sub_table(4,2,5,3)=2

        sub_table(8,1,1,1)=1
        sub_table(8,1,1,2)=1
        sub_table(8,1,2,1)=1
        sub_table(8,1,2,2)=2
        sub_table(8,1,3,1)=1
        sub_table(8,1,3,2)=1
        sub_table(8,1,4,1)=1
        sub_table(8,1,4,2)=2
        sub_table(8,1,5,1)=2
        sub_table(8,1,5,2)=3
        sub_table(8,1,5,3)=4

        sub_table(12,1,1,1)=1
        sub_table(12,1,1,2)=1
        sub_table(12,1,2,1)=1
        sub_table(12,1,2,2)=2
        sub_table(12,1,3,1)=1
        sub_table(12,1,3,2)=2
        sub_table(12,1,4,1)=1
        sub_table(12,1,4,2)=1
        sub_table(12,1,5,1)=2
        sub_table(12,1,5,2)=3
        sub_table(12,1,5,3)=4

        sub_table(26,1,1,1)=1
        sub_table(26,1,1,2)=1
        sub_table(26,1,2,1)=1
        sub_table(26,1,2,2)=1
        sub_table(26,1,3,1)=1
        sub_table(26,1,3,2)=2
        sub_table(26,1,4,1)=1
        sub_table(26,1,4,2)=2
        sub_table(26,1,5,1)=2
        sub_table(26,1,5,2)=3
        sub_table(26,1,5,3)=4
 
     CASE(25)
!
!  D_3d
!
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=2
        sub_table(2,1,3,2)=1
        sub_table(2,1,3,3)=1
        sub_table(2,1,4,1)=1
        sub_table(2,1,4,2)=2
        sub_table(2,1,5,1)=1
        sub_table(2,1,5,2)=2
        sub_table(2,1,6,1)=2
        sub_table(2,1,6,2)=2
        sub_table(2,1,6,3)=2

        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=2
        sub_table(3,1,3,2)=1
        sub_table(3,1,3,3)=2
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=2
        sub_table(3,1,5,1)=1
        sub_table(3,1,5,2)=1
        sub_table(3,1,6,1)=2
        sub_table(3,1,6,2)=1
        sub_table(3,1,6,3)=2

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=2
        sub_table(4,1,3,1)=2
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=1
        sub_table(4,1,5,1)=1
        sub_table(4,1,5,2)=2
        sub_table(4,1,6,1)=2
        sub_table(4,1,6,2)=1
        sub_table(4,1,6,3)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=1
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=2
        sub_table(5,1,3,3)=3
        sub_table(5,1,4,1)=1
        sub_table(5,1,4,2)=1
        sub_table(5,1,5,1)=1
        sub_table(5,1,5,2)=1
        sub_table(5,1,6,1)=2
        sub_table(5,1,6,2)=2
        sub_table(5,1,6,3)=3

        sub_table(9,1,1,1)=1
        sub_table(9,1,1,2)=1
        sub_table(9,1,2,1)=1
        sub_table(9,1,2,2)=2
        sub_table(9,1,3,1)=2
        sub_table(9,1,3,2)=3
        sub_table(9,1,3,3)=3
        sub_table(9,1,4,1)=1
        sub_table(9,1,4,2)=1
        sub_table(9,1,5,1)=1
        sub_table(9,1,5,2)=2
        sub_table(9,1,6,1)=2
        sub_table(9,1,6,2)=3
        sub_table(9,1,6,3)=3

        sub_table(13,1,1,1)=1
        sub_table(13,1,1,2)=1
        sub_table(13,1,2,1)=1
        sub_table(13,1,2,2)=2
        sub_table(13,1,3,1)=2
        sub_table(13,1,3,2)=3
        sub_table(13,1,3,3)=3
        sub_table(13,1,4,1)=1
        sub_table(13,1,4,2)=2
        sub_table(13,1,5,1)=1
        sub_table(13,1,5,2)=1
        sub_table(13,1,6,1)=2
        sub_table(13,1,6,2)=3
        sub_table(13,1,6,3)=3

        sub_table(16,1,1,1)=1
        sub_table(16,1,1,2)=1
        sub_table(16,1,2,1)=1
        sub_table(16,1,2,2)=2
        sub_table(16,1,3,1)=2
        sub_table(16,1,3,2)=1
        sub_table(16,1,3,3)=2
        sub_table(16,1,4,1)=1
        sub_table(16,1,4,2)=3
        sub_table(16,1,5,1)=1
        sub_table(16,1,5,2)=4
        sub_table(16,1,6,1)=2
        sub_table(16,1,6,2)=3
        sub_table(16,1,6,3)=4

        sub_table(27,1,1,1)=1
        sub_table(27,1,1,2)=1
        sub_table(27,1,2,1)=1
        sub_table(27,1,2,2)=1
        sub_table(27,1,3,1)=2
        sub_table(27,1,3,2)=2
        sub_table(27,1,3,3)=3
        sub_table(27,1,4,1)=1
        sub_table(27,1,4,2)=4
        sub_table(27,1,5,1)=1
        sub_table(27,1,5,2)=4
        sub_table(27,1,6,1)=2
        sub_table(27,1,6,2)=5
        sub_table(27,1,6,3)=6

     CASE(26)
!
!  S_4
!
        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2
     CASE(27)
!
!  S_6
!
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=1
        sub_table(2,1,3,2)=1
        sub_table(2,1,4,1)=1
        sub_table(2,1,4,2)=2
        sub_table(2,1,5,1)=1
        sub_table(2,1,5,2)=2
        sub_table(2,1,6,1)=1
        sub_table(2,1,6,2)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=2
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=3
        sub_table(5,1,4,1)=1
        sub_table(5,1,4,2)=1
        sub_table(5,1,5,1)=1
        sub_table(5,1,5,2)=2
        sub_table(5,1,6,1)=1
        sub_table(5,1,6,2)=3

     CASE(28)
!
!  T
!
        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=1
        sub_table(4,1,4,1)=3
        sub_table(4,1,4,2)=1
        sub_table(4,1,4,3)=2
        sub_table(4,1,4,4)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=3
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=2
        sub_table(5,1,4,1)=3
        sub_table(5,1,4,2)=1
        sub_table(5,1,4,3)=2
        sub_table(5,1,4,4)=3

        sub_table(8,1,1,1)=1
        sub_table(8,1,1,2)=1
        sub_table(8,1,2,1)=1
        sub_table(8,1,2,2)=1
        sub_table(8,1,3,1)=1
        sub_table(8,1,3,2)=1
        sub_table(8,1,4,1)=3
        sub_table(8,1,4,2)=2
        sub_table(8,1,4,3)=3
        sub_table(8,1,4,4)=4

     CASE(29)
!
!  T_h
!
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=1
        sub_table(2,1,3,2)=1
        sub_table(2,1,4,1)=3
        sub_table(2,1,4,2)=1
        sub_table(2,1,4,3)=1
        sub_table(2,1,4,4)=1
        sub_table(2,1,5,1)=1
        sub_table(2,1,5,2)=2
        sub_table(2,1,6,1)=1
        sub_table(2,1,6,2)=2
        sub_table(2,1,7,1)=1
        sub_table(2,1,7,2)=2
        sub_table(2,1,8,1)=3
        sub_table(2,1,8,2)=2
        sub_table(2,1,8,3)=2
        sub_table(2,1,8,4)=2

        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=1
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=1
        sub_table(3,1,4,1)=3
        sub_table(3,1,4,2)=1
        sub_table(3,1,4,3)=2
        sub_table(3,1,4,4)=2
        sub_table(3,1,5,1)=1
        sub_table(3,1,5,2)=2
        sub_table(3,1,6,1)=1
        sub_table(3,1,6,2)=2
        sub_table(3,1,7,1)=1
        sub_table(3,1,7,2)=2
        sub_table(3,1,8,1)=3
        sub_table(3,1,8,2)=1
        sub_table(3,1,8,3)=1
        sub_table(3,1,8,4)=2

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=1
        sub_table(4,1,4,1)=3
        sub_table(4,1,4,2)=1
        sub_table(4,1,4,3)=2
        sub_table(4,1,4,4)=2
        sub_table(4,1,5,1)=1
        sub_table(4,1,5,2)=1
        sub_table(4,1,6,1)=1
        sub_table(4,1,6,2)=1
        sub_table(4,1,7,1)=1
        sub_table(4,1,7,2)=1
        sub_table(4,1,8,1)=3
        sub_table(4,1,8,2)=1
        sub_table(4,1,8,3)=2
        sub_table(4,1,8,4)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=3
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=2
        sub_table(5,1,4,1)=3
        sub_table(5,1,4,2)=1
        sub_table(5,1,4,3)=2
        sub_table(5,1,4,4)=3
        sub_table(5,1,5,1)=1
        sub_table(5,1,5,2)=1
        sub_table(5,1,6,1)=1
        sub_table(5,1,6,2)=3
        sub_table(5,1,7,1)=1
        sub_table(5,1,7,2)=2
        sub_table(5,1,8,1)=3
        sub_table(5,1,8,2)=1
        sub_table(5,1,8,3)=2
        sub_table(5,1,8,4)=3

        sub_table(8,1,1,1)=1
        sub_table(8,1,1,2)=1
        sub_table(8,1,2,1)=1
        sub_table(8,1,2,2)=1
        sub_table(8,1,3,1)=1
        sub_table(8,1,3,2)=1
        sub_table(8,1,4,1)=3
        sub_table(8,1,4,2)=2
        sub_table(8,1,4,3)=3
        sub_table(8,1,4,4)=4
        sub_table(8,1,5,1)=1
        sub_table(8,1,5,2)=1
        sub_table(8,1,6,1)=1
        sub_table(8,1,6,2)=1
        sub_table(8,1,7,1)=1
        sub_table(8,1,7,2)=1
        sub_table(8,1,8,1)=3
        sub_table(8,1,8,2)=2
        sub_table(8,1,8,3)=3
        sub_table(8,1,8,4)=4

        sub_table(12,1,1,1)=1
        sub_table(12,1,1,2)=1
        sub_table(12,1,2,1)=1
        sub_table(12,1,2,2)=1
        sub_table(12,1,3,1)=1
        sub_table(12,1,3,2)=1
        sub_table(12,1,4,1)=3
        sub_table(12,1,4,2)=2
        sub_table(12,1,4,3)=3
        sub_table(12,1,4,4)=4
        sub_table(12,1,5,1)=1
        sub_table(12,1,5,2)=2
        sub_table(12,1,6,1)=1
        sub_table(12,1,6,2)=2
        sub_table(12,1,7,1)=1
        sub_table(12,1,7,2)=2
        sub_table(12,1,8,1)=3
        sub_table(12,1,8,2)=1
        sub_table(12,1,8,3)=3
        sub_table(12,1,8,4)=4
 
        sub_table(16,1,1,1)=1
        sub_table(16,1,1,2)=1
        sub_table(16,1,2,1)=1
        sub_table(16,1,2,2)=1
        sub_table(16,1,3,1)=1
        sub_table(16,1,3,2)=1
        sub_table(16,1,4,1)=3
        sub_table(16,1,4,2)=1
        sub_table(16,1,4,3)=2
        sub_table(16,1,4,4)=2
        sub_table(16,1,5,1)=1
        sub_table(16,1,5,2)=3
        sub_table(16,1,6,1)=1
        sub_table(16,1,6,2)=3
        sub_table(16,1,7,1)=1
        sub_table(16,1,7,2)=3
        sub_table(16,1,8,1)=3
        sub_table(16,1,8,2)=3
        sub_table(16,1,8,3)=4
        sub_table(16,1,8,4)=4

        sub_table(20,1,1,1)=1
        sub_table(20,1,1,2)=1
        sub_table(20,1,2,1)=1
        sub_table(20,1,2,2)=1
        sub_table(20,1,3,1)=1
        sub_table(20,1,3,2)=1
        sub_table(20,1,4,1)=3
        sub_table(20,1,4,2)=2
        sub_table(20,1,4,3)=3
        sub_table(20,1,4,4)=4
        sub_table(20,1,5,1)=1
        sub_table(20,1,5,2)=5
        sub_table(20,1,6,1)=1
        sub_table(20,1,6,2)=5
        sub_table(20,1,7,1)=1
        sub_table(20,1,7,2)=5
        sub_table(20,1,8,1)=3
        sub_table(20,1,8,2)=6
        sub_table(20,1,8,3)=7
        sub_table(20,1,8,4)=8

        sub_table(27,1,1,1)=1
        sub_table(27,1,1,2)=1
        sub_table(27,1,2,1)=1
        sub_table(27,1,2,2)=3
        sub_table(27,1,3,1)=1
        sub_table(27,1,3,2)=2
        sub_table(27,1,4,1)=3
        sub_table(27,1,4,2)=1
        sub_table(27,1,4,3)=2
        sub_table(27,1,4,4)=3
        sub_table(27,1,5,1)=1
        sub_table(27,1,5,2)=4
        sub_table(27,1,6,1)=1
        sub_table(27,1,6,2)=6
        sub_table(27,1,7,1)=1
        sub_table(27,1,7,2)=5
        sub_table(27,1,8,1)=3
        sub_table(27,1,8,2)=4
        sub_table(27,1,8,3)=5
        sub_table(27,1,8,4)=6

        sub_table(28,1,1,1)=1
        sub_table(28,1,1,2)=1
        sub_table(28,1,2,1)=1
        sub_table(28,1,2,2)=2
        sub_table(28,1,3,1)=1
        sub_table(28,1,3,2)=3
        sub_table(28,1,4,1)=3
        sub_table(28,1,4,2)=4
        sub_table(28,1,4,3)=4
        sub_table(28,1,4,4)=4
        sub_table(28,1,5,1)=1
        sub_table(28,1,5,2)=1
        sub_table(28,1,6,1)=1
        sub_table(28,1,6,2)=2
        sub_table(28,1,7,1)=1
        sub_table(28,1,7,2)=3
        sub_table(28,1,8,1)=3
        sub_table(28,1,8,2)=4
        sub_table(28,1,8,3)=4
        sub_table(28,1,8,4)=4
 
     CASE(30)
!
!  T_d
!
        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=2
        sub_table(3,1,3,2)=1
        sub_table(3,1,3,3)=2
        sub_table(3,1,4,1)=3
        sub_table(3,1,4,2)=1
        sub_table(3,1,4,3)=1
        sub_table(3,1,4,4)=2
        sub_table(3,1,5,1)=3
        sub_table(3,1,5,2)=1
        sub_table(3,1,5,3)=2
        sub_table(3,1,5,4)=2

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=2
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=1
        sub_table(4,1,4,1)=3
        sub_table(4,1,4,2)=1
        sub_table(4,1,4,3)=2
        sub_table(4,1,4,4)=2
        sub_table(4,1,5,1)=3
        sub_table(4,1,5,2)=1
        sub_table(4,1,5,3)=2
        sub_table(4,1,5,4)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=1
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=2
        sub_table(5,1,3,3)=3
        sub_table(5,1,4,1)=3
        sub_table(5,1,4,2)=1
        sub_table(5,1,4,3)=2
        sub_table(5,1,4,4)=3
        sub_table(5,1,5,1)=3
        sub_table(5,1,5,2)=1
        sub_table(5,1,5,3)=2
        sub_table(5,1,5,4)=3

        sub_table(8,1,1,1)=1
        sub_table(8,1,1,2)=1
        sub_table(8,1,2,1)=1
        sub_table(8,1,2,2)=1
        sub_table(8,1,3,1)=2
        sub_table(8,1,3,2)=1
        sub_table(8,1,3,3)=1
        sub_table(8,1,4,1)=3
        sub_table(8,1,4,2)=2
        sub_table(8,1,4,3)=3
        sub_table(8,1,4,4)=4
        sub_table(8,1,5,1)=3
        sub_table(8,1,5,2)=2
        sub_table(8,1,5,3)=3
        sub_table(8,1,5,4)=4

        sub_table(12,1,1,1)=1
        sub_table(12,1,1,2)=1
        sub_table(12,1,2,1)=1
        sub_table(12,1,2,2)=2
        sub_table(12,1,3,1)=2
        sub_table(12,1,3,2)=1
        sub_table(12,1,3,3)=2
        sub_table(12,1,4,1)=3
        sub_table(12,1,4,2)=1
        sub_table(12,1,4,3)=2
        sub_table(12,1,4,4)=3
        sub_table(12,1,5,1)=3
        sub_table(12,1,5,2)=2
        sub_table(12,1,5,3)=3
        sub_table(12,1,5,4)=4

        sub_table(13,1,1,1)=1
        sub_table(13,1,1,2)=1
        sub_table(13,1,2,1)=1
        sub_table(13,1,2,2)=2
        sub_table(13,1,3,1)=2
        sub_table(13,1,3,2)=3
        sub_table(13,1,3,3)=3
        sub_table(13,1,4,1)=3
        sub_table(13,1,4,2)=1
        sub_table(13,1,4,3)=3
        sub_table(13,1,4,4)=3
        sub_table(13,1,5,1)=3
        sub_table(13,1,5,2)=2
        sub_table(13,1,5,3)=3
        sub_table(13,1,5,4)=3

        sub_table(24,1,1,1)=1
        sub_table(24,1,1,2)=1
        sub_table(24,1,2,1)=1
        sub_table(24,1,2,2)=3
        sub_table(24,1,3,1)=2
        sub_table(24,1,3,2)=1
        sub_table(24,1,3,3)=3
        sub_table(24,1,4,1)=3
        sub_table(24,1,4,2)=4
        sub_table(24,1,4,3)=5
        sub_table(24,1,4,4)=5
        sub_table(24,1,5,1)=3
        sub_table(24,1,5,2)=2
        sub_table(24,1,5,3)=5
        sub_table(24,1,5,4)=5

        sub_table(26,1,1,1)=1
        sub_table(26,1,1,2)=1
        sub_table(26,1,2,1)=1
        sub_table(26,1,2,2)=2
        sub_table(26,1,3,1)=2
        sub_table(26,1,3,2)=1
        sub_table(26,1,3,3)=2
        sub_table(26,1,4,1)=3
        sub_table(26,1,4,2)=2
        sub_table(26,1,4,3)=3
        sub_table(26,1,4,4)=4
        sub_table(26,1,5,1)=3
        sub_table(26,1,5,2)=1
        sub_table(26,1,5,3)=3
        sub_table(26,1,5,4)=4

        sub_table(28,1,1,1)=1
        sub_table(28,1,1,2)=1
        sub_table(28,1,2,1)=1
        sub_table(28,1,2,2)=1
        sub_table(28,1,3,1)=2
        sub_table(28,1,3,2)=2
        sub_table(28,1,3,3)=3
        sub_table(28,1,4,1)=3
        sub_table(28,1,4,2)=4
        sub_table(28,1,4,3)=4
        sub_table(28,1,4,4)=4
        sub_table(28,1,5,1)=3
        sub_table(28,1,5,2)=4
        sub_table(28,1,5,3)=4
        sub_table(28,1,5,4)=4

     CASE(31)
!
!  O
!
        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=2
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=1
        sub_table(4,1,4,1)=3
        sub_table(4,1,4,2)=1
        sub_table(4,1,4,3)=2
        sub_table(4,1,4,4)=2
        sub_table(4,1,5,1)=3
        sub_table(4,1,5,2)=1
        sub_table(4,1,5,3)=2
        sub_table(4,1,5,4)=2

        sub_table(4,2,1,1)=1
        sub_table(4,2,1,2)=1
        sub_table(4,2,2,1)=1
        sub_table(4,2,2,2)=2
        sub_table(4,2,3,1)=2
        sub_table(4,2,3,2)=1
        sub_table(4,2,3,3)=2
        sub_table(4,2,4,1)=3
        sub_table(4,2,4,2)=1
        sub_table(4,2,4,3)=2
        sub_table(4,2,4,4)=2
        sub_table(4,2,5,1)=3
        sub_table(4,2,5,2)=1
        sub_table(4,2,5,3)=1
        sub_table(4,2,5,4)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=1
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=2
        sub_table(5,1,3,3)=3
        sub_table(5,1,4,1)=3
        sub_table(5,1,4,2)=1
        sub_table(5,1,4,3)=2
        sub_table(5,1,4,4)=3
        sub_table(5,1,5,1)=3
        sub_table(5,1,5,2)=1
        sub_table(5,1,5,3)=2
        sub_table(5,1,5,4)=3

        sub_table(6,1,1,1)=1
        sub_table(6,1,1,2)=1
        sub_table(6,1,2,1)=1
        sub_table(6,1,2,2)=2
        sub_table(6,1,3,1)=2
        sub_table(6,1,3,2)=1
        sub_table(6,1,3,3)=2
        sub_table(6,1,4,1)=3
        sub_table(6,1,4,2)=1
        sub_table(6,1,4,3)=3
        sub_table(6,1,4,4)=4
        sub_table(6,1,5,1)=3
        sub_table(6,1,5,2)=2
        sub_table(6,1,5,3)=3
        sub_table(6,1,5,4)=4

        sub_table(8,1,1,1)=1
        sub_table(8,1,1,2)=1
        sub_table(8,1,2,1)=1
        sub_table(8,1,2,2)=1
        sub_table(8,1,3,1)=2
        sub_table(8,1,3,2)=1
        sub_table(8,1,3,3)=1
        sub_table(8,1,4,1)=3
        sub_table(8,1,4,2)=2
        sub_table(8,1,4,3)=3
        sub_table(8,1,4,4)=4
        sub_table(8,1,5,1)=3
        sub_table(8,1,5,2)=2
        sub_table(8,1,5,3)=3
        sub_table(8,1,5,4)=4

        sub_table(8,2,1,1)=1
        sub_table(8,2,1,2)=1
        sub_table(8,2,2,1)=1
        sub_table(8,2,2,2)=2
        sub_table(8,2,3,1)=2
        sub_table(8,2,3,2)=1
        sub_table(8,2,3,3)=2
        sub_table(8,2,4,1)=3
        sub_table(8,2,4,2)=2
        sub_table(8,2,4,3)=3
        sub_table(8,2,4,4)=4
        sub_table(8,2,5,1)=3
        sub_table(8,2,5,2)=1
        sub_table(8,2,5,3)=3
        sub_table(8,2,5,4)=4

        sub_table(9,1,1,1)=1
        sub_table(9,1,1,2)=1
        sub_table(9,1,2,1)=1
        sub_table(9,1,2,2)=2
        sub_table(9,1,3,1)=2
        sub_table(9,1,3,2)=3
        sub_table(9,1,3,3)=3
        sub_table(9,1,4,1)=3
        sub_table(9,1,4,2)=2
        sub_table(9,1,4,3)=3
        sub_table(9,1,4,4)=3
        sub_table(9,1,5,1)=3
        sub_table(9,1,5,2)=1
        sub_table(9,1,5,3)=3
        sub_table(9,1,5,4)=3

        sub_table(10,1,1,1)=1
        sub_table(10,1,1,2)=1
        sub_table(10,1,2,1)=1
        sub_table(10,1,2,2)=3
        sub_table(10,1,3,1)=2
        sub_table(10,1,3,2)=1
        sub_table(10,1,3,3)=3
        sub_table(10,1,4,1)=3
        sub_table(10,1,4,2)=2
        sub_table(10,1,4,3)=4
        sub_table(10,1,4,4)=4
        sub_table(10,1,5,1)=3
        sub_table(10,1,5,2)=3
        sub_table(10,1,5,3)=4
        sub_table(10,1,5,4)=4

        sub_table(28,1,1,1)=1
        sub_table(28,1,1,2)=1
        sub_table(28,1,2,1)=1
        sub_table(28,1,2,2)=1
        sub_table(28,1,3,1)=2
        sub_table(28,1,3,2)=2
        sub_table(28,1,3,3)=3
        sub_table(28,1,4,1)=3
        sub_table(28,1,4,2)=4
        sub_table(28,1,4,3)=4
        sub_table(28,1,4,4)=4
        sub_table(28,1,5,1)=3
        sub_table(28,1,5,2)=4
        sub_table(28,1,5,3)=4
        sub_table(28,1,5,4)=4

     CASE(32)
!
!  O_h
!
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=2
        sub_table(2,1,3,2)=1
        sub_table(2,1,3,3)=1
        sub_table(2,1,4,1)=3
        sub_table(2,1,4,2)=1
        sub_table(2,1,4,3)=1
        sub_table(2,1,4,4)=1
        sub_table(2,1,5,1)=3
        sub_table(2,1,5,2)=1
        sub_table(2,1,5,3)=1
        sub_table(2,1,5,4)=1
        sub_table(2,1,6,1)=1
        sub_table(2,1,6,2)=2
        sub_table(2,1,7,1)=1
        sub_table(2,1,7,2)=2
        sub_table(2,1,8,1)=2
        sub_table(2,1,8,2)=2
        sub_table(2,1,8,3)=2
        sub_table(2,1,9,1)=3
        sub_table(2,1,9,2)=2
        sub_table(2,1,9,3)=2
        sub_table(2,1,9,4)=2
        sub_table(2,1,10,1)=3
        sub_table(2,1,10,2)=2
        sub_table(2,1,10,3)=2
        sub_table(2,1,10,4)=2

        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=1
        sub_table(3,1,3,1)=2
        sub_table(3,1,3,2)=1
        sub_table(3,1,3,3)=1
        sub_table(3,1,4,1)=3
        sub_table(3,1,4,2)=1
        sub_table(3,1,4,3)=2
        sub_table(3,1,4,4)=2
        sub_table(3,1,5,1)=3
        sub_table(3,1,5,2)=1
        sub_table(3,1,5,3)=2
        sub_table(3,1,5,4)=2
        sub_table(3,1,6,1)=1
        sub_table(3,1,6,2)=2
        sub_table(3,1,7,1)=1
        sub_table(3,1,7,2)=2
        sub_table(3,1,8,1)=2
        sub_table(3,1,8,2)=2
        sub_table(3,1,8,3)=2
        sub_table(3,1,9,1)=3
        sub_table(3,1,9,2)=1
        sub_table(3,1,9,3)=1
        sub_table(3,1,9,4)=2
        sub_table(3,1,10,1)=3
        sub_table(3,1,10,2)=1
        sub_table(3,1,10,3)=1
        sub_table(3,1,10,4)=2

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=2
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=1
        sub_table(4,1,4,1)=3
        sub_table(4,1,4,2)=1
        sub_table(4,1,4,3)=2
        sub_table(4,1,4,4)=2
        sub_table(4,1,5,1)=3
        sub_table(4,1,5,2)=1
        sub_table(4,1,5,3)=2
        sub_table(4,1,5,4)=2
        sub_table(4,1,6,1)=1
        sub_table(4,1,6,2)=1
        sub_table(4,1,7,1)=1
        sub_table(4,1,7,2)=1
        sub_table(4,1,8,1)=2
        sub_table(4,1,8,2)=1
        sub_table(4,1,8,3)=1
        sub_table(4,1,9,1)=3
        sub_table(4,1,9,2)=1
        sub_table(4,1,9,3)=2
        sub_table(4,1,9,4)=2
        sub_table(4,1,10,1)=3
        sub_table(4,1,10,2)=1
        sub_table(4,1,10,3)=2
        sub_table(4,1,10,4)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=1
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=2
        sub_table(5,1,3,3)=3
        sub_table(5,1,4,1)=3
        sub_table(5,1,4,2)=1
        sub_table(5,1,4,3)=2
        sub_table(5,1,4,4)=3
        sub_table(5,1,5,1)=3
        sub_table(5,1,5,2)=1
        sub_table(5,1,5,3)=2
        sub_table(5,1,5,4)=3
        sub_table(5,1,6,1)=1
        sub_table(5,1,6,2)=1
        sub_table(5,1,7,1)=1
        sub_table(5,1,7,2)=1
        sub_table(5,1,8,1)=2
        sub_table(5,1,8,2)=2
        sub_table(5,1,8,3)=3
        sub_table(5,1,9,1)=3
        sub_table(5,1,9,2)=1
        sub_table(5,1,9,3)=2
        sub_table(5,1,9,4)=3
        sub_table(5,1,10,1)=3
        sub_table(5,1,10,2)=1
        sub_table(5,1,10,3)=2
        sub_table(5,1,10,4)=3

        sub_table(6,1,1,1)=1
        sub_table(6,1,1,2)=1
        sub_table(6,1,2,1)=1
        sub_table(6,1,2,2)=2
        sub_table(6,1,3,1)=2
        sub_table(6,1,3,2)=1
        sub_table(6,1,3,3)=2
        sub_table(6,1,4,1)=3
        sub_table(6,1,4,2)=1
        sub_table(6,1,4,3)=3
        sub_table(6,1,4,4)=4
        sub_table(6,1,5,1)=3
        sub_table(6,1,5,2)=2
        sub_table(6,1,5,3)=3
        sub_table(6,1,5,4)=4
        sub_table(6,1,6,1)=1
        sub_table(6,1,6,2)=1
        sub_table(6,1,7,1)=1
        sub_table(6,1,7,2)=2
        sub_table(6,1,8,1)=2
        sub_table(6,1,8,2)=1
        sub_table(6,1,8,3)=2
        sub_table(6,1,9,1)=3
        sub_table(6,1,9,2)=1
        sub_table(6,1,9,3)=3
        sub_table(6,1,9,4)=4
        sub_table(6,1,10,1)=3
        sub_table(6,1,10,2)=2
        sub_table(6,1,10,3)=3
        sub_table(6,1,10,4)=4

        sub_table(8,1,1,1)=1
        sub_table(8,1,1,2)=1
        sub_table(8,1,2,1)=1
        sub_table(8,1,2,2)=1
        sub_table(8,1,3,1)=2
        sub_table(8,1,3,2)=1
        sub_table(8,1,3,3)=1
        sub_table(8,1,4,1)=3
        sub_table(8,1,4,2)=2
        sub_table(8,1,4,3)=3
        sub_table(8,1,4,4)=4
        sub_table(8,1,5,1)=3
        sub_table(8,1,5,2)=2
        sub_table(8,1,5,3)=3
        sub_table(8,1,5,4)=4
        sub_table(8,1,6,1)=1
        sub_table(8,1,6,2)=1
        sub_table(8,1,7,1)=1
        sub_table(8,1,7,2)=1
        sub_table(8,1,8,1)=2
        sub_table(8,1,8,2)=1
        sub_table(8,1,8,3)=1
        sub_table(8,1,9,1)=3
        sub_table(8,1,9,2)=2
        sub_table(8,1,9,3)=3
        sub_table(8,1,9,4)=4
        sub_table(8,1,10,1)=3
        sub_table(8,1,10,2)=2
        sub_table(8,1,10,3)=3
        sub_table(8,1,10,4)=4

        sub_table(9,1,1,1)=1
        sub_table(9,1,1,2)=1
        sub_table(9,1,2,1)=1
        sub_table(9,1,2,2)=2
        sub_table(9,1,3,1)=2
        sub_table(9,1,3,2)=3
        sub_table(9,1,3,3)=3
        sub_table(9,1,4,1)=3
        sub_table(9,1,4,2)=2
        sub_table(9,1,4,3)=3
        sub_table(9,1,4,4)=3
        sub_table(9,1,5,1)=3
        sub_table(9,1,5,2)=1
        sub_table(9,1,5,3)=3
        sub_table(9,1,5,4)=3
        sub_table(9,1,6,1)=1
        sub_table(9,1,6,2)=1
        sub_table(9,1,7,1)=1
        sub_table(9,1,7,2)=2
        sub_table(9,1,8,1)=2
        sub_table(9,1,8,2)=3
        sub_table(9,1,8,3)=3
        sub_table(9,1,9,1)=3
        sub_table(9,1,9,2)=2
        sub_table(9,1,9,3)=3
        sub_table(9,1,9,4)=3
        sub_table(9,1,10,1)=3
        sub_table(9,1,10,2)=1
        sub_table(9,1,10,3)=3
        sub_table(9,1,10,4)=3

        sub_table(10,1,1,1)=1
        sub_table(10,1,1,2)=1
        sub_table(10,1,2,1)=1
        sub_table(10,1,2,2)=3
        sub_table(10,1,3,1)=2
        sub_table(10,1,3,2)=1
        sub_table(10,1,3,3)=3
        sub_table(10,1,4,1)=3
        sub_table(10,1,4,2)=2
        sub_table(10,1,4,3)=5
        sub_table(10,1,4,4)=5
        sub_table(10,1,5,1)=3
        sub_table(10,1,5,2)=4
        sub_table(10,1,5,3)=5
        sub_table(10,1,5,4)=5
        sub_table(10,1,6,1)=1
        sub_table(10,1,6,2)=1
        sub_table(10,1,7,1)=1
        sub_table(10,1,7,2)=3
        sub_table(10,1,8,1)=2
        sub_table(10,1,8,2)=1
        sub_table(10,1,8,3)=3
        sub_table(10,1,9,1)=3
        sub_table(10,1,9,2)=2
        sub_table(10,1,9,3)=5
        sub_table(10,1,9,4)=5
        sub_table(10,1,10,1)=3
        sub_table(10,1,10,2)=4
        sub_table(10,1,10,3)=5
        sub_table(10,1,10,4)=5

        sub_table(12,1,1,1)=1
        sub_table(12,1,1,2)=1
        sub_table(12,1,2,1)=1
        sub_table(12,1,2,2)=1
        sub_table(12,1,3,1)=2
        sub_table(12,1,3,2)=1
        sub_table(12,1,3,3)=1
        sub_table(12,1,4,1)=3
        sub_table(12,1,4,2)=2
        sub_table(12,1,4,3)=3
        sub_table(12,1,4,4)=4
        sub_table(12,1,5,1)=3
        sub_table(12,1,5,2)=2
        sub_table(12,1,5,3)=3
        sub_table(12,1,5,4)=4
        sub_table(12,1,6,1)=1
        sub_table(12,1,6,2)=2
        sub_table(12,1,7,1)=1
        sub_table(12,1,7,2)=2
        sub_table(12,1,8,1)=2
        sub_table(12,1,8,2)=2
        sub_table(12,1,8,3)=2
        sub_table(12,1,9,1)=3
        sub_table(12,1,9,2)=1
        sub_table(12,1,9,3)=3
        sub_table(12,1,9,4)=4
        sub_table(12,1,10,1)=3
        sub_table(12,1,10,2)=1
        sub_table(12,1,10,3)=3
        sub_table(12,1,10,4)=4

        sub_table(12,2,1,1)=1
        sub_table(12,2,1,2)=1
        sub_table(12,2,2,1)=1
        sub_table(12,2,2,2)=2
        sub_table(12,2,3,1)=2
        sub_table(12,2,3,2)=1
        sub_table(12,2,3,3)=2
        sub_table(12,2,4,1)=3
        sub_table(12,2,4,2)=2
        sub_table(12,2,4,3)=3
        sub_table(12,2,4,4)=4
        sub_table(12,2,5,1)=3
        sub_table(12,2,5,2)=2
        sub_table(12,2,5,3)=3
        sub_table(12,2,5,4)=4
        sub_table(12,2,6,1)=1
        sub_table(12,2,6,2)=2
        sub_table(12,2,7,1)=1
        sub_table(12,2,7,2)=1
        sub_table(12,2,8,1)=2
        sub_table(12,2,8,2)=1
        sub_table(12,2,8,3)=2
        sub_table(12,2,9,1)=3
        sub_table(12,2,9,2)=1
        sub_table(12,2,9,3)=3
        sub_table(12,2,9,4)=4
        sub_table(12,2,10,1)=3
        sub_table(12,2,10,2)=1
        sub_table(12,2,10,3)=3
        sub_table(12,2,10,4)=4

        sub_table(12,3,1,1)=1
        sub_table(12,3,1,2)=1
        sub_table(12,3,2,1)=1
        sub_table(12,3,2,2)=3
        sub_table(12,3,3,1)=2
        sub_table(12,3,3,2)=1
        sub_table(12,3,3,3)=3
        sub_table(12,3,4,1)=3
        sub_table(12,3,4,2)=2
        sub_table(12,3,4,3)=3
        sub_table(12,3,4,4)=4
        sub_table(12,3,5,1)=3
        sub_table(12,3,5,2)=1
        sub_table(12,3,5,3)=2
        sub_table(12,3,5,4)=4
        sub_table(12,3,6,1)=1
        sub_table(12,3,6,2)=2
        sub_table(12,3,7,1)=1
        sub_table(12,3,7,2)=4
        sub_table(12,3,8,1)=2
        sub_table(12,3,8,2)=2
        sub_table(12,3,8,3)=4
        sub_table(12,3,9,1)=3
        sub_table(12,3,9,2)=1
        sub_table(12,3,9,3)=3
        sub_table(12,3,9,4)=4
        sub_table(12,3,10,1)=3
        sub_table(12,3,10,2)=1
        sub_table(12,3,10,3)=2
        sub_table(12,3,10,4)=3

        sub_table(13,1,1,1)=1
        sub_table(13,1,1,2)=1
        sub_table(13,1,2,1)=1
        sub_table(13,1,2,2)=2
        sub_table(13,1,3,1)=2
        sub_table(13,1,3,2)=3
        sub_table(13,1,3,3)=3
        sub_table(13,1,4,1)=3
        sub_table(13,1,4,2)=2
        sub_table(13,1,4,3)=3
        sub_table(13,1,4,4)=3
        sub_table(13,1,5,1)=3
        sub_table(13,1,5,2)=1
        sub_table(13,1,5,3)=3
        sub_table(13,1,5,4)=3
        sub_table(13,1,6,1)=1
        sub_table(13,1,6,2)=2
        sub_table(13,1,7,1)=1
        sub_table(13,1,7,2)=1
        sub_table(13,1,8,1)=2
        sub_table(13,1,8,2)=3
        sub_table(13,1,8,3)=3
        sub_table(13,1,9,1)=3
        sub_table(13,1,9,2)=1
        sub_table(13,1,9,3)=3
        sub_table(13,1,9,4)=3
        sub_table(13,1,10,1)=3
        sub_table(13,1,10,2)=2
        sub_table(13,1,10,3)=3
        sub_table(13,1,10,4)=3

        sub_table(14,1,1,1)=1
        sub_table(14,1,1,2)=1
        sub_table(14,1,2,1)=1
        sub_table(14,1,2,2)=3
        sub_table(14,1,3,1)=2
        sub_table(14,1,3,2)=1
        sub_table(14,1,3,3)=3
        sub_table(14,1,4,1)=3
        sub_table(14,1,4,2)=2
        sub_table(14,1,4,3)=5
        sub_table(14,1,4,4)=5
        sub_table(14,1,5,1)=3
        sub_table(14,1,5,2)=4
        sub_table(14,1,5,3)=5
        sub_table(14,1,5,4)=5
        sub_table(14,1,6,1)=1
        sub_table(14,1,6,2)=2
        sub_table(14,1,7,1)=1
        sub_table(14,1,7,2)=4
        sub_table(14,1,8,1)=2
        sub_table(14,1,8,2)=2
        sub_table(14,1,8,3)=4
        sub_table(14,1,9,1)=3
        sub_table(14,1,9,2)=1
        sub_table(14,1,9,3)=5
        sub_table(14,1,9,4)=5
        sub_table(14,1,10,1)=3
        sub_table(14,1,10,2)=3
        sub_table(14,1,10,3)=5
        sub_table(14,1,10,4)=5

        sub_table(16,1,1,1)=1
        sub_table(16,1,1,2)=1
        sub_table(16,1,2,1)=1
        sub_table(16,1,2,2)=1
        sub_table(16,1,3,1)=2
        sub_table(16,1,3,2)=1
        sub_table(16,1,3,3)=1
        sub_table(16,1,4,1)=3
        sub_table(16,1,4,2)=1
        sub_table(16,1,4,3)=2
        sub_table(16,1,4,4)=2
        sub_table(16,1,5,1)=3
        sub_table(16,1,5,2)=1
        sub_table(16,1,5,3)=2
        sub_table(16,1,5,4)=2
        sub_table(16,1,6,1)=1
        sub_table(16,1,6,2)=3
        sub_table(16,1,7,1)=1
        sub_table(16,1,7,2)=3
        sub_table(16,1,8,1)=2
        sub_table(16,1,8,2)=3
        sub_table(16,1,8,3)=3
        sub_table(16,1,9,1)=3
        sub_table(16,1,9,2)=3
        sub_table(16,1,9,3)=4
        sub_table(16,1,9,4)=4
        sub_table(16,1,10,1)=3
        sub_table(16,1,10,2)=3
        sub_table(16,1,10,3)=4
        sub_table(16,1,10,4)=4

        sub_table(16,2,1,1)=1
        sub_table(16,2,1,2)=1
        sub_table(16,2,2,1)=1
        sub_table(16,2,2,2)=2
        sub_table(16,2,3,1)=2
        sub_table(16,2,3,2)=1
        sub_table(16,2,3,3)=2
        sub_table(16,2,4,1)=3
        sub_table(16,2,4,2)=1
        sub_table(16,2,4,3)=2
        sub_table(16,2,4,4)=2
        sub_table(16,2,5,1)=3
        sub_table(16,2,5,2)=1
        sub_table(16,2,5,3)=1
        sub_table(16,2,5,4)=2
        sub_table(16,2,6,1)=1
        sub_table(16,2,6,2)=3
        sub_table(16,2,7,1)=1
        sub_table(16,2,7,2)=4
        sub_table(16,2,8,1)=2
        sub_table(16,2,8,2)=3
        sub_table(16,2,8,3)=4
        sub_table(16,2,9,1)=3
        sub_table(16,2,9,2)=3
        sub_table(16,2,9,3)=4
        sub_table(16,2,9,4)=4
        sub_table(16,2,10,1)=3
        sub_table(16,2,10,2)=3
        sub_table(16,2,10,3)=3
        sub_table(16,2,10,4)=4

        sub_table(18,1,1,1)=1
        sub_table(18,1,1,2)=1
        sub_table(18,1,2,1)=1
        sub_table(18,1,2,2)=2
        sub_table(18,1,3,1)=2
        sub_table(18,1,3,2)=1
        sub_table(18,1,3,3)=2
        sub_table(18,1,4,1)=3
        sub_table(18,1,4,2)=1
        sub_table(18,1,4,3)=3
        sub_table(18,1,4,4)=4
        sub_table(18,1,5,1)=3
        sub_table(18,1,5,2)=2
        sub_table(18,1,5,3)=3
        sub_table(18,1,5,4)=4
        sub_table(18,1,6,1)=1
        sub_table(18,1,6,2)=5
        sub_table(18,1,7,1)=1
        sub_table(18,1,7,2)=6
        sub_table(18,1,8,1)=2
        sub_table(18,1,8,2)=5
        sub_table(18,1,8,3)=6
        sub_table(18,1,9,1)=3
        sub_table(18,1,9,2)=5
        sub_table(18,1,9,3)=7
        sub_table(18,1,9,4)=8
        sub_table(18,1,10,1)=3
        sub_table(18,1,10,2)=6
        sub_table(18,1,10,3)=7
        sub_table(18,1,10,4)=8

        sub_table(20,1,1,1)=1
        sub_table(20,1,1,2)=1
        sub_table(20,1,2,1)=1
        sub_table(20,1,2,2)=1
        sub_table(20,1,3,1)=2
        sub_table(20,1,3,2)=1
        sub_table(20,1,3,3)=1
        sub_table(20,1,4,1)=3
        sub_table(20,1,4,2)=2
        sub_table(20,1,4,3)=3
        sub_table(20,1,4,4)=4
        sub_table(20,1,5,1)=3
        sub_table(20,1,5,2)=2
        sub_table(20,1,5,3)=3
        sub_table(20,1,5,4)=4
        sub_table(20,1,6,1)=1
        sub_table(20,1,6,2)=5
        sub_table(20,1,7,1)=1
        sub_table(20,1,7,2)=5
        sub_table(20,1,8,1)=2
        sub_table(20,1,8,2)=5
        sub_table(20,1,8,3)=5
        sub_table(20,1,9,1)=3
        sub_table(20,1,9,2)=6
        sub_table(20,1,9,3)=7
        sub_table(20,1,9,4)=8
        sub_table(20,1,10,1)=3
        sub_table(20,1,10,2)=6
        sub_table(20,1,10,3)=7
        sub_table(20,1,10,4)=8

        sub_table(20,2,1,1)=1
        sub_table(20,2,1,2)=1
        sub_table(20,2,2,1)=1
        sub_table(20,2,2,2)=2
        sub_table(20,2,3,1)=2
        sub_table(20,2,3,2)=1
        sub_table(20,2,3,3)=2
        sub_table(20,2,4,1)=3
        sub_table(20,2,4,2)=3
        sub_table(20,2,4,3)=4
        sub_table(20,2,4,4)=5
        sub_table(20,2,5,1)=3
        sub_table(20,2,5,2)=3
        sub_table(20,2,5,3)=4
        sub_table(20,2,5,4)=5
        sub_table(20,2,6,1)=1
        sub_table(20,2,6,2)=6
        sub_table(20,2,7,1)=1
        sub_table(20,2,7,2)=7
        sub_table(20,2,8,1)=2
        sub_table(20,2,8,2)=6
        sub_table(20,2,8,3)=7
        sub_table(20,2,9,1)=3
        sub_table(20,2,9,2)=8
        sub_table(20,2,9,3)=9
        sub_table(20,2,9,4)=10
        sub_table(20,2,10,1)=3
        sub_table(20,2,10,2)=8
        sub_table(20,2,10,3)=9
        sub_table(20,2,10,4)=10

        sub_table(22,1,1,1)=1
        sub_table(22,1,1,2)=1
        sub_table(22,1,2,1)=1
        sub_table(22,1,2,2)=3
        sub_table(22,1,3,1)=2
        sub_table(22,1,3,2)=1
        sub_table(22,1,3,3)=3
        sub_table(22,1,4,1)=3
        sub_table(22,1,4,2)=2
        sub_table(22,1,4,3)=5
        sub_table(22,1,4,4)=5
        sub_table(22,1,5,1)=3
        sub_table(22,1,5,2)=4
        sub_table(22,1,5,3)=5
        sub_table(22,1,5,4)=5
        sub_table(22,1,6,1)=1
        sub_table(22,1,6,2)=6
        sub_table(22,1,7,1)=1
        sub_table(22,1,7,2)=8
        sub_table(22,1,8,1)=2
        sub_table(22,1,8,2)=6
        sub_table(22,1,8,3)=8
        sub_table(22,1,9,1)=3
        sub_table(22,1,9,2)=7
        sub_table(22,1,9,3)=10
        sub_table(22,1,9,4)=10
        sub_table(22,1,10,1)=3
        sub_table(22,1,10,2)=9
        sub_table(22,1,10,3)=10
        sub_table(22,1,10,4)=10

        sub_table(24,1,1,1)=1
        sub_table(24,1,1,2)=1
        sub_table(24,1,2,1)=1
        sub_table(24,1,2,2)=3
        sub_table(24,1,3,1)=2
        sub_table(24,1,3,2)=1
        sub_table(24,1,3,3)=3
        sub_table(24,1,4,1)=3
        sub_table(24,1,4,2)=2
        sub_table(24,1,4,3)=5
        sub_table(24,1,4,4)=5
        sub_table(24,1,5,1)=3
        sub_table(24,1,5,2)=4
        sub_table(24,1,5,3)=5
        sub_table(24,1,5,4)=5
        sub_table(24,1,6,1)=1
        sub_table(24,1,6,2)=3
        sub_table(24,1,7,1)=1
        sub_table(24,1,7,2)=1
        sub_table(24,1,8,1)=2
        sub_table(24,1,8,2)=1
        sub_table(24,1,8,3)=3
        sub_table(24,1,9,1)=3
        sub_table(24,1,9,2)=4
        sub_table(24,1,9,3)=5
        sub_table(24,1,9,4)=5
        sub_table(24,1,10,1)=3
        sub_table(24,1,10,2)=2
        sub_table(24,1,10,3)=5
        sub_table(24,1,10,4)=5

        sub_table(25,1,1,1)=1
        sub_table(25,1,1,2)=1
        sub_table(25,1,2,1)=1
        sub_table(25,1,2,2)=2
        sub_table(25,1,3,1)=2
        sub_table(25,1,3,2)=3
        sub_table(25,1,3,3)=3
        sub_table(25,1,4,1)=3
        sub_table(25,1,4,2)=2
        sub_table(25,1,4,3)=3
        sub_table(25,1,4,4)=3
        sub_table(25,1,5,1)=3
        sub_table(25,1,5,2)=1
        sub_table(25,1,5,3)=3
        sub_table(25,1,5,4)=3
        sub_table(25,1,6,1)=1
        sub_table(25,1,6,2)=4
        sub_table(25,1,7,1)=1
        sub_table(25,1,7,2)=5
        sub_table(25,1,8,1)=2
        sub_table(25,1,8,2)=6
        sub_table(25,1,8,3)=6
        sub_table(25,1,9,1)=3
        sub_table(25,1,9,2)=5
        sub_table(25,1,9,3)=6
        sub_table(25,1,9,4)=6
        sub_table(25,1,10,1)=3
        sub_table(25,1,10,2)=4
        sub_table(25,1,10,3)=6
        sub_table(25,1,10,4)=6

        sub_table(26,1,1,1)=1
        sub_table(26,1,1,2)=1
        sub_table(26,1,2,1)=1
        sub_table(26,1,2,2)=2
        sub_table(26,1,3,1)=2
        sub_table(26,1,3,2)=1
        sub_table(26,1,3,3)=2
        sub_table(26,1,4,1)=3
        sub_table(26,1,4,2)=1
        sub_table(26,1,4,3)=3
        sub_table(26,1,4,4)=4
        sub_table(26,1,5,1)=3
        sub_table(26,1,5,2)=2
        sub_table(26,1,5,3)=3
        sub_table(26,1,5,4)=4
        sub_table(26,1,6,1)=1
        sub_table(26,1,6,2)=2
        sub_table(26,1,7,1)=1
        sub_table(26,1,7,2)=1
        sub_table(26,1,8,1)=2
        sub_table(26,1,8,2)=1
        sub_table(26,1,8,3)=2
        sub_table(26,1,9,1)=3
        sub_table(26,1,9,2)=2
        sub_table(26,1,9,3)=3
        sub_table(26,1,9,4)=4
        sub_table(26,1,10,1)=3
        sub_table(26,1,10,2)=1
        sub_table(26,1,10,3)=3
        sub_table(26,1,10,4)=4

        sub_table(27,1,1,1)=1
        sub_table(27,1,1,2)=1
        sub_table(27,1,2,1)=1
        sub_table(27,1,2,2)=1
        sub_table(27,1,3,1)=2
        sub_table(27,1,3,2)=2
        sub_table(27,1,3,3)=3
        sub_table(27,1,4,1)=3
        sub_table(27,1,4,2)=1
        sub_table(27,1,4,3)=2
        sub_table(27,1,4,4)=3
        sub_table(27,1,5,1)=3
        sub_table(27,1,5,2)=1
        sub_table(27,1,5,3)=2
        sub_table(27,1,5,4)=3
        sub_table(27,1,6,1)=1
        sub_table(27,1,6,2)=4
        sub_table(27,1,7,1)=1
        sub_table(27,1,7,2)=4
        sub_table(27,1,8,1)=2
        sub_table(27,1,8,2)=5
        sub_table(27,1,8,3)=6
        sub_table(27,1,9,1)=3
        sub_table(27,1,9,2)=4
        sub_table(27,1,9,3)=5
        sub_table(27,1,9,4)=6
        sub_table(27,1,10,1)=3
        sub_table(27,1,10,2)=4
        sub_table(27,1,10,3)=5
        sub_table(27,1,10,4)=6

        sub_table(28,1,1,1)=1
        sub_table(28,1,1,2)=1
        sub_table(28,1,2,1)=1
        sub_table(28,1,2,2)=1
        sub_table(28,1,3,1)=2
        sub_table(28,1,3,2)=2
        sub_table(28,1,3,3)=3
        sub_table(28,1,4,1)=3
        sub_table(28,1,4,2)=4
        sub_table(28,1,4,3)=4
        sub_table(28,1,4,4)=4
        sub_table(28,1,5,1)=3
        sub_table(28,1,5,2)=4
        sub_table(28,1,5,3)=4
        sub_table(28,1,5,4)=4
        sub_table(28,1,6,1)=1
        sub_table(28,1,6,2)=1
        sub_table(28,1,7,1)=1
        sub_table(28,1,7,2)=1
        sub_table(28,1,8,1)=2
        sub_table(28,1,8,2)=2
        sub_table(28,1,8,3)=3
        sub_table(28,1,9,1)=3
        sub_table(28,1,9,2)=4
        sub_table(28,1,9,3)=4
        sub_table(28,1,9,4)=4
        sub_table(28,1,10,1)=3
        sub_table(28,1,10,2)=4
        sub_table(28,1,10,3)=4
        sub_table(28,1,10,4)=4

        sub_table(29,1,1,1)=1
        sub_table(29,1,1,2)=1
        sub_table(29,1,2,1)=1
        sub_table(29,1,2,2)=1
        sub_table(29,1,3,1)=2
        sub_table(29,1,3,2)=2
        sub_table(29,1,3,3)=3
        sub_table(29,1,4,1)=3
        sub_table(29,1,4,2)=4
        sub_table(29,1,4,3)=4
        sub_table(29,1,4,4)=4
        sub_table(29,1,5,1)=3
        sub_table(29,1,5,2)=4
        sub_table(29,1,5,3)=4
        sub_table(29,1,5,4)=4
        sub_table(29,1,6,1)=1
        sub_table(29,1,6,2)=5
        sub_table(29,1,7,1)=1
        sub_table(29,1,7,2)=5
        sub_table(29,1,8,1)=2
        sub_table(29,1,8,2)=6
        sub_table(29,1,8,3)=7
        sub_table(29,1,9,1)=3
        sub_table(29,1,9,2)=8
        sub_table(29,1,9,3)=8
        sub_table(29,1,9,4)=8
        sub_table(29,1,10,1)=3
        sub_table(29,1,10,2)=8
        sub_table(29,1,10,3)=8
        sub_table(29,1,10,4)=8

        sub_table(30,1,1,1)=1
        sub_table(30,1,1,2)=1
        sub_table(30,1,2,1)=1
        sub_table(30,1,2,2)=2
        sub_table(30,1,3,1)=2
        sub_table(30,1,3,2)=2
        sub_table(30,1,3,3)=3
        sub_table(30,1,4,1)=3
        sub_table(30,1,4,2)=4
        sub_table(30,1,4,3)=4
        sub_table(30,1,4,4)=4
        sub_table(30,1,5,1)=3
        sub_table(30,1,5,2)=5
        sub_table(30,1,5,3)=5
        sub_table(30,1,5,4)=5
        sub_table(30,1,6,1)=1
        sub_table(30,1,6,2)=2
        sub_table(30,1,7,1)=1
        sub_table(30,1,7,2)=1
        sub_table(30,1,8,1)=2
        sub_table(30,1,8,2)=3
        sub_table(30,1,8,3)=3
        sub_table(30,1,9,1)=3
        sub_table(30,1,9,2)=5
        sub_table(30,1,9,3)=5
        sub_table(30,1,9,4)=5
        sub_table(30,1,10,1)=3
        sub_table(30,1,10,2)=4
        sub_table(30,1,10,3)=4
        sub_table(30,1,10,4)=4

        sub_table(31,1,1,1)=1
        sub_table(31,1,1,2)=1
        sub_table(31,1,2,1)=1
        sub_table(31,1,2,2)=2
        sub_table(31,1,3,1)=2
        sub_table(31,1,3,2)=3
        sub_table(31,1,3,3)=3
        sub_table(31,1,4,1)=3
        sub_table(31,1,4,2)=4
        sub_table(31,1,4,3)=4
        sub_table(31,1,4,4)=4
        sub_table(31,1,5,1)=3
        sub_table(31,1,5,2)=5
        sub_table(31,1,5,3)=5
        sub_table(31,1,5,4)=5
        sub_table(31,1,6,1)=1
        sub_table(31,1,6,2)=1
        sub_table(31,1,7,1)=1
        sub_table(31,1,7,2)=2
        sub_table(31,1,8,1)=2
        sub_table(31,1,8,2)=3
        sub_table(31,1,8,3)=3
        sub_table(31,1,9,1)=3
        sub_table(31,1,9,2)=4
        sub_table(31,1,9,3)=4
        sub_table(31,1,9,4)=4
        sub_table(31,1,10,1)=3
        sub_table(31,1,10,2)=5
        sub_table(31,1,10,3)=5
        sub_table(31,1,10,4)=5

     CASE DEFAULT 
        CALL errore('convert_one_rap','Input point group uncorrect',1)
  END SELECT

  ndeg = sub_table(group_out, aux_ind, rap, 1)
  IF (ndeg==0) THEN
     WRITE(stdout,'("group_in, group_out, representation",4i5)') group_in, &
                                                group_out, rap, aux_ind
     CALL errore('convert_one_rap','problem representation not found',1)
  END IF
  DO ideg=1, ndeg
     rap_list(ideg) = sub_table(group_out, aux_ind, rap, ideg+1)
  ENDDO

  RETURN

  END SUBROUTINE convert_one_rap
!
!  
  SUBROUTINE convert_one_rap_so(rap, ndeg, rap_list, group_in, group_out)
!
!  This routine sets the subduction table for the group subgroup relationship.
!  This subduction table is organized like this. It is set for each group_in,
!  for all the possibile subgroups. The first index is group_out. In the
!  relativistic case aux_ind (see the scalar relativistic routine) is not
!  necessary. For compatibility with the scalar relativistic routine we keep
!  the index, that however is always one.
!  The third index is rap, and the fourth index contains in the first position
!  the degeneracy of the rappresentation (1, 2, 3, 4) and in the four
!  following positions the indices of the representations.
!  The first part of the routine set the information for all
!  the representations of group_in,
!  and the final instructions copy in ndeg and rap_list only the
!  information for the required representation.
!  The representations numbers are syncronized with those defined in
!  the routine set_irr_so (in PW/src/divide_class_so.f90).
!
!  This routine must be used in the fully relativistic case with spin-orbit,
!  for the corresponding scalar relativistic routine use convert_one_rap.
!


  USE io_global, ONLY : stdout

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: group_in, group_out, rap
  INTEGER, INTENT(OUT) :: ndeg, rap_list(4)
  INTEGER :: sub_table(32, 1, 12, 5)
  INTEGER :: ideg

  sub_table=0
!
!  C_1  has only representation 1
!
  sub_table(1,:,:,1)=1
  sub_table(1,:,:,2)=1

  SELECT CASE (group_in) 
!
!    
!
     CASE(1,2,3,4,5)

     CASE(6)
!
!  C_4
!
        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=2
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=1
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2

     CASE(7)
!
!  C_6
!
        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=2
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=1
        sub_table(4,1,5,1)=1
        sub_table(4,1,5,2)=1
        sub_table(4,1,6,1)=1
        sub_table(4,1,6,2)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=2
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=1
        sub_table(5,1,4,1)=1
        sub_table(5,1,4,2)=2
        sub_table(5,1,5,1)=1
        sub_table(5,1,5,2)=3
        sub_table(5,1,6,1)=1
        sub_table(5,1,6,2)=3

     CASE(8)
!
! D_2
!
        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2

     CASE(9)
!
! D_3
!
        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2

        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=3
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=3


     CASE(10)
!
!  D_4
!
        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2

        sub_table(6,1,1,1)=2
        sub_table(6,1,1,2)=1
        sub_table(6,1,1,3)=2
        sub_table(6,1,2,1)=2
        sub_table(6,1,2,2)=3
        sub_table(6,1,2,3)=4

        sub_table(8,1,1,1)=2
        sub_table(8,1,1,2)=1
        sub_table(8,1,1,3)=1
        sub_table(8,1,2,1)=2
        sub_table(8,1,2,2)=1
        sub_table(8,1,2,3)=1


     CASE(11)
!
!  D_6
!
        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2
        sub_table(4,1,3,1)=2
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=2

        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=2
        sub_table(5,1,2,2)=1
        sub_table(5,1,2,3)=2
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=3
        sub_table(5,1,3,3)=3

        sub_table(7,1,1,1)=2
        sub_table(7,1,1,2)=1
        sub_table(7,1,1,3)=2
        sub_table(7,1,2,1)=2
        sub_table(7,1,2,2)=5
        sub_table(7,1,2,3)=6
        sub_table(7,1,3,1)=2
        sub_table(7,1,3,2)=3
        sub_table(7,1,3,3)=4

        sub_table(9,1,1,1)=2
        sub_table(9,1,1,2)=1
        sub_table(9,1,1,3)=2
        sub_table(9,1,2,1)=2
        sub_table(9,1,2,2)=1
        sub_table(9,1,2,3)=2
        sub_table(9,1,3,1)=2
        sub_table(9,1,3,2)=3
        sub_table(9,1,3,3)=3

     CASE(12)
!
!   C_2v
!
        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        
        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        
     CASE(13)
!
!   C_3v
!
        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=1
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=2

        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=3
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=3


     CASE(14)
!
!   C_4v
!
        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=2
        sub_table(3,1,2,2)=1
        sub_table(3,1,2,3)=2

        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2

        sub_table(6,1,1,1)=2
        sub_table(6,1,1,2)=1
        sub_table(6,1,1,3)=2
        sub_table(6,1,2,1)=2
        sub_table(6,1,2,2)=3
        sub_table(6,1,2,3)=4

        sub_table(12,1,1,1)=2
        sub_table(12,1,1,2)=1
        sub_table(12,1,1,3)=1
        sub_table(12,1,2,1)=2
        sub_table(12,1,2,2)=1
        sub_table(12,1,2,3)=1

     CASE(15)
!
!   C_6v
!
        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=2
        sub_table(3,1,2,2)=1
        sub_table(3,1,2,3)=2
        sub_table(3,1,3,1)=2
        sub_table(3,1,3,2)=1
        sub_table(3,1,3,3)=2

        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2
        sub_table(4,1,3,1)=2
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=2

        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=2
        sub_table(5,1,2,2)=1
        sub_table(5,1,2,3)=2
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=3
        sub_table(5,1,3,3)=3
 
        sub_table(7,1,1,1)=2
        sub_table(7,1,1,2)=1
        sub_table(7,1,1,3)=2
        sub_table(7,1,2,1)=2
        sub_table(7,1,2,2)=5
        sub_table(7,1,2,3)=6
        sub_table(7,1,3,1)=2
        sub_table(7,1,3,2)=3
        sub_table(7,1,3,3)=4

        sub_table(12,1,1,1)=2
        sub_table(12,1,1,2)=1
        sub_table(12,1,1,3)=1
        sub_table(12,1,2,1)=2
        sub_table(12,1,2,2)=1
        sub_table(12,1,2,3)=1
        sub_table(12,1,3,1)=2
        sub_table(12,1,3,2)=1
        sub_table(12,1,3,3)=1

        sub_table(13,1,1,1)=2
        sub_table(13,1,1,2)=1
        sub_table(13,1,1,3)=1
        sub_table(13,1,2,1)=2
        sub_table(13,1,2,2)=1
        sub_table(13,1,2,3)=1
        sub_table(13,1,3,1)=2
        sub_table(13,1,3,2)=2
        sub_table(13,1,3,3)=3

     CASE(16)
!
! C_2h
!
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=1
        sub_table(2,1,3,2)=2
        sub_table(2,1,4,1)=1
        sub_table(2,1,4,2)=2

        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=2
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=1

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=2
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=1
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2

     CASE(17)
!
! C_3h
!
        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=2
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=1
        sub_table(3,1,5,1)=1
        sub_table(3,1,5,2)=1
        sub_table(3,1,6,1)=1
        sub_table(3,1,6,2)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=2
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=1
        sub_table(5,1,4,1)=1
        sub_table(5,1,4,2)=2
        sub_table(5,1,5,1)=1
        sub_table(5,1,5,2)=3
        sub_table(5,1,6,1)=1
        sub_table(5,1,6,2)=3

     CASE(18)
!
!   C_4h
!
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=1
        sub_table(2,1,3,2)=1
        sub_table(2,1,4,1)=1
        sub_table(2,1,4,2)=1
        sub_table(2,1,5,1)=1
        sub_table(2,1,5,2)=2
        sub_table(2,1,6,1)=1
        sub_table(2,1,6,2)=2
        sub_table(2,1,7,1)=1
        sub_table(2,1,7,2)=2
        sub_table(2,1,8,1)=1
        sub_table(2,1,8,2)=2

        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=1
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=2
        sub_table(3,1,5,1)=1
        sub_table(3,1,5,2)=2
        sub_table(3,1,6,1)=1
        sub_table(3,1,6,2)=1
        sub_table(3,1,7,1)=1
        sub_table(3,1,7,2)=2
        sub_table(3,1,8,1)=1
        sub_table(3,1,8,2)=1

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=2
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=1
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2
        sub_table(4,1,5,1)=1
        sub_table(4,1,5,2)=1
        sub_table(4,1,6,1)=1
        sub_table(4,1,6,2)=2
        sub_table(4,1,7,1)=1
        sub_table(4,1,7,2)=1
        sub_table(4,1,8,1)=1
        sub_table(4,1,8,2)=2

        sub_table(6,1,1,1)=1
        sub_table(6,1,1,2)=1
        sub_table(6,1,2,1)=1
        sub_table(6,1,2,2)=2
        sub_table(6,1,3,1)=1
        sub_table(6,1,3,2)=3
        sub_table(6,1,4,1)=1
        sub_table(6,1,4,2)=4
        sub_table(6,1,5,1)=1
        sub_table(6,1,5,2)=1
        sub_table(6,1,6,1)=1
        sub_table(6,1,6,2)=2
        sub_table(6,1,7,1)=1
        sub_table(6,1,7,2)=3
        sub_table(6,1,8,1)=1
        sub_table(6,1,8,2)=4

        sub_table(16,1,1,1)=1
        sub_table(16,1,1,2)=1
        sub_table(16,1,2,1)=1
        sub_table(16,1,2,2)=2
        sub_table(16,1,3,1)=1
        sub_table(16,1,3,2)=1
        sub_table(16,1,4,1)=1
        sub_table(16,1,4,2)=2
        sub_table(16,1,5,1)=1
        sub_table(16,1,5,2)=3
        sub_table(16,1,6,1)=1
        sub_table(16,1,6,2)=3
        sub_table(16,1,7,1)=1
        sub_table(16,1,7,2)=4
        sub_table(16,1,8,1)=1
        sub_table(16,1,8,2)=4

        sub_table(26,1,1,1)=1
        sub_table(26,1,1,2)=1
        sub_table(26,1,2,1)=1
        sub_table(26,1,2,2)=2
        sub_table(26,1,3,1)=1
        sub_table(26,1,3,2)=3
        sub_table(26,1,4,1)=1
        sub_table(26,1,4,2)=4
        sub_table(26,1,5,1)=1
        sub_table(26,1,5,2)=3
        sub_table(26,1,6,1)=1
        sub_table(26,1,6,2)=4
        sub_table(26,1,7,1)=1
        sub_table(26,1,7,2)=1
        sub_table(26,1,8,1)=1
        sub_table(26,1,8,2)=2

     CASE(19)
!
!  C_6h
!
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=1
        sub_table(2,1,3,2)=1
        sub_table(2,1,4,1)=1
        sub_table(2,1,4,2)=1
        sub_table(2,1,5,1)=1
        sub_table(2,1,5,2)=1
        sub_table(2,1,6,1)=1
        sub_table(2,1,6,2)=1
        sub_table(2,1,7,1)=1
        sub_table(2,1,7,2)=2
        sub_table(2,1,8,1)=1
        sub_table(2,1,8,2)=2
        sub_table(2,1,9,1)=1
        sub_table(2,1,9,2)=2
        sub_table(2,1,10,1)=1
        sub_table(2,1,10,2)=2
        sub_table(2,1,11,1)=1
        sub_table(2,1,11,2)=2
        sub_table(2,1,12,1)=1
        sub_table(2,1,12,2)=2

        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=2
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=1
        sub_table(3,1,5,1)=1
        sub_table(3,1,5,2)=1
        sub_table(3,1,6,1)=1
        sub_table(3,1,6,2)=2
        sub_table(3,1,7,1)=1
        sub_table(3,1,7,2)=2
        sub_table(3,1,8,1)=1
        sub_table(3,1,8,2)=1
        sub_table(3,1,9,1)=1
        sub_table(3,1,9,2)=1
        sub_table(3,1,10,1)=1
        sub_table(3,1,10,2)=2
        sub_table(3,1,11,1)=1
        sub_table(3,1,11,2)=2
        sub_table(3,1,12,1)=1
        sub_table(3,1,12,2)=1

        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=2
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=1
        sub_table(4,1,5,1)=1
        sub_table(4,1,5,2)=1
        sub_table(4,1,6,1)=1
        sub_table(4,1,6,2)=2
        sub_table(4,1,7,1)=1
        sub_table(4,1,7,2)=1
        sub_table(4,1,8,1)=1
        sub_table(4,1,8,2)=2
        sub_table(4,1,9,1)=1
        sub_table(4,1,9,2)=2
        sub_table(4,1,10,1)=1
        sub_table(4,1,10,2)=1
        sub_table(4,1,11,1)=1
        sub_table(4,1,11,2)=1
        sub_table(4,1,12,1)=1
        sub_table(4,1,12,2)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=2
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=1
        sub_table(5,1,4,1)=1
        sub_table(5,1,4,2)=2
        sub_table(5,1,5,1)=1
        sub_table(5,1,5,2)=3
        sub_table(5,1,6,1)=1
        sub_table(5,1,6,2)=3
        sub_table(5,1,7,1)=1
        sub_table(5,1,7,2)=1
        sub_table(5,1,8,1)=1
        sub_table(5,1,8,2)=2
        sub_table(5,1,9,1)=1
        sub_table(5,1,9,2)=1
        sub_table(5,1,10,1)=1
        sub_table(5,1,10,2)=2
        sub_table(5,1,11,1)=1
        sub_table(5,1,11,2)=3
        sub_table(5,1,12,1)=1
        sub_table(5,1,12,2)=3

        sub_table(7,1,1,1)=1
        sub_table(7,1,1,2)=1
        sub_table(7,1,2,1)=1
        sub_table(7,1,2,2)=2
        sub_table(7,1,3,1)=1
        sub_table(7,1,3,2)=5
        sub_table(7,1,4,1)=1
        sub_table(7,1,4,2)=6
        sub_table(7,1,5,1)=1
        sub_table(7,1,5,2)=4
        sub_table(7,1,6,1)=1
        sub_table(7,1,6,2)=3
        sub_table(7,1,7,1)=1
        sub_table(7,1,7,2)=1
        sub_table(7,1,8,1)=1
        sub_table(7,1,8,2)=2
        sub_table(7,1,9,1)=1
        sub_table(7,1,9,2)=5
        sub_table(7,1,10,1)=1
        sub_table(7,1,10,2)=6
        sub_table(7,1,11,1)=1
        sub_table(7,1,11,2)=4
        sub_table(7,1,12,1)=1
        sub_table(7,1,12,2)=3

        sub_table(16,1,1,1)=1
        sub_table(16,1,1,2)=1
        sub_table(16,1,2,1)=1
        sub_table(16,1,2,2)=2
        sub_table(16,1,3,1)=1
        sub_table(16,1,3,2)=2
        sub_table(16,1,4,1)=1
        sub_table(16,1,4,2)=1
        sub_table(16,1,5,1)=1
        sub_table(16,1,5,2)=1
        sub_table(16,1,6,1)=1
        sub_table(16,1,6,2)=2
        sub_table(16,1,7,1)=1
        sub_table(16,1,7,2)=3
        sub_table(16,1,8,1)=1
        sub_table(16,1,8,2)=4
        sub_table(16,1,9,1)=1
        sub_table(16,1,9,2)=4
        sub_table(16,1,10,1)=1
        sub_table(16,1,10,2)=3
        sub_table(16,1,11,1)=1
        sub_table(16,1,11,2)=3
        sub_table(16,1,12,1)=1
        sub_table(16,1,12,2)=4

        sub_table(17,1,1,1)=1
        sub_table(17,1,1,2)=1
        sub_table(17,1,2,1)=1
        sub_table(17,1,2,2)=2
        sub_table(17,1,3,1)=1
        sub_table(17,1,3,2)=3
        sub_table(17,1,4,1)=1
        sub_table(17,1,4,2)=4
        sub_table(17,1,5,1)=1
        sub_table(17,1,5,2)=5
        sub_table(17,1,6,1)=1
        sub_table(17,1,6,2)=6
        sub_table(17,1,7,1)=1
        sub_table(17,1,7,2)=3
        sub_table(17,1,8,1)=1
        sub_table(17,1,8,2)=4
        sub_table(17,1,9,1)=1
        sub_table(17,1,9,2)=1
        sub_table(17,1,10,1)=1
        sub_table(17,1,10,2)=2
        sub_table(17,1,11,1)=1
        sub_table(17,1,11,2)=6
        sub_table(17,1,12,1)=1
        sub_table(17,1,12,2)=5

        sub_table(27,1,1,1)=1
        sub_table(27,1,1,2)=1
        sub_table(27,1,2,1)=1
        sub_table(27,1,2,2)=2
        sub_table(27,1,3,1)=1
        sub_table(27,1,3,2)=1
        sub_table(27,1,4,1)=1
        sub_table(27,1,4,2)=2
        sub_table(27,1,5,1)=1
        sub_table(27,1,5,2)=3
        sub_table(27,1,6,1)=1
        sub_table(27,1,6,2)=3
        sub_table(27,1,7,1)=1
        sub_table(27,1,7,2)=4
        sub_table(27,1,8,1)=1
        sub_table(27,1,8,2)=5
        sub_table(27,1,9,1)=1
        sub_table(27,1,9,2)=4
        sub_table(27,1,10,1)=1
        sub_table(27,1,10,2)=5
        sub_table(27,1,11,1)=1
        sub_table(27,1,11,2)=6
        sub_table(27,1,12,1)=1
        sub_table(27,1,12,2)=6
 
     CASE(20)
!
!  D_2h
!
        sub_table(2,1,1,1)=2
        sub_table(2,1,1,2)=1
        sub_table(2,1,1,3)=1
        sub_table(2,1,2,1)=2
        sub_table(2,1,2,2)=2
        sub_table(2,1,2,3)=2

        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=2
        sub_table(3,1,2,2)=1
        sub_table(3,1,2,3)=2

        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2

        sub_table(8,1,1,1)=2
        sub_table(8,1,1,2)=1
        sub_table(8,1,1,3)=1
        sub_table(8,1,2,1)=2
        sub_table(8,1,2,2)=1
        sub_table(8,1,2,3)=1

        sub_table(12,1,1,1)=2
        sub_table(12,1,1,2)=1
        sub_table(12,1,1,3)=1
        sub_table(12,1,2,1)=2
        sub_table(12,1,2,2)=1
        sub_table(12,1,2,3)=1

        sub_table(16,1,1,1)=2
        sub_table(16,1,1,2)=1
        sub_table(16,1,1,3)=2
        sub_table(16,1,2,1)=2
        sub_table(16,1,2,2)=3
        sub_table(16,1,2,3)=4

     CASE(21)
!
!  D_3h
!
        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=2
        sub_table(3,1,2,2)=1
        sub_table(3,1,2,3)=2
        sub_table(3,1,3,1)=2
        sub_table(3,1,3,2)=1
        sub_table(3,1,3,3)=2

        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2
        sub_table(4,1,3,1)=2
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=2

        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=2
        sub_table(5,1,2,2)=1
        sub_table(5,1,2,3)=2
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=3
        sub_table(5,1,3,3)=3

        sub_table(9,1,1,1)=2
        sub_table(9,1,1,2)=1
        sub_table(9,1,1,3)=1
        sub_table(9,1,2,1)=2
        sub_table(9,1,2,2)=1
        sub_table(9,1,2,3)=1
        sub_table(9,1,3,1)=2
        sub_table(9,1,3,2)=2
        sub_table(9,1,3,3)=3

        sub_table(12,1,1,1)=2
        sub_table(12,1,1,2)=1
        sub_table(12,1,1,3)=1
        sub_table(12,1,2,1)=2
        sub_table(12,1,2,2)=1
        sub_table(12,1,2,3)=1
        sub_table(12,1,3,1)=2
        sub_table(12,1,3,2)=1
        sub_table(12,1,3,3)=1

        sub_table(13,1,1,1)=2
        sub_table(13,1,1,2)=1
        sub_table(13,1,1,3)=1
        sub_table(13,1,2,1)=2
        sub_table(13,1,2,2)=1
        sub_table(13,1,2,3)=1
        sub_table(13,1,3,1)=2
        sub_table(13,1,3,2)=2
        sub_table(13,1,3,3)=3

        sub_table(17,1,1,1)=2
        sub_table(17,1,1,2)=1
        sub_table(17,1,1,3)=2
        sub_table(17,1,2,1)=2
        sub_table(17,1,2,2)=3
        sub_table(17,1,2,3)=4
        sub_table(17,1,3,1)=2
        sub_table(17,1,3,2)=5
        sub_table(17,1,3,3)=6

 
     CASE(22)
!
!  D_4h
!
        sub_table(2,1,1,1)=2
        sub_table(2,1,1,2)=1
        sub_table(2,1,1,3)=1
        sub_table(2,1,2,1)=2
        sub_table(2,1,2,2)=1
        sub_table(2,1,2,3)=1
        sub_table(2,1,3,1)=2
        sub_table(2,1,3,2)=2
        sub_table(2,1,3,3)=2
        sub_table(2,1,4,1)=2
        sub_table(2,1,4,2)=2
        sub_table(2,1,4,3)=2

        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=2
        sub_table(3,1,2,2)=1
        sub_table(3,1,2,3)=2
        sub_table(3,1,3,1)=2
        sub_table(3,1,3,2)=1
        sub_table(3,1,3,3)=2
        sub_table(3,1,4,1)=2
        sub_table(3,1,4,2)=1
        sub_table(3,1,4,3)=2

        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2
        sub_table(4,1,3,1)=2
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=2
        sub_table(4,1,4,1)=2
        sub_table(4,1,4,2)=1
        sub_table(4,1,4,3)=2

        sub_table(6,1,1,1)=2
        sub_table(6,1,1,2)=1
        sub_table(6,1,1,3)=2
        sub_table(6,1,2,1)=2
        sub_table(6,1,2,2)=3
        sub_table(6,1,2,3)=4
        sub_table(6,1,3,1)=2
        sub_table(6,1,3,2)=1
        sub_table(6,1,3,3)=2
        sub_table(6,1,4,1)=2
        sub_table(6,1,4,2)=3
        sub_table(6,1,4,3)=4

        sub_table(8,1,1,1)=2
        sub_table(8,1,1,2)=1
        sub_table(8,1,1,3)=1
        sub_table(8,1,2,1)=2
        sub_table(8,1,2,2)=1
        sub_table(8,1,2,3)=1
        sub_table(8,1,3,1)=2
        sub_table(8,1,3,2)=1
        sub_table(8,1,3,3)=1
        sub_table(8,1,4,1)=2
        sub_table(8,1,4,2)=1
        sub_table(8,1,4,3)=1

        sub_table(10,1,1,1)=2
        sub_table(10,1,1,2)=1
        sub_table(10,1,1,3)=1
        sub_table(10,1,2,1)=2
        sub_table(10,1,2,2)=2
        sub_table(10,1,2,3)=2
        sub_table(10,1,3,1)=2
        sub_table(10,1,3,2)=1
        sub_table(10,1,3,3)=1
        sub_table(10,1,4,1)=2
        sub_table(10,1,4,2)=2
        sub_table(10,1,4,3)=2

        sub_table(12,1,1,1)=2
        sub_table(12,1,1,2)=1
        sub_table(12,1,1,3)=1
        sub_table(12,1,2,1)=2
        sub_table(12,1,2,2)=1
        sub_table(12,1,2,3)=1
        sub_table(12,1,3,1)=2
        sub_table(12,1,3,2)=1
        sub_table(12,1,3,3)=1
        sub_table(12,1,4,1)=2
        sub_table(12,1,4,2)=1
        sub_table(12,1,4,3)=1

        sub_table(14,1,1,1)=2
        sub_table(14,1,1,2)=1
        sub_table(14,1,1,3)=1
        sub_table(14,1,2,1)=2
        sub_table(14,1,2,2)=2
        sub_table(14,1,2,3)=2
        sub_table(14,1,3,1)=2
        sub_table(14,1,3,2)=1
        sub_table(14,1,3,3)=1
        sub_table(14,1,4,1)=2
        sub_table(14,1,4,2)=2
        sub_table(14,1,4,3)=2

        sub_table(16,1,1,1)=2
        sub_table(16,1,1,2)=1
        sub_table(16,1,1,3)=2
        sub_table(16,1,2,1)=2
        sub_table(16,1,2,2)=1
        sub_table(16,1,2,3)=2
        sub_table(16,1,3,1)=2
        sub_table(16,1,3,2)=3
        sub_table(16,1,3,3)=4
        sub_table(16,1,4,1)=2
        sub_table(16,1,4,2)=3
        sub_table(16,1,4,3)=4

        sub_table(18,1,1,1)=2
        sub_table(18,1,1,2)=1
        sub_table(18,1,1,3)=2
        sub_table(18,1,2,1)=2
        sub_table(18,1,2,2)=3
        sub_table(18,1,2,3)=4
        sub_table(18,1,3,1)=2
        sub_table(18,1,3,2)=5
        sub_table(18,1,3,3)=6
        sub_table(18,1,4,1)=2
        sub_table(18,1,4,2)=7
        sub_table(18,1,4,3)=8

        sub_table(20,1,1,1)=2
        sub_table(20,1,1,2)=1
        sub_table(20,1,1,3)=1
        sub_table(20,1,2,1)=2
        sub_table(20,1,2,2)=1
        sub_table(20,1,2,3)=1
        sub_table(20,1,3,1)=2
        sub_table(20,1,3,2)=2
        sub_table(20,1,3,3)=2
        sub_table(20,1,4,1)=2
        sub_table(20,1,4,2)=2
        sub_table(20,1,4,3)=2

        sub_table(24,1,1,1)=2
        sub_table(24,1,1,2)=1
        sub_table(24,1,1,3)=1
        sub_table(24,1,2,1)=2
        sub_table(24,1,2,2)=2
        sub_table(24,1,2,3)=2
        sub_table(24,1,3,1)=2
        sub_table(24,1,3,2)=2
        sub_table(24,1,3,3)=2
        sub_table(24,1,4,1)=2
        sub_table(24,1,4,2)=1
        sub_table(24,1,4,3)=1

        sub_table(26,1,1,1)=2
        sub_table(26,1,1,2)=1
        sub_table(26,1,1,3)=2
        sub_table(26,1,2,1)=2
        sub_table(26,1,2,2)=3
        sub_table(26,1,2,3)=4
        sub_table(26,1,3,1)=2
        sub_table(26,1,3,2)=3
        sub_table(26,1,3,3)=4
        sub_table(26,1,4,1)=2
        sub_table(26,1,4,2)=1
        sub_table(26,1,4,3)=2

     CASE(23)
!
! D_6h
!
        sub_table(2,1,1,1)=2
        sub_table(2,1,1,2)=1
        sub_table(2,1,1,3)=1
        sub_table(2,1,2,1)=2
        sub_table(2,1,2,2)=1
        sub_table(2,1,2,3)=1
        sub_table(2,1,3,1)=2
        sub_table(2,1,3,2)=1
        sub_table(2,1,3,3)=1
        sub_table(2,1,4,1)=2
        sub_table(2,1,4,2)=2
        sub_table(2,1,4,3)=2
        sub_table(2,1,5,1)=2
        sub_table(2,1,5,2)=2
        sub_table(2,1,5,3)=2
        sub_table(2,1,6,1)=2
        sub_table(2,1,6,2)=2
        sub_table(2,1,6,3)=2

        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=2
        sub_table(3,1,2,2)=1
        sub_table(3,1,2,3)=2
        sub_table(3,1,3,1)=2
        sub_table(3,1,3,2)=1
        sub_table(3,1,3,3)=2
        sub_table(3,1,4,1)=2
        sub_table(3,1,4,2)=1
        sub_table(3,1,4,3)=2
        sub_table(3,1,5,1)=2
        sub_table(3,1,5,2)=1
        sub_table(3,1,5,3)=2
        sub_table(3,1,6,1)=2
        sub_table(3,1,6,2)=1
        sub_table(3,1,6,3)=2

        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2
        sub_table(4,1,3,1)=2
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=2
        sub_table(4,1,4,1)=2
        sub_table(4,1,4,2)=1
        sub_table(4,1,4,3)=2
        sub_table(4,1,5,1)=2
        sub_table(4,1,5,2)=1
        sub_table(4,1,5,3)=2
        sub_table(4,1,6,1)=2
        sub_table(4,1,6,2)=1
        sub_table(4,1,6,3)=2

        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=2
        sub_table(5,1,2,2)=1
        sub_table(5,1,2,3)=2
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=3
        sub_table(5,1,3,3)=3
        sub_table(5,1,4,1)=2
        sub_table(5,1,4,2)=1
        sub_table(5,1,4,3)=2
        sub_table(5,1,5,1)=2
        sub_table(5,1,5,2)=1
        sub_table(5,1,5,3)=2
        sub_table(5,1,6,1)=2
        sub_table(5,1,6,2)=3
        sub_table(5,1,6,3)=3

        sub_table(7,1,1,1)=2
        sub_table(7,1,1,2)=1
        sub_table(7,1,1,3)=2
        sub_table(7,1,2,1)=2
        sub_table(7,1,2,2)=5
        sub_table(7,1,2,3)=6
        sub_table(7,1,3,1)=2
        sub_table(7,1,3,2)=4
        sub_table(7,1,3,3)=3
        sub_table(7,1,4,1)=2
        sub_table(7,1,4,2)=1
        sub_table(7,1,4,3)=2
        sub_table(7,1,5,1)=2
        sub_table(7,1,5,2)=5
        sub_table(7,1,5,3)=6
        sub_table(7,1,6,1)=2
        sub_table(7,1,6,2)=3
        sub_table(7,1,6,3)=4

        sub_table(8,1,1,1)=2
        sub_table(8,1,1,2)=1
        sub_table(8,1,1,3)=1
        sub_table(8,1,2,1)=2
        sub_table(8,1,2,2)=1
        sub_table(8,1,2,3)=1
        sub_table(8,1,3,1)=2
        sub_table(8,1,3,2)=1
        sub_table(8,1,3,3)=1
        sub_table(8,1,4,1)=2
        sub_table(8,1,4,2)=1
        sub_table(8,1,4,3)=1
        sub_table(8,1,5,1)=2
        sub_table(8,1,5,2)=1
        sub_table(8,1,5,3)=1
        sub_table(8,1,6,1)=2
        sub_table(8,1,6,2)=1
        sub_table(8,1,6,3)=1

        sub_table(9,1,1,1)=2
        sub_table(9,1,1,2)=1
        sub_table(9,1,1,3)=1
        sub_table(9,1,2,1)=2
        sub_table(9,1,2,2)=1
        sub_table(9,1,2,3)=1
        sub_table(9,1,3,1)=2
        sub_table(9,1,3,2)=2
        sub_table(9,1,3,3)=3
        sub_table(9,1,4,1)=2
        sub_table(9,1,4,2)=1
        sub_table(9,1,4,3)=1
        sub_table(9,1,5,1)=2
        sub_table(9,1,5,2)=1
        sub_table(9,1,5,3)=1
        sub_table(9,1,6,1)=2
        sub_table(9,1,6,2)=2
        sub_table(9,1,6,3)=3

        sub_table(11,1,1,1)=2
        sub_table(11,1,1,2)=1
        sub_table(11,1,1,3)=1
        sub_table(11,1,2,1)=2
        sub_table(11,1,2,2)=2
        sub_table(11,1,2,3)=2
        sub_table(11,1,3,1)=2
        sub_table(11,1,3,2)=3
        sub_table(11,1,3,3)=3
        sub_table(11,1,4,1)=2
        sub_table(11,1,4,2)=1
        sub_table(11,1,4,3)=1
        sub_table(11,1,5,1)=2
        sub_table(11,1,5,2)=2
        sub_table(11,1,5,3)=2
        sub_table(11,1,6,1)=2
        sub_table(11,1,6,2)=3
        sub_table(11,1,6,3)=3

        sub_table(12,1,1,1)=2
        sub_table(12,1,1,2)=1
        sub_table(12,1,1,3)=1
        sub_table(12,1,2,1)=2
        sub_table(12,1,2,2)=1
        sub_table(12,1,2,3)=1
        sub_table(12,1,3,1)=2
        sub_table(12,1,3,2)=1
        sub_table(12,1,3,3)=1
        sub_table(12,1,4,1)=2
        sub_table(12,1,4,2)=1
        sub_table(12,1,4,3)=1
        sub_table(12,1,5,1)=2
        sub_table(12,1,5,2)=1
        sub_table(12,1,5,3)=1
        sub_table(12,1,6,1)=2
        sub_table(12,1,6,2)=1
        sub_table(12,1,6,3)=1

        sub_table(13,1,1,1)=2
        sub_table(13,1,1,2)=1
        sub_table(13,1,1,3)=1
        sub_table(13,1,2,1)=2
        sub_table(13,1,2,2)=1
        sub_table(13,1,2,3)=1
        sub_table(13,1,3,1)=2
        sub_table(13,1,3,2)=2
        sub_table(13,1,3,3)=3
        sub_table(13,1,4,1)=2
        sub_table(13,1,4,2)=1
        sub_table(13,1,4,3)=1
        sub_table(13,1,5,1)=2
        sub_table(13,1,5,2)=1
        sub_table(13,1,5,3)=1
        sub_table(13,1,6,1)=2
        sub_table(13,1,6,2)=2
        sub_table(13,1,6,3)=3

        sub_table(15,1,1,1)=2
        sub_table(15,1,1,2)=1
        sub_table(15,1,1,3)=1
        sub_table(15,1,2,1)=2
        sub_table(15,1,2,2)=2
        sub_table(15,1,2,3)=2
        sub_table(15,1,3,1)=2
        sub_table(15,1,3,2)=3
        sub_table(15,1,3,3)=3
        sub_table(15,1,4,1)=2
        sub_table(15,1,4,2)=1
        sub_table(15,1,4,3)=1
        sub_table(15,1,5,1)=2
        sub_table(15,1,5,2)=2
        sub_table(15,1,5,3)=2
        sub_table(15,1,6,1)=2
        sub_table(15,1,6,2)=3
        sub_table(15,1,6,3)=3

        sub_table(16,1,1,1)=2
        sub_table(16,1,1,2)=1
        sub_table(16,1,1,3)=2
        sub_table(16,1,2,1)=2
        sub_table(16,1,2,2)=1
        sub_table(16,1,2,3)=2
        sub_table(16,1,3,1)=2
        sub_table(16,1,3,2)=1
        sub_table(16,1,3,3)=2
        sub_table(16,1,4,1)=2
        sub_table(16,1,4,2)=3
        sub_table(16,1,4,3)=4
        sub_table(16,1,5,1)=2
        sub_table(16,1,5,2)=3
        sub_table(16,1,5,3)=4
        sub_table(16,1,6,1)=2
        sub_table(16,1,6,2)=3
        sub_table(16,1,6,3)=4

        sub_table(17,1,1,1)=2
        sub_table(17,1,1,2)=1
        sub_table(17,1,1,3)=2
        sub_table(17,1,2,1)=2
        sub_table(17,1,2,2)=3
        sub_table(17,1,2,3)=4
        sub_table(17,1,3,1)=2
        sub_table(17,1,3,2)=5
        sub_table(17,1,3,3)=6
        sub_table(17,1,4,1)=2
        sub_table(17,1,4,2)=3
        sub_table(17,1,4,3)=4
        sub_table(17,1,5,1)=2
        sub_table(17,1,5,2)=1
        sub_table(17,1,5,3)=2
        sub_table(17,1,6,1)=2
        sub_table(17,1,6,2)=5
        sub_table(17,1,6,3)=6

        sub_table(19,1,1,1)=2
        sub_table(19,1,1,2)=1
        sub_table(19,1,1,3)=2
        sub_table(19,1,2,1)=2
        sub_table(19,1,2,2)=5
        sub_table(19,1,2,3)=6
        sub_table(19,1,3,1)=2
        sub_table(19,1,3,2)=3
        sub_table(19,1,3,3)=4
        sub_table(19,1,4,1)=2
        sub_table(19,1,4,2)=7
        sub_table(19,1,4,3)=8
        sub_table(19,1,5,1)=2
        sub_table(19,1,5,2)=11
        sub_table(19,1,5,3)=12
        sub_table(19,1,6,1)=2
        sub_table(19,1,6,2)=9
        sub_table(19,1,6,3)=10

        sub_table(20,1,1,1)=2
        sub_table(20,1,1,2)=1
        sub_table(20,1,1,3)=1
        sub_table(20,1,2,1)=2
        sub_table(20,1,2,2)=1
        sub_table(20,1,2,3)=1
        sub_table(20,1,3,1)=2
        sub_table(20,1,3,2)=1
        sub_table(20,1,3,3)=1
        sub_table(20,1,4,1)=2
        sub_table(20,1,4,2)=2
        sub_table(20,1,4,3)=2
        sub_table(20,1,5,1)=2
        sub_table(20,1,5,2)=2
        sub_table(20,1,5,3)=2
        sub_table(20,1,6,1)=2
        sub_table(20,1,6,2)=2
        sub_table(20,1,6,3)=2

        sub_table(21,1,1,1)=2
        sub_table(21,1,1,2)=1
        sub_table(21,1,1,3)=1
        sub_table(21,1,2,1)=2
        sub_table(21,1,2,2)=2
        sub_table(21,1,2,3)=2
        sub_table(21,1,3,1)=2
        sub_table(21,1,3,2)=3
        sub_table(21,1,3,3)=3
        sub_table(21,1,4,1)=2
        sub_table(21,1,4,2)=1
        sub_table(21,1,4,3)=1
        sub_table(21,1,5,1)=2
        sub_table(21,1,5,2)=2
        sub_table(21,1,5,3)=2
        sub_table(21,1,6,1)=2
        sub_table(21,1,6,2)=3
        sub_table(21,1,6,3)=3

        sub_table(25,1,1,1)=2
        sub_table(25,1,1,2)=1
        sub_table(25,1,1,3)=2
        sub_table(25,1,2,1)=2
        sub_table(25,1,2,2)=1
        sub_table(25,1,2,3)=2
        sub_table(25,1,3,1)=2
        sub_table(25,1,3,2)=3
        sub_table(25,1,3,3)=3
        sub_table(25,1,4,1)=2
        sub_table(25,1,4,2)=4
        sub_table(25,1,4,3)=5
        sub_table(25,1,5,1)=2
        sub_table(25,1,5,2)=4
        sub_table(25,1,5,3)=5
        sub_table(25,1,6,1)=2
        sub_table(25,1,6,2)=6
        sub_table(25,1,6,3)=6


     CASE(24)
!
! D_2d
!
        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=2
        sub_table(3,1,2,2)=1
        sub_table(3,1,2,3)=2

        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2

        sub_table(8,1,1,1)=2
        sub_table(8,1,1,2)=1
        sub_table(8,1,1,3)=1
        sub_table(8,1,2,1)=2
        sub_table(8,1,2,2)=1
        sub_table(8,1,2,3)=1

        sub_table(12,1,1,1)=2
        sub_table(12,1,1,2)=1
        sub_table(12,1,1,3)=1
        sub_table(12,1,2,1)=2
        sub_table(12,1,2,2)=1
        sub_table(12,1,2,3)=1

        sub_table(26,1,1,1)=2
        sub_table(26,1,1,2)=1
        sub_table(26,1,1,3)=2
        sub_table(26,1,2,1)=2
        sub_table(26,1,2,2)=3
        sub_table(26,1,2,3)=4

     CASE(25)
!
!  D_3d
!
        sub_table(2,1,1,1)=2
        sub_table(2,1,1,2)=1
        sub_table(2,1,1,3)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=1
        sub_table(2,1,3,2)=1
        sub_table(2,1,4,1)=2
        sub_table(2,1,4,2)=2
        sub_table(2,1,4,3)=2
        sub_table(2,1,5,1)=1
        sub_table(2,1,5,2)=2
        sub_table(2,1,6,1)=1
        sub_table(2,1,6,2)=2

        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=1
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=2
        sub_table(3,1,4,1)=2
        sub_table(3,1,4,2)=1
        sub_table(3,1,4,3)=2
        sub_table(3,1,5,1)=1
        sub_table(3,1,5,2)=2
        sub_table(3,1,6,1)=1
        sub_table(3,1,6,2)=1

        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=2
        sub_table(4,1,4,2)=1
        sub_table(4,1,4,3)=2
        sub_table(4,1,5,1)=1
        sub_table(4,1,5,2)=1
        sub_table(4,1,6,1)=1
        sub_table(4,1,6,2)=2

        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=3
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=3
        sub_table(5,1,4,1)=2
        sub_table(5,1,4,2)=1
        sub_table(5,1,4,3)=2
        sub_table(5,1,5,1)=1
        sub_table(5,1,5,2)=3
        sub_table(5,1,6,1)=1
        sub_table(5,1,6,2)=3

        sub_table(9,1,1,1)=2
        sub_table(9,1,1,2)=1
        sub_table(9,1,1,3)=1
        sub_table(9,1,2,1)=1
        sub_table(9,1,2,2)=2
        sub_table(9,1,3,1)=1
        sub_table(9,1,3,2)=3
        sub_table(9,1,4,1)=2
        sub_table(9,1,4,2)=1
        sub_table(9,1,4,3)=1
        sub_table(9,1,5,1)=1
        sub_table(9,1,5,2)=2
        sub_table(9,1,6,1)=1
        sub_table(9,1,6,2)=3

        sub_table(13,1,1,1)=2
        sub_table(13,1,1,2)=1
        sub_table(13,1,1,3)=1
        sub_table(13,1,2,1)=1
        sub_table(13,1,2,2)=2
        sub_table(13,1,3,1)=1
        sub_table(13,1,3,2)=3
        sub_table(13,1,4,1)=2
        sub_table(13,1,4,2)=1
        sub_table(13,1,4,3)=1
        sub_table(13,1,5,1)=1
        sub_table(13,1,5,2)=3
        sub_table(13,1,6,1)=1
        sub_table(13,1,6,2)=2

        sub_table(16,1,1,1)=2
        sub_table(16,1,1,2)=1
        sub_table(16,1,1,3)=2
        sub_table(16,1,2,1)=1
        sub_table(16,1,2,2)=1
        sub_table(16,1,3,1)=1
        sub_table(16,1,3,2)=2
        sub_table(16,1,4,1)=2
        sub_table(16,1,4,2)=3
        sub_table(16,1,4,3)=4
        sub_table(16,1,5,1)=1
        sub_table(16,1,5,2)=3
        sub_table(16,1,6,1)=1
        sub_table(16,1,6,2)=4

        sub_table(27,1,1,1)=2
        sub_table(27,1,1,2)=1
        sub_table(27,1,1,3)=2
        sub_table(27,1,2,1)=1
        sub_table(27,1,2,2)=3
        sub_table(27,1,3,1)=1
        sub_table(27,1,3,2)=3
        sub_table(27,1,4,1)=2
        sub_table(27,1,4,2)=4
        sub_table(27,1,4,3)=5
        sub_table(27,1,5,1)=1
        sub_table(27,1,5,2)=6
        sub_table(27,1,6,1)=1
        sub_table(27,1,6,2)=6

     CASE(26)
!
!  S_4
!
        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=2
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=1
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2

     CASE(27)
!
!  S_6
!
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=1
        sub_table(2,1,3,2)=1
        sub_table(2,1,4,1)=1
        sub_table(2,1,4,2)=2
        sub_table(2,1,5,1)=1
        sub_table(2,1,5,2)=2
        sub_table(2,1,6,1)=1
        sub_table(2,1,6,2)=2

        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=2
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=3
        sub_table(5,1,4,1)=1
        sub_table(5,1,4,2)=1
        sub_table(5,1,5,1)=1
        sub_table(5,1,5,2)=2
        sub_table(5,1,6,1)=1
        sub_table(5,1,6,2)=3

     CASE(28)
!
!  T
!
        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2
        sub_table(4,1,3,1)=2
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=2

        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=2
        sub_table(5,1,2,2)=3
        sub_table(5,1,2,3)=1
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=3
        sub_table(5,1,3,3)=2


        sub_table(8,1,1,1)=2
        sub_table(8,1,1,2)=1
        sub_table(8,1,1,3)=1
        sub_table(8,1,2,1)=2
        sub_table(8,1,2,2)=1
        sub_table(8,1,2,3)=1
        sub_table(8,1,3,1)=2
        sub_table(8,1,3,2)=1
        sub_table(8,1,3,3)=1


     CASE(29)
!
!  T_h
!

        sub_table(2,1,1,1)=2
        sub_table(2,1,1,2)=1
        sub_table(2,1,1,3)=1
        sub_table(2,1,2,1)=2
        sub_table(2,1,2,2)=1
        sub_table(2,1,2,3)=1
        sub_table(2,1,3,1)=2
        sub_table(2,1,3,2)=1
        sub_table(2,1,3,3)=1
        sub_table(2,1,4,1)=2
        sub_table(2,1,4,2)=2
        sub_table(2,1,4,3)=2
        sub_table(2,1,5,1)=2
        sub_table(2,1,5,2)=2
        sub_table(2,1,5,3)=2
        sub_table(2,1,6,1)=2
        sub_table(2,1,6,2)=2
        sub_table(2,1,6,3)=2

        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=2
        sub_table(3,1,2,2)=1
        sub_table(3,1,2,3)=2
        sub_table(3,1,3,1)=2
        sub_table(3,1,3,2)=1
        sub_table(3,1,3,3)=2
        sub_table(3,1,4,1)=2
        sub_table(3,1,4,2)=1
        sub_table(3,1,4,3)=2
        sub_table(3,1,5,1)=2
        sub_table(3,1,5,2)=1
        sub_table(3,1,5,3)=2
        sub_table(3,1,6,1)=2
        sub_table(3,1,6,2)=1
        sub_table(3,1,6,3)=2

        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2
        sub_table(4,1,3,1)=2
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=2
        sub_table(4,1,4,1)=2
        sub_table(4,1,4,2)=1
        sub_table(4,1,4,3)=2
        sub_table(4,1,5,1)=2
        sub_table(4,1,5,2)=1
        sub_table(4,1,5,3)=2
        sub_table(4,1,6,1)=2
        sub_table(4,1,6,2)=1
        sub_table(4,1,6,3)=2

        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=2
        sub_table(5,1,2,2)=1
        sub_table(5,1,2,3)=3
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=1
        sub_table(5,1,3,3)=3
        sub_table(5,1,4,1)=2
        sub_table(5,1,4,2)=1
        sub_table(5,1,4,3)=2
        sub_table(5,1,5,1)=2
        sub_table(5,1,5,2)=1
        sub_table(5,1,5,3)=3
        sub_table(5,1,6,1)=2
        sub_table(5,1,6,2)=3
        sub_table(5,1,6,3)=2

        sub_table(8,1,1,1)=2
        sub_table(8,1,1,2)=1
        sub_table(8,1,1,3)=1
        sub_table(8,1,2,1)=2
        sub_table(8,1,2,2)=1
        sub_table(8,1,2,3)=1
        sub_table(8,1,3,1)=2
        sub_table(8,1,3,2)=1
        sub_table(8,1,3,3)=1
        sub_table(8,1,4,1)=2
        sub_table(8,1,4,2)=1
        sub_table(8,1,4,3)=1
        sub_table(8,1,5,1)=2
        sub_table(8,1,5,2)=1
        sub_table(8,1,5,3)=1
        sub_table(8,1,6,1)=2
        sub_table(8,1,6,2)=1
        sub_table(8,1,6,3)=1

        sub_table(12,1,1,1)=2
        sub_table(12,1,1,2)=1
        sub_table(12,1,1,3)=1
        sub_table(12,1,2,1)=2
        sub_table(12,1,2,2)=1
        sub_table(12,1,2,3)=1
        sub_table(12,1,3,1)=2
        sub_table(12,1,3,2)=1
        sub_table(12,1,3,3)=1
        sub_table(12,1,4,1)=2
        sub_table(12,1,4,2)=1
        sub_table(12,1,4,3)=1
        sub_table(12,1,5,1)=2
        sub_table(12,1,5,2)=1
        sub_table(12,1,5,3)=1
        sub_table(12,1,6,1)=2
        sub_table(12,1,6,2)=1
        sub_table(12,1,6,3)=1

        sub_table(16,1,1,1)=2
        sub_table(16,1,1,2)=1
        sub_table(16,1,1,3)=2
        sub_table(16,1,2,1)=2
        sub_table(16,1,2,2)=1
        sub_table(16,1,2,3)=2
        sub_table(16,1,3,1)=2
        sub_table(16,1,3,2)=1
        sub_table(16,1,3,3)=2
        sub_table(16,1,4,1)=2
        sub_table(16,1,4,2)=3
        sub_table(16,1,4,3)=4
        sub_table(16,1,5,1)=2
        sub_table(16,1,5,2)=3
        sub_table(16,1,5,3)=4
        sub_table(16,1,6,1)=2
        sub_table(16,1,6,2)=3
        sub_table(16,1,6,3)=4

        sub_table(20,1,1,1)=2
        sub_table(20,1,1,2)=1
        sub_table(20,1,1,3)=1
        sub_table(20,1,2,1)=2
        sub_table(20,1,2,2)=1
        sub_table(20,1,2,3)=1
        sub_table(20,1,3,1)=2
        sub_table(20,1,3,2)=1
        sub_table(20,1,3,3)=1
        sub_table(20,1,4,1)=2
        sub_table(20,1,4,2)=2
        sub_table(20,1,4,3)=2
        sub_table(20,1,5,1)=2
        sub_table(20,1,5,2)=2
        sub_table(20,1,5,3)=2
        sub_table(20,1,6,1)=2
        sub_table(20,1,6,2)=2
        sub_table(20,1,6,3)=2

        sub_table(27,1,1,1)=2
        sub_table(27,1,1,2)=1
        sub_table(27,1,1,3)=2
        sub_table(27,1,2,1)=2
        sub_table(27,1,2,2)=2
        sub_table(27,1,2,3)=3
        sub_table(27,1,3,1)=2
        sub_table(27,1,3,2)=1
        sub_table(27,1,3,3)=3
        sub_table(27,1,4,1)=2
        sub_table(27,1,4,2)=4
        sub_table(27,1,4,3)=5
        sub_table(27,1,5,1)=2
        sub_table(27,1,5,2)=5
        sub_table(27,1,5,3)=6
        sub_table(27,1,6,1)=2
        sub_table(27,1,6,2)=4
        sub_table(27,1,6,3)=6

        sub_table(28,1,1,1)=2
        sub_table(28,1,1,2)=1
        sub_table(28,1,1,3)=1
        sub_table(28,1,2,1)=2
        sub_table(28,1,2,2)=2
        sub_table(28,1,2,3)=2
        sub_table(28,1,3,1)=2
        sub_table(28,1,3,2)=3
        sub_table(28,1,3,3)=3
        sub_table(28,1,4,1)=2
        sub_table(28,1,4,2)=1
        sub_table(28,1,4,3)=1
        sub_table(28,1,5,1)=2
        sub_table(28,1,5,2)=2
        sub_table(28,1,5,3)=2
        sub_table(28,1,6,1)=2
        sub_table(28,1,6,2)=3
        sub_table(28,1,6,3)=3

 
     CASE(30)
!
!  T_d
!
        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=2
        sub_table(3,1,2,2)=1
        sub_table(3,1,2,3)=2
        sub_table(3,1,3,1)=4
        sub_table(3,1,3,2)=1
        sub_table(3,1,3,3)=1
        sub_table(3,1,3,4)=2
        sub_table(3,1,3,5)=2

        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2
        sub_table(4,1,3,1)=4
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=1
        sub_table(4,1,3,4)=2
        sub_table(4,1,3,5)=2

        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=2
        sub_table(5,1,2,2)=1
        sub_table(5,1,2,3)=2
        sub_table(5,1,3,1)=4
        sub_table(5,1,3,2)=1
        sub_table(5,1,3,3)=2
        sub_table(5,1,3,4)=3
        sub_table(5,1,3,5)=3

        sub_table(8,1,1,1)=2
        sub_table(8,1,1,2)=1
        sub_table(8,1,1,3)=1
        sub_table(8,1,2,1)=2
        sub_table(8,1,2,2)=1
        sub_table(8,1,2,3)=1
        sub_table(8,1,3,1)=4
        sub_table(8,1,3,2)=1
        sub_table(8,1,3,3)=1
        sub_table(8,1,3,4)=1
        sub_table(8,1,3,5)=1


        sub_table(12,1,1,1)=2
        sub_table(12,1,1,2)=1
        sub_table(12,1,1,3)=1
        sub_table(12,1,2,1)=2
        sub_table(12,1,2,2)=1
        sub_table(12,1,2,3)=1
        sub_table(12,1,3,1)=4
        sub_table(12,1,3,2)=1
        sub_table(12,1,3,3)=1
        sub_table(12,1,3,4)=1
        sub_table(12,1,3,5)=1

        sub_table(13,1,1,1)=2
        sub_table(13,1,1,2)=1
        sub_table(13,1,1,3)=1
        sub_table(13,1,2,1)=2
        sub_table(13,1,2,2)=1
        sub_table(13,1,2,3)=1
        sub_table(13,1,3,1)=4
        sub_table(13,1,3,2)=1
        sub_table(13,1,3,3)=1
        sub_table(13,1,3,4)=2
        sub_table(13,1,3,5)=3

        sub_table(24,1,1,1)=2
        sub_table(24,1,1,2)=1
        sub_table(24,1,1,3)=1
        sub_table(24,1,2,1)=2
        sub_table(24,1,2,2)=2
        sub_table(24,1,2,3)=2
        sub_table(24,1,3,1)=4
        sub_table(24,1,3,2)=1
        sub_table(24,1,3,3)=1
        sub_table(24,1,3,4)=2
        sub_table(24,1,3,5)=2

        sub_table(26,1,1,1)=2
        sub_table(26,1,1,2)=1
        sub_table(26,1,1,3)=2
        sub_table(26,1,2,1)=2
        sub_table(26,1,2,2)=3
        sub_table(26,1,2,3)=4
        sub_table(26,1,3,1)=4
        sub_table(26,1,3,2)=1
        sub_table(26,1,3,3)=2
        sub_table(26,1,3,4)=3
        sub_table(26,1,3,5)=4

        sub_table(28,1,1,1)=2
        sub_table(28,1,1,2)=1
        sub_table(28,1,1,3)=1
        sub_table(28,1,2,1)=2
        sub_table(28,1,2,2)=1
        sub_table(28,1,2,3)=1
        sub_table(28,1,3,1)=4
        sub_table(28,1,3,2)=2
        sub_table(28,1,3,3)=2
        sub_table(28,1,3,4)=3
        sub_table(28,1,3,5)=3

     CASE(31)
!
!  O
!
        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2
        sub_table(4,1,3,1)=2
        sub_table(4,1,3,2)=4
        sub_table(4,1,3,3)=1
        sub_table(4,1,3,4)=2
        sub_table(4,1,3,5)=2

        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=2
        sub_table(5,1,2,2)=1
        sub_table(5,1,2,3)=2
        sub_table(5,1,3,1)=4
        sub_table(5,1,3,2)=1
        sub_table(5,1,3,3)=2
        sub_table(5,1,3,4)=3
        sub_table(5,1,3,5)=4

        sub_table(6,1,1,1)=2
        sub_table(6,1,1,2)=1
        sub_table(6,1,1,3)=2
        sub_table(6,1,2,1)=2
        sub_table(6,1,2,2)=3
        sub_table(6,1,2,3)=4
        sub_table(6,1,3,1)=4
        sub_table(6,1,3,2)=1
        sub_table(6,1,3,3)=2
        sub_table(6,1,3,4)=3
        sub_table(6,1,3,5)=4

        sub_table(8,1,1,1)=2
        sub_table(8,1,1,2)=1
        sub_table(8,1,1,3)=1
        sub_table(8,1,2,1)=2
        sub_table(8,1,2,2)=1
        sub_table(8,1,2,3)=1
        sub_table(8,1,3,1)=4
        sub_table(8,1,3,2)=1
        sub_table(8,1,3,3)=1
        sub_table(8,1,3,4)=1
        sub_table(8,1,3,5)=1

        sub_table(9,1,1,1)=2
        sub_table(9,1,1,2)=1
        sub_table(9,1,1,3)=1
        sub_table(9,1,2,1)=2
        sub_table(9,1,2,2)=1
        sub_table(9,1,2,3)=1
        sub_table(9,1,3,1)=4
        sub_table(9,1,3,2)=1
        sub_table(9,1,3,3)=1
        sub_table(9,1,3,4)=2
        sub_table(9,1,3,5)=3

        sub_table(10,1,1,1)=2
        sub_table(10,1,1,2)=1
        sub_table(10,1,1,3)=1
        sub_table(10,1,2,1)=2
        sub_table(10,1,2,2)=2
        sub_table(10,1,2,3)=2
        sub_table(10,1,3,1)=4
        sub_table(10,1,3,2)=1
        sub_table(10,1,3,3)=1
        sub_table(10,1,3,4)=2
        sub_table(10,1,3,5)=2

        sub_table(28,1,1,1)=2
        sub_table(28,1,1,2)=1
        sub_table(28,1,1,3)=1
        sub_table(28,1,2,1)=2
        sub_table(28,1,2,2)=1
        sub_table(28,1,2,3)=1
        sub_table(28,1,3,1)=4
        sub_table(28,1,3,2)=2
        sub_table(28,1,3,3)=2
        sub_table(28,1,3,4)=3
        sub_table(28,1,3,5)=3


     CASE(32)
!
!  O_h
!
        sub_table(2,1,1,1)=2
        sub_table(2,1,1,2)=1
        sub_table(2,1,1,3)=1
        sub_table(2,1,2,1)=2
        sub_table(2,1,2,2)=1
        sub_table(2,1,2,3)=1
        sub_table(2,1,3,1)=4
        sub_table(2,1,3,2)=1
        sub_table(2,1,3,3)=1
        sub_table(2,1,3,4)=1
        sub_table(2,1,3,5)=1
        sub_table(2,1,4,1)=2
        sub_table(2,1,4,2)=2
        sub_table(2,1,4,3)=2
        sub_table(2,1,5,1)=2
        sub_table(2,1,5,2)=2
        sub_table(2,1,5,3)=2
        sub_table(2,1,6,1)=4
        sub_table(2,1,6,2)=2
        sub_table(2,1,6,3)=2
        sub_table(2,1,6,4)=2
        sub_table(2,1,6,5)=2

        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=2
        sub_table(3,1,2,2)=1
        sub_table(3,1,2,3)=2
        sub_table(3,1,3,1)=4
        sub_table(3,1,3,2)=1
        sub_table(3,1,3,3)=1
        sub_table(3,1,3,4)=2
        sub_table(3,1,3,5)=2
        sub_table(3,1,4,1)=2
        sub_table(3,1,4,2)=1
        sub_table(3,1,4,3)=2
        sub_table(3,1,5,1)=2
        sub_table(3,1,5,2)=1
        sub_table(3,1,5,3)=2
        sub_table(3,1,6,1)=4
        sub_table(3,1,6,2)=1
        sub_table(3,1,6,3)=1
        sub_table(3,1,6,4)=2
        sub_table(3,1,6,5)=2

        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2
        sub_table(4,1,3,1)=4
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=1
        sub_table(4,1,3,4)=2
        sub_table(4,1,3,5)=2
        sub_table(4,1,4,1)=2
        sub_table(4,1,4,2)=1
        sub_table(4,1,4,3)=2
        sub_table(4,1,5,1)=2
        sub_table(4,1,5,2)=1
        sub_table(4,1,5,3)=2
        sub_table(4,1,6,1)=4
        sub_table(4,1,6,2)=1
        sub_table(4,1,6,3)=1
        sub_table(4,1,6,4)=2
        sub_table(4,1,6,5)=2

        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=2
        sub_table(5,1,2,2)=1
        sub_table(5,1,2,3)=2
        sub_table(5,1,3,1)=4
        sub_table(5,1,3,2)=1
        sub_table(5,1,3,3)=2
        sub_table(5,1,3,4)=3
        sub_table(5,1,3,5)=4
        sub_table(5,1,4,1)=2
        sub_table(5,1,4,2)=1
        sub_table(5,1,4,3)=2
        sub_table(5,1,5,1)=2
        sub_table(5,1,5,2)=1
        sub_table(5,1,5,3)=2
        sub_table(5,1,6,1)=4
        sub_table(5,1,6,2)=1
        sub_table(5,1,6,3)=2
        sub_table(5,1,6,4)=3
        sub_table(5,1,6,5)=4

        sub_table(6,1,1,1)=2
        sub_table(6,1,1,2)=1
        sub_table(6,1,1,3)=2
        sub_table(6,1,2,1)=2
        sub_table(6,1,2,2)=3
        sub_table(6,1,2,3)=4
        sub_table(6,1,3,1)=4
        sub_table(6,1,3,2)=1
        sub_table(6,1,3,3)=2
        sub_table(6,1,3,4)=3
        sub_table(6,1,3,5)=4
        sub_table(6,1,4,1)=2
        sub_table(6,1,4,2)=1
        sub_table(6,1,4,3)=2
        sub_table(6,1,5,1)=2
        sub_table(6,1,5,2)=3
        sub_table(6,1,5,3)=4
        sub_table(6,1,6,1)=4
        sub_table(6,1,6,2)=1
        sub_table(6,1,6,3)=2
        sub_table(6,1,6,4)=3
        sub_table(6,1,6,5)=4

        sub_table(8,1,1,1)=2
        sub_table(8,1,1,2)=1
        sub_table(8,1,1,3)=1
        sub_table(8,1,2,1)=2
        sub_table(8,1,2,2)=1
        sub_table(8,1,2,3)=1
        sub_table(8,1,3,1)=4
        sub_table(8,1,3,2)=1
        sub_table(8,1,3,3)=1
        sub_table(8,1,3,4)=1
        sub_table(8,1,3,5)=1
        sub_table(8,1,4,1)=2
        sub_table(8,1,4,2)=1
        sub_table(8,1,4,3)=1
        sub_table(8,1,5,1)=2
        sub_table(8,1,5,2)=1
        sub_table(8,1,5,3)=1
        sub_table(8,1,6,1)=4
        sub_table(8,1,6,2)=1
        sub_table(8,1,6,3)=1
        sub_table(8,1,6,4)=1
        sub_table(8,1,6,5)=1

        sub_table(9,1,1,1)=2
        sub_table(9,1,1,2)=1
        sub_table(9,1,1,3)=1
        sub_table(9,1,2,1)=2
        sub_table(9,1,2,2)=1
        sub_table(9,1,2,3)=1
        sub_table(9,1,3,1)=4
        sub_table(9,1,3,2)=1
        sub_table(9,1,3,3)=1
        sub_table(9,1,3,4)=2
        sub_table(9,1,3,5)=3
        sub_table(9,1,4,1)=2
        sub_table(9,1,4,2)=1
        sub_table(9,1,4,3)=1
        sub_table(9,1,5,1)=2
        sub_table(9,1,5,2)=1
        sub_table(9,1,5,3)=1
        sub_table(9,1,6,1)=4
        sub_table(9,1,6,2)=1
        sub_table(9,1,6,3)=1
        sub_table(9,1,6,4)=2
        sub_table(9,1,6,5)=3

        sub_table(10,1,1,1)=2
        sub_table(10,1,1,2)=1
        sub_table(10,1,1,3)=1
        sub_table(10,1,2,1)=2
        sub_table(10,1,2,2)=2
        sub_table(10,1,2,3)=2
        sub_table(10,1,3,1)=4
        sub_table(10,1,3,2)=1
        sub_table(10,1,3,3)=1
        sub_table(10,1,3,4)=2
        sub_table(10,1,3,5)=2
        sub_table(10,1,4,1)=2
        sub_table(10,1,4,2)=1
        sub_table(10,1,4,3)=1
        sub_table(10,1,5,1)=2
        sub_table(10,1,5,2)=2
        sub_table(10,1,5,3)=2
        sub_table(10,1,6,1)=4
        sub_table(10,1,6,2)=1
        sub_table(10,1,6,3)=1
        sub_table(10,1,6,4)=2
        sub_table(10,1,6,5)=2

        sub_table(12,1,1,1)=2
        sub_table(12,1,1,2)=1
        sub_table(12,1,1,3)=1
        sub_table(12,1,2,1)=2
        sub_table(12,1,2,2)=1
        sub_table(12,1,2,3)=1
        sub_table(12,1,3,1)=4
        sub_table(12,1,3,2)=1
        sub_table(12,1,3,3)=1
        sub_table(12,1,3,4)=1
        sub_table(12,1,3,5)=1
        sub_table(12,1,4,1)=2
        sub_table(12,1,4,2)=1
        sub_table(12,1,4,3)=1
        sub_table(12,1,5,1)=2
        sub_table(12,1,5,2)=1
        sub_table(12,1,5,3)=1
        sub_table(12,1,6,1)=4
        sub_table(12,1,6,2)=1
        sub_table(12,1,6,3)=1
        sub_table(12,1,6,4)=1
        sub_table(12,1,6,5)=1
 
        sub_table(13,1,1,1)=2
        sub_table(13,1,1,2)=1
        sub_table(13,1,1,3)=1
        sub_table(13,1,2,1)=2
        sub_table(13,1,2,2)=1
        sub_table(13,1,2,3)=1
        sub_table(13,1,3,1)=4
        sub_table(13,1,3,2)=1
        sub_table(13,1,3,3)=1
        sub_table(13,1,3,4)=2
        sub_table(13,1,3,5)=3
        sub_table(13,1,4,1)=2
        sub_table(13,1,4,2)=1
        sub_table(13,1,4,3)=1
        sub_table(13,1,5,1)=2
        sub_table(13,1,5,2)=1
        sub_table(13,1,5,3)=1
        sub_table(13,1,6,1)=4
        sub_table(13,1,6,2)=1
        sub_table(13,1,6,3)=1
        sub_table(13,1,6,4)=2
        sub_table(13,1,6,5)=3

        sub_table(14,1,1,1)=2
        sub_table(14,1,1,2)=1
        sub_table(14,1,1,3)=1
        sub_table(14,1,2,1)=2
        sub_table(14,1,2,2)=2
        sub_table(14,1,2,3)=2
        sub_table(14,1,3,1)=4
        sub_table(14,1,3,2)=1
        sub_table(14,1,3,3)=1
        sub_table(14,1,3,4)=2
        sub_table(14,1,3,5)=2
        sub_table(14,1,4,1)=2
        sub_table(14,1,4,2)=1
        sub_table(14,1,4,3)=1
        sub_table(14,1,5,1)=2
        sub_table(14,1,5,2)=2
        sub_table(14,1,5,3)=2
        sub_table(14,1,6,1)=4
        sub_table(14,1,6,2)=1
        sub_table(14,1,6,3)=1
        sub_table(14,1,6,4)=2
        sub_table(14,1,6,5)=2

        sub_table(16,1,1,1)=2
        sub_table(16,1,1,2)=1
        sub_table(16,1,1,3)=2
        sub_table(16,1,2,1)=2
        sub_table(16,1,2,2)=1
        sub_table(16,1,2,3)=2
        sub_table(16,1,3,1)=4
        sub_table(16,1,3,2)=1
        sub_table(16,1,3,3)=1
        sub_table(16,1,3,4)=2
        sub_table(16,1,3,5)=2
        sub_table(16,1,4,1)=2
        sub_table(16,1,4,2)=3
        sub_table(16,1,4,3)=4
        sub_table(16,1,5,1)=2
        sub_table(16,1,5,2)=3
        sub_table(16,1,5,3)=4
        sub_table(16,1,6,1)=4
        sub_table(16,1,6,2)=3
        sub_table(16,1,6,3)=3
        sub_table(16,1,6,4)=4
        sub_table(16,1,6,5)=4
 
        sub_table(18,1,1,1)=2
        sub_table(18,1,1,2)=1
        sub_table(18,1,1,3)=2
        sub_table(18,1,2,1)=2
        sub_table(18,1,2,2)=3
        sub_table(18,1,2,3)=4
        sub_table(18,1,3,1)=4
        sub_table(18,1,3,2)=1
        sub_table(18,1,3,3)=2
        sub_table(18,1,3,4)=3
        sub_table(18,1,3,5)=4
        sub_table(18,1,4,1)=2
        sub_table(18,1,4,2)=5
        sub_table(18,1,4,3)=6
        sub_table(18,1,5,1)=2
        sub_table(18,1,5,2)=7
        sub_table(18,1,5,3)=8
        sub_table(18,1,6,1)=4
        sub_table(18,1,6,2)=5
        sub_table(18,1,6,3)=6
        sub_table(18,1,6,4)=7
        sub_table(18,1,6,5)=8

        sub_table(20,1,1,1)=2
        sub_table(20,1,1,2)=1
        sub_table(20,1,1,3)=1
        sub_table(20,1,2,1)=2
        sub_table(20,1,2,2)=1
        sub_table(20,1,2,3)=1
        sub_table(20,1,3,1)=4
        sub_table(20,1,3,2)=1
        sub_table(20,1,3,3)=1
        sub_table(20,1,3,4)=1
        sub_table(20,1,3,5)=1
        sub_table(20,1,4,1)=2
        sub_table(20,1,4,2)=2
        sub_table(20,1,4,3)=2
        sub_table(20,1,5,1)=2
        sub_table(20,1,5,2)=2
        sub_table(20,1,5,3)=2
        sub_table(20,1,6,1)=4
        sub_table(20,1,6,2)=2
        sub_table(20,1,6,3)=2
        sub_table(20,1,6,4)=2
        sub_table(20,1,6,5)=2

        sub_table(22,1,1,1)=2
        sub_table(22,1,1,2)=1
        sub_table(22,1,1,3)=1
        sub_table(22,1,2,1)=2
        sub_table(22,1,2,2)=2
        sub_table(22,1,2,3)=2
        sub_table(22,1,3,1)=4
        sub_table(22,1,3,2)=1
        sub_table(22,1,3,3)=1
        sub_table(22,1,3,4)=2
        sub_table(22,1,3,5)=2
        sub_table(22,1,4,1)=2
        sub_table(22,1,4,2)=3
        sub_table(22,1,4,3)=3
        sub_table(22,1,5,1)=2
        sub_table(22,1,5,2)=4
        sub_table(22,1,5,3)=4
        sub_table(22,1,6,1)=4
        sub_table(22,1,6,2)=3
        sub_table(22,1,6,3)=3
        sub_table(22,1,6,4)=4
        sub_table(22,1,6,5)=4

        sub_table(24,1,1,1)=2
        sub_table(24,1,1,2)=1
        sub_table(24,1,1,3)=1
        sub_table(24,1,2,1)=2
        sub_table(24,1,2,2)=2
        sub_table(24,1,2,3)=2
        sub_table(24,1,3,1)=4
        sub_table(24,1,3,2)=1
        sub_table(24,1,3,3)=1
        sub_table(24,1,3,4)=2
        sub_table(24,1,3,5)=2
        sub_table(24,1,4,1)=2
        sub_table(24,1,4,2)=1
        sub_table(24,1,4,3)=1
        sub_table(24,1,5,1)=2
        sub_table(24,1,5,2)=2
        sub_table(24,1,5,3)=2
        sub_table(24,1,6,1)=4
        sub_table(24,1,6,2)=1
        sub_table(24,1,6,3)=1
        sub_table(24,1,6,4)=2
        sub_table(24,1,6,5)=2

        sub_table(25,1,1,1)=2
        sub_table(25,1,1,2)=1
        sub_table(25,1,1,3)=1
        sub_table(25,1,2,1)=2
        sub_table(25,1,2,2)=1
        sub_table(25,1,2,3)=1
        sub_table(25,1,3,1)=4
        sub_table(25,1,3,2)=1
        sub_table(25,1,3,3)=1
        sub_table(25,1,3,4)=2
        sub_table(25,1,3,5)=3
        sub_table(25,1,4,1)=2
        sub_table(25,1,4,2)=4
        sub_table(25,1,4,3)=4
        sub_table(25,1,5,1)=2
        sub_table(25,1,5,2)=4
        sub_table(25,1,5,3)=4
        sub_table(25,1,6,1)=4
        sub_table(25,1,6,2)=4
        sub_table(25,1,6,3)=4
        sub_table(25,1,6,4)=5
        sub_table(25,1,6,5)=6

        sub_table(26,1,1,1)=2
        sub_table(26,1,1,2)=1
        sub_table(26,1,1,3)=2
        sub_table(26,1,2,1)=2
        sub_table(26,1,2,2)=3
        sub_table(26,1,2,3)=4
        sub_table(26,1,3,1)=4
        sub_table(26,1,3,2)=1
        sub_table(26,1,3,3)=2
        sub_table(26,1,3,4)=3
        sub_table(26,1,3,5)=4
        sub_table(26,1,4,1)=2
        sub_table(26,1,4,2)=1
        sub_table(26,1,4,3)=2
        sub_table(26,1,5,1)=2
        sub_table(26,1,5,2)=3
        sub_table(26,1,5,3)=4
        sub_table(26,1,6,1)=4
        sub_table(26,1,6,2)=1
        sub_table(26,1,6,3)=2
        sub_table(26,1,6,4)=3
        sub_table(26,1,6,5)=4

        sub_table(27,1,1,1)=2
        sub_table(27,1,1,2)=1
        sub_table(27,1,1,3)=2
        sub_table(27,1,2,1)=2
        sub_table(27,1,2,2)=1
        sub_table(27,1,2,3)=2
        sub_table(27,1,3,1)=4
        sub_table(27,1,3,2)=1
        sub_table(27,1,3,3)=2
        sub_table(27,1,3,4)=3
        sub_table(27,1,3,5)=3
        sub_table(27,1,4,1)=2
        sub_table(27,1,4,2)=4
        sub_table(27,1,4,3)=5
        sub_table(27,1,5,1)=2
        sub_table(27,1,5,2)=4
        sub_table(27,1,5,3)=5
        sub_table(27,1,6,1)=4
        sub_table(27,1,6,2)=4
        sub_table(27,1,6,3)=5
        sub_table(27,1,6,4)=6
        sub_table(27,1,6,5)=6

        sub_table(28,1,1,1)=2
        sub_table(28,1,1,2)=1
        sub_table(28,1,1,3)=1
        sub_table(28,1,2,1)=2
        sub_table(28,1,2,2)=1
        sub_table(28,1,2,3)=1
        sub_table(28,1,3,1)=4
        sub_table(28,1,3,2)=2
        sub_table(28,1,3,3)=2
        sub_table(28,1,3,4)=3
        sub_table(28,1,3,5)=3
        sub_table(28,1,4,1)=2
        sub_table(28,1,4,2)=1
        sub_table(28,1,4,3)=1
        sub_table(28,1,5,1)=2
        sub_table(28,1,5,2)=2
        sub_table(28,1,5,3)=2
        sub_table(28,1,6,1)=4
        sub_table(28,1,6,2)=3
        sub_table(28,1,6,3)=3
        sub_table(28,1,6,4)=4
        sub_table(28,1,6,5)=4

        sub_table(29,1,1,1)=2
        sub_table(29,1,1,2)=1
        sub_table(29,1,1,3)=1
        sub_table(29,1,2,1)=2
        sub_table(29,1,2,2)=1
        sub_table(29,1,2,3)=1
        sub_table(29,1,3,1)=4
        sub_table(29,1,3,2)=2
        sub_table(29,1,3,3)=2
        sub_table(29,1,3,4)=3
        sub_table(29,1,3,5)=3
        sub_table(29,1,4,1)=2
        sub_table(29,1,4,2)=4
        sub_table(29,1,4,3)=4
        sub_table(29,1,5,1)=2
        sub_table(29,1,5,2)=4
        sub_table(29,1,5,3)=4
        sub_table(29,1,6,1)=4
        sub_table(29,1,6,2)=5
        sub_table(29,1,6,3)=5
        sub_table(29,1,6,4)=6
        sub_table(29,1,6,5)=6

        sub_table(30,1,1,1)=2
        sub_table(30,1,1,2)=1
        sub_table(30,1,1,3)=1
        sub_table(30,1,2,1)=2
        sub_table(30,1,2,2)=2
        sub_table(30,1,2,3)=2
        sub_table(30,1,3,1)=4
        sub_table(30,1,3,2)=3
        sub_table(30,1,3,3)=3
        sub_table(30,1,3,4)=3
        sub_table(30,1,3,5)=3
        sub_table(30,1,4,1)=2
        sub_table(30,1,4,2)=2
        sub_table(30,1,4,3)=2
        sub_table(30,1,5,1)=2
        sub_table(30,1,5,2)=1
        sub_table(30,1,5,3)=1
        sub_table(30,1,6,1)=4
        sub_table(30,1,6,2)=3
        sub_table(30,1,6,3)=3
        sub_table(30,1,6,4)=3
        sub_table(30,1,6,5)=3

        sub_table(31,1,1,1)=2
        sub_table(31,1,1,2)=1
        sub_table(31,1,1,3)=1
        sub_table(31,1,2,1)=2
        sub_table(31,1,2,2)=2
        sub_table(31,1,2,3)=2
        sub_table(31,1,3,1)=4
        sub_table(31,1,3,2)=3
        sub_table(31,1,3,3)=3
        sub_table(31,1,3,4)=3
        sub_table(31,1,3,5)=3
        sub_table(31,1,4,1)=2
        sub_table(31,1,4,2)=1
        sub_table(31,1,4,3)=1
        sub_table(31,1,5,1)=2
        sub_table(31,1,5,2)=2
        sub_table(31,1,5,3)=2
        sub_table(31,1,6,1)=4
        sub_table(31,1,6,2)=3
        sub_table(31,1,6,3)=3
        sub_table(31,1,6,4)=3
        sub_table(31,1,6,5)=3

     CASE DEFAULT 
        CALL errore('convert_one_rap_so','Input point group uncorrect',1)
  END SELECT

  ndeg = sub_table(group_out, 1, rap, 1)
  IF (ndeg==0) THEN
     WRITE(stdout,'("group_in, group_out, representation",3i8)') group_in, &
                                           group_out, rap
     CALL errore('convert_one_rap_so','problem representation not found',1)
  END IF
  DO ideg=1, ndeg
     rap_list(ideg) = sub_table(group_out, 1, rap, ideg+1)
  ENDDO

  RETURN

  END SUBROUTINE convert_one_rap_so


  SUBROUTINE find_aux_ind_two_groups(nsym_a, nsym_b, sk_a, sk_b, at, bg, &
                                     group_a, group_b, aux_ind)
!
!  This routine assumes that the point group_b is a subgroup of point
!  group_a and find the auxiliary index that tells which type of 
!  subgroup it is. It receives as input the rotation matrices of 
!  both groups and uses them when necessary to distinguish the different 
!  options. The codes of the point group are:
!
!   1  "C_1 "     11 "D_6 "     21 "D_3h"     31 "O   " 
!   2  "C_i "     12 "C_2v"     22 "D_4h"     32 "O_h "  
!   3  "C_s "     13 "C_3v"     23 "D_6h" 
!   4  "C_2 "     14 "C_4v"     24 "D_2d" 
!   5  "C_3 "     15 "C_6v"     25 "D_3d" 
!   6  "C_4 "     16 "C_2h"     26 "S_4 " 
!   7  "C_6 "     17 "C_3h"     27 "S_6 " 
!   8  "D_2 "     18 "C_4h"     28 "T   " 
!   9  "D_3 "     19 "C_6h"     29 "T_h " 
!   10 "D_4 "     20 "D_2h"     30 "T_d "
!
!  The possible subgroups are the following:
!
!  1 C_1, 2 C_i, 3 C_s, 4 C_2, 5 C_3 : (1)
!  C_1
!
!  6 C_4 : (2)
!  C_1, C_2 
!
!  7 C_6 : (3)
!  C_1, C_2, C_3
!
!  8 D_2 : (4)
!  C_1, C_2_1, C_2_2, C_2_3
!
!  9 D_3 : (3)
!  C_1, C_2, C_3
!
!  10 D_4 : (7)
!  C_1, C_2_1, C_2_2, C_2_3, D_2_1, D_2_2, C_4
!
!  11 D_6 : (8)
!  C_1, C_2_1, C_2_2, C_2_3, C_3, C_6, D_3_1, D_3_2
!
!  12 C_2v : (4)
!  C_1, C_s_1, C_s_2, C_2
!
!  13 C_3v : (3)
!  C_1, C_s, C_3
!
!  14 C_4v : (7)
!  C_1, C_s_1, C_s_2, C_2, C_4, C_2v_1, C_2v_2
!
!  15 C_6v : (9)
!  C_1, C_s_1, C_s_2, C_2, C_3, C_6, C_2v, C_3v_1, C_3v_2
!
!  16 C_2h : (4)
!  C_1, C_i, C_s, C_2
!
!  17 C_3h : (3)
!  C_1, C_s, C_3
!
!  18 C_4h : (7)
!  C_1, C_i, C_s, C_2, C_4, C_2h, S_4
! 
!  19 C_6h : (9)
!  C_1, C_i, C_2, C_s, C_3, C_6, C_2h, C_3h, S_6
!
!  20 D_2h : (15)
!  C_1, C_i, C_s_1, C_s_2, C_s_3, C_2_1, C_2_2, C_2_3, D_2, C_2h_1, 
!  C_2h_2, C_2h_3, C_2v_1, C_2v_2, C_2v_3 
!
!  21 D_3h : (9)
!  C_1, C_s_1, C_s_2, C_2, C_3, D_3, C_2v, C_3v, C_3h
!
!  22 D_4h : (26)
!  C_1, C_i, C_s_1, C_s_2, C_s_3, C_2_1, C_2_2, C_2_3, C_4, D_2_1, D_2_2, D_4,
!  C_2v_1, C_2v_2, C_2v_3, C_2v_4, C_4v, C_2h_1, C_2h_2, C_2h_3, C_4h, D_2h_1, 
!  D_2h_2, D_2d_1, D_2d_2, S_4
!
!  23 D_6h : (31)
!  C_1, C_i, C_s_1, C_s_2, C_s_3, C_2_1, C_2_2, C_2_3, C_3, C_6, D_2, D_3_1, 
!  D_3_2, D_6, C_2v_1, C_2v_2, C_2v_3, C_3v_1, C_3v_2, C_6v, C_2h_1, C_2h_2, 
!  C_2h_3, C_3h, C_6h, D_2h, D_3h_1, D_3h_2, D_3d_1, D_3d_2, S_6 
!
!  24 D_2d : (7)
!  C_1, C_s, C_2_1, C_2_2, D_2, C_2v, S_4
!  
!  25 D_3d : (9)
!  C_1, C_i, C_s, C_2, C_3, D_3, C_3v, C_2h, S_6
!
!  26 S_4 : (2)
!  C_1, C_2
!
!  27 S_6 : (3)
!  C_1, C_i, C_3
!
!  28 T : (4)
!  C_1, C_2, C_3, D_2
!
!  29 T_h : (11)
!  C_1, C_i, C_s, C_2, C_3, D_2, C_2v, C_2h, D_2h, S_6, T
!
!  30 T_d : (10)
!  C_1, C_s, C_2, C_3, D_2, C_2v, C_3v, D_2d, S_4, T
!
!  31 O : (10)
!  C_1, C_2_1, C_2_2, C_3, C_4, D_2_1, D_2_2, D_3, D_4, T
!
!  32 O_h : (28)
!  C_1, C_i, C_s, C_2, C_3, C_4, D_2, D_3, D_4, C_2v_1, C_2v_2, C_2v_3, 
!  C_3v, C_4v, C_2h_1, C_2h_2, C_4h, D_2h_1, D_2h_2, D_4h, D_2d, D_3d, 
!  S_4, S_6, T, T_h, T_d, O
!

  USE kinds, ONLY : DP
  USE constants, ONLY : pi
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nsym_a, nsym_b, group_a, group_b
  INTEGER, INTENT(IN) :: sk_a(3,3,nsym_a), sk_b(3,3,nsym_b)
  INTEGER, INTENT(OUT) :: aux_ind
  LOGICAL :: is_axis
  REAL(DP), INTENT(IN) :: at(3,3), bg(3,3)
  REAL(DP) :: angle, prod

  REAL(DP) :: sr_a(3,3,nsym_a), sr_b(3,3,nsym_b), ax(3), bx(3), saxis(3,3)
  LOGICAL :: equal
  INTEGER :: isym, ipol, jpol, imirror, iaxis, four_axis
  INTEGER :: xaxis, yaxis, zaxis, naxis, isave
  LOGICAL :: is_parallel
  INTEGER :: tipo_sym

  IF (group_b==1) THEN
!
!  C_1 is a subgroup of any groups
!
     aux_ind=1
     RETURN
  ENDIF

  CALL transform_s_to_cart(sk_a, sr_a, nsym_a, at, bg)
  CALL transform_s_to_cart(sk_b, sr_b, nsym_b, at, bg)

  SELECT CASE (group_a) 
     CASE(1,2,3,4,5)
!
!   C_1, C_i, C_s, C_2, C_3  have no subgroups except C_1
!
        SELECT CASE (group_b)
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(6)
!
!   C_4 
!
        SELECT CASE (group_b)
           CASE (2)
              aux_ind=1
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(7)
!
!   C_6
!
        SELECT CASE (group_b)
           CASE(4,5)
              aux_ind=1
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(8)
!
!  D_2
!
        SELECT CASE (group_b)
           CASE(2)
      !
      !    C_2
      !
              DO isym=1,nsym_b
!
!   find the axis of order 2 and check to which axis it is parallel
!
                IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                   CALL versor(sr_b(1,1,isym),ax)
                   IF (is_axis(ax,1)) THEN
                      aux_ind=3
                   ELSEIF (is_axis(ax,2)) THEN
                      aux_ind=2
                   ELSEIF (is_axis(ax,3)) THEN
                      aux_ind=1
                   ENDIF
                ENDIF
             ENDDO

             
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   

     CASE(9)
!
!  D_3
!
        SELECT CASE (group_b)
           CASE(4,5)
      !
      !    C_2, C_3
      !
              aux_ind=1
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(10)
!
!  D_4
!
        SELECT CASE (group_b)
           CASE(2)
              DO isym=1,nsym_b
!
!   find the axis of order 2 and check to which axis it is parallel
!
                 IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                    CALL versor(sr_b(1,1,isym),ax)
                    IF (is_axis(ax,1).OR.is_axis(ax,2)) THEN
                       aux_ind=2
                    ELSEIF (is_axis(ax,3)) THEN
                       aux_ind=1
                    ELSE
                       aux_ind=3
                    ENDIF
                 ENDIF
              ENDDO
           CASE(6)
      !
      !    C_4
      !
              aux_ind=1
           CASE(8)
      !
      !  D_2
      !
              aux_ind = 1
              DO isym = 1, nsym_b
      !
      !   check if there is an axis not parallel to x, y, or z, in this case
      !   it is of type 2
      !
                 IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                    CALL versor(sr_b(1,1,isym),ax)
                    IF (.NOT.( is_axis(ax,1) .OR. is_axis(ax,2) &
                                             .OR. is_axis(ax,3) )) aux_ind=2
                 END IF
              END DO
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   

     CASE(11)
!
!  D_6
!
        SELECT CASE (group_b)
           CASE(2)
              DO isym=1,nsym_b
!
!   find the axis of order 2 and check to which axis it is parallel
!
                 IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                    CALL versor(sr_b(1,1,isym),ax)
                    IF (is_axis(ax,3)) THEN
                       aux_ind=1
                    ELSE
                       angle=ACOS(ax(1))*180.0_DP / pi
                       IF (MOD(NINT(angle), 60)==0) THEN
                          aux_ind=2
                       ELSEIF (MOD(NINT(angle), 30)==0) THEN
                          aux_ind=3
                       ELSE
                          CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
                       ENDIF
                    ENDIF
                 ENDIF
              ENDDO
           CASE(5,7)
      !
      !    C_3, C_6
      !
              aux_ind=1
           CASE(9)
      !
      !    D_3
      !
              DO isym=1,nsym_b
!
!   find the axis of order 2 and check to which axis it is parallel
!
                 IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                    CALL versor(sr_b(1,1,isym),ax)
                    angle=ACOS(ax(1))*180.0_DP / pi
                    IF (MOD(NINT(angle), 60)==0) THEN
                       aux_ind=1
                    ELSEIF (MOD(NINT(angle), 30)==0) THEN
                       aux_ind=2
                    ELSE
                       CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
                    ENDIF
                 ENDIF
              ENDDO
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(12)
!
!  C_2v
!
        SELECT CASE (group_b)
           CASE (3)
      !
      !  C_s, there is no way to distinguish the two mirrors. We assume
      !  aux_ind=1 if the first mirror of group a is the same as the C_s mirror.
      !
              DO isym=1, nsym_a
                 IF (tipo_sym(sr_a(1,1,isym))==5) THEN
                    imirror = isym
                    EXIT
                 ENDIF
              ENDDO
              equal=.TRUE.
              DO ipol=1,3
                 DO jpol=1,3
                    equal=equal.AND.sk_a(ipol,jpol,imirror)==sk_b(ipol,jpol,2)
                 ENDDO
              ENDDO
              IF (equal) THEN
                 aux_ind=1
              ELSE
                 aux_ind=2
              ENDIF
           CASE (4)
      !
      !    C_2
      !
              aux_ind=1
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(13)
!
!  C_3v
!
        SELECT CASE (group_b)
           CASE(3,5)
      !
      !    C_s, C_3
      !
              aux_ind=1
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(14)
!
!  C_4v
!
        SELECT CASE (group_b)
           CASE(3)
        !
        !  C_s
        !
        !  The second operation is the mirror
        !
             CALL mirror_axis(sr_b(1,1,2),ax)

             IF (is_axis(ax,1).OR.is_axis(ax,2).OR.is_axis(ax,3)) THEN
                aux_ind=1
             ELSE
                aux_ind=2
             END IF

           CASE(4,6)
      !
      !    C_2, C_4
      !
              aux_ind=1
           CASE(12)
      !
      !   C_2v
      !   There are two cases, mirror perpendicular to the axis x, y or
      !   to x=y and x=-y
      !
              DO isym=1,nsym_b
!
!   find one of the mirrors
!
                 IF (tipo_sym(sr_b(1,1,isym))==5) imirror=isym
              ENDDO
              CALL mirror_axis(sr_b(1,1,imirror),bx)
              IF (is_axis(bx,1).OR.is_axis(bx,2).OR.is_axis(bx,3)) THEN
                 aux_ind=1
              ELSE
                 aux_ind=2
              ENDIF
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(15)
!
!  C_6v
!
        SELECT CASE (group_b)
           CASE(3)
      !
      !    C_s
      !
              CALL mirror_axis(sr_b(1,1,2),bx)
              angle=ACOS(bx(1))*180.0_DP / pi
              IF (MOD(NINT(angle), 60)==0) THEN
                 aux_ind=1
              ELSEIF (MOD(NINT(angle), 30)==0) THEN
                 aux_ind=2
              ELSE
                 CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
              ENDIF
           CASE(4,5,7,12)
      !
      !    C_2, C_3, C_6, C_2v
      !
              aux_ind=1
           CASE (13)
      !
      !  C_3v
      !  find one mirror and check the angle of its axis with the x axis
      !
              DO isym=1,nsym_b
                 IF (tipo_sym(sr_b(1,1,isym))==5) imirror=isym
              ENDDO
              CALL mirror_axis(sr_b(1,1,imirror),bx)
              angle=ACOS(bx(1))*180.0_DP / pi
              IF (MOD(NINT(angle), 60)==0) THEN
                 aux_ind=1
              ELSEIF (MOD(NINT(angle), 30)==0) THEN
                 aux_ind=2
              ELSE
                 CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
              ENDIF
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(16)
!
!  C_2h
!
        SELECT CASE (group_b)
           CASE(2,3,4)
      !
      !    C_i, C_s, C_2
      !
              aux_ind=1
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(17)
!
!  C_3h
!
        SELECT CASE (group_b)
           CASE(3,5)
      !
      !    C_s, C_3
      !
              aux_ind=1
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(18)
!
!  C_4h
!
        SELECT CASE (group_b)
           CASE(2,3,4,6,16,26)
      !
      !    C_i, C_s, C_2, C_4, C_2h, S_4
      !
              aux_ind=1
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(19)
!
!  C_6h
!
        SELECT CASE (group_b)
           CASE(2,3,4,5,7,16,17,27)
      !
      !    C_i, C_s, C_2, C_3, C_6, C_2h, C_3h, S_6
      !
              aux_ind=1

           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(20)
!
!  D_2h
!
       SELECT CASE (group_b)
          CASE(2)
             aux_ind=1
          CASE(3)
!
!   find the mirror normal
!
             CALL mirror_axis(sr_b(1,1,2),bx)
             IF (is_axis(bx,1)) THEN
                aux_ind=1
             ELSEIF (is_axis(bx,2)) THEN
                aux_ind=2
             ELSEIF (is_axis(bx,3)) THEN
                aux_ind=3
             ENDIF
          CASE(4)
!
!   find the C_2 axis normal
!
             CALL versor(sr_b(1,1,2),ax)
             IF (is_axis(ax,1)) THEN
                aux_ind=1
             ELSEIF (is_axis(ax,2)) THEN
                aux_ind=2
             ELSEIF (is_axis(ax,3)) THEN
                aux_ind=3
             ENDIF
          CASE (8)
             aux_ind=1
          CASE (12)
       !
       !  C_2v
       !
       !
       !  Find three axis of order 2 of D_2h
       !
              naxis=0
              xaxis=0
              yaxis=0
              zaxis=0
              DO isym=1,nsym_a
                 IF (tipo_sym(sr_a(1,1,isym))==4) THEN
                    CALL versor(sr_a(1,1,isym),ax)
                    naxis=naxis+1
                    IF (naxis > 3) CALL errore('find_aux_ind_two_groups',&
                                               'two many axis',1)
                    saxis(:,naxis)=ax(:)
                    IF (is_axis(ax,1)) xaxis=naxis
                    IF (is_axis(ax,2)) yaxis=naxis
                    IF (is_axis(ax,3)) zaxis=naxis
                 ENDIF
              ENDDO
              IF (naxis /=3 ) CALL errore('find_aux_ind_two_groups',&
                                               'missing axis',1)
              IF (zaxis ==0 ) CALL errore('find_aux_ind_two_groups',&
                                               'missing z axis',1)
              
              write(6,*) 'xaxis,yaxis,zaxis', xaxis, yaxis, zaxis
              IF (xaxis==0.OR.yaxis==0) THEN
                 IF (zaxis==1) THEN
                    xaxis=2
                    yaxis=3
                 ELSEIF (zaxis==2) THEN
                    xaxis=3
                    yaxis=1
                 ELSEIF (zaxis==3) THEN
                    xaxis=1
                    yaxis=2
                 ENDIF
                 IF (is_right_oriented(saxis(1,xaxis), saxis(1,yaxis), &
                                             saxis(1,zaxis) ) ) THEN
                    isave=xaxis
                    xaxis=yaxis
                    yaxis=isave
                 ENDIF
              ENDIF
              write(6,*) 'xaxis,yaxis,zaxis', xaxis, yaxis, zaxis
              write(6,*) 'xaxis,', saxis(:,xaxis)
              write(6,*) 'yaxis,', saxis(:,yaxis)
              write(6,*) 'zaxis,', saxis(:,zaxis)
                 
              DO isym=1,nsym_b
!
!   find the axis of order 2 and check to which axis it is parallel
!
                IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                   CALL versor(sr_b(1,1,isym),ax)
                   WRITE(6,*) ax(:)
                   IF (is_parallel(ax,saxis(1,xaxis))) THEN
                      aux_ind=3
                   ELSEIF (is_parallel(ax,saxis(1,yaxis))) THEN
                      aux_ind=2
                   ELSEIF (is_parallel(ax,saxis(1,zaxis))) THEN
                      aux_ind=1
                   ENDIF 
                ENDIF
             ENDDO
         CASE(16)
!
!   find the C_2 axis and check with which axis it is parallel
!
             DO isym=1,nsym_b
                IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                   CALL versor(sr_b(1,1,isym),ax)
                   IF (is_axis(ax,1)) THEN
                      aux_ind=1
                   ELSEIF (is_axis(ax,2)) THEN
                      aux_ind=2
                   ELSEIF (is_axis(ax,3)) THEN
                      aux_ind=3
                   ENDIF
                ENDIF
             ENDDO
          CASE DEFAULT
             CALL errore('find_aux_ind_two_groups','This is not a subgroup',1)
       END SELECT 

     CASE(21)
!
!  D_3h
!
        SELECT CASE (group_b)
           CASE(3)
      !
      !   C_s
      !
              CALL mirror_axis(sr_b(1,1,2),bx)
              IF (is_axis(bx,3)) THEN
                 aux_ind=1
              ELSE
                 aux_ind=2
              ENDIF
           CASE(4,5,9,12,13,17)
      !
      !    C_2, C_3, D_3, C_2v, C_3v, C_3h
      !
              aux_ind=1
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(22)
!
!   D_4h
!
       SELECT CASE (group_b)
          CASE(3)
          !
          !  C_s
          !
             CALL mirror_axis(sr_b(1,1,2),bx)
             IF (is_axis(bx,3)) THEN
                aux_ind=1
             ELSEIF (is_axis(bx,1).OR.is_axis(bx,2)) THEN
                aux_ind=2
             ELSE
                aux_ind=3
             ENDIF
          CASE(4)
          !
          !  C_2
          !
             CALL versor(sr_b(1,1,2),ax)
             IF (is_axis(ax,3)) THEN
                aux_ind=1
             ELSEIF (is_axis(ax,1).OR.is_axis(ax,2)) THEN
                aux_ind=2
             ELSE
                aux_ind=3
             ENDIF
          CASE (2,6)
             aux_ind=1
          CASE (8)
             aux_ind=1
             DO isym=1,nsym_b
                IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                   CALL versor(sr_b(1,1,isym),ax)
                   IF (.NOT.(is_axis(ax,1).OR.is_axis(ax,2).OR.is_axis(ax,3))) &
                      aux_ind=2
                END IF
             ENDDO
          CASE (10)
             aux_ind=1
          CASE (12)
       !
       !  C_2v
       !
             DO isym=1,nsym_a
!
!   find the direction of the axis of order 4 in D_4h
!
                IF (tipo_sym(sr_a(1,1,isym))==3) iaxis=isym
             END DO
             CALL versor(sr_a(1,1,iaxis),ax)
             four_axis=0
             IF (is_axis(ax,1)) four_axis=1
             IF (is_axis(ax,2)) four_axis=2
             IF (is_axis(ax,3)) four_axis=3
             WRITE(6,*) 'four_axis', four_axis
             IF (four_axis==0) CALL errore('find_aux_ind_two_groups',  &
                                       'problem with fourfold axis',1)
             

             DO isym=1,nsym_b
!
!   find the axis of order 2 and one of the mirrors.
!
                IF (tipo_sym(sr_b(1,1,isym))==4) iaxis=isym
                IF (tipo_sym(sr_b(1,1,isym))==5) imirror=isym
             ENDDO
             CALL versor(sr_b(1,1,iaxis),ax)
             CALL mirror_axis(sr_b(1,1,imirror),bx)


             IF (is_axis(ax,four_axis)) THEN
!
!   aux_num 1 and 2 the twofold axis is the z axis and in 1 the mirror
!   are perpendicular to the x and y axis, in 2 the mirror are perpendicular
!   to the 110 and 1-10 directions.
!

                IF (is_axis(bx,1) .OR. is_axis(bx,2)) THEN
                   aux_ind=1
                ELSE
                   aux_ind=2
                ENDIF
             ELSEIF (is_axis(ax,1).OR.is_axis(ax,2).OR.is_axis(ax,3)) THEN
!
!  aux_num 3 or 5 when the axis is parallel to x or y
!
                IF (four_axis==1) THEN
                   IF (is_axis(ax,2)) THEN
                      aux_ind=3
                   ELSE
                      aux_ind=5
                   ENDIF
                ELSEIF (four_axis==2) THEN
                   IF (is_axis(ax,1)) THEN
                      aux_ind=5
                   ELSE
                      aux_ind=3
                   ENDIF
                ELSE
                   IF (is_axis(ax,1)) THEN
                      aux_ind=3
                   ELSE
                      aux_ind=3
                   ENDIF
                ENDIF
             ELSE
!
!  aux_num 4 or 6 when the twofold axis is parallel to 110 or 1-10 directions
!
                IF (four_axis==1) THEN
                   IF (ABS(ax(3)-ax(2))<1.d-6) THEN
                      aux_ind=4
                   ELSE 
                      aux_ind=6
                   ENDIF
                ELSEIF (four_axis==2) THEN
                   IF (ABS(ax(1)-ax(3))<1.d-6) THEN
                      aux_ind=6
                   ELSE 
                      aux_ind=4
                   ENDIF
                ELSE
                   IF (ABS(ax(1)-ax(2))<1.d-6) THEN
                      aux_ind=4
                   ELSE 
                      aux_ind=6
                   ENDIF
                ENDIF
             END IF
          CASE (16)
             DO isym=1,nsym_b
!
!   find the axis of order 2.
!
                IF (tipo_sym(sr_b(1,1,isym))==4) iaxis=isym
             ENDDO
             CALL versor(sr_b(1,1,iaxis),ax)
             IF (is_axis(ax,3)) THEN
                aux_ind=1
             ELSEIF(is_axis(ax,1).OR.is_axis(ax,2)) THEN
                aux_ind=2
             ELSE
                aux_ind=3
             ENDIF
          CASE (14,18)
       !
       !  C_4v, C_4h
       !
             aux_ind=1
          CASE(20,24)
       !
       !  D_2h, D_2d
       !
             aux_ind=1
             DO isym=1,nsym_b
                IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                   CALL versor(sr_b(1,1,isym),ax)
                   IF (.NOT.(is_axis(ax,1).OR.is_axis(ax,2).OR.is_axis(ax,3))) &
                      aux_ind=2
                END IF
             ENDDO
          CASE(26)
       !
       !  S_4
       ! 
             aux_ind=1
          CASE DEFAULT
             CALL errore('find_aux_ind_two_groups','This is not a subgroup',1)
       END SELECT 
     CASE(23)
!
!  D_6h
!
        SELECT CASE (group_b)
           CASE(2)
        !
        !  C_i
        !
              aux_ind=1
           CASE(3)
        !
        !  C_s
        ! 
              CALL mirror_axis(sr_b(1,1,2),bx)
              IF (is_axis(bx,3)) THEN
                 aux_ind=1
              ELSE
                 angle=ACOS(bx(1))*180.0_DP / pi
                 IF (MOD(NINT(angle), 60)==0) THEN
                    aux_ind=2
                 ELSEIF (MOD(NINT(angle), 30)==0) THEN
                    aux_ind=3
                 ELSE
                    CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
                 ENDIF
              END IF
           CASE(4)
        !
        !  C_2
        ! 
              CALL versor(sr_b(1,1,2),bx)
              IF (is_axis(bx,3)) THEN
                 aux_ind=1
              ELSE
                 angle=ACOS(bx(1))*180.0_DP / pi
                 IF (MOD(NINT(angle), 60)==0) THEN
                    aux_ind=2
                 ELSEIF (MOD(NINT(angle), 30)==0) THEN
                    aux_ind=3
                 ELSE
                    CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
                 ENDIF
              END IF
           CASE(5,7)
      !
      !  C_3, C_6 
      !
               aux_ind=1
           CASE(9)
      !
      !    D_3
      !
              DO isym=1,nsym_b
                 IF (tipo_sym(sr_b(1,1,isym))==4) iaxis=isym
              END DO
              CALL versor(sr_b(1,1,isym),ax)
              angle=ACOS(ax(1))*180.0_DP / pi
              IF (MOD(NINT(angle), 60)==0) THEN
                 aux_ind=1
              ELSEIF (MOD(NINT(angle), 30)==0) THEN
                 aux_ind=2
              ELSE
                    CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
              ENDIF
           CASE(12,16)
      !
      !    C_2v, C_2h
      !
              DO isym=1,nsym_b
                 IF (tipo_sym(sr_b(1,1,isym))==4) iaxis=isym
              END DO
              CALL versor(sr_b(1,1,iaxis),ax)
              IF (is_axis(ax,3)) THEN
                 aux_ind=1
              ELSE
                 angle=ACOS(ax(1))*180.0_DP / pi
                 IF (MOD(NINT(angle), 60)==0) THEN
                    aux_ind=2
                 ELSEIF (MOD(NINT(angle), 30)==0) THEN
                    aux_ind=3
                 ELSE
                       CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
                 ENDIF
              ENDIF
           CASE(13)
      !
      !     C_3v
      !
              DO isym=1,nsym_b
                 IF (tipo_sym(sr_b(1,1,isym))==5) imirror=isym
              END DO
              CALL versor(sr_b(1,1,imirror),bx)
              angle=ACOS(bx(1))*180.0_DP / pi
              IF (MOD(NINT(angle), 60)==0) THEN
                 aux_ind=1
              ELSEIF (MOD(NINT(angle), 30)==0) THEN
                 aux_ind=2
              ELSE
                 CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
              ENDIF
           CASE(21)
      !
      ! D_3h
      !
              DO isym=1,nsym_b
                 IF (tipo_sym(sr_b(1,1,isym))==4) iaxis=isym
              END DO
              CALL versor(sr_b(1,1,isym),ax)
              angle=ACOS(ax(1))*180.0_DP / pi
              IF (MOD(NINT(angle), 60)==0) THEN
                 aux_ind=1
              ELSEIF (MOD(NINT(angle), 30)==0) THEN
                 aux_ind=2
              ELSE
                 CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
              ENDIF
           CASE(25)
      !
      !     D_3d
      !
              DO isym=1,nsym_b
                 IF (tipo_sym(sr_b(1,1,isym))==4) iaxis=isym
              END DO
              CALL versor(sr_b(1,1,iaxis),ax)
              angle=ACOS(ax(1))*180.0_DP / pi
              IF (MOD(NINT(angle), 60)==0) THEN
                 aux_ind=1
              ELSEIF (MOD(NINT(angle), 30)==0) THEN
                 aux_ind=2
              ELSE
                 CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
              ENDIF
           CASE(8,11,15,17,19,20)
      !
      !     D_2, D_6, C_6v, C_3h, C_6h, D_2h
      !
              aux_ind=1

           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(24)
!
!  D_2d
!
        SELECT CASE (group_b)
           CASE(3)
              aux_ind=1
           CASE(4)
!
!   Compare the direction of the versor of C_2 and of the rotation -4 in D_2d
!   If they are perpendicular the subgroup C_2 is of aux_ind=2
!
              CALL versor(sr_b(1,1,2),ax)
              DO isym=1,nsym_a
                 IF (tipo_sym(sr_a(1,1,isym))==6)  &
                    CALL versor(sr_a(1,1,isym),bx)
              ENDDO
              prod = ax(1)*bx(1) + ax(2)*bx(2) + ax(3)*bx(3)
              IF (prod > 1.d-6) THEN
                 aux_ind=1
              ELSE
                 aux_ind=2
              ENDIF

           CASE(8,12,26)
      !
      !     C_s, D_2, C_2v, S_4
      !
              aux_ind=1

           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(25)
!
!  D_3d
!
        SELECT CASE (group_b)
           CASE(2,3,4,5,9,13,16,27)
      !
      !     C_i, C_s, C_2, C_3, D_3, C_3v, C_2h, S_6
      !
              aux_ind=1

           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(26)
!
!  S_4
!
        SELECT CASE(group_b)
           CASE (4)
       !
       !  C_2
       !
              aux_ind=1

           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(27)
!
! S_6
!
        SELECT CASE(group_b)
           CASE (2,5)
       !
       !  C_i, C_3
       !
              aux_ind=1

           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(28)
        SELECT CASE(group_b)
           CASE (4,5,8)
       !
       !  C_2, C_3, D_2
       !
              aux_ind=1

           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(29)
        SELECT CASE(group_b)
           CASE (2,3,4,5,8,12,16,20,27,28)
       !
       !  C_i, C_s, C_2, C_3, D_2, C_2v, C_2h, D_2h, S_6, T
       !
              aux_ind=1

           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(30)
        SELECT CASE(group_b)
           CASE (3,4,5,8,12,13,24,26,28)
       !
       !  C_s, C_2, C_3, D_2, C_2v, C_3v, D_2d, S_4, T
       !
              aux_ind=1

           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(31)
        SELECT CASE(group_b)
           CASE (4)
       !
       !  C_2
       !
              CALL versor(sr_b(1,1,2),ax)
              IF (is_axis(ax,1).OR.is_axis(ax,2).OR.is_axis(ax,3)) THEN
                 aux_ind=1
              ELSE
                 aux_ind=2
              ENDIF
           CASE (8)
       !
       !  D_2
       !
              aux_ind=1
              DO isym=1,nsym_b
                 IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                    CALL versor(sr_b(1,1,isym),ax)
                    IF (.NOT.(is_axis(ax,1).OR.is_axis(ax,2).OR.is_axis(ax,3))) &
                       aux_ind=2
                 END IF
              ENDDO
           CASE (5,6,9,10,28)
       !
       !  C_3, C_4, D_3, D_4, T
       !
              aux_ind=1

           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(32)
!
!  O_h
!
        SELECT CASE (group_b)
           CASE(2,3,4,5,6,8,9,10)
      !
      !   C_i, C_s, C_2, C_3, C_4, D_2, D_3, D_4
      !
               aux_ind = 1 
           CASE(12)
      !
      !     C_2v  to choose between 1 2 or 3
      !
              DO isym=1,nsym_b
!
!   find the axis of order 2 and one of the mirrors.
!
                 IF (tipo_sym(sr_b(1,1,isym))==4) iaxis=isym
                 IF (tipo_sym(sr_b(1,1,isym))==5) imirror=isym
              ENDDO
              CALL versor(sr_b(1,1,iaxis),ax)
              CALL mirror_axis(sr_b(1,1,imirror),bx)
              IF (is_axis(ax,1).OR.is_axis(ax,2).OR.is_axis(ax,3)) THEN
                 IF (is_axis(bx,1).OR.is_axis(bx,2).OR.is_axis(bx,3)) THEN
                    aux_ind=1
                 ELSE
                    aux_ind=2
                 ENDIF
              ELSE
                 aux_ind=3
              ENDIF
           CASE (16)
       !
       !  C_2h
       !
              DO isym=1,nsym_b
!
!   find the axis of order 2 and check if it is parallel to x, y, or z 
!
                 IF (tipo_sym(sr_b(1,1,isym))==4) iaxis=isym
              ENDDO
              CALL versor(sr_b(1,1,iaxis),ax)
              IF (is_axis(ax,1).OR.is_axis(ax,2).OR.is_axis(ax,3)) THEN
                 aux_ind=1
              ELSE
                 aux_ind=2
              ENDIF
          CASE (20)
       !
       !  D_2h
       !
              aux_ind=1
              DO isym=1,nsym_b
!
!   find if one axis of order 2 is not parallel to x, y, or z 
!
                 IF (tipo_sym(sr_b(1,1,isym))==4) iaxis=isym
                 CALL versor(sr_b(1,1,iaxis),ax)
                 IF (.NOT.(is_axis(ax,1).OR.is_axis(ax,2).OR.is_axis(ax,3))) &
                    aux_ind=2
              END DO
           CASE (13,14,18,22,24,25,26,27,28,29,30,31)
       !
       !  C_3v, C_4v, C_2h, D_4h, D_2d, D_3d, S_4, S_6, T, T_h, T_d, O
       !
              aux_ind=1

           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE DEFAULT
       CALL errore('find_aux_ind_two_groups',' group not available',1)
  END SELECT

  RETURN

  END SUBROUTINE find_aux_ind_two_groups

SUBROUTINE  transform_s_to_cart( sk, sr, nsym, at, bg )
  !----------------------------------------------------------------------
  !
  !     This routine transforms symmetry matrices expressed in the
  !     basis of the crystal axis into rotations in cartesian axis
  !
  USE kinds
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nsym
  INTEGER, INTENT(IN) :: sk(3,3,nsym)
  REAL(DP), INTENT(IN) :: at(3,3), bg(3,3)
  REAL(DP), INTENT(OUT) :: sr(3,3,nsym)
  INTEGER :: isym
  REAL(DP):: sa(3,3), sb(3,3)
  !
  DO isym = 1,nsym
     sa (:,:) = dble ( sk(:,:,isym) )
     sb = matmul ( bg, sa )
     sr (:,:, isym) = matmul ( at, transpose (sb) )
  ENDDO

  RETURN

  END SUBROUTINE transform_s_to_cart

  LOGICAL FUNCTION has_sigma_h(code_group)

  INTEGER, INTENT(IN) :: code_group

  SELECT CASE (code_group)
     CASE (16,17,18,19,20,21,22,23,32)
        has_sigma_h=.TRUE.
     CASE DEFAULT
        has_sigma_h=.FALSE.
  END SELECT
  RETURN
  END FUNCTION has_sigma_h

  FUNCTION is_right_oriented (a,b,c)
  !
  !  This functions receives 3 vectors a,b,c and gives as output true
  !  if they are oriented as x,y,z, .false. if they have the opposite
  !  orientation. It calculates the determinant of the the matrix of
  !  three vectors and gives true if it is positive. If it is zero
  !  the result is undermined
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: a(3), b(3), c(3)
  LOGICAL :: is_right_oriented

  REAL(DP) :: vect(3), det

  vect(1) = b(2) * c(3) - b(3) * c(1)
  vect(2) = b(3) * c(1) - b(1) * c(3)
  vect(3) = b(1) * c(2) - b(2) * c(1)

  det= vect(1) * a(1) + vect(2) * a(2) + vect(3) * a(3)

  is_right_oriented = (det > 0.0_DP)

  RETURN
  END FUNCTION is_right_oriented

  END MODULE point_group
