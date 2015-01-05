!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE point_group
!
!  This module contains variables and routines to deal with the crystallographic
!  point group symmetry. It complements the routines in find_group.f90,
!  divide_class.f90 and divide_class_so.f90 in the PW/src directory of the
!  QE package.
!  The conventions, such as the code group, the symmetry operation types,
!  the irreducible representations etc. are the same.
!  Presently it has routines to perform the following task:
!  Given two point groups, the second a subgroup of the first,
!  an a list of representations of the first point group, it
!  transform it in a list of representations of the second group.
!  Given two point groups, the second a subgroup of the first, 
!  find which type it is. The different cases are discussed in the point-group
!  manual in the thermo_pw/Doc directory. The routine find_aux_ind_two_groups,
!  receives the rotation matrices of the two group and gives an index that
!  correspond to the case.
!  Double groups are supported, however in this case the distinction between
!  different subgroup cases is irrelevant and not used.
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
  INTEGER :: ideg, iax, ibx

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
        !
        ! and C_2
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
        !
        ! and C_2
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
        !
        ! and C_3
        !
        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=1
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=2
        sub_table(5,1,4,1)=1
        sub_table(5,1,4,2)=3
        sub_table(5,1,5,1)=1
        sub_table(5,1,5,2)=3
        sub_table(5,1,6,1)=1
        sub_table(5,1,6,2)=2

     CASE(8)
!
! D_2
!
        !
        ! and C_2
        !
        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2

        sub_table(4,2,1,1)=1
        sub_table(4,2,1,2)=1
        sub_table(4,2,2,1)=1
        sub_table(4,2,2,2)=2
        sub_table(4,2,3,1)=1
        sub_table(4,2,3,2)=1
        sub_table(4,2,4,1)=1
        sub_table(4,2,4,2)=2

        sub_table(4,3,1,1)=1
        sub_table(4,3,1,2)=1
        sub_table(4,3,2,1)=1
        sub_table(4,3,2,2)=2
        sub_table(4,3,3,1)=1
        sub_table(4,3,3,2)=2
        sub_table(4,3,4,1)=1
        sub_table(4,3,4,2)=1

     CASE(9)
!
! D_3
!
        !
        !  with C_2
        !
        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=2
        sub_table(4,1,3,1)=2
        sub_table(4,1,3,2)=1
        sub_table(4,1,3,3)=2
        !
        !  with C_3
        !
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
        !
        !  with C_2
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
        !
        !  with C_4
        !
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
        !
        !  with D_2
        !
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

        sub_table(8,3,1,1)=1
        sub_table(8,3,1,2)=1
        sub_table(8,3,2,1)=1
        sub_table(8,3,2,2)=3
        sub_table(8,3,3,1)=1
        sub_table(8,3,3,2)=1
        sub_table(8,3,4,1)=1
        sub_table(8,3,4,2)=3
        sub_table(8,3,5,1)=2
        sub_table(8,3,5,2)=2
        sub_table(8,3,5,3)=4

        sub_table(8,4,1,1)=1
        sub_table(8,4,1,2)=1
        sub_table(8,4,2,1)=1
        sub_table(8,4,2,2)=4
        sub_table(8,4,3,1)=1
        sub_table(8,4,3,2)=1
        sub_table(8,4,4,1)=1
        sub_table(8,4,4,2)=4
        sub_table(8,4,5,1)=2
        sub_table(8,4,5,2)=2
        sub_table(8,4,5,3)=3

     CASE(11)
!
!  D_6
!
        !
        !  with C_2
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
        !
        !  with C_3
        !
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
        !
        !  with C_6
        !
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
        !
        !  with D_2
        !
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
        !
        !  with D_3
        !
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
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2

     CASE(13)
!
!   C_3v
!
        !
        !  with C_s
        !
        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=2
        sub_table(3,1,3,2)=1
        sub_table(3,1,3,3)=2
        !
        !  with C_3v
        !
        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=1
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=2
        sub_table(5,1,3,3)=3

     CASE(14)
!
!   C_4v
!
        !
        !  with C_s
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
        !
        !  with C_2
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
        !
        !  with C_4
        !
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
        !
        !  with C_2v
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
        !
        !  with C_s
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
        !
        !  with C_2
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
        sub_table(4,1,5,2)=2
        sub_table(4,1,5,3)=2
        sub_table(4,1,6,1)=2
        sub_table(4,1,6,2)=1
        sub_table(4,1,6,3)=1
        !
        !  with C_3
        !
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
        !
        !  with C_6
        !
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
        !
        !  with C_2v
        !
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
        !
        !  with C_3v
        !
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
        !
        !  with C_i
        !
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=1
        sub_table(2,1,3,2)=2
        sub_table(2,1,4,1)=1
        sub_table(2,1,4,2)=2
        !
        !  with C_s
        !
        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=2
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=1
        !
        !  with C_2
        !
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
        !
        !  with C_s
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
        !
        !  with C_3
        !
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
        !
        !  with C_i
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
        !
        !  with C_s
        !
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
        !
        !  with C_2
        !
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
        !
        !  with C_4
        !
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
        !
        !  with C_2h
        !
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
        !
        !  with S_4
        !
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
        !
        !  with C_i
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
        !
        !  with C_s
        !
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
        !
        !  with C_2
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
        !
        !  with C_3
        !
        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=1
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=2
        sub_table(5,1,4,1)=1
        sub_table(5,1,4,2)=3
        sub_table(5,1,5,1)=1
        sub_table(5,1,5,2)=3
        sub_table(5,1,6,1)=1
        sub_table(5,1,6,2)=2
        sub_table(5,1,7,1)=1
        sub_table(5,1,7,2)=1
        sub_table(5,1,8,1)=1
        sub_table(5,1,8,2)=1
        sub_table(5,1,9,1)=1
        sub_table(5,1,9,2)=2
        sub_table(5,1,10,1)=1
        sub_table(5,1,10,2)=3
        sub_table(5,1,11,1)=1
        sub_table(5,1,11,2)=3
        sub_table(5,1,12,1)=1
        sub_table(5,1,12,2)=2
        !
        !  with C_6
        !
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
        !
        !  with C_2h
        !
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
        !
        !  with C_3h
        !
        sub_table(17,1,1,1)=1
        sub_table(17,1,1,2)=1
        sub_table(17,1,2,1)=1
        sub_table(17,1,2,2)=4
        sub_table(17,1,3,1)=1
        sub_table(17,1,3,2)=5
        sub_table(17,1,4,1)=1
        sub_table(17,1,4,2)=6
        sub_table(17,1,5,1)=1
        sub_table(17,1,5,2)=3
        sub_table(17,1,6,1)=1
        sub_table(17,1,6,2)=2
        sub_table(17,1,7,1)=1
        sub_table(17,1,7,2)=4
        sub_table(17,1,8,1)=1
        sub_table(17,1,8,2)=1
        sub_table(17,1,9,1)=1
        sub_table(17,1,9,2)=2
        sub_table(17,1,10,1)=1
        sub_table(17,1,10,2)=3
        sub_table(17,1,11,1)=1
        sub_table(17,1,11,2)=6
        sub_table(17,1,12,1)=1
        sub_table(17,1,12,2)=5
        !
        !  with S_6
        !
        sub_table(27,1,1,1)=1
        sub_table(27,1,1,2)=1
        sub_table(27,1,2,1)=1
        sub_table(27,1,2,2)=1
        sub_table(27,1,3,1)=1
        sub_table(27,1,3,2)=2
        sub_table(27,1,4,1)=1
        sub_table(27,1,4,2)=3
        sub_table(27,1,5,1)=1
        sub_table(27,1,5,2)=3
        sub_table(27,1,6,1)=1
        sub_table(27,1,6,2)=2
        sub_table(27,1,7,1)=1
        sub_table(27,1,7,2)=4
        sub_table(27,1,8,1)=1
        sub_table(27,1,8,2)=4
        sub_table(27,1,9,1)=1
        sub_table(27,1,9,2)=5
        sub_table(27,1,10,1)=1
        sub_table(27,1,10,2)=6
        sub_table(27,1,11,1)=1
        sub_table(27,1,11,2)=6
        sub_table(27,1,12,1)=1
        sub_table(27,1,12,2)=5
 
     CASE(20)
!
!  D_2h
!
        !
        !  with C_i
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
        !
        !  with C_s
        !
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
        sub_table(3,3,2,2)=2
        sub_table(3,3,3,1)=1
        sub_table(3,3,3,2)=2
        sub_table(3,3,4,1)=1
        sub_table(3,3,4,2)=1
        sub_table(3,3,5,1)=1
        sub_table(3,3,5,2)=2
        sub_table(3,3,6,1)=1
        sub_table(3,3,6,2)=1
        sub_table(3,3,7,1)=1
        sub_table(3,3,7,2)=1
        sub_table(3,3,8,1)=1
        sub_table(3,3,8,2)=2
        !
        !  with C_2
        !
        sub_table(4,1,1,1)=1
        sub_table(4,1,1,2)=1
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        sub_table(4,1,4,1)=1
        sub_table(4,1,4,2)=2
        sub_table(4,1,5,1)=1
        sub_table(4,1,5,2)=2
        sub_table(4,1,6,1)=1
        sub_table(4,1,6,2)=2
        sub_table(4,1,7,1)=1
        sub_table(4,1,7,2)=1
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
        sub_table(4,2,5,2)=2
        sub_table(4,2,6,1)=1
        sub_table(4,2,6,2)=1
        sub_table(4,2,7,1)=1
        sub_table(4,2,7,2)=2
        sub_table(4,2,8,1)=1
        sub_table(4,2,8,2)=1

        sub_table(4,3,1,1)=1
        sub_table(4,3,1,2)=1
        sub_table(4,3,2,1)=1
        sub_table(4,3,2,2)=2
        sub_table(4,3,3,1)=1
        sub_table(4,3,3,2)=2
        sub_table(4,3,4,1)=1
        sub_table(4,3,4,2)=1
        sub_table(4,3,5,1)=1
        sub_table(4,3,5,2)=1
        sub_table(4,3,6,1)=1
        sub_table(4,3,6,2)=2
        sub_table(4,3,7,1)=1
        sub_table(4,3,7,2)=2
        sub_table(4,3,8,1)=1
        sub_table(4,3,8,2)=1
        !
        !  with D_2
        !
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
        !
        !  with C_2v
        !
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

        sub_table(12,4,1,1)=1
        sub_table(12,4,1,2)=1
        sub_table(12,4,2,1)=1
        sub_table(12,4,2,2)=3
        sub_table(12,4,3,1)=1
        sub_table(12,4,3,2)=4
        sub_table(12,4,4,1)=1
        sub_table(12,4,4,2)=2
        sub_table(12,4,5,1)=1
        sub_table(12,4,5,2)=2
        sub_table(12,4,6,1)=1
        sub_table(12,4,6,2)=4
        sub_table(12,4,7,1)=1
        sub_table(12,4,7,2)=3
        sub_table(12,4,8,1)=1
        sub_table(12,4,8,2)=1

        sub_table(12,5,1,1)=1
        sub_table(12,5,1,2)=1
        sub_table(12,5,2,1)=1
        sub_table(12,5,2,2)=4
        sub_table(12,5,3,1)=1
        sub_table(12,5,3,2)=2
        sub_table(12,5,4,1)=1
        sub_table(12,5,4,2)=3
        sub_table(12,5,5,1)=1
        sub_table(12,5,5,2)=2
        sub_table(12,5,6,1)=1
        sub_table(12,5,6,2)=3
        sub_table(12,5,7,1)=1
        sub_table(12,5,7,2)=1
        sub_table(12,5,8,1)=1
        sub_table(12,5,8,2)=4

        sub_table(12,6,1,1)=1
        sub_table(12,6,1,2)=1
        sub_table(12,6,2,1)=1
        sub_table(12,6,2,2)=4
        sub_table(12,6,3,1)=1
        sub_table(12,6,3,2)=3
        sub_table(12,6,4,1)=1
        sub_table(12,6,4,2)=2
        sub_table(12,6,5,1)=1
        sub_table(12,6,5,2)=2
        sub_table(12,6,6,1)=1
        sub_table(12,6,6,2)=3
        sub_table(12,6,7,1)=1
        sub_table(12,6,7,2)=4
        sub_table(12,6,8,1)=1
        sub_table(12,6,8,2)=1
        !
        !  with C_2h
        !
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
        sub_table(16,3,2,2)=2
        sub_table(16,3,3,1)=1
        sub_table(16,3,3,2)=2
        sub_table(16,3,4,1)=1
        sub_table(16,3,4,2)=1
        sub_table(16,3,5,1)=1
        sub_table(16,3,5,2)=3
        sub_table(16,3,6,1)=1
        sub_table(16,3,6,2)=4
        sub_table(16,3,7,1)=1
        sub_table(16,3,7,2)=4
        sub_table(16,3,8,1)=1
        sub_table(16,3,8,2)=3

     CASE(21)
!
!  D_3h
!
        !
        !  with C_s
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
        !
        !  with C_2
        !
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
        !
        !  with C_3
        !
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
        !
        !  with D_3
        !
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
        !
        !  with C_2v
        !
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

        sub_table(12,2,1,1)=1
        sub_table(12,2,1,2)=1
        sub_table(12,2,2,1)=1
        sub_table(12,2,2,2)=4
        sub_table(12,2,3,1)=2
        sub_table(12,2,3,2)=1
        sub_table(12,2,3,3)=4
        sub_table(12,2,4,1)=1
        sub_table(12,2,4,2)=2
        sub_table(12,2,5,1)=1
        sub_table(12,2,5,2)=3
        sub_table(12,2,6,1)=2
        sub_table(12,2,6,2)=2
        sub_table(12,2,6,3)=3
        !
        !  with C_3v
        !
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
        !
        !  with C_3h
        !
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
        !
        !  with C_i
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
        !
        !  with C_s
        !
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
        !
        !  with C_2
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
        !
        !  with C_4
        !
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
        !
        !  with D_2
        !
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
        !
        !  with D_4
        !
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
        !  with C_2v 
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
        sub_table(12,4,2,2)=4
        sub_table(12,4,3,1)=1
        sub_table(12,4,3,2)=1
        sub_table(12,4,4,1)=1
        sub_table(12,4,4,2)=4
        sub_table(12,4,5,1)=2
        sub_table(12,4,5,2)=2
        sub_table(12,4,5,3)=3
        sub_table(12,4,6,1)=1
        sub_table(12,4,6,2)=2
        sub_table(12,4,7,1)=1
        sub_table(12,4,7,2)=3
        sub_table(12,4,8,1)=1
        sub_table(12,4,8,2)=2
        sub_table(12,4,9,1)=1
        sub_table(12,4,9,2)=3
        sub_table(12,4,10,1)=2
        sub_table(12,4,10,2)=1
        sub_table(12,4,10,3)=4

        sub_table(12,5,1,1)=1
        sub_table(12,5,1,2)=1
        sub_table(12,5,2,1)=1
        sub_table(12,5,2,2)=3
        sub_table(12,5,3,1)=1
        sub_table(12,5,3,2)=3
        sub_table(12,5,4,1)=1
        sub_table(12,5,4,2)=1
        sub_table(12,5,5,1)=2
        sub_table(12,5,5,2)=2
        sub_table(12,5,5,3)=4
        sub_table(12,5,6,1)=1
        sub_table(12,5,6,2)=2
        sub_table(12,5,7,1)=1
        sub_table(12,5,7,2)=4
        sub_table(12,5,8,1)=1
        sub_table(12,5,8,2)=4
        sub_table(12,5,9,1)=1
        sub_table(12,5,9,2)=2
        sub_table(12,5,10,1)=2
        sub_table(12,5,10,2)=1
        sub_table(12,5,10,3)=3
        !
        !  with C_4v 
        !
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
        !
        !  with C_2h 
        !
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
        !
        !  with C_4h 
        !
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
        !
        !  with D_2h 
        !
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
        !
        !  with D_2d 
        !
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
        !
        !  with S_4 
        !
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
        !
        !  with C_i
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
        !
        !  with C_s
        !
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
        !
        !  with C_2
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
        !
        !  with C_3
        !
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
        !
        !  with C_6
        !
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
        !
        !  with D_2
        !
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
        !
        !  with D_3
        !
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
        !
        !  with D_6
        !
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
        !
        !  with C_2v
        !
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
        sub_table(12,1,9,2)=3
        sub_table(12,1,10,1)=1
        sub_table(12,1,10,2)=4
        sub_table(12,1,11,1)=2
        sub_table(12,1,11,2)=3
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
        sub_table(12,2,11,2)=1
        sub_table(12,2,11,3)=3
        sub_table(12,2,12,1)=2
        sub_table(12,2,12,2)=2
        sub_table(12,2,12,3)=4

        sub_table(12,3,1,1)=1
        sub_table(12,3,1,2)=1
        sub_table(12,3,2,1)=1
        sub_table(12,3,2,2)=4
        sub_table(12,3,3,1)=1
        sub_table(12,3,3,2)=2
        sub_table(12,3,4,1)=1
        sub_table(12,3,4,2)=3
        sub_table(12,3,5,1)=2
        sub_table(12,3,5,2)=2
        sub_table(12,3,5,3)=3
        sub_table(12,3,6,1)=2
        sub_table(12,3,6,2)=1
        sub_table(12,3,6,3)=4
        sub_table(12,3,7,1)=1
        sub_table(12,3,7,2)=2
        sub_table(12,3,8,1)=1
        sub_table(12,3,8,2)=3
        sub_table(12,3,9,1)=1
        sub_table(12,3,9,2)=1
        sub_table(12,3,10,1)=1
        sub_table(12,3,10,2)=4
        sub_table(12,3,11,1)=2
        sub_table(12,3,11,2)=1
        sub_table(12,3,11,3)=4
        sub_table(12,3,12,1)=2
        sub_table(12,3,12,2)=2
        sub_table(12,3,12,3)=3

        sub_table(12,4,1,1)=1
        sub_table(12,4,1,2)=1
        sub_table(12,4,2,1)=1
        sub_table(12,4,2,2)=3
        sub_table(12,4,3,1)=1
        sub_table(12,4,3,2)=4
        sub_table(12,4,4,1)=1
        sub_table(12,4,4,2)=2
        sub_table(12,4,5,1)=2
        sub_table(12,4,5,2)=2
        sub_table(12,4,5,3)=4
        sub_table(12,4,6,1)=2
        sub_table(12,4,6,2)=1
        sub_table(12,4,6,3)=3
        sub_table(12,4,7,1)=1
        sub_table(12,4,7,2)=2
        sub_table(12,4,8,1)=1
        sub_table(12,4,8,2)=4
        sub_table(12,4,9,1)=1
        sub_table(12,4,9,2)=3
        sub_table(12,4,10,1)=1
        sub_table(12,4,10,2)=1
        sub_table(12,4,11,1)=2
        sub_table(12,4,11,2)=1
        sub_table(12,4,11,3)=3
        sub_table(12,4,12,1)=2
        sub_table(12,4,12,2)=2
        sub_table(12,4,12,3)=4
        !
        !  with C_3v
        !
        sub_table(13,1,1,1)=1
        sub_table(13,1,1,2)=1
        sub_table(13,1,2,1)=1
        sub_table(13,1,2,2)=2
        sub_table(13,1,3,1)=1
        sub_table(13,1,3,2)=2
        sub_table(13,1,4,1)=1
        sub_table(13,1,4,2)=1
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
        sub_table(13,1,9,2)=1
        sub_table(13,1,10,1)=1
        sub_table(13,1,10,2)=2
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
        sub_table(13,2,3,2)=1
        sub_table(13,2,4,1)=1
        sub_table(13,2,4,2)=2
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
        sub_table(13,2,9,2)=2
        sub_table(13,2,10,1)=1
        sub_table(13,2,10,2)=1
        sub_table(13,2,11,1)=2
        sub_table(13,2,11,2)=3
        sub_table(13,2,11,3)=3
        sub_table(13,2,12,1)=2
        sub_table(13,2,12,2)=3
        sub_table(13,2,12,3)=3
        !
        !  with C_6v
        !
        sub_table(15,1,1,1)=1
        sub_table(15,1,1,2)=1
        sub_table(15,1,2,1)=1
        sub_table(15,1,2,2)=2
        sub_table(15,1,3,1)=1
        sub_table(15,1,3,2)=4
        sub_table(15,1,4,1)=1
        sub_table(15,1,4,2)=3
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
        sub_table(15,1,9,2)=3
        sub_table(15,1,10,1)=1
        sub_table(15,1,10,2)=4
        sub_table(15,1,11,1)=2
        sub_table(15,1,11,2)=5
        sub_table(15,1,11,3)=5
        sub_table(15,1,12,1)=2
        sub_table(15,1,12,2)=6
        sub_table(15,1,12,3)=6
        !
        !  with C_2h
        !
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
        sub_table(16,3,9,2)=4
        sub_table(16,3,10,1)=1
        sub_table(16,3,10,2)=3
        sub_table(16,3,11,1)=2
        sub_table(16,3,11,2)=3
        sub_table(16,3,11,3)=4
        sub_table(16,3,12,1)=2
        sub_table(16,3,12,2)=3
        sub_table(16,3,12,3)=4
        !
        !  with C_3h
        !
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
        !
        !  with C_6h
        !
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
        !
        !  with D_2h
        !
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
        !
        !  with D_3h
        !
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
        !
        !  with D_3d
        !
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

        sub_table(27,1,1,1)=1
        sub_table(27,1,1,2)=1
        sub_table(27,1,2,1)=1
        sub_table(27,1,2,2)=1
        sub_table(27,1,3,1)=1
        sub_table(27,1,3,2)=1
        sub_table(27,1,4,1)=1
        sub_table(27,1,4,2)=1
        sub_table(27,1,5,1)=2
        sub_table(27,1,5,2)=2
        sub_table(27,1,5,3)=3
        sub_table(27,1,6,1)=2
        sub_table(27,1,6,2)=2
        sub_table(27,1,6,3)=3
        sub_table(27,1,7,1)=1
        sub_table(27,1,7,2)=4
        sub_table(27,1,8,1)=1
        sub_table(27,1,8,2)=4
        sub_table(27,1,9,1)=1
        sub_table(27,1,9,2)=4
        sub_table(27,1,10,1)=1
        sub_table(27,1,10,2)=4
        sub_table(27,1,11,1)=2
        sub_table(27,1,11,2)=5
        sub_table(27,1,11,3)=6
        sub_table(27,1,12,1)=2
        sub_table(27,1,12,2)=5
        sub_table(27,1,12,3)=6

     CASE(24)
!
! D_2d
!
        !
        !  with C_s
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
        !
        !  with C_2
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
        !
        !  with D_2
        !
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
        sub_table(8,2,2,2)=3
        sub_table(8,2,3,1)=1
        sub_table(8,2,3,2)=1
        sub_table(8,2,4,1)=1
        sub_table(8,2,4,2)=3
        sub_table(8,2,5,1)=2
        sub_table(8,2,5,2)=1
        sub_table(8,2,5,3)=4

        sub_table(8,3,1,1)=1
        sub_table(8,3,1,2)=1
        sub_table(8,3,2,1)=1
        sub_table(8,3,2,2)=4
        sub_table(8,3,3,1)=1
        sub_table(8,3,3,2)=1
        sub_table(8,3,4,1)=1
        sub_table(8,3,4,2)=4
        sub_table(8,3,5,1)=2
        sub_table(8,3,5,2)=2
        sub_table(8,3,5,3)=3
        !
        !  with C_2v
        !
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
        !
        !  with S_4
        !
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
        !
        !  with C_i
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
        !
        !  with C_s
        !
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
        !
        !  with C_2
        !
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
        !
        !  with C_3
        !
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
        !
        !  with D_3
        !
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
        !
        !  with C_3v
        !
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
        !
        !  with C_2h
        !
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
        !
        !  with S_6
        !
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
        !
        !  with C_2
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
        !
        !  with C_i
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
        !
        !  with C_3
        !
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
        !
        !  with C_2
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
        !
        !  with C_3
        !
        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=2
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=3
        sub_table(5,1,4,1)=3
        sub_table(5,1,4,2)=1
        sub_table(5,1,4,3)=2
        sub_table(5,1,4,4)=3
        !
        !  with D_2
        !
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
        !
        !  with C_i
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
        !
        !  with C_s
        !
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
        !
        !  with C_2
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
        !
        !  with C_3
        !
        sub_table(5,1,1,1)=1
        sub_table(5,1,1,2)=1
        sub_table(5,1,2,1)=1
        sub_table(5,1,2,2)=2
        sub_table(5,1,3,1)=1
        sub_table(5,1,3,2)=3
        sub_table(5,1,4,1)=3
        sub_table(5,1,4,2)=1
        sub_table(5,1,4,3)=2
        sub_table(5,1,4,4)=3
        sub_table(5,1,5,1)=1
        sub_table(5,1,5,2)=1
        sub_table(5,1,6,1)=1
        sub_table(5,1,6,2)=2
        sub_table(5,1,7,1)=1
        sub_table(5,1,7,2)=3
        sub_table(5,1,8,1)=3
        sub_table(5,1,8,2)=1
        sub_table(5,1,8,3)=2
        sub_table(5,1,8,4)=3
        !
        !  with D_2
        !
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
        !
        !  with C_2v
        !
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
        !
        !  with C_2h
        !
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
        !
        !  with D_2h
        !
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
        !
        !  with S_6
        !
        sub_table(27,1,1,1)=1
        sub_table(27,1,1,2)=1
        sub_table(27,1,2,1)=1
        sub_table(27,1,2,2)=2
        sub_table(27,1,3,1)=1
        sub_table(27,1,3,2)=3
        sub_table(27,1,4,1)=3
        sub_table(27,1,4,2)=1
        sub_table(27,1,4,3)=2
        sub_table(27,1,4,4)=3
        sub_table(27,1,5,1)=1
        sub_table(27,1,5,2)=4
        sub_table(27,1,6,1)=1
        sub_table(27,1,6,2)=5
        sub_table(27,1,7,1)=1
        sub_table(27,1,7,2)=6
        sub_table(27,1,8,1)=3
        sub_table(27,1,8,2)=4
        sub_table(27,1,8,3)=5
        sub_table(27,1,8,4)=6
        !
        !  with T
        !
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
        !
        !  with C_s
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
        sub_table(3,1,4,3)=2
        sub_table(3,1,4,4)=2
        sub_table(3,1,5,1)=3
        sub_table(3,1,5,2)=1
        sub_table(3,1,5,3)=1
        sub_table(3,1,5,4)=2
        !
        !  with C_2
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
        !
        !  with C_3
        !
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
        !
        !  with D_2
        !
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
        !
        !  with C_2v
        !
        sub_table(12,1,1,1)=1
        sub_table(12,1,1,2)=1
        sub_table(12,1,2,1)=1
        sub_table(12,1,2,2)=2
        sub_table(12,1,3,1)=2
        sub_table(12,1,3,2)=1
        sub_table(12,1,3,3)=2
        sub_table(12,1,4,1)=3
        sub_table(12,1,4,2)=2
        sub_table(12,1,4,3)=3
        sub_table(12,1,4,4)=4
        sub_table(12,1,5,1)=3
        sub_table(12,1,5,2)=1
        sub_table(12,1,5,3)=3
        sub_table(12,1,5,4)=4
        !
        !  with C_3v
        !
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
        !
        !  with D_2d
        !
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
        !
        !  with S_4
        !
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
        !
        !  with T
        !
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
        !
        !  with C_2
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
        !
        !  with C_3
        !
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
        !
        !  with C_4
        !
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
        !
        !  with D_2
        !
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
        !
        !  with D_3
        !
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
        !
        !  with D_4
        !
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
        !
        !  with T
        !
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
        !
        !  with C_i
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
        !
        !  with C_s
        !
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

        sub_table(3,2,1,1)=1
        sub_table(3,2,1,2)=1
        sub_table(3,2,2,1)=1
        sub_table(3,2,2,2)=2
        sub_table(3,2,3,1)=2
        sub_table(3,2,3,2)=1
        sub_table(3,2,3,3)=2
        sub_table(3,2,4,1)=3
        sub_table(3,2,4,2)=1
        sub_table(3,2,4,3)=2
        sub_table(3,2,4,4)=2
        sub_table(3,2,5,1)=3
        sub_table(3,2,5,2)=1
        sub_table(3,2,5,3)=1
        sub_table(3,2,5,4)=2
        sub_table(3,2,6,1)=1
        sub_table(3,2,6,2)=2
        sub_table(3,2,7,1)=1
        sub_table(3,2,7,2)=1
        sub_table(3,2,8,1)=2
        sub_table(3,2,8,2)=1
        sub_table(3,2,8,3)=2
        sub_table(3,2,9,1)=3
        sub_table(3,2,9,2)=1
        sub_table(3,2,9,3)=1
        sub_table(3,2,9,4)=2
        sub_table(3,2,10,1)=3
        sub_table(3,2,10,2)=1
        sub_table(3,2,10,3)=2
        sub_table(3,2,10,4)=2
        !
        !  with C_2
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
        sub_table(4,2,6,1)=1
        sub_table(4,2,6,2)=1
        sub_table(4,2,7,1)=1
        sub_table(4,2,7,2)=2
        sub_table(4,2,8,1)=2
        sub_table(4,2,8,2)=1
        sub_table(4,2,8,3)=2
        sub_table(4,2,9,1)=3
        sub_table(4,2,9,2)=1
        sub_table(4,2,9,3)=2
        sub_table(4,2,9,4)=2
        sub_table(4,2,10,1)=3
        sub_table(4,2,10,2)=1
        sub_table(4,2,10,3)=1
        sub_table(4,2,10,4)=2
        !
        !  with C_3
        !
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
        !
        !  with C_4
        !
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
        !
        !  with D_2
        !
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
        sub_table(8,2,6,1)=1
        sub_table(8,2,6,2)=1
        sub_table(8,2,7,1)=1
        sub_table(8,2,7,2)=2
        sub_table(8,2,8,1)=2
        sub_table(8,2,8,2)=1
        sub_table(8,2,8,3)=2
        sub_table(8,2,9,1)=3
        sub_table(8,2,9,2)=2
        sub_table(8,2,9,3)=3
        sub_table(8,2,9,4)=4
        sub_table(8,2,10,1)=3
        sub_table(8,2,10,2)=1
        sub_table(8,2,10,3)=3
        sub_table(8,2,10,4)=4
        !
        !  with D_3
        !
        sub_table(9,1,1,1)=1
        sub_table(9,1,1,2)=1
        sub_table(9,1,2,1)=1
        sub_table(9,1,2,2)=2
        sub_table(9,1,3,1)=2
        sub_table(9,1,3,2)=3
        sub_table(9,1,3,3)=3
        sub_table(9,1,4,1)=3
        sub_table(9,1,4,2)=1
        sub_table(9,1,4,3)=3
        sub_table(9,1,4,4)=3
        sub_table(9,1,5,1)=3
        sub_table(9,1,5,2)=2
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
        sub_table(9,1,9,2)=1
        sub_table(9,1,9,3)=3
        sub_table(9,1,9,4)=3
        sub_table(9,1,10,1)=3
        sub_table(9,1,10,2)=2
        sub_table(9,1,10,3)=3
        sub_table(9,1,10,4)=3
        !
        !  with D_4
        !
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
        !
        !  with C_2v
        !
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
        sub_table(12,2,5,2)=1
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
        sub_table(12,2,10,2)=2
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
        !
        !  with C_3v
        !
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
        !
        !  with C_4v
        !
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
        !
        !  with C_2h
        !
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
        !
        !  with C_4h
        !
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
        !
        !  with D_2h
        !
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
        sub_table(20,2,4,2)=2
        sub_table(20,2,4,3)=3
        sub_table(20,2,4,4)=4
        sub_table(20,2,5,1)=3
        sub_table(20,2,5,2)=1
        sub_table(20,2,5,3)=3
        sub_table(20,2,5,4)=4
        sub_table(20,2,6,1)=1
        sub_table(20,2,6,2)=5
        sub_table(20,2,7,1)=1
        sub_table(20,2,7,2)=6
        sub_table(20,2,8,1)=2
        sub_table(20,2,8,2)=5
        sub_table(20,2,8,3)=6
        sub_table(20,2,9,1)=3
        sub_table(20,2,9,2)=6
        sub_table(20,2,9,3)=7
        sub_table(20,2,9,4)=8
        sub_table(20,2,10,1)=3
        sub_table(20,2,10,2)=5
        sub_table(20,2,10,3)=7
        sub_table(20,2,10,4)=8
        !
        !  with D_4h
        !
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
        !
        !  with D_2d
        !
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

        sub_table(24,2,1,1)=1
        sub_table(24,2,1,2)=1
        sub_table(24,2,2,1)=1
        sub_table(24,2,2,2)=4
        sub_table(24,2,3,1)=2
        sub_table(24,2,3,2)=1
        sub_table(24,2,3,3)=4
        sub_table(24,2,4,1)=3
        sub_table(24,2,4,2)=2
        sub_table(24,2,4,3)=5
        sub_table(24,2,4,4)=5
        sub_table(24,2,5,1)=3
        sub_table(24,2,5,2)=3
        sub_table(24,2,5,3)=5
        sub_table(24,2,5,4)=5
        sub_table(24,2,6,1)=1
        sub_table(24,2,6,2)=3
        sub_table(24,2,7,1)=1
        sub_table(24,2,7,2)=2
        sub_table(24,2,8,1)=2
        sub_table(24,2,8,2)=2
        sub_table(24,2,8,3)=3
        sub_table(24,2,9,1)=3
        sub_table(24,2,9,2)=4
        sub_table(24,2,9,3)=5
        sub_table(24,2,9,4)=5
        sub_table(24,2,10,1)=3
        sub_table(24,2,10,2)=1
        sub_table(24,2,10,3)=5
        sub_table(24,2,10,4)=5
        !
        !  with D_3d
        !
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
        !
        !  with S_4
        !
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
        !
        !  with S_6
        !
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
        !
        !  with T
        !
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
        !
        !  with T_h
        !
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
        !
        !  with T_d
        !
        sub_table(30,1,1,1)=1
        sub_table(30,1,1,2)=1
        sub_table(30,1,2,1)=1
        sub_table(30,1,2,2)=2
        sub_table(30,1,3,1)=2
        sub_table(30,1,3,2)=3
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
        !
        !  with O
        !
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
!  The first part of the routine sets the information for all
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
        !
        !  with C_2
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
        !
        !  with C_2
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
        sub_table(4,1,5,2)=2
        sub_table(4,1,6,1)=1
        sub_table(4,1,6,2)=1

        !
        !  with C_3
        !
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
        !
        !  with C_2
        !
        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2

     CASE(9)
!
! D_3
!
        !
        !  with C_2
        !
        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=1
        sub_table(4,1,2,2)=1
        sub_table(4,1,3,1)=1
        sub_table(4,1,3,2)=2
        !
        !  with C_3
        !
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
        !
        !  with C_2
        !
        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2
        !
        !  with C_4
        !
        sub_table(6,1,1,1)=2
        sub_table(6,1,1,2)=1
        sub_table(6,1,1,3)=2
        sub_table(6,1,2,1)=2
        sub_table(6,1,2,2)=3
        sub_table(6,1,2,3)=4
        !
        !  with D_2
        !
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
        !
        !  with C_2
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
        !
        !  with C_3
        !
        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=2
        sub_table(5,1,2,2)=1
        sub_table(5,1,2,3)=2
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=3
        sub_table(5,1,3,3)=3
        !
        !  with C_6
        !
        sub_table(7,1,1,1)=2
        sub_table(7,1,1,2)=1
        sub_table(7,1,1,3)=2
        sub_table(7,1,2,1)=2
        sub_table(7,1,2,2)=3
        sub_table(7,1,2,3)=4
        sub_table(7,1,3,1)=2
        sub_table(7,1,3,2)=5
        sub_table(7,1,3,3)=6
        !
        !  with D_2
        !
        sub_table(8,1,1,1)=2
        sub_table(8,1,1,2)=1
        sub_table(8,1,1,3)=1
        sub_table(8,1,2,1)=2
        sub_table(8,1,2,2)=1
        sub_table(8,1,2,3)=1
        sub_table(8,1,3,1)=2
        sub_table(8,1,3,2)=1
        sub_table(8,1,3,3)=1
        !
        !  with D_3
        !
        sub_table(9,1,1,1)=2
        sub_table(9,1,1,2)=1
        sub_table(9,1,1,3)=1
        sub_table(9,1,2,1)=2
        sub_table(9,1,2,2)=1
        sub_table(9,1,2,3)=1
        sub_table(9,1,3,1)=2
        sub_table(9,1,3,2)=2
        sub_table(9,1,3,3)=3

     CASE(12)
!
!   C_2v
!
        !
        !  with C_s
        !
        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        !
        !  with C_2
        !
        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        
     CASE(13)
!
!   C_3v
!
        !
        !  with C_s
        !
        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=1
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=2
        !
        !  with C_3
        !
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
        !
        !  with C_3
        !
        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=2
        sub_table(3,1,2,2)=1
        sub_table(3,1,2,3)=2
        !
        !  with C_2
        !
        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2
        !
        !  with C_4
        !
        sub_table(6,1,1,1)=2
        sub_table(6,1,1,2)=1
        sub_table(6,1,1,3)=2
        sub_table(6,1,2,1)=2
        sub_table(6,1,2,2)=3
        sub_table(6,1,2,3)=4
        !
        !  with C_2v
        !
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
        !
        !  with C_s
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
        !
        !  with C_2
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
        !
        !  with C_3
        !
        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=2
        sub_table(5,1,2,2)=1
        sub_table(5,1,2,3)=2
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=3
        sub_table(5,1,3,3)=3
        !
        !  with C_6
        !
        sub_table(7,1,1,1)=2
        sub_table(7,1,1,2)=1
        sub_table(7,1,1,3)=2
        sub_table(7,1,2,1)=2
        sub_table(7,1,2,2)=3
        sub_table(7,1,2,3)=4
        sub_table(7,1,3,1)=2
        sub_table(7,1,3,2)=5
        sub_table(7,1,3,3)=6
        !
        !  with C_2v
        !
        sub_table(12,1,1,1)=2
        sub_table(12,1,1,2)=1
        sub_table(12,1,1,3)=1
        sub_table(12,1,2,1)=2
        sub_table(12,1,2,2)=1
        sub_table(12,1,2,3)=1
        sub_table(12,1,3,1)=2
        sub_table(12,1,3,2)=1
        sub_table(12,1,3,3)=1
        !
        !  with C_3v
        !
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
        !
        !  with C_i
        !
        sub_table(2,1,1,1)=1
        sub_table(2,1,1,2)=1
        sub_table(2,1,2,1)=1
        sub_table(2,1,2,2)=1
        sub_table(2,1,3,1)=1
        sub_table(2,1,3,2)=2
        sub_table(2,1,4,1)=1
        sub_table(2,1,4,2)=2
        !
        !  with C_s
        !
        sub_table(3,1,1,1)=1
        sub_table(3,1,1,2)=1
        sub_table(3,1,2,1)=1
        sub_table(3,1,2,2)=2
        sub_table(3,1,3,1)=1
        sub_table(3,1,3,2)=2
        sub_table(3,1,4,1)=1
        sub_table(3,1,4,2)=1
        !
        !  with C_2
        !
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
        !
        !  with C_s
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
        sub_table(3,1,5,2)=2
        sub_table(3,1,6,1)=1
        sub_table(3,1,6,2)=1
        !
        !  with C_3
        !
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
        !
        !  with C_i
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
        !
        !  with C_s
        !
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
        !
        !  with C_2
        !
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
        !
        !  with C_4
        !
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
        !
        !  with C_2h
        !
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
        sub_table(16,1,6,2)=4
        sub_table(16,1,7,1)=1
        sub_table(16,1,7,2)=3
        sub_table(16,1,8,1)=1
        sub_table(16,1,8,2)=4
        !
        !  with S_4
        !
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
        !
        !  with C_i
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
        !
        !  with C_s
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
        sub_table(3,1,5,2)=2
        sub_table(3,1,6,1)=1
        sub_table(3,1,6,2)=1
        sub_table(3,1,7,1)=1
        sub_table(3,1,7,2)=2
        sub_table(3,1,8,1)=1
        sub_table(3,1,8,2)=1
        sub_table(3,1,9,1)=1
        sub_table(3,1,9,2)=1
        sub_table(3,1,10,1)=1
        sub_table(3,1,10,2)=2
        sub_table(3,1,11,1)=1
        sub_table(3,1,11,2)=1
        sub_table(3,1,12,1)=1
        sub_table(3,1,12,2)=2
        !
        !  with C_2
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
        sub_table(4,1,5,2)=2
        sub_table(4,1,6,1)=1
        sub_table(4,1,6,2)=1
        sub_table(4,1,7,1)=1
        sub_table(4,1,7,2)=1
        sub_table(4,1,8,1)=1
        sub_table(4,1,8,2)=2
        sub_table(4,1,9,1)=1
        sub_table(4,1,9,2)=2
        sub_table(4,1,10,1)=1
        sub_table(4,1,10,2)=1
        sub_table(4,1,11,1)=1
        sub_table(4,1,11,2)=2
        sub_table(4,1,12,1)=1
        sub_table(4,1,12,2)=1
        !
        !  with C_3
        !
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
        !
        !  with C_6
        !
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
        !
        !  with C_2h
        !
        sub_table(16,1,1,1)=1
        sub_table(16,1,1,2)=1
        sub_table(16,1,2,1)=1
        sub_table(16,1,2,2)=2
        sub_table(16,1,3,1)=1
        sub_table(16,1,3,2)=2
        sub_table(16,1,4,1)=1
        sub_table(16,1,4,2)=1
        sub_table(16,1,5,1)=1
        sub_table(16,1,5,2)=2
        sub_table(16,1,6,1)=1
        sub_table(16,1,6,2)=1
        sub_table(16,1,7,1)=1
        sub_table(16,1,7,2)=3
        sub_table(16,1,8,1)=1
        sub_table(16,1,8,2)=4
        sub_table(16,1,9,1)=1
        sub_table(16,1,9,2)=4
        sub_table(16,1,10,1)=1
        sub_table(16,1,10,2)=3
        sub_table(16,1,11,1)=1
        sub_table(16,1,11,2)=4
        sub_table(16,1,12,1)=1
        sub_table(16,1,12,2)=3
        !
        !  with C_3h
        !
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
        !
        !  with S_6
        !
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
        !
        !  with C_i
        !
        sub_table(2,1,1,1)=2
        sub_table(2,1,1,2)=1
        sub_table(2,1,1,3)=1
        sub_table(2,1,2,1)=2
        sub_table(2,1,2,2)=2
        sub_table(2,1,2,3)=2
        !
        !  with C_s
        !
        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=2
        sub_table(3,1,2,2)=1
        sub_table(3,1,2,3)=2
        !
        !  with C_2
        !
        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2
        !
        !  with D_2
        !
        sub_table(8,1,1,1)=2
        sub_table(8,1,1,2)=1
        sub_table(8,1,1,3)=1
        sub_table(8,1,2,1)=2
        sub_table(8,1,2,2)=1
        sub_table(8,1,2,3)=1
        !
        !  with C_2v
        !
        sub_table(12,1,1,1)=2
        sub_table(12,1,1,2)=1
        sub_table(12,1,1,3)=1
        sub_table(12,1,2,1)=2
        sub_table(12,1,2,2)=1
        sub_table(12,1,2,3)=1
        !
        !  with C_2h
        !
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
        !
        !  with C_s
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
        !
        !  with C_2
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
        !
        !  with C_3
        !
        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=2
        sub_table(5,1,2,2)=1
        sub_table(5,1,2,3)=2
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=3
        sub_table(5,1,3,3)=3
        !
        !  with D_3
        !
        sub_table(9,1,1,1)=2
        sub_table(9,1,1,2)=1
        sub_table(9,1,1,3)=1
        sub_table(9,1,2,1)=2
        sub_table(9,1,2,2)=1
        sub_table(9,1,2,3)=1
        sub_table(9,1,3,1)=2
        sub_table(9,1,3,2)=2
        sub_table(9,1,3,3)=3
        !
        !  with C_2v
        !
        sub_table(12,1,1,1)=2
        sub_table(12,1,1,2)=1
        sub_table(12,1,1,3)=1
        sub_table(12,1,2,1)=2
        sub_table(12,1,2,2)=1
        sub_table(12,1,2,3)=1
        sub_table(12,1,3,1)=2
        sub_table(12,1,3,2)=1
        sub_table(12,1,3,3)=1
        !
        !  with C_3v
        !
        sub_table(13,1,1,1)=2
        sub_table(13,1,1,2)=1
        sub_table(13,1,1,3)=1
        sub_table(13,1,2,1)=2
        sub_table(13,1,2,2)=1
        sub_table(13,1,2,3)=1
        sub_table(13,1,3,1)=2
        sub_table(13,1,3,2)=2
        sub_table(13,1,3,3)=3
        !
        !  with C_3h
        !
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
        !
        !  with C_i
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
        !
        !  with C_s
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
        sub_table(3,1,4,1)=2
        sub_table(3,1,4,2)=1
        sub_table(3,1,4,3)=2
        !
        !  with C_2
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
        sub_table(4,1,4,1)=2
        sub_table(4,1,4,2)=1
        sub_table(4,1,4,3)=2
        !
        !  with C_4
        !
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
        !
        !  with D_2
        !
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
        !
        !  with D_4
        !
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
        !
        !  with C_2v
        !
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
        !
        !  with C_4v
        !
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
        !
        !  with C_2h
        !
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
        !
        !  with C_4h
        !
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
        !
        !  with D_2h
        !
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
        !
        !  with D_2d
        !
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
        !
        !  with S_4
        !
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
        !
        !  with C_i
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
        !
        !  with C_s
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
        sub_table(3,1,4,1)=2
        sub_table(3,1,4,2)=1
        sub_table(3,1,4,3)=2
        sub_table(3,1,5,1)=2
        sub_table(3,1,5,2)=1
        sub_table(3,1,5,3)=2
        sub_table(3,1,6,1)=2
        sub_table(3,1,6,2)=1
        sub_table(3,1,6,3)=2
        !
        !  with C_2
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
        sub_table(4,1,4,1)=2
        sub_table(4,1,4,2)=1
        sub_table(4,1,4,3)=2
        sub_table(4,1,5,1)=2
        sub_table(4,1,5,2)=1
        sub_table(4,1,5,3)=2
        sub_table(4,1,6,1)=2
        sub_table(4,1,6,2)=1
        sub_table(4,1,6,3)=2
        !
        !  with C_3
        !
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
        !
        !  with C_6
        !
        sub_table(7,1,1,1)=2
        sub_table(7,1,1,2)=1
        sub_table(7,1,1,3)=2
        sub_table(7,1,2,1)=2
        sub_table(7,1,2,2)=3
        sub_table(7,1,2,3)=4
        sub_table(7,1,3,1)=2
        sub_table(7,1,3,2)=5
        sub_table(7,1,3,3)=6
        sub_table(7,1,4,1)=2
        sub_table(7,1,4,2)=1
        sub_table(7,1,4,3)=2
        sub_table(7,1,5,1)=2
        sub_table(7,1,5,2)=3
        sub_table(7,1,5,3)=4
        sub_table(7,1,6,1)=2
        sub_table(7,1,6,2)=5
        sub_table(7,1,6,3)=6
        !
        !  with D_2
        !
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
        !
        !  with D_3
        !
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
        !
        !  with D_6
        !
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
        !
        !  with C_2v
        !
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
        !
        !  with C_3v
        !
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
        !
        !  with C_6v
        !
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
        !
        !  with C_2h
        !
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
        !
        !  with C_3h
        !
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
        !
        !  with C_6h
        !
        sub_table(19,1,1,1)=2
        sub_table(19,1,1,2)=1
        sub_table(19,1,1,3)=2
        sub_table(19,1,2,1)=2
        sub_table(19,1,2,2)=3
        sub_table(19,1,2,3)=4
        sub_table(19,1,3,1)=2
        sub_table(19,1,3,2)=5
        sub_table(19,1,3,3)=6
        sub_table(19,1,4,1)=2
        sub_table(19,1,4,2)=7
        sub_table(19,1,4,3)=8
        sub_table(19,1,5,1)=2
        sub_table(19,1,5,2)=9
        sub_table(19,1,5,3)=10
        sub_table(19,1,6,1)=2
        sub_table(19,1,6,2)=11
        sub_table(19,1,6,3)=12
        !
        !  with D_2h
        !
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
        !
        !  with D_3h
        !
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
        sub_table(21,1,4,2)=2
        sub_table(21,1,4,3)=2
        sub_table(21,1,5,1)=2
        sub_table(21,1,5,2)=1
        sub_table(21,1,5,3)=1
        sub_table(21,1,6,1)=2
        sub_table(21,1,6,2)=3
        sub_table(21,1,6,3)=3
        !
        !  with D_3d
        !
        sub_table(25,1,1,1)=2
        sub_table(25,1,1,2)=1
        sub_table(25,1,1,3)=1
        sub_table(25,1,2,1)=2
        sub_table(25,1,2,2)=1
        sub_table(25,1,2,3)=1
        sub_table(25,1,3,1)=2
        sub_table(25,1,3,2)=2
        sub_table(25,1,3,3)=3
        sub_table(25,1,4,1)=2
        sub_table(25,1,4,2)=4
        sub_table(25,1,4,3)=4
        sub_table(25,1,5,1)=2
        sub_table(25,1,5,2)=4
        sub_table(25,1,5,3)=4
        sub_table(25,1,6,1)=2
        sub_table(25,1,6,2)=5
        sub_table(25,1,6,3)=6
        !
        !  with S_6
        !
        sub_table(27,1,1,1)=2
        sub_table(27,1,1,2)=1
        sub_table(27,1,1,3)=2
        sub_table(27,1,2,1)=2
        sub_table(27,1,2,2)=1
        sub_table(27,1,2,3)=2
        sub_table(27,1,3,1)=2
        sub_table(27,1,3,2)=3
        sub_table(27,1,3,3)=3
        sub_table(27,1,4,1)=2
        sub_table(27,1,4,2)=4
        sub_table(27,1,4,3)=5
        sub_table(27,1,5,1)=2
        sub_table(27,1,5,2)=4
        sub_table(27,1,5,3)=5
        sub_table(27,1,6,1)=2
        sub_table(27,1,6,2)=6
        sub_table(27,1,6,3)=6

     CASE(24)
!
! D_2d
!
        !
        !  with C_s
        !
        sub_table(3,1,1,1)=2
        sub_table(3,1,1,2)=1
        sub_table(3,1,1,3)=2
        sub_table(3,1,2,1)=2
        sub_table(3,1,2,2)=1
        sub_table(3,1,2,3)=2
        !
        !  with C_2
        !
        sub_table(4,1,1,1)=2
        sub_table(4,1,1,2)=1
        sub_table(4,1,1,3)=2
        sub_table(4,1,2,1)=2
        sub_table(4,1,2,2)=1
        sub_table(4,1,2,3)=2
        !
        !  with D_2
        !
        sub_table(8,1,1,1)=2
        sub_table(8,1,1,2)=1
        sub_table(8,1,1,3)=1
        sub_table(8,1,2,1)=2
        sub_table(8,1,2,2)=1
        sub_table(8,1,2,3)=1
        !
        !  with C_2v
        !
        sub_table(12,1,1,1)=2
        sub_table(12,1,1,2)=1
        sub_table(12,1,1,3)=1
        sub_table(12,1,2,1)=2
        sub_table(12,1,2,2)=1
        sub_table(12,1,2,3)=1
        !
        !  with S_4
        !
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
        !
        !  with C_i
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
        !
        !  with C_s
        !
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
        !
        !  with C_2
        !
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
        !
        !  with C_3
        !
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
        !
        !  with D_3
        !
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
        !
        !  with C_3v
        !
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
        !
        !  with C_2h
        !
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
        !
        !  with S_6
        !
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
        !
        !  with C_2
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
        !
        !  with C_i
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
        !
        !  with C_3
        !
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
        !
        !  with C_2
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
        !
        !  with C_3
        !
        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=2
        sub_table(5,1,2,2)=1
        sub_table(5,1,2,3)=3
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=2
        sub_table(5,1,3,3)=3
        !
        !  with D_2
        !
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
        !
        !  with C_i
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
        !
        !  with C_s
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
        sub_table(3,1,4,1)=2
        sub_table(3,1,4,2)=1
        sub_table(3,1,4,3)=2
        sub_table(3,1,5,1)=2
        sub_table(3,1,5,2)=1
        sub_table(3,1,5,3)=2
        sub_table(3,1,6,1)=2
        sub_table(3,1,6,2)=1
        sub_table(3,1,6,3)=2
        !
        !  with C_2
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
        sub_table(4,1,4,1)=2
        sub_table(4,1,4,2)=1
        sub_table(4,1,4,3)=2
        sub_table(4,1,5,1)=2
        sub_table(4,1,5,2)=1
        sub_table(4,1,5,3)=2
        sub_table(4,1,6,1)=2
        sub_table(4,1,6,2)=1
        sub_table(4,1,6,3)=2
        !
        !  with C_3
        !
        sub_table(5,1,1,1)=2
        sub_table(5,1,1,2)=1
        sub_table(5,1,1,3)=2
        sub_table(5,1,2,1)=2
        sub_table(5,1,2,2)=1
        sub_table(5,1,2,3)=3
        sub_table(5,1,3,1)=2
        sub_table(5,1,3,2)=2
        sub_table(5,1,3,3)=3
        sub_table(5,1,4,1)=2
        sub_table(5,1,4,2)=1
        sub_table(5,1,4,3)=2
        sub_table(5,1,5,1)=2
        sub_table(5,1,5,2)=1
        sub_table(5,1,5,3)=3
        sub_table(5,1,6,1)=2
        sub_table(5,1,6,2)=2
        sub_table(5,1,6,3)=3
        !
        !  with D_2
        !
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
        !
        !  with C_2v
        !
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
        !
        !  with C_2h
        !
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
        !
        !  with D_2h
        !
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
        !
        !  with S_6
        !
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
        !
        !  with T
        !
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
        !
        !  with C_s
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
        !
        !  with C_2
        !
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
        !
        !  with C_3
        !
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
        !
        !  with D_2
        !
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
        !
        !  with C_2v
        !
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
        !
        !  with C_3v
        !
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
        !
        !  with D_2d
        !
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
        !
        !  with S_4
        !
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
        !
        !  with T
        !
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
        !
        !  with C_2
        !
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
        !
        !  with C_3
        !
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
        !
        !  with C_4
        !
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
        !
        !  with D_2
        !
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
        !
        !  with D_3
        !
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
        !
        !  with D_4
        !
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
        !
        !  with T
        !
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
        !
        !  with C_i
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
        !
        !  with C_s
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
        !
        !  with C_2
        !
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
        !
        !  with C_3
        !
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
        sub_table(5,1,6,5)=3
        !
        !  with C_4
        !
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
        !
        !  with D_2
        !
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
        !
        !  with D_3
        !
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
        !
        !  with D_4
        !
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
        !
        !  with C_2v
        !
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
        !
        !  with C_3v
        !
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
        !
        !  with C_4v
        !
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
        !
        !  with C_2h
        !
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
        !
        !  with C_4h
        !
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
        !
        !  with D_2h
        !
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
        !
        !  with D_4h
        !
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
        !
        !  with D_2d
        !
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
        !
        !  with D_3d
        !
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
        !
        !  with S_4
        !
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
        !
        !  with S_6
        !
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
        !
        !  with T
        !
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
        sub_table(28,1,5,2)=1
        sub_table(28,1,5,3)=1
        sub_table(28,1,6,1)=4
        sub_table(28,1,6,2)=2
        sub_table(28,1,6,3)=2
        sub_table(28,1,6,4)=3
        sub_table(28,1,6,5)=3
        !
        !  with T_h
        !
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
        !
        !  with T_d
        !
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
        !
        !  with O
        !
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
!  10 D_4 : (9)
!  C_1, C_2_1, C_2_2, C_2_3, C_4, D_2_1, D_2_2, D_2_3, D_2_4
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

  REAL(DP) :: sr_a(3,3,nsym_a), sr_b(3,3,nsym_b), ax(3), bx(3), cx(3), saxis(3,3)
  LOGICAL :: equal
  INTEGER :: isym, ipol, jpol, imirror, iaxis, four_axis, iax, ibx, icx, id2
  INTEGER :: imax, imbx, imcx, ic4, is4
  INTEGER :: ic2, ic21, ic211, isv, isv1, ts, ind2(3)
  INTEGER :: xaxis, yaxis, zaxis, naxis, isave
  LOGICAL :: is_parallel, isok, isok1
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
!
!   first determine the type of D_2: 1 axes parallel to x,y,z; 2 some axis parallel
!   to some other C_2 axis
!
        id2=1
        DO isym=1,nsym_a
           IF (id2==1.AND. tipo_sym(sr_a(1,1,isym))==4) THEN
              CALL versor(sr_a(1,1,isym),ax)
              CALL which_c2(ax, iax)
              IF (iax > 3) id2=2
           ENDIF
        ENDDO 

        SELECT CASE (group_b)
           CASE(2)
      !
      !    C_2
      !
              DO isym=1,nsym_b
!
!   find the axis of order 2 of C_2 and check to which axis it is parallel
!
                IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                   CALL versor(sr_b(1,1,isym),ax)
                   CALL which_c2(ax, iax)
                   SELECT CASE (iax)
!
!   See point group manual for the different possibilities.
!
                      CASE (1)
                         aux_ind=1
                      CASE (2)
                         IF (id2==1) THEN
                           aux_ind=2
                         ELSE
                           aux_ind=1
                         ENDIF
                      CASE (3)
                         IF (id2==1) THEN
                           aux_ind=3
                         ELSE
                           aux_ind=1
                         ENDIF
                      CASE (4,6,8,10,11) 
                         aux_ind=2
                      CASE (5,7,9,12,13) 
                         aux_ind=3
                   END SELECT
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
!
!   find the axis C_4 of D_4
!
        DO isym=1,nsym_a
           IF (tipo_sym(sr_a(1,1,isym))==3) THEN
              CALL versor(sr_a(1,1,isym),ax)
              CALL which_c2(ax, iax)
              IF (iax > 3) CALL errore('find_aux_ind_two_groups','wrong C_4',1)
           ENDIF
        ENDDO
             
        SELECT CASE (group_b)
           CASE(4)
              DO isym=1,nsym_b
!
!   find the axis of order 2 and check if it is parallel to the C_4 axis of D_4
!   (case 1), to x,y, or z axes (case 2) or to another axis (case 3)
!
                 IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                    CALL versor(sr_b(1,1,isym),ax)
                    IF (is_axis(ax,iax)) THEN
                       aux_ind=1
                    ELSEIF (is_axis(ax,1).OR.is_axis(ax,2).OR.is_axis(ax,3)) THEN
                       aux_ind=2
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
!
!   first determine the type of D_2: 1 all axes parallel to x,y,z; 
!   2 some axis parallel to some other C_2 axis
!
              id2=1
              DO isym=1,nsym_b
                 IF (id2==1.AND. tipo_sym(sr_b(1,1,isym))==4) THEN
                   CALL versor(sr_b(1,1,isym),bx)
                   CALL which_c2(bx, ibx)
                   IF (ibx > 3) id2=2
                 ENDIF
              ENDDO
              IF (id2==1) THEN
!
!  In this case we have only to check which axis of D_2 coincides
!  with C_4 of D_4. If C_4 is x it is C_2'', if C_4 is y it is C_2', if C_4 is
!  z it is C_2.
!
                 IF (iax==1) THEN
                    aux_ind=4
                 ELSEIF (iax==2) THEN
                    aux_ind=3
                 ELSEIF (iax==3) THEN
                    aux_ind=1
                 ENDIF
              ELSE
!
!   In this case C_2 of D_2 coincides with C_4 of D_4 and C_2' and C_2'' of D_2
!   coincide with C_2'' of D_4
!
                 aux_ind=2
              ENDIF
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   

     CASE(11)
!
!  D_6
!
        SELECT CASE (group_b)
           CASE(4)
      !
      !    C_2
      !
              DO isym=1,nsym_b
!
!   find the axis of order 2 and check to which axis it is parallel
!
                 IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                    CALL versor(sr_b(1,1,isym),ax)
                    CALL which_c2(ax, iax)
                    IF (iax==3) THEN
                       aux_ind=1
                    ELSEIF (iax==1.OR.iax==10.OR.iax==11) THEN
                       aux_ind=2
                    ELSEIF (iax==2.OR.iax==12.OR.iax==13) THEN
                       aux_ind=3
                    ELSE
                       CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
                    ENDIF
                 ENDIF
              ENDDO
           CASE(5,7,8)
      !
      !    C_3, C_6, D_2
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
                    CALL which_c2(ax, iax)
                    IF (iax==1.OR.iax==10.OR.iax==11) THEN
                       aux_ind=1
                    ELSEIF (iax==2.OR.iax==12.OR.iax==13) THEN
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
              !  with C_s
              !
              !  In ibx and icx the two perpendiculars of the two mirrors of C_2v
              !
              ibx=0
              DO isym=1, nsym_a
                 IF (tipo_sym(sr_a(1,1,isym))==5) THEN
                    CALL mirror_axis(sr_a(1,1,isym),ax)
                    IF (ibx==0) THEN
                       CALL which_c2(ax, ibx)
                    ELSE
                       CALL which_c2(ax, icx)
                    ENDIF
                 ENDIF
              ENDDO
              !
              !  In iax the perpendicular to the C_s mirror
              !
              DO isym=1, nsym_b
                 IF (tipo_sym(sr_b(1,1,isym))==5) THEN
                    CALL mirror_axis(sr_b(1,1,isym),ax)
                    CALL which_c2(ax, iax)
                 ENDIF
              ENDDO
              imirror=0
              IF (ibx==iax) imirror=icx 
              IF (icx==iax) imirror=ibx 
              IF (imirror==0) CALL errore('find_aux_ind_two_groups',&
                                           'C_2v and C_s have no common mirror',1)

              IF (iax==1 .OR. iax==5 .OR. iax==7 .OR. iax==9 .OR. iax==10 .OR. &
                                                                  iax==11) THEN
                  aux_ind=2
              ELSEIF (iax==2) THEN
                  aux_ind=1
              ELSEIF (iax==3) THEN
                  IF (imirror==2) THEN
                     aux_ind=2
                  ELSE
                     aux_ind=1
                  ENDIF
              ELSEIF (iax==4) THEN
                  IF (imirror==1) THEN
                     aux_ind=2
                  ELSE
                     aux_ind=1
                  ENDIF
              ELSEIF (iax==6) THEN
                  IF (imirror==2) THEN
                     aux_ind=2
                  ELSE
                     aux_ind=1
                  ENDIF
              ELSEIF (iax==8.OR.iax==12.OR.iax==13) THEN
                  IF (imirror==3) THEN
                     aux_ind=2
                  ELSE
                     aux_ind=1
                  ENDIF
              ELSE
                 CALL errore('find_aux_ind_two_groups','C_2v/C_s problem with axis',1)
              ENDIF 
           CASE (4)
      !
      !   with C_2
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
              CALL which_c2(bx, ibx)
              IF ( ibx==2 .OR. ibx==12 .OR. ibx==13 ) THEN
                 aux_ind=1
              ELSEIF ( ibx==1 .OR. ibx==10 .OR. ibx==11 ) THEN
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
              CALL which_c2(bx, ibx)
              IF ( ibx==2 .OR. ibx==12 .OR. ibx==13 ) THEN
                 aux_ind=1
              ELSEIF ( ibx==1 .OR. ibx==10 .OR. ibx==11 ) THEN
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
!
!      First analyze D_2h and determine C_2, C_2', and C_2''
!
        iax=0
        ibx=0
        icx=0
        imax=0
        imbx=0
        imcx=0
        DO isym=2,nsym_a
           ts=tipo_sym(sr_a(1,1,isym))
           IF (ts==4) THEN
              CALL versor(sr_a(1,1,isym),ax)
              IF (iax==0) THEN
                 CALL which_c2(ax, iax)
              ELSEIF (ibx==0) THEN
                 CALL which_c2(ax, ibx)
              ELSEIF (icx==0) THEN
                 CALL which_c2(ax, icx)
              ELSE
                 CALL errore('find_aux_ind_two_groups','D_2h too many C_2 axis',1)
              ENDIF
           ELSEIF (ts==5) THEN
              CALL mirror_axis(sr_a(1,1,isym),ax)
              IF (imax==0) THEN
                 CALL which_c2(ax, imax)
              ELSEIF (imbx==0) THEN
                 CALL which_c2(ax, imbx)
              ELSEIF (imcx==0) THEN
                 CALL which_c2(ax, imcx)
              ELSE
                 CALL errore('find_aux_ind_two_groups','D_2h too many mirrors',1)
              ENDIF
           ELSEIF (ts /= 2) THEN
              CALL errore('find_aux_ind_two_groups','D_2h operation not recognized',1)
           ENDIF
        ENDDO
        CALL is_d2(iax, ibx, icx, ind2)

        IF (ind2(1)==1) ic2=iax
        IF (ind2(1)==2) ic2=ibx
        IF (ind2(1)==3) ic2=icx
        IF (ind2(2)==1) ic21=iax
        IF (ind2(2)==2) ic21=ibx
        IF (ind2(2)==3) ic21=icx
        IF (ind2(3)==1) ic211=iax
        IF (ind2(3)==2) ic211=ibx
        IF (ind2(3)==3) ic211=icx

       SELECT CASE (group_b)
          CASE(2)
             aux_ind=1
          CASE(3)
!
!   find the mirror normal
!
             CALL mirror_axis(sr_b(1,1,2),bx)
             CALL which_c2(bx, ibx)
             IF (ibx==ic2) THEN
                aux_ind=1
             ELSEIF (ibx==ic21) THEN
                aux_ind=2
             ELSEIF (ibx==ic211) THEN
                aux_ind=3
             ENDIF
          CASE(4)
!
!   find the C_2 axis normal
!
             CALL versor(sr_b(1,1,2),bx)
             CALL which_c2(bx, ibx)
             IF (ibx==ic2) THEN
                aux_ind=1
             ELSEIF (ibx==ic21) THEN
                aux_ind=2
             ELSEIF (ibx==ic211) THEN
                aux_ind=3
             ENDIF
          CASE (8)
             aux_ind=1
          CASE (12)
       !
       !  C_2v
       !
              iax=0
              ibx=0
              icx=0
              DO isym=1,nsym_b
!
!   find the axis of order 2 and the two mirrors of C_2v
!
                 IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                    CALL versor(sr_b(1,1,isym),ax)
                    CALL which_c2(ax, iax)
                 ELSEIF (tipo_sym(sr_b(1,1,isym))==5) THEN
                    CALL mirror_axis(sr_b(1,1,isym),ax)
                    IF (ibx==0) THEN
                       CALL which_c2(ax, ibx)
                    ELSEIF (icx==0) THEN
                       CALL which_c2(ax, icx)
                    ENDIF
                 ENDIF
              ENDDO
              CALL is_c2v(iax, ibx, icx, isok)
              IF (isok) THEN
                 isv=ibx
                 isv1=icx
              ELSE
                 CALL is_c2v(iax, icx, ibx, isok1)
                 IF (.NOT. isok1) CALL errore('find_aux_ind_two_groups',&
                                               'problem D_2h C_2v',1)
                 isv=icx
                 isv1=ibx
              ENDIF

              IF (iax==ic2) THEN
                 IF (isv==ic21 .AND. isv1==ic211) THEN
                    aux_ind=1
                 ELSEIF (isv==ic211 .AND. isv1==ic21) THEN
                    aux_ind=4
                 ELSE
                    CALL errore('find_aux_ind_two_groups','problem D_2h C_2v',2)
                 ENDIF
              ELSEIF (iax==ic21) THEN
                 IF (isv==ic2 .AND. isv1==ic211) THEN
                    aux_ind=2
                 ELSEIF (isv==ic211 .AND. isv1==ic2) THEN
                    aux_ind=5
                 ELSE
                    CALL errore('find_aux_ind_two_groups','problem D_2h C_2v',3)
                 ENDIF
              ELSEIF (iax==ic211) THEN
                 IF (isv==ic2 .AND. isv1==ic21) THEN
                    aux_ind=3
                 ELSEIF (isv==ic21 .AND. isv1==ic2) THEN
                    aux_ind=6
                 ELSE
                    CALL errore('find_aux_ind_two_groups','problem D_2h C_2v',4)
                 ENDIF
              ELSE
                 CALL errore('find_aux_ind_two_groups','problem of C_2 D_2h C_2v',5)
              ENDIF
         CASE(16)
!
!   find the C_2 axis and check with which axis it is parallel
!
             DO isym=1,nsym_b
                IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                   CALL versor(sr_b(1,1,isym),bx)
                   CALL which_c2(bx,ibx)
                   IF ( ibx==ic2 ) THEN
                      aux_ind=1
                   ELSEIF ( ibx==ic21 ) THEN
                      aux_ind=2
                   ELSEIF ( ibx==ic211 ) THEN
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
           CASE(4,5,9,13,17)
      !
      !    C_2, C_3, D_3, C_2v, C_3v, C_3h
      !
              aux_ind=1
           CASE(12)
             iax=0
             ibx=0
             icx=0
             DO isym=1,nsym_b
!
!   find the axis of order 2 and the two mirrors of C_2v
!
                IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                   CALL versor(sr_b(1,1,isym),ax)
                   CALL which_c2(ax, iax)
                ELSEIF (tipo_sym(sr_b(1,1,isym))==5) THEN
                   CALL mirror_axis(sr_b(1,1,isym),ax)
                   IF (ibx==0) THEN
                      CALL which_c2(ax, ibx)
                   ELSEIF (icx==0) THEN
                      CALL which_c2(ax, icx)
                   ENDIF
                ENDIF
             ENDDO
             CALL is_c2v(iax, ibx, icx, isok)
             IF (isok) THEN
                isv=ibx
                isv1=icx
             ELSE
                CALL is_c2v(iax, icx, ibx, isok1)
                IF (.NOT. isok1) CALL errore('find_aux_ind_two_groups',&
                                              'problem D_2h C_2v',1)
                isv=icx
                isv1=ibx
             ENDIF
             IF (isv==3) THEN
                aux_ind=1
             ELSEIF (isv1==3) THEN
                aux_ind=2
             ELSE
                CALL errore('find_aux_ind_two_groups',&
                                              'problem D_2h C_2v',2)
             ENDIF
           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(22)
!
!   D_4h
!
!
!   first find the axis C_4 of D_4h
!
       iax=0
       DO isym=2,nsym_a
          ts=tipo_sym(sr_a(1,1,isym))
          IF (ts==3) THEN
             CALL versor(sr_a(1,1,isym),ax)
             CALL which_c2(ax, iax)
          ENDIF
       ENDDO
       IF (iax > 3) CALL errore('find_aux_ind_two_groups','problem with D_4h',1)
       ic4=iax

       SELECT CASE (group_b)
          CASE(3)
          !
          !  C_s
          !
             CALL mirror_axis(sr_b(1,1,2),bx)
             CALL which_c2(bx, ibx)
             IF (ibx==ic4) THEN
                aux_ind=1
             ELSEIF ( ibx < 4 ) THEN
                aux_ind=2
             ELSE
                aux_ind=3
             ENDIF
          CASE(4)
          !
          !  C_2
          !
             CALL versor(sr_b(1,1,2),bx)
             CALL which_c2(bx, ibx)
             IF (ibx==ic4) THEN
                aux_ind=1
             ELSEIF (ibx < 4) THEN
                aux_ind=2
             ELSE
                aux_ind=3
             ENDIF
          CASE (2,6,10)
          !
          ! C_i, C_4, D_4
          !
             aux_ind=1
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
          CASE (12)
       !
       !  C_2v
       !
             iax=0
             ibx=0
             icx=0
             DO isym=1,nsym_b
!
!   find the axis of order 2 and the two mirrors of C_2v
!
                IF (tipo_sym(sr_b(1,1,isym))==4) THEN
                   CALL versor(sr_b(1,1,isym),ax)
                   CALL which_c2(ax, iax)
                ELSEIF (tipo_sym(sr_b(1,1,isym))==5) THEN
                   CALL mirror_axis(sr_b(1,1,isym),ax)
                   IF (ibx==0) THEN
                      CALL which_c2(ax, ibx)
                   ELSEIF (icx==0) THEN
                      CALL which_c2(ax, icx)
                   ENDIF
                ENDIF
             ENDDO
             CALL is_c2v(iax, ibx, icx, isok)
             IF (isok) THEN
                isv=ibx
                isv1=icx
             ELSE
                CALL is_c2v(iax, icx, ibx, isok1)
                IF (.NOT. isok1) CALL errore('find_aux_ind_two_groups',&
                                              'problem D_2h C_2v',1)
                isv=icx
                isv1=ibx
             ENDIF
             IF (iax==ic4) THEN
!
!   aux_num 1 and 2 the twofold axis is the z axis and in 1 the mirror
!   are perpendicular to the x and y axis, in 2 the mirror are perpendicular
!   to the 110 and 1-10 directions.
!

                IF ( isv < 4 ) THEN
                   aux_ind=1
                ELSE
                   aux_ind=2
                ENDIF
             ELSEIF (iax < 4) THEN
!
!  aux_num 3 when the axis is parallel to x,y or z
!
                IF (isv==ic4) THEN
                   aux_ind=3
                ELSEIF (isv1==ic4) THEN
                   aux_ind=4
                ELSE
                   CALL errore('find_aux_ind_two_groups', &
                               'D_4h problem with sigma_h',1)
                ENDIF
             ELSE
!
!  aux_num 5 when the axis of C_2v is not parallel to x,y or z
!
                aux_ind=5
             END IF

          CASE (16)
             DO isym=1,nsym_b
!
!   find the axis of order 2.
!
                IF (tipo_sym(sr_b(1,1,isym))==4) iaxis=isym
             ENDDO
             CALL versor(sr_b(1,1,iaxis),bx)
             CALL which_c2(bx, ibx)
             IF ( ibx==ic4 ) THEN
                aux_ind=1
             ELSEIF(ibx < 4) THEN
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
             iaxis=0
             DO isym=1,nsym_b
                IF (iaxis==0 .AND. tipo_sym(sr_b(1,1,isym))==4) THEN
                   CALL versor(sr_b(1,1,iaxis),bx)
                   CALL which_c2(bx, ibx)
                   IF (ibx /= ic4) iaxis=ibx
                ENDIF
             ENDDO
             IF (iaxis < 4) THEN
                aux_ind=1
             ELSE
                aux_ind=2
             ENDIF
          CASE(26)
       !
       !  S_4
       ! 
             aux_ind=1
          CASE DEFAULT
             CALL errore('find_aux_ind_two_groups','D_4h this is not a subgroup',1)
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
              CALL which_c2(bx, ibx)
              IF ( ibx==3 ) THEN
                 aux_ind=1
              ELSEIF (ibx==1 .OR. ibx==10 .OR. ibx==11) THEN
                 aux_ind=2
              ELSEIF (ibx==2 .OR. ibx==12 .OR. ibx==13) THEN
                 aux_ind=3
              ELSE
                 CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
              END IF
           CASE(4)
        !
        !  C_2
        ! 
              CALL versor(sr_b(1,1,2),bx)
              CALL which_c2(bx, ibx)
              IF ( ibx==3 ) THEN
                 aux_ind=1
              ELSEIF (ibx==1 .OR. ibx==10 .OR. ibx==11) THEN
                 aux_ind=2
              ELSEIF (ibx==2 .OR. ibx==12 .OR. ibx==13) THEN
                 aux_ind=3
              ELSE
                 CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
              END IF
           CASE(5,7,8)
      !
      !  C_3, C_6, D_2
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
           CASE(12)
      !
      !    C_2v
      !
              DO isym=1,nsym_b
                 IF (tipo_sym(sr_b(1,1,isym))==4) iaxis=isym
              END DO
              CALL versor(sr_b(1,1,iaxis),ax)
              CALL which_c2(ax, iax)
              IF (iax==3) THEN
                 aux_ind=1
              ELSEIF (iax==1) THEN
                 aux_ind=3
              ELSEIF (iax==10 .OR. iax==11) THEN
                 aux_ind=2
              ELSEIF (iax==2 .OR. iax==12 .OR. iax==13) THEN
                 aux_ind=4
              ELSE
                       CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
              ENDIF
           CASE(13)
      !
      !     C_3v
      !
              DO isym=1,nsym_b
                 IF (tipo_sym(sr_b(1,1,isym))==5) imirror=isym
              END DO
              CALL mirror_axis(sr_b(1,1,imirror),ax)
              CALL which_c2(ax, iax) 
              IF (iax==1 .OR. iax==10 .OR. iax==11) THEN
                 aux_ind=2
              ELSEIF (iax==2 .OR. iax==12 .OR. iax==13) THEN
                 aux_ind=1
              ELSE
                 CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
              ENDIF
          CASE(16)
      !
      !   C_2h
      !
              DO isym=1,nsym_b
                 IF (tipo_sym(sr_b(1,1,isym))==4) iaxis=isym
              END DO
              CALL versor(sr_b(1,1,iaxis),ax)
              CALL which_c2(ax, iax)
              IF (iax==3) THEN
                 aux_ind=1
              ELSEIF (iax==1 .OR. iax==10 .OR. iax==11) THEN
                 aux_ind=2
              ELSEIF (iax==2 .OR. iax==12 .OR. iax==13) THEN
                 aux_ind=3
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
           CASE(11,15,17,19,20,27)
      !
      !      D_6, C_6v, C_3h, C_6h, D_2h, S_6
      !
              aux_ind=1

           CASE DEFAULT
                CALL errore('find_aux_ind_two_groups',' group not available',1)
        END SELECT   
     CASE(24)
!
!  D_2d
!
        DO isym=1,nsym_a
           IF (tipo_sym(sr_a(1,1,isym))==6)  &
              CALL versor(sr_a(1,1,isym),bx)
        ENDDO
        CALL which_c2(bx, is4)

        SELECT CASE (group_b)
           CASE(3)
              aux_ind=1
           CASE(4)
!
!   Compare the direction of the versor of C_2 and of the rotation -4 in D_2d
!   If they are perpendicular the subgroup C_2 is of aux_ind=2
!
              CALL versor(sr_b(1,1,2),ax)
              prod = ax(1)*bx(1) + ax(2)*bx(2) + ax(3)*bx(3)
              IF (prod > 1.d-6) THEN
                 aux_ind=1
              ELSE
                 aux_ind=2
              ENDIF
           CASE (8)
        !
        !    D_2
        !
!
!   determine the three axes of D_2
!

              iax=0
              ibx=0
              icx=0
              DO isym=2,nsym_b
                 ts=tipo_sym(sr_b(1,1,isym))
                 IF (ts==4) THEN
                    CALL versor(sr_b(1,1,isym),ax)
                    IF (iax==0) THEN
                       CALL which_c2(ax, iax)
                    ELSEIF (ibx==0) THEN
                       CALL which_c2(ax, ibx)
                    ELSEIF (icx==0) THEN
                       CALL which_c2(ax, icx)
                    ELSE
                      CALL errore('find_aux_ind_two_groups','D_2 problem C_2 axis',1)
                    ENDIF
                 ENDIF
              ENDDO
              CALL is_d2(iax,ibx,icx,ind2)

              IF (ind2(1)==1) ic2=iax
              IF (ind2(1)==2) ic21=iax
              IF (ind2(1)==3) ic211=iax
              IF (ind2(2)==1) ic2=ibx
              IF (ind2(2)==2) ic21=ibx
              IF (ind2(2)==3) ic211=ibx
              IF (ind2(3)==1) ic2=icx
              IF (ind2(3)==2) ic21=icx
              IF (ind2(3)==3) ic211=icx
              IF (ic2==is4) THEN
                 aux_ind=1
              ELSEIF (ic21==is4) THEN 
                 aux_ind=2
              ELSEIF (ic211==is4) THEN
                 aux_ind=3
              ELSE
                 CALL errore('find_aux_ind_two_groups','problem D_2d and D_2',1)
              ENDIF

           CASE(12,26)
      !
      !      C_2v, S_4
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
!
!   T
!
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
!
!   T_h
!
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
!
!   T_d
!
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
!
!   O
!
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
           CASE(2,5,6,9,10)
      !
      !   C_i, C_3, C_4, D_3, D_4
      !
               aux_ind = 1 
           CASE(3)
      !
      !  C_s
      !
              CALL mirror_axis(sr_b(1,1,2),ax)
              CALL which_c2(ax, iax)
              IF (iax < 4) THEN
                 aux_ind=1
              ELSE
                 aux_ind=2
              ENDIF
           CASE(4)
      !
      !  C_2
      !
              CALL versor(sr_b(1,1,2),ax)
              CALL which_c2(ax, iax)
              IF (iax < 4) THEN
                 aux_ind=1
              ELSE
                 aux_ind=2
              ENDIF

           CASE(8)
!
!   Find the three axis of D_2
!
              iax=0
              ibx=0
              icx=0
              DO isym=2,nsym_b
                 ts=tipo_sym(sr_b(1,1,isym))
                 IF (ts==4) THEN
                    CALL versor(sr_b(1,1,isym),ax)
                    IF (iax==0) THEN
                       CALL which_c2(ax, iax)
                    ELSEIF (ibx==0) THEN
                       CALL which_c2(ax, ibx)
                    ELSEIF (icx==0) THEN
                       CALL which_c2(ax, icx)
                    ELSE
                      CALL errore('find_aux_ind_two_groups','D_2 problem C_2 axis',1)
                    ENDIF
                 ENDIF
              ENDDO

              IF (iax < 4 .AND. ibx < 4 .AND. icx < 4) THEN
                 aux_ind=1
              ELSE
                 aux_ind=2
              ENDIF

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
