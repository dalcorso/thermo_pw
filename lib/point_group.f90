!
! Copyright (C) 2014-2016 Andrea Dal Corso 
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

  CHARACTER(LEN=30) :: color_rap(12)

  DATA color_rap / 'color_red', 'color_green', 'color_blue','color_cyan', &
              'color_magenta', 'color_gold', 'color_pink', 'color_black', &
              'color_olive', 'color_brown', 'color_light_blue', 'color_orange' /

  CHARACTER(LEN=8) :: sym_label(64)
  DATA sym_label / 'E',   '2z', '2y',  '2x',   '2xy', '2x-y', '4-z', '4z',     &
              '2xz',   '2x-z', '4y', '4-y',   '2yz', '2y-z', '4-x', '4x',     &
              '3-x-y-z', '3-xyz',  '3xy-z', '3x-yz', '3xyz', '3-xy-z',        &
              '3x-y-z', '3-x-yz',  '6z', '6-z', '3z', '3-z', '21-10', '2210', & 
              '2010', '2110', &
              'i',   'i2z',   'i2y', 'i2x',  'i2xy', 'i2x-y', 'i4-z', 'i4z',  &
              'i2xz', 'i2x-z', 'i4y', 'i4-y', 'i2yz', 'i2y-z', 'i4-x', 'i4x', &
              'i3-x-y-z', 'i3-xyz', 'i3xy-z', 'i3x-yz', 'i3xyz', 'i3-xy-z',   &
              'i3x-y-z', 'i3-x-yz', 'i6z', 'i6-z', 'i3z', 'i3-z', 'i21-10', &
              'i2210', 'i2010', 'i2110' /

  PUBLIC convert_rap, find_aux_ind_two_groups, has_sigma_h, is_right_oriented,&
         color_rap, find_group_info_ext, sym_label

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
           CALL errore('convert_rap','Problem with degeneracy',1)
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
        sub_table(12,4,2,2)=2
        sub_table(12,4,3,1)=1
        sub_table(12,4,3,2)=4
        sub_table(12,4,4,1)=1
        sub_table(12,4,4,2)=3
        sub_table(12,4,5,1)=1
        sub_table(12,4,5,2)=2
        sub_table(12,4,6,1)=1
        sub_table(12,4,6,2)=1
        sub_table(12,4,7,1)=1
        sub_table(12,4,7,2)=3
        sub_table(12,4,8,1)=1
        sub_table(12,4,8,2)=4


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
  USE io_global, ONLY : stdout
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups',' group not available',2)
        END SELECT   
     CASE(7)
!
!   C_6
!
        SELECT CASE (group_b)
           CASE(4,5)
              aux_ind=1
           CASE DEFAULT
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups',' group not available',3)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups',' group not available',4)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups',' group not available',5)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups',' group not available',6)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups',' group not available',7)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups',' group not available',8)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups',' group not available',9)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups',' group not available',10)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups',' group not available',11)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups',' group not available',12)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups',' group not available',13)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups',' group not available',14)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups',' group not available',15)
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
             WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                              group_b
             CALL errore('find_aux_ind_two_groups','Group not available',16)
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
             WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                              group_b
             CALL errore('find_aux_ind_two_groups',' group not available',17)
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
             WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                              group_b
             CALL errore('find_aux_ind_two_groups','Group not available',17)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                              group_b
              CALL errore('find_aux_ind_two_groups',' Group not available',18)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                              group_b
              CALL errore('find_aux_ind_two_groups','Group not available',19)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups','Group not available',20)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups','Group not available',21)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups','Group not available',22)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups','Group not available',23)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups','Group not available',24)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups',' group not available',25)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups','group not available',26)
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
              WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                               group_b
              CALL errore('find_aux_ind_two_groups','group not available',27)
        END SELECT   
     CASE DEFAULT
       WRITE(stdout,'(5x,"Group a",i3," Group b", i3)') group_a, &
                                                        group_b
       CALL errore('find_aux_ind_two_groups',' group not available',1)
  END SELECT

  RETURN

  END SUBROUTINE find_aux_ind_two_groups


  SUBROUTINE find_group_info_ext(nsym, smat, code_group, code_group_ext, &
                                                  which_elem, group_desc)
!
!  This subroutine extends the find_group_info routine of QE.
!  In addition to the code of the group 1-32, it provides an extended code
!  of the group that identifies unambigously which point group it is 
!  among all the possible orientations of the axis. Accounting for the
!  orientation, there are a total of 136 different point groups. Furthermore 
!  it orders the symmetry elements, so that
!  it corresponds to the standard ordering.
!
!  The following table gives the extended code of each group and the
!  sequence of operations in each group. The numbers are the
!  same as in the symm_base routine of QE that lists 64 symmetry operations.
!  The extended point group number is identified from the sum of the square 
!  of the number of the symmetry operations that gives a unique signature 
!  to the group.
!  See also the point_group manual of thermo_pw for a table of the rotation
!  matrices that correspond to each element.
!
!
!  1)   C_1  E           1          = 1
!  2)   C_2  E  2z       1  2       = 5
!  3)   C_2  E  2y       1  3       = 10
!  4)   C_2  E  2x       1  4       = 17
!  5)   C_2  E  2xy      1  5       = 26
!  6)   C_2  E  2x-y     1  6       = 37
!  7)   C_2  E  2xz      1  9       = 82
!  8)   C_2  E  2x-z     1  10      = 101
!  9)   C_2  E  2yz      1  13      = 170
! 10)   C_2  E  2y-z     1  14      = 197
! 11)   C_2  E  21-10    1  29      = 842
! 12)   C_2  E  2210     1  30      = 901
! 13)   C_2  E  2010     1  31      = 962
! 14)   C_2  E  2110     1  32      = 1025
! 15)   C_s  E  i2z      1  34      = 1157
! 16)   C_s  E  i2y      1  35      = 1226
! 17)   C_s  E  i2x      1  36      = 1297
! 18)   C_s  E  i2xy     1  37      = 1370
! 19)   C_s  E  i2x-y    1  38      = 1445
! 20)   C_s  E  i2xz     1  41      = 1682
! 21)   C_s  E  i2x-z    1  42      = 1765
! 22)   C_s  E  i2yz     1  45      = 2026
! 23)   C_s  E  i2y-z    1  46      = 2117
! 24)   C_s  E  i21-10   1  61      = 3722
! 25)   C_s  E  i2210    1  62      = 3845
! 26)   C_s  E  i2010    1  63      = 3970
! 27)   C_s  E  i2110    1  64      = 4097
! 28)   C_i  E  i         1  33      = 1090
! 29)   C_3  E  3xyz    3-x-y-z 1  21  17   = 731
! 30)   C_3  E  3-xyz   3x-y-z  1  18  23   = 854
! 31)   C_3  E  3-x-yz  3xy-z   1  24  19   = 938
! 32)   C_3  E  3x-yz   3-xy-z  1  20  22   = 885
! 33)   C_3  E  3z      3-z     1  27  28   = 1514
! 34)   C_4  E  4z  2z  4-z     1  8   2  7  = 118
! 35)   C_4  E  4y  2y  4-y     1  11  3  12 = 275
! 36)   C_4  E  4x  2x  4-x     1  16  4  15 = 498
! 37)   C_6  E  6z  3z  2z  3-z 6-z  1  25  27  2  28  26 = 2819
! 38)   D_2  E  2z  2y  2x      1  2   3   4   = 30
! 39)   D_2  E  2z  2xy 2-xy    1  2   5   6   = 66
! 40)   D_2  E  2y  2xz 2-xz    1  3   9   10  = 191
! 41)   D_2  E  2x  2yz 2y-z    1  4   13  14  = 382
! 42)   D_2  E  2z  21-10 2110 1  2   29  32  = 1870
! 43)   D_2  E  2z  2210  2010 1  2   30  31  = 1866
! 44)   D_3  E  3z  3-z 2x 2110 2010   1  27 28 4 32 31 = 3515
! 45)   D_3  E  3z  3-z 2210  2y 21-10  1  27 28 30 3 29 = 3264
! 46)   D_3  E  3xyz   3-x-y-z 2y-z 2x-y 2x-z  1  21 17 14 6  10 = 1063
! 47)   D_3  E  3x-yz  3-xy-z 2yz 2x-z 2xy 1  20 22 13 10 5 = 1179 
! 48)   D_3  E  3-xyz  3x-y-z 2xz 2xy 2y-z   1  18 23 9 5  14 = 1156
! 49)   D_3  E  3-x-yz 3xy-z  2xz 2yz 2x-y 1  24 19 9  13 6 = 1224
! 50)   D_4  E  4z 2z  4-z 2x 2xy 2y 2x-y 1 8 2 7 4 5 3 6 = 204
! 51)   D_4  E  4y 2y  4-y 2z 2xz 2x 2x-z 1 11 3 12 2 9 4 10 = 476
! 52)   D_4  E  4x 2x  4-x 2y 2yz 2z 2y-z 1 15 4 16 3 13 2 14 =876
! 53)   D_6  E  6z 3z 2z 3-z 6-z 2x 2210 2110 2y 2010 21-10 
!                                   1 25 27 2 28 26 4 30 32 3 31 29 = 6570
! 54)   C_2v E  2x   i2y    i2z     1   4   35  34  = 2398
! 55)   C_2v E  2x   i2yz   i2y-z   1   4   45  46  = 4158
! 56)   C_2v E  2y   i2z    i2x     1   3   34  36  = 2462
! 57)   C_2v E  2y   i2xz   i2x-z   1   3   41  42  = 3455
! 58)   C_2v E  2z   i2y    i2x     1   2   35  36  = 2526
! 59)   C_2v E  2z   i2xy   i2x-y   1   2   37  38  = 2818
! 60)   C_2v E  2xy  i2z    i2x-y   1   5   34  38  = 2626
! 61)   C_2v E  2x-y i2z    i2xy    1   6   34  37  = 2562
! 62)   C_2v E  2xz  i2y    i2x-z   1   9   35  42  = 3071
! 63)   C_2v E  2x-z i2y    i2xz    1   10  35  41  = 3007
! 64)   C_2v E  2yz  i2x    i2y-z   1   13  36  46  = 3582
! 65)   C_2v E  2y-z i2x    i2yz    1   14  36  45  = 3518
! 66)   C_2v E  2z   i2210  i2010   1   2   62  63  = 7818
! 67)   C_2v E  2z   i21-10 i2110   1   2   61  64  = 7822
! 68)   C_2v E  2210  i2z   i2010   1   30  34  63  = 6026
! 69)   C_2v E  21-10 i2z   i2110   1   29  34  64  = 6094
! 70)   C_2v E  2110  i2z   i21-10  1   32  34  61  = 5902
! 71)   C_2v E  2010  i2z   i2210   1   31  34  62  = 5962
! 72)   C_3v E  3z 3-z i2x i2110 i2010 1 27 28 36 64 63 = 10875
! 73)   C_3v E  3z 3-z i2210 i2y i21-10  1 27 28 62 35 61 = 10304
! 74)   C_3v E  3xyz 3-x-y-z i2y-z i2x-y i2x-z 1 21 17 46 38 42 = 6055
! 75)   C_3v E  3x-yz 3-xy-z i2yz  i2x-z i2xy 1 20 22 45 42 37 = 6043
! 76)   C_3v E  3-xyz 3x-y-z i2xz i2xy  i2y-z  1 18 23 41 37 46 = 6020
! 77)   C_3v E  3-x-yz 3xy-z i2xz  i2yz i2x-y  1 24 19 41 45 38 = 6088
! 78)   C_4v E  4z 2z 4-z  i2x  i2xy  i2y i2x-y  1 8 2 7 36 37 35 38 = 5452
! 79)   C_4v E  4y 2y 4-y  i2z  i2xz  i2x i2x-z  1 11 3 12 34 41 36 42 = 6172
! 80)   C_4v E  4x 2x 4-x  i2y  i2yz  i2z i2y-z  1 16 4 15 35 45 34 46 = 7020
! 81)   C_6v E  6z 3z 2z   3-z  6-z   i2x i2210  i2110  i2y i2010 i21-10 
!                1 25 27 2 28 26 36 62 64 35 63 61 = 20970
! 82)   C_2h E  2z    i  i2z    1  2  33  34 = 2250
! 83)   C_2h E  2y    i  i2y    1  3  33  35 = 2324
! 84)   C_2h E  2x    i  i2x    1  4  33  36 = 2402
! 85)   C_2h E  2xy   i  i2xy   1  5  33  37 = 2484
! 86)   C_2h E  2x-y  i  i2x-y  1  6  33  38 = 2570
! 87)   C_2h E  2xz   i  i2xz   1  9  33  41 = 2852
! 88)   C_2h E  2x-z  i  i2x-z  1  10 33  42 = 2954
! 89)   C_2h E  2yz   i  i2yz   1  13 33  45 = 3284
! 90)   C_2h E  2y-z  i  i2y-z  1  14 33  46 = 3402
! 91)   C_2h E  2110  i  i2110  1  32 33  64 = 6210
! 92)   C_2h E  2010  i  i2010  1  31 33  63 = 6020 
! 93)   C_2h E  2210  i  i2210  1  30 33  62 = 5834 
! 94)   C_2h E  21-10 i  i21-10 1 29 33  61 = 5652
! 95)   C_3h E  6z  3z i2z 3-z i6-z 1 57 27 34 28 58  = 9283

! 96)   C_4h E  4z 2z 4-z i i24z i2z i24-z 1 8 2 7  33 40 34 39 = 5484
! 97)   C_4h E  4y 2y 4-y i i24y i2y i24-y 1 11 3 12 33 43 35 44 = 6374
! 98)   C_4h E  4x 2x 4-x i i24x i2x i24-x 1 16 4 15 33 48 36 47 = 7396
! 99)   C_6h E  6z 3z 2z 3-z 6-z i i6z i3z i2z i3-z i6-z 
!               1  25 27 2 28 26 33 57 59 34 60 58 = 18758 
! 100)  D_2h E  2z 2y   2x    i i2z i2y   i2x  1 2 3 4 33 34 35 36 = 4796
! 101)  D_2h E  2x 2yz  2y-z  i i2x i2yz  i2y-z  1 4 13 14 33 36 45 46 = 6908
! 102)  D_2h E  2y 2xz  2x-z  i i2y i2xz  i2x-z  1 3 9 10 33 35 41 42 = 5950
! 103)  D_2h E  2z 2xy  2x-y  i i2z i2xy  i2x-y  1 2 5 6 33 34 37 38 = 5124
! 104)  D_2h E  2z 2110 21-10 i i2z i2110 i21-10 1 2 32 29 33 34 64 61 = 11932
! 105)  D_2h E  2z 2010 2210  i i2z i2010 i2210  1 2 31 30 33 34 63 62 = 11924
! 106)  D_3h E  i6z 3z i2z 3-z i6-z 2x i2210 2110 i2y 2010 i21-10   
!                      1  57  27  34  28  58  4  62  32  35 31  61 = 20074
! 107)  D_3h E  i6z 3z i2z 3-z i6-z 2y i2010  21-10 i2x 2210 i2110   
!                      1  57  27  34  28  58  3  63  29  36  30  64 = 20394 
! 108)  D_4h E  4z 2z 4-z 2x 2xy 2y 2x-y i i4z i2z i4-z i2x 
!            s_xy s_y s_x-y 1 8 2 7 4  5  3  6  33 40 34 39 36 37 35 38 = 10904
! 109)  D_4h E  4y 2y 4-y 2z 2xz 2x 2x-z  i i4y i2y i4-y i2z i2xz i2x i2x-z 
!                           1 11 3 12  2  9 4 10 33 43 35 44 34 41 36 42 =12472
! 110)  D_4h E 4x 2x 4-x 2y 2yz 2z 2y-z i i4x i2x i4-x i2y i2yz i2z i2y-z 
!                          1 16 4 15 3 13 2 14 33 48 36 47 35 45 34 46 =  14296
! 111)  D_6h E 6z 3z 2z 3-z 6-z 2x 2210 2110 2y 2010 21-10
!            i i6z i3z i2z i3-z i6-z i2x  i2210  i2110  i2y i2010 i21-10
!            1 25 27 2 28 26 4 30 32 3 31 29 33 57 59 34 60 58 36 62 64 35
!            63 61 = 40660
! 112)  D_2d E i4-z 2z i4z 2x i2x-y 2y i2xy 1 39 2 40 4 38 3 37 = 5964

! 113)  D_2d E i4-z 2z i4z 2xy i2x   2x-y i2y   1 39 2 40 5 36 6 35   = 5708
! 114)  D_2d E i4-y 2y i4y 2x  i2xz  2z   i2x-z 1 44 3 43 4 41 2 42   = 7260
! 115)  D_2d E i4-y 2y i4y 2xz i2z   2x-z i2x   1 44 3 43 9 34 10 36  = 6428
! 116)  D_2d E i4-x 2x i4x 2y  i2y-z 2z   i2yz  1 47 4 48 3 46 2 45   = 8684
! 117)  D_2d E i4-x 2x i4x 2yz i2y   2y-z i2z   1 47 4 48 13 35 14 34 = 7276

! 118)  D_3d E 3z 3-z 2x 2110 2010 i i3z i3-z i2x i2110 i2010
!                    1 27 28 4 32 31 33 59 60 36 64 63 = 21046
! 119)  D_3d E 3z 3-z 2210 2y 21-10 i i3z i3-z i2210  i2y i21-10  
!                    1 27 28 30 3 29 33 59 60 62 35 61 = 20224
! 120)  D_3d E 3xyz 3-x-y-z 2y-z 2x-y 2x-z i i3xyz i3-x-y-z i2x-y 
!            i2x-z i2y-z  1 21 17 14 6 10 33 53 49 46 38 42  = 12686
! 121)  D_3d E 3x-yz 3-xy-z 2yz 2x-z 2xy i i3x-yz i3-xy-z i2yz 
!            i2x-z i2xy  1 20 22 13 10 5 33 52 54 45 42 37 = 13046
! 122)  D_3d E 3-xyz 3x-y-z 2xz 2xy 2y-z  i i3-xyz i3x-y-z i2xz i2xy 
!            i2y-z  1 18 23 9 5 14 33 50 55 41 37 46 = 12936
! 123)  D_3d E 3-x-yz 3xy-z 2xz 2yz 2x-y i i3-x-yz i3xy-z i2xz 
!            i2yz i2x-y 1 24 19 9 13 6 33 56 51 41 45 38 = 13200

! 124)  S_4  E i4z 2z i4-z 1 40 2 39 = 3126  
! 125)  S_4  E i4y 2y i4-y 1 43 3 44 = 3795
! 126)  S_4  E i4x 2x i4-x 1 48 4 47 = 4530

! 127)  S_6  E 3z 3-z i i3z i3-z 1 27 28 33 59 60 = 9684
! 128)  S_6  E 3xyz 3-x-y-z i i3xyz i3-x-y-z  1 21 17 33 53 49 = 7030  
! 129)  S_6  E 3-xyz 3x-y-z i i3-xyz i3x-y-z  1 18 23 33 50 55 = 7468
! 130)  S_6  E 3-x-yz 3xy-z i i3-x-yz i3xy-z  1 24 19 33 56 51 = 7764
! 131)  S_6  E 3x-yz 3-xy-z i i3x-yz i3-xy-z  1 20 22 33 52 54 = 7594
! 132)  T    E 2z  2y 2x 3xyz 3x-y-z 3-x-yz 3x-yz 3-x-y-z 
!              3-xyz  3xy-z 3x-yz 1 2 3 4 21 23 24 22 17 18 19 20 = 3434
! 133)  T_h  E 2z 2y 2x 3xyz 3x-y-z 3-x-yz 3x-yz 3-x-y-z 
!            3-xyz  3xy-z 3x-yz i i2z i2y i2x i3xyz i3x-y-z i3-x-yz 
!            i3x-yz i3-x-y-z  i3-xyz  i3xy-z i3x-yz 1 2 3 4 21 23 24 
!            22 17 18 19 20 33 34 35 26 53 55 56 54 49 50 51 52 = 30292
! 134)  T_d  E i4z i4-z i4y i4-y i4x i4-x 2z 2y 2x i2xy i2x-y 
!            i2yz i2y-z i2xz i2x-z 3xyz 3x-y-z 3-z-yz 3-xy-z 
!            3-x-y-z 3-xyz 3xy-z 3x-yz  1 39 40 43 44 47 48 2 3 4  
!            37 38 45 46 41 42  21 23 24 22 17 18 19 20  = 25252
! 135)  O    E 4z 4-z 4y 4-y 4x 4-x 2z 2y 2x 2xy 2x-y 
!            2yz 2y-z 2xz 2x-z 3xyz 3x-y-z 3-z-yz 3-xy-z 
!            3-x-y-z 3-xyz 3xy-z 3x-yz 1 7 8 11 12 15 16 2 3 4 5 6 
!            13 14 9 10 21 23 24 22 17 18 19 20  = 4900
! 136)  O_h  E 4z 4-z 4y 4-y 4x 4-x 2z 2y 2x 2xy 2x-y 
!            2yz 2y-z 2xz 2x-z 3xyz 3x-y-z 3-z-yz 3-xy-z 
!            3-x-y-z 3-xyz 3xy-z 3x-yz i i4z i4-z i4y i4-y i4x 
!            i4-x i2z i2y i2x i2xy i2x-y i2yz i2y-z i2xz i2x-z i3xyz 
!            i3x-y-z i3-z-yz i3-xy-z i3-x-y-z i3-xyz i3xy-z i3x-yz 
!            1 7 8 11 12 15 16 2 3 4 5 6  13 14 9 10 21 23 24 22 17 18 19 20  
!            33 39 40 43 44 47 48 34 35 36 37 38 45 46 41 42 53 55 56 54 49 
!            50 51 52  = 53576

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nsym
  REAL(DP), INTENT(IN) :: smat(3,3,nsym)
  INTEGER, INTENT(OUT) :: code_group, code_group_ext 
  INTEGER, INTENT(OUT) :: which_elem(nsym)
  INTEGER, INTENT(OUT) :: group_desc(48)

  INTEGER, PARAMETER :: npg = 136
  INTEGER :: hash_tab (npg)
  DATA hash_tab /   1,    5,   10,   17,   26,   37,   82,  101,  170,  197, &
                  842,  901,  962, 1025, 1157, 1226, 1297, 1370, 1445, 1682, &
                 1765, 2026, 2117, 3722, 3845, 3970, 4097, 1090,  731,  854, &
                  938,  885, 1514,  118,  275,  498, 2819,   30,   66,  191, &
                  382, 1870, 1866, 3515, 3264, 1063, 1179, 1156, 1224,  204, &
                  476,  876, 6570, 2398, 4158, 2462, 3455, 2526, 2818, 2626, &
                 2562, 3071, 3007, 3582, 3518, 7818, 7822, 6026, 6094, 5902, &
                 5962,10875,10304, 6055, 6043, 6020, 6088, 5452, 6172, 7020, &
                20970, 2250, 2324, 2402, 2484, 2570, 2852, 2954, 3284, 3402, &
                 6210, 6020, 5834, 5652, 9283, 5484, 6374, 7396,18758, 4796, &
                 6908, 5950, 5124,11932,11924,20074,20394,10904,12472,14296, & 
                40660, 5964, 5708, 7260, 6428, 8684, 7276,21046,20224,12686, &
                13046,12936,13200, 3126, 3795, 4530, 9684, 7030, 7468, 7764, &
                 7594, 3434,30292,25252, 4900,53576 /  

  INTEGER :: group_tags(nsym), group_tag

  INTEGER :: isym, jsym, i, nsym_

  CALL find_group_tags(nsym, smat, group_tags)

  group_tag=0
  DO isym=1,nsym
     group_tag = group_tag + group_tags(isym)**2
  ENDDO

  code_group_ext=0
  DO i=1,npg
     IF (group_tag==hash_tab(i)) code_group_ext=i
  ENDDO
  IF (code_group_ext==0) &
      CALL errore('find_group_info_ext','input group unknown',1)

  code_group = group_index_from_ext(code_group_ext)
!
!  set the description of the point group. This is the order of the
!  elements used in the table of projective representations
!
  CALL set_group_desc(group_desc, nsym_, code_group_ext)
  IF (nsym_ /= nsym) &
     CALL errore('find_group_info_ext','problem with code_group_ext',1)
 
!
!  for each element of the input symmetry group find its position
!  in the group description
!
  which_elem=0
  DO isym=1, nsym
     DO jsym=1,nsym
        IF (group_tags(isym)==group_desc(jsym)) THEN
           which_elem(isym)=jsym
           EXIT
        ENDIF
     ENDDO
     IF (which_elem(isym)==0) &
        CALL errore('find_group_info_ext','element not found',1)
  ENDDO

  RETURN
  END SUBROUTINE find_group_info_ext

  FUNCTION group_index_from_ext(code_group_ext)

  IMPLICIT NONE
  INTEGER :: group_index_from_ext, code_group_ext
  INTEGER, PARAMETER :: npg = 136
  INTEGER :: code_group_tab(npg)
  DATA code_group_tab / 1,  4,  4,  4,  4,  4,  4,  4,  4,  4,   &
                        4,  4,  4,  4,  3,  3,  3,  3,  3,  3,   &
                        3,  3,  3,  3,  3,  3,  3,  2,  5,  5,   &
                        5,  5,  5,  6,  6,  6,  7,  8,  8,  8,   &
                        8,  8,  8,  9,  9,  9,  9,  9,  9, 10,   &
                       10, 10, 11, 12, 12, 12, 12, 12, 12, 12,   &
                       12, 12, 12, 12, 12, 12, 12, 12, 12, 12,   &
                       12, 13, 13, 13, 13, 13, 13, 14, 14, 14,   &
                       15, 16, 16, 16, 16, 16, 16, 16, 16, 16,   &
                       16, 16, 16, 16, 17, 18, 18, 18, 19, 20,   &
                       20, 20, 20, 20, 20, 21, 21, 22, 22, 22,   &
                       23, 24, 24, 24, 24, 24, 24, 25, 25, 25,   &
                       25, 25, 25, 26, 26, 26, 27, 27, 27, 27,   &
                       27, 28, 29, 30, 31, 32 /  

  group_index_from_ext=code_group_tab(code_group_ext) 

  RETURN
  END FUNCTION group_index_from_ext


  SUBROUTINE set_group_desc(group_desc, nsym, code_group_ext)

  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: group_desc(48)
  INTEGER, INTENT(IN) :: code_group_ext
  INTEGER, INTENT(OUT) :: nsym

  INTEGER :: isym

  group_desc=0
  group_desc(1)=1
  SELECT CASE (code_group_ext)
     CASE (1)
!
!  C_2
!
     CASE (2)
       group_desc(2)=2
     CASE (3)
       group_desc(2)=3
     CASE (4)
       group_desc(2)=4
     CASE (5)
       group_desc(2)=5
     CASE (6)
       group_desc(2)=6
     CASE (7)
       group_desc(2)=9
     CASE (8)
       group_desc(2)=10
     CASE (9)
       group_desc(2)=13
     CASE (10)
       group_desc(2)=14
     CASE (11)
       group_desc(2)=29
     CASE (12)
       group_desc(2)=30
     CASE (13)
       group_desc(2)=31
     CASE (14)
       group_desc(2)=32
!
!   C_s
!
     CASE (15)
       group_desc(2)=34
     CASE (16)
       group_desc(2)=35
     CASE (17)
       group_desc(2)=36
     CASE (18)
       group_desc(2)=37
     CASE (19)
       group_desc(2)=38
     CASE (20)
       group_desc(2)=41
     CASE (21)
       group_desc(2)=42
     CASE (22)
       group_desc(2)=45
     CASE (23)
       group_desc(2)=46
     CASE (24)
       group_desc(2)=61
     CASE (25)
       group_desc(2)=62
     CASE (26)
       group_desc(2)=63
     CASE (27)
       group_desc(2)=64
!
!  C_i
!
     CASE (28)
       group_desc(2)=33
!
!  C_3    E   C_3   C_3^-1
!
     CASE (29)
       group_desc(2)=21
       group_desc(3)=17
     CASE (30)
       group_desc(2)=18
       group_desc(3)=23
     CASE (31)
       group_desc(2)=24
       group_desc(3)=19
     CASE (32)
       group_desc(2)=20
       group_desc(3)=22
     CASE (33)
       group_desc(2)=27
       group_desc(3)=28
!
!  C_4   E  C_4  C_2  C_4^-1
!
     CASE (34)
       group_desc(2)=8
       group_desc(3)=2
       group_desc(4)=7
     CASE (35)
       group_desc(2)=11
       group_desc(3)=3
       group_desc(4)=12
     CASE (36)
       group_desc(2)=16
       group_desc(3)=4
       group_desc(4)=15
!
!  C_6   C_6z  C_3z  C_2z  C_3-z  C_6-z
!
     CASE (37)
       group_desc(2)=25
       group_desc(3)=27
       group_desc(4)=2
       group_desc(5)=28
       group_desc(6)=26
!
!  D_2  order written above
!
     CASE (38)
       group_desc(2)=2
       group_desc(3)=3
       group_desc(4)=4
     CASE (39)
       group_desc(2)=2
       group_desc(3)=5
       group_desc(4)=6
     CASE (40)
       group_desc(2)=3
       group_desc(3)=9
       group_desc(4)=10
     CASE (41)
       group_desc(2)=4
       group_desc(3)=13
       group_desc(4)=14
     CASE (42)
       group_desc(2)=2
       group_desc(3)=29
       group_desc(4)=32
     CASE (43)
       group_desc(2)=2
       group_desc(3)=30
       group_desc(4)=31
!
!  D_3
!
     CASE (44)
       group_desc(2)=27
       group_desc(3)=28
       group_desc(4)=4
       group_desc(5)=32
       group_desc(6)=31
     CASE (45)
       group_desc(2)=27
       group_desc(3)=28
       group_desc(4)=30
       group_desc(5)=3
       group_desc(6)=29
     CASE (46)
       group_desc(2)=21
       group_desc(3)=17
       group_desc(4)=14
       group_desc(5)=6
       group_desc(6)=10
     CASE (47)
       group_desc(2)=20
       group_desc(3)=22
       group_desc(4)=13
       group_desc(5)=10
       group_desc(6)=5
     CASE (48)
       group_desc(2)=18
       group_desc(3)=23
       group_desc(4)=9
       group_desc(5)=5
       group_desc(6)=14
     CASE (49)
       group_desc(2)=24
       group_desc(3)=19
       group_desc(4)=9
       group_desc(5)=13
       group_desc(6)=6
!
!  D_4
!
     CASE (50)
       group_desc(2)=8
       group_desc(3)=2
       group_desc(4)=7
       group_desc(5)=4
       group_desc(6)=5
       group_desc(7)=3
       group_desc(8)=6
     CASE (51)
       group_desc(2)=11
       group_desc(3)=3
       group_desc(4)=12
       group_desc(5)=2
       group_desc(6)=9
       group_desc(7)=4
       group_desc(8)=10
     CASE (52)
       group_desc(2)=16
       group_desc(3)=4
       group_desc(4)=15
       group_desc(5)=3
       group_desc(6)=13
       group_desc(7)=2
       group_desc(8)=14
!
!  D_6
!
     CASE (53)
       group_desc(2)=25
       group_desc(3)=27
       group_desc(4)=2
       group_desc(5)=28
       group_desc(6)=26
       group_desc(7)=4
       group_desc(8)=30
       group_desc(9)=32
       group_desc(10)=3
       group_desc(11)=31
       group_desc(12)=29

!  C_2v
!
     CASE (54)
       group_desc(2)=4
       group_desc(3)=35
       group_desc(4)=34
     CASE (55)
       group_desc(2)=4
       group_desc(3)=45
       group_desc(4)=46
     CASE (56)
       group_desc(2)=3
       group_desc(3)=34
       group_desc(4)=36
     CASE (57)
       group_desc(2)=3
       group_desc(3)=41
       group_desc(4)=42
     CASE (58)
       group_desc(2)=2
       group_desc(3)=35
       group_desc(4)=36
     CASE (59)
       group_desc(2)=2
       group_desc(3)=37
       group_desc(4)=38
     CASE (60)
       group_desc(2)=5
       group_desc(3)=34
       group_desc(4)=38
     CASE (61)
       group_desc(2)=6
       group_desc(3)=34
       group_desc(4)=37
     CASE (62)
       group_desc(2)=9
       group_desc(3)=35
       group_desc(4)=42
     CASE (63)
       group_desc(2)=10
       group_desc(3)=35
       group_desc(4)=41
     CASE (64)
       group_desc(2)=13
       group_desc(3)=36
       group_desc(4)=46
     CASE (65)
       group_desc(2)=14
       group_desc(3)=36
       group_desc(4)=45
     CASE (66)
       group_desc(2)=2
       group_desc(3)=62
       group_desc(4)=63
     CASE (67)
       group_desc(2)=2
       group_desc(3)=61
       group_desc(4)=64
     CASE (68)
       group_desc(2)=30
       group_desc(3)=34
       group_desc(4)=63
     CASE (69)
       group_desc(2)=29
       group_desc(3)=34
       group_desc(4)=64
     CASE (70)
       group_desc(2)=32
       group_desc(3)=34
       group_desc(4)=61
     CASE (71)
       group_desc(2)=31
       group_desc(3)=34
       group_desc(4)=62
!
!  C_3v
!
     CASE (72)
       group_desc(2)=27
       group_desc(3)=28
       group_desc(4)=36
       group_desc(5)=64
       group_desc(6)=63
     CASE (73)
       group_desc(2)=27
       group_desc(3)=28
       group_desc(4)=62
       group_desc(5)=35
       group_desc(6)=61
     CASE (74)
       group_desc(2)=21
       group_desc(3)=17
       group_desc(4)=46
       group_desc(5)=38
       group_desc(6)=42
     CASE (75)
       group_desc(2)=20
       group_desc(3)=22
       group_desc(4)=45
       group_desc(5)=42
       group_desc(6)=37
     CASE (76)
       group_desc(2)=18
       group_desc(3)=23
       group_desc(4)=41
       group_desc(5)=37
       group_desc(6)=46
     CASE (77)
       group_desc(2)=24
       group_desc(3)=19
       group_desc(4)=41
       group_desc(5)=45
       group_desc(6)=38
!
!  C_4v
!
     CASE (78)
       group_desc(2)=8
       group_desc(3)=2
       group_desc(4)=7
       group_desc(5)=36
       group_desc(6)=37
       group_desc(7)=35
       group_desc(8)=38
     CASE (79)
       group_desc(2)=11
       group_desc(3)=3
       group_desc(4)=12
       group_desc(5)=34
       group_desc(6)=41
       group_desc(7)=36
       group_desc(8)=42

     CASE (80)
       group_desc(2)=16
       group_desc(3)=4
       group_desc(4)=15
       group_desc(5)=35
       group_desc(6)=45
       group_desc(7)=34
       group_desc(8)=46
!
!   C_6v
!
     CASE (81)
       group_desc(2)=25
       group_desc(3)=27
       group_desc(4)=2
       group_desc(5)=28
       group_desc(6)=26
       group_desc(7)=36
       group_desc(8)=62
       group_desc(9)=64
       group_desc(10)=35
       group_desc(11)=63
       group_desc(12)=61
!
!   C_2h
!
     CASE (82)
       group_desc(2)=2
       group_desc(3)=33
       group_desc(4)=34
     CASE (83)
       group_desc(2)=3
       group_desc(3)=33
       group_desc(4)=35
     CASE (84)
       group_desc(2)=4
       group_desc(3)=33
       group_desc(4)=36
     CASE (85)
       group_desc(2)=5
       group_desc(3)=33
       group_desc(4)=37
     CASE (86)
       group_desc(2)=6
       group_desc(3)=33
       group_desc(4)=38
     CASE (87)
       group_desc(2)=9
       group_desc(3)=33
       group_desc(4)=41
     CASE (88)
       group_desc(2)=10
       group_desc(3)=33
       group_desc(4)=42
     CASE (89)
       group_desc(2)=13
       group_desc(3)=33
       group_desc(4)=45
     CASE (90)
       group_desc(2)=14
       group_desc(3)=33
       group_desc(4)=46
     CASE (91)
       group_desc(2)=32
       group_desc(3)=33
       group_desc(4)=64
     CASE (92)
       group_desc(2)=31
       group_desc(3)=33
       group_desc(4)=63
     CASE (93)
       group_desc(2)=30
       group_desc(3)=33
       group_desc(4)=62
     CASE (94)
       group_desc(2)=29
       group_desc(3)=33
       group_desc(4)=61
!
!  C_3h
!
     CASE (95)
       group_desc(2)=57
       group_desc(3)=27
       group_desc(4)=34
       group_desc(5)=28
       group_desc(6)=58
!
!  C_4h
!
     CASE (96)
       group_desc(2)=8
       group_desc(3)=2
       group_desc(4)=7
       group_desc(5)=33
       group_desc(6)=40
       group_desc(7)=34
       group_desc(8)=39
     CASE (97)
       group_desc(2)=11
       group_desc(3)=3
       group_desc(4)=12
       group_desc(5)=33
       group_desc(6)=43
       group_desc(7)=35
       group_desc(8)=44
     CASE (98)
       group_desc(2)=16
       group_desc(3)=4
       group_desc(4)=15
       group_desc(5)=33
       group_desc(6)=48
       group_desc(7)=36
       group_desc(8)=47
!
!  C_6h
!
     CASE (99)
       group_desc(2)=25
       group_desc(3)=27
       group_desc(4)=2
       group_desc(5)=28
       group_desc(6)=26
       group_desc(7)=33
       group_desc(8)=57
       group_desc(9)=59
       group_desc(10)=34
       group_desc(11)=60
       group_desc(12)=58
!
!  D_2h
!
     CASE (100)
       group_desc(2)=2
       group_desc(3)=3
       group_desc(4)=4
       group_desc(5)=33
       group_desc(6)=34
       group_desc(7)=35
       group_desc(8)=36
     CASE (101)
       group_desc(2)=4
       group_desc(3)=13
       group_desc(4)=14
       group_desc(5)=33
       group_desc(6)=36
       group_desc(7)=45
       group_desc(8)=46
     CASE (102)
       group_desc(2)=3
       group_desc(3)=9
       group_desc(4)=10
       group_desc(5)=33
       group_desc(6)=35
       group_desc(7)=41
       group_desc(8)=42
     CASE (103)
       group_desc(2)=2
       group_desc(3)=5
       group_desc(4)=6
       group_desc(5)=33
       group_desc(6)=34
       group_desc(7)=37
       group_desc(8)=38
     CASE (104)
       group_desc(2)=2
       group_desc(3)=32
       group_desc(4)=29
       group_desc(5)=33
       group_desc(6)=34
       group_desc(7)=64
       group_desc(8)=61
     CASE (105)
       group_desc(2)=2
       group_desc(3)=31
       group_desc(4)=30
       group_desc(5)=33
       group_desc(6)=34
       group_desc(7)=63
       group_desc(8)=62
!
!  D_3h
!
     CASE (106)
!       group_desc(2)=58
!       group_desc(3)=28
!       group_desc(4)=34
!       group_desc(5)=27
!       group_desc(6)=57
!       group_desc(7)=4
!       group_desc(8)=61
!       group_desc(9)=31
!       group_desc(10)=35
!       group_desc(11)=32
!       group_desc(12)=62
       group_desc(2)=57
       group_desc(3)=27
       group_desc(4)=34
       group_desc(5)=28
       group_desc(6)=58
       group_desc(7)=4
       group_desc(8)=62
       group_desc(9)=32
       group_desc(10)=35
       group_desc(11)=31
       group_desc(12)=61
     CASE (107)
!       group_desc(2)=58
!       group_desc(3)=28
!       group_desc(4)=34
!       group_desc(5)=27
!       group_desc(6)=57
!       group_desc(7)=3
!       group_desc(8)=63
!       group_desc(9)=30
!       group_desc(10)=36
!       group_desc(11)=29
!       group_desc(12)=64
       group_desc(2)=57
       group_desc(3)=27
       group_desc(4)=34
       group_desc(5)=28
       group_desc(6)=58
       group_desc(7)=3
       group_desc(8)=63
       group_desc(9)=29
       group_desc(10)=36
       group_desc(11)=30
       group_desc(12)=64
!
!  D_4h
!
     CASE (108)
       group_desc(2)=8
       group_desc(3)=2
       group_desc(4)=7
       group_desc(5)=4
       group_desc(6)=5
       group_desc(7)=3
       group_desc(8)=6
       group_desc(9)=33
       group_desc(10)=40
       group_desc(11)=34
       group_desc(12)=39
       group_desc(13)=36
       group_desc(14)=37
       group_desc(15)=35
       group_desc(16)=38

     CASE (109)
       group_desc(2)=11
       group_desc(3)=3
       group_desc(4)=12
       group_desc(5)=2
       group_desc(6)=9
       group_desc(7)=4
       group_desc(8)=10
       group_desc(9)=33
       group_desc(10)=43
       group_desc(11)=35
       group_desc(12)=44
       group_desc(13)=34
       group_desc(14)=41
       group_desc(15)=36
       group_desc(16)=42

     CASE (110)
       group_desc(2)=16
       group_desc(3)=4
       group_desc(4)=15
       group_desc(5)=3
       group_desc(6)=13
       group_desc(7)=2
       group_desc(8)=14
       group_desc(9)=33
       group_desc(10)=48
       group_desc(11)=36
       group_desc(12)=47
       group_desc(13)=35
       group_desc(14)=45
       group_desc(15)=34
       group_desc(16)=46

!
!  D_6h
!
     CASE (111)
       group_desc(2)=25
       group_desc(3)=27
       group_desc(4)=2
       group_desc(5)=28
       group_desc(6)=26
       group_desc(7)=4
       group_desc(8)=30
       group_desc(9)=32
       group_desc(10)=3
       group_desc(11)=31
       group_desc(12)=29
       group_desc(13)=33
       group_desc(14)=57
       group_desc(15)=59
       group_desc(16)=34
       group_desc(17)=60
       group_desc(18)=58
       group_desc(19)=36
       group_desc(20)=62
       group_desc(21)=64
       group_desc(22)=35
       group_desc(23)=63
       group_desc(24)=61
!
!  D_2d
!
     CASE (112)
       group_desc(2)=39
       group_desc(3)=2
       group_desc(4)=40
       group_desc(5)=4
       group_desc(6)=38
       group_desc(7)=3
       group_desc(8)=37
     CASE (113)
       group_desc(2)=39
       group_desc(3)=2
       group_desc(4)=40
       group_desc(5)=5
       group_desc(6)=36
       group_desc(7)=6
       group_desc(8)=35
     CASE (114)
       group_desc(2)=44
       group_desc(3)=3
       group_desc(4)=43
       group_desc(5)=4
       group_desc(6)=41
       group_desc(7)=2
       group_desc(8)=42
     CASE (115)
       group_desc(2)=44
       group_desc(3)=3
       group_desc(4)=43
       group_desc(5)=9
       group_desc(6)=34
       group_desc(7)=10
       group_desc(8)=36
     CASE (116)
       group_desc(2)=47
       group_desc(3)=4
       group_desc(4)=48
       group_desc(5)=3
       group_desc(6)=46
       group_desc(7)=2
       group_desc(8)=45
     CASE (117)
       group_desc(2)=47
       group_desc(3)=4
       group_desc(4)=48
       group_desc(5)=13
       group_desc(6)=35
       group_desc(7)=14
       group_desc(8)=34
!
!  D_3d
!
     CASE (118)
       group_desc(2)=27
       group_desc(3)=28
       group_desc(4)=4
       group_desc(5)=32
       group_desc(6)=31
       group_desc(7)=33
       group_desc(8)=59
       group_desc(9)=60
       group_desc(10)=36
       group_desc(11)=64
       group_desc(12)=63

     CASE (119)
       group_desc(2)=27
       group_desc(3)=28
       group_desc(4)=30
       group_desc(5)=3
       group_desc(6)=29
       group_desc(7)=33
       group_desc(8)=59
       group_desc(9)=60
       group_desc(10)=62
       group_desc(11)=35
       group_desc(12)=61

     CASE (120)
       group_desc(2)=21
       group_desc(3)=17
       group_desc(4)=14
       group_desc(5)=6
       group_desc(6)=10
       group_desc(7)=33
       group_desc(8)=53
       group_desc(9)=49
       group_desc(10)=46
       group_desc(11)=38
       group_desc(12)=42

     CASE (121)
       group_desc(2)=20
       group_desc(3)=22
       group_desc(4)=13
       group_desc(5)=10
       group_desc(6)=5
       group_desc(7)=33
       group_desc(8)=52
       group_desc(9)=54
       group_desc(10)=45
       group_desc(11)=42
       group_desc(12)=37
     CASE (122)
       group_desc(2)=18
       group_desc(3)=23
       group_desc(4)=9
       group_desc(5)=5
       group_desc(6)=14
       group_desc(7)=33
       group_desc(8)=50
       group_desc(9)=55
       group_desc(10)=41
       group_desc(11)=37
       group_desc(12)=46
      CASE (123)
       group_desc(2)=24
       group_desc(3)=19
       group_desc(4)=9
       group_desc(5)=13
       group_desc(6)=6
       group_desc(7)=33
       group_desc(8)=56
       group_desc(9)=51
       group_desc(10)=41
       group_desc(11)=45
       group_desc(12)=38
!
!  S_4
!
     CASE (124)
       group_desc(2)=40
       group_desc(3)=2
       group_desc(4)=39
     CASE (125)
       group_desc(2)=43
       group_desc(3)=3
       group_desc(4)=44
     CASE (126)
       group_desc(2)=48
       group_desc(3)=4
       group_desc(4)=47
!
!  S_6
!
     CASE (127)
       group_desc(2)=27
       group_desc(3)=28
       group_desc(4)=33
       group_desc(5)=59
       group_desc(6)=60
     CASE (128)
       group_desc(2)=21
       group_desc(3)=17
       group_desc(4)=33
       group_desc(5)=53
       group_desc(6)=49
     CASE (129)
       group_desc(2)=18
       group_desc(3)=23
       group_desc(4)=33
       group_desc(5)=50
       group_desc(6)=55
     CASE (130)
       group_desc(2)=24
       group_desc(3)=19
       group_desc(4)=33
       group_desc(5)=56
       group_desc(6)=51
     CASE (131)
       group_desc(2)=20
       group_desc(3)=22
       group_desc(4)=33
       group_desc(5)=52
       group_desc(6)=54
!
!  T
!
     CASE (132)
       group_desc(2)=2
       group_desc(3)=3
       group_desc(4)=4
       group_desc(5)=21
       group_desc(6)=23
       group_desc(7)=24
       group_desc(8)=22
       group_desc(9)=17
       group_desc(10)=18
       group_desc(11)=19
       group_desc(12)=20
!
!  T_h
!
     CASE (133)

       group_desc(2)=2
       group_desc(3)=3
       group_desc(4)=4
       group_desc(5)=21
       group_desc(6)=23
       group_desc(7)=24
       group_desc(8)=22
       group_desc(9)=17
       group_desc(10)=18
       group_desc(11)=19
       group_desc(12)=20
       DO isym=1,12
          group_desc(12+isym)=group_desc(isym)+32
       ENDDO
!
!  T_d
!
     CASE (134)
       group_desc(2)=39
       group_desc(3)=40
       group_desc(4)=43
       group_desc(5)=44
       group_desc(6)=47
       group_desc(7)=48
       group_desc(8)=2
       group_desc(9)=3
       group_desc(10)=4
       group_desc(11)=37
       group_desc(12)=38
       group_desc(13)=45
       group_desc(14)=46
       group_desc(15)=41
       group_desc(16)=42
       group_desc(17)=21
       group_desc(18)=23
       group_desc(19)=24
       group_desc(20)=22
       group_desc(21)=17
       group_desc(22)=18
       group_desc(23)=19
       group_desc(24)=20
!
!   O
!
     CASE (135)
       group_desc(2)=7
       group_desc(3)=8
       group_desc(4)=11
       group_desc(5)=12
       group_desc(6)=15
       group_desc(7)=16
       group_desc(8)=2
       group_desc(9)=3
       group_desc(10)=4
       group_desc(11)=5
       group_desc(12)=6
       group_desc(13)=13
       group_desc(14)=14
       group_desc(15)=9
       group_desc(16)=10
       group_desc(17)=21
       group_desc(18)=23
       group_desc(19)=24
       group_desc(20)=22
       group_desc(21)=17
       group_desc(22)=18
       group_desc(23)=19
       group_desc(24)=20
!
!   O_h
!
     CASE (136)
       group_desc(2)=7
       group_desc(3)=8
       group_desc(4)=11
       group_desc(5)=12
       group_desc(6)=15
       group_desc(7)=16
       group_desc(8)=2
       group_desc(9)=3
       group_desc(10)=4
       group_desc(11)=5
       group_desc(12)=6
       group_desc(13)=13
       group_desc(14)=14
       group_desc(15)=9
       group_desc(16)=10
       group_desc(17)=21
       group_desc(18)=23
       group_desc(19)=24
       group_desc(20)=22
       group_desc(21)=17
       group_desc(22)=18
       group_desc(23)=19
       group_desc(24)=20
       DO isym=1,24
          group_desc(isym+24) = 32+group_desc(isym)
       ENDDO
     CASE DEFAULT
       CALL errore('set_group_desc','group index not found',1)
  END SELECT

  nsym=0
  DO isym=1,48
     IF (group_desc(isym) > 0) nsym=isym
  END DO

  RETURN
  END SUBROUTINE set_group_desc



SUBROUTINE find_group_tags(nsym, smat, group_tags)
!
!  This routine finds for all the symmetries of a point group the
!  number of the symmetry operation in the list of symmetries
!
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nsym
  REAL(DP), INTENT(IN) :: smat(3,3,nsym)
  INTEGER, INTENT(OUT) :: group_tags(nsym)

  INTEGER :: isym, i
  INTEGER :: tipo_sym, ts  
  REAL(DP), PARAMETER :: sqr2=SQRT(2.0_DP), sqr3=SQRT(3.0_DP)
  REAL(DP) :: a(3,32), b(3,32), angl(32), angle_rot, ang, s(3,3), ax(3)

  group_tags=0
  a=0.0_DP
  a(3,2)=1.0_DP
  a(2,3)=1.0_DP
  a(1,4)=1.0_DP
  a(1,5)=1.0_DP/sqr2
  a(2,5)=1.0_DP/sqr2
  a(1,6)=1.0_DP/sqr2
  a(2,6)=-1.0_DP/sqr2
  a(1,9)=1.0_DP/sqr2
  a(3,9)=1.0_DP/sqr2
  a(1,10)=1.0_DP/sqr2
  a(3,10)=-1.0_DP/sqr2
  a(2,13)=1.0_DP/sqr2
  a(3,13)=1.0_DP/sqr2
  a(2,14)=1.0_DP/sqr2
  a(3,14)=-1.0_DP/sqr2
  a(1,29)=sqr3/2.0_DP
  a(2,29)=-0.5_DP
  a(1,30)=sqr3/2.0_DP
  a(2,30)=0.5_DP
  a(1,31)=-0.5_DP
  a(2,31)=sqr3/2.0_DP
  a(1,32)=0.5_DP
  a(2,32)=sqr3/2.0_DP
  b=0.0_DP
  angl=0.0_DP
  angl(7)=270._DP
  b(3,7)=1.0_DP
  angl(8)=90._DP
  b(3,8)=1.0_DP
  angl(11)=90._DP
  b(2,11)=1.0_DP
  angl(12)=270._DP
  b(2,12)=1.0_DP
  angl(15)=270._DP
  b(1,15)=1.0_DP
  angl(16)=90._DP
  b(1,16)=1.0_DP
  angl(17)=240._DP
  b(1,17)=1.0_DP/sqr3
  b(2,17)=1.0_DP/sqr3
  b(3,17)=1.0_DP/sqr3
  angl(18)=120._DP
  b(1,18)=-1.0_DP/sqr3
  b(2,18)=1.0_DP/sqr3
  b(3,18)=1.0_DP/sqr3
  angl(19)=240._DP
  b(1,19)=-1.0_DP/sqr3
  b(2,19)=-1.0_DP/sqr3
  b(3,19)=1.0_DP/sqr3
  angl(20)=120._DP
  b(1,20)=1.0_DP/sqr3
  b(2,20)=-1.0_DP/sqr3
  b(3,20)=1.0_DP/sqr3
  angl(21)=120._DP
  b(:,21)=b(:,17)
  angl(22)=240._DP
  b(:,22)=b(:,20)
  angl(23)=240._DP
  b(:,23)=b(:,18)
  angl(24)=120._DP
  b(:,24)=b(:,19)
  angl(25)=60._DP
  b(3,25)=1.0_DP
  angl(26)=300._DP
  b(3,26)=1.0_DP
  angl(27)=120._DP
  b(3,27)=1.0_DP
  angl(28)=240._DP
  b(3,28)=1.0_DP

  DO isym=1,nsym
     ts=tipo_sym(smat(1,1,isym))
     IF (ts==1) THEN
        group_tags(isym) = 1
     ELSEIF (ts==2) THEN
        group_tags(isym) = 33
     ELSEIF (ts==3) THEN
        CALL versor(smat(1,1,isym),ax)
        ang=angle_rot(smat(1,1,isym))
        DO i=2,32
           IF (ABS(ax(1)*b(1,i)+ax(2)*b(2,i)+ax(3)*b(3,i)-1.0_DP)<1.D-8 &
              .AND. ABS(ang-angl(i))<1.D-8 ) THEN
              group_tags(isym)=i
            END IF
        END DO
     ELSEIF (ts==4) THEN
        CALL versor(smat(1,1,isym),ax)
        DO i=2,32
           IF (ABS(ax(1)*a(1,i)+ax(2)*a(2,i)+ax(3)*a(3,i)-1.0_DP)<1.D-8 &
              .OR. ABS(ax(1)*a(1,i)+ax(2)*a(2,i)+ax(3)*a(3,i)+1.0_DP)<1.D-8) &
                                                                          THEN
              group_tags(isym)=i
           ENDIF
        ENDDO
     ELSEIF (ts==5) THEN
        s=-smat(:,:,isym)
        CALL versor(s,ax)
        DO i=2,32
           IF (ABS(ax(1)*a(1,i)+ax(2)*a(2,i)+ax(3)*a(3,i)-1.0_DP)<1.D-8 &
              .OR. ABS(ax(1)*a(1,i)+ax(2)*a(2,i)+ax(3)*a(3,i)+1.0_DP)<1.D-8) &
                                                                          THEN
              group_tags(isym)=i+32
           ENDIF
        ENDDO
     ELSEIF (ts==6) THEN
        s=-smat(:,:,isym)
        CALL versor(s,ax)
        ang=angle_rot(s)
        DO i=2,32
           IF (ABS(ax(1)*b(1,i)+ax(2)*b(2,i)+ax(3)*b(3,i)-1.0_DP)<1.D-8 &
              .AND. ABS(ang-angl(i))<1.D-8 ) THEN
              group_tags(isym)=i+32
            END IF
        END DO
     ENDIF   
     IF (group_tags(isym)==0) &
        CALL errore('find_group_tags','problem identifying symmetry', isym)   
  END DO


RETURN
END SUBROUTINE find_group_tags

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


  FUNCTION is_subgroup(group_in_ext, group_out_ext)
  !
  !  The function returns .TRUE. if group_out is a subgroup of group_in
  !
  IMPLICIT NONE
  LOGICAL :: is_subgroup
  INTEGER, INTENT(IN) :: group_in_ext, group_out_ext

  INTEGER :: nsym_in, nsym_out
  INTEGER :: group_desc_in(48), group_desc_out(48)
  INTEGER :: isym, jsym, irot

  CALL set_group_desc(group_desc_in,nsym_in,group_in_ext)
  CALL set_group_desc(group_desc_out,nsym_out,group_out_ext)
  
  is_subgroup=.TRUE.
  DO isym=1, nsym_out
     irot=group_desc_out(isym)
     DO jsym=1, nsym_in
        IF (group_desc_in(jsym)==irot) GOTO 100
     ENDDO
     is_subgroup=.FALSE.
     RETURN
100  CONTINUE
  ENDDO
  RETURN
  END FUNCTION is_subgroup

  SUBROUTINE set_sym_matrix(sym_num, smat, sinv)
  !
  !  This routine uses the Cayley-Klein parameters to set the rotations
  !  matrix for the 32 proper rotations. See sym_label for the name of each
  !  rotation. The integer sinv is set to 1 if the operation has not inversion,
  !  to -1 if the operation is to be multiplied by inversion. Note that
  !  inversion does not act on spin, so the matrices remain the same.
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: sym_num
  INTEGER, INTENT(OUT) :: sinv
  COMPLEX(DP), INTENT(OUT) :: smat(2,2)

  INTEGER :: snum
  COMPLEX(DP) :: a(32), b(32)
  REAL(DP) :: sqrt2=SQRT(2.0_DP), sqrt3=SQRT(3.0_DP)

  IF (sym_num < 1 .OR. sym_num > 64) CALL errore('set_sym_matrix', &
              'problem with symmetry number',1)

  a=(0.0_DP,0.0_DP)
  b=(0.0_DP,0.0_DP)

  a(1)=(1.0_DP,0.0_DP)
  a(2)=(0.0_DP,-1.0_DP)
  b(3)=(-1.0_DP,0.0_DP)
  b(4)=(0.0_DP,-1.0_DP)
  b(5)=(-1.0_DP,-1.0_DP) / sqrt2
  b(6)=(1.0_DP,-1.0_DP) / sqrt2
  a(7)=(1.0_DP,1.0_DP) / sqrt2
  a(8)=(1.0_DP,-1.0_DP) / sqrt2
  a(9)=(0.0_DP,-1.0_DP) / sqrt2
  b(9)=(0.0_DP,-1.0_DP) / sqrt2
  a(10)=(0.0_DP,-1.0_DP) / sqrt2
  b(10)=(0.0_DP,1.0_DP) / sqrt2
  a(11)=(1.0_DP,0.0_DP) / sqrt2
  b(11)=(-1.0_DP,0.0_DP) / sqrt2
  a(12)=(1.0_DP,0.0_DP) / sqrt2
  b(12)=(1.0_DP,0.0_DP) / sqrt2
  a(13)=(0.0_DP,-1.0_DP) / sqrt2
  b(13)=(-1.0_DP,0.0_DP) / sqrt2
  a(14)=(0.0_DP,-1.0_DP) / sqrt2
  b(14)=(1.0_DP,0.0_DP) / sqrt2
  a(15)=(1.0_DP,0.0_DP) / sqrt2
  b(15)=(0.0_DP,1.0_DP) / sqrt2
  a(16)=(1.0_DP,0.0_DP) / sqrt2
  b(16)=(0.0_DP,-1.0_DP) / sqrt2
  a(17)=(1.0_DP,1.0_DP) / 2.0_DP
  b(17)=(1.0_DP,1.0_DP) / 2.0_DP
  a(18)=(1.0_DP,-1.0_DP) / 2.0_DP
  b(18)=(-1.0_DP,1.0_DP) / 2.0_DP
  a(19)=(1.0_DP,1.0_DP) / 2.0_DP
  b(19)=(-1.0_DP,-1.0_DP) / 2.0_DP
  a(20)=(1.0_DP,-1.0_DP) / 2.0_DP
  b(20)=(1.0_DP,-1.0_DP) / 2.0_DP
  a(21)=(1.0_DP,-1.0_DP) / 2.0_DP
  b(21)=(-1.0_DP,-1.0_DP) / 2.0_DP
  a(22)=(1.0_DP,1.0_DP) / 2.0_DP
  b(22)=(-1.0_DP,1.0_DP) / 2.0_DP
  a(23)=(1.0_DP,1.0_DP) / 2.0_DP
  b(23)=(1.0_DP,-1.0_DP) / 2.0_DP
  a(24)=(1.0_DP,-1.0_DP) / 2.0_DP
  b(24)=(1.0_DP,1.0_DP) / 2.0_DP
  a(25)=CMPLX(sqrt3,-1.0_DP) / 2.0_DP
  a(26)=CMPLX(sqrt3,1.0_DP) / 2.0_DP
  a(27)=CMPLX(1.0_DP,-sqrt3) / 2.0_DP
  a(28)=CMPLX(1.0_DP,sqrt3) / 2.0_DP
  b(29)=CMPLX(1.0_DP,-sqrt3) / 2.0_DP
  b(30)=CMPLX(-1.0_DP,-sqrt3) / 2.0_DP
  b(31)=CMPLX(sqrt3,-1.0_DP) / 2.0_DP
  b(32)=CMPLX(-sqrt3,-1.0_DP) / 2.0_DP

  sinv=1
  snum=sym_num
  IF (sym_num > 32) THEN
     sinv=-1
     snum=sym_num-32
  ENDIF
  smat(1,1) = a(snum)
  smat(2,2) = CONJG(a(snum))
  smat(1,2) = b(snum)
  smat(2,1) = -CONJG(b(snum))

  RETURN
  END SUBROUTINE set_sym_matrix

  SUBROUTINE find_product_table(prd, code_group_ext)

  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: prd(48,48)
  INTEGER, INTENT(IN) :: code_group_ext

  REAL(DP) :: s0(3,3,32), group_mat(3,3,48), a_mat(3,3), b_mat(3,3), &
              c_mat(3,3) 
  real(DP), PARAMETER :: sin3 = 0.866025403784438597d0, cos3 = 0.5d0, &
                        msin3 =-0.866025403784438597d0, mcos3 = -0.5d0
  INTEGER :: group_desc(48)
  INTEGER :: isym, jsym, ksym, nsym

  ! full name of the rotational part of each symmetry operation

  data s0/ 1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
          -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
          -1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
           1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
           0.d0,  1.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
           0.d0, -1.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
           0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
           0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
           0.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  0.d0, &
           0.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0,  0.d0, &
           0.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0,  0.d0, &
           0.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  0.d0, &
          -1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0, &
          -1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0, &
           1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0, &
           1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0, &
           0.d0,  0.d0,  1.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, &
           0.d0,  0.d0, -1.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, &
           0.d0,  0.d0, -1.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, &
           0.d0,  0.d0,  1.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, &
           0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  1.d0,  0.d0,  0.d0, &
           0.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  1.d0,  0.d0,  0.d0, &
           0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0, -1.d0,  0.d0,  0.d0, &
           0.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0, -1.d0,  0.d0,  0.d0, &
           cos3,  sin3, 0.d0, msin3,  cos3, 0.d0, 0.d0, 0.d0,  1.d0, &
           cos3, msin3, 0.d0,  sin3,  cos3, 0.d0, 0.d0, 0.d0,  1.d0, &
          mcos3,  sin3, 0.d0, msin3, mcos3, 0.d0, 0.d0, 0.d0,  1.d0, &
          mcos3, msin3, 0.d0,  sin3, mcos3, 0.d0, 0.d0, 0.d0,  1.d0, &
           cos3, msin3, 0.d0, msin3, mcos3, 0.d0, 0.d0, 0.d0, -1.d0, &
           cos3,  sin3, 0.d0,  sin3, mcos3, 0.d0, 0.d0, 0.d0, -1.d0, &
          mcos3, msin3, 0.d0, msin3,  cos3, 0.d0, 0.d0, 0.d0, -1.d0, &
          mcos3,  sin3, 0.d0,  sin3,  cos3, 0.d0, 0.d0, 0.d0, -1.d0 /


  CALL set_group_desc(group_desc,nsym,code_group_ext)

  DO isym=1, nsym
     IF (group_desc(isym) < 33 ) THEN
        group_mat(:,:,isym) = s0(:,:,group_desc(isym))
     ELSE
        group_mat(:,:,isym) = -s0(:,:,group_desc(isym)-32)
     ENDIF
  ENDDO

  prd=0
  DO isym=1, nsym
     a_mat(:,:) = group_mat(:,:,isym)
     DO jsym=1, nsym
        b_mat(:,:) = group_mat(:,:,jsym)
        c_mat=MATMUL(a_mat, b_mat)
        DO ksym = 1, nsym
           IF (SUM(ABS(c_mat(:,:) - group_mat(:,:,ksym))) < 1.D-8) &
                                 prd(isym,jsym)=ksym
        END DO
        IF (prd(isym,jsym)==0) &
           CALL errore('find_product_table','problem with matrices',1)
     END DO
  END DO

  RETURN
  END SUBROUTINE find_product_table

 SUBROUTINE find_double_product_table(prd, epos, code_group_ext)
!
!  This routine provides the product table prd of a given double group.
!  The symmetry operations are only those without E'. If necessary
!  the table can be extended to the entire group using EE=E, EE'=E'E=E', E'E'=E
!
!  The routine provides also the factor system epos that one should use 
!  if the character tables of the double group are used as projective
!  representations of the point group
!
 USE kinds, ONLY : DP
 IMPLICIT NONE
 INTEGER, INTENT(IN)  :: code_group_ext
 INTEGER, INTENT(OUT) :: prd(48,48), epos(48,48)

 INTEGER :: group_desc(48), sinv(48), isym, jsym, ksym, sinv1, sinv2, sinvp, &
            pcount, nsym
 COMPLEX(DP) :: a_mat(2,2), b_mat(2,2), c_mat(2,2), group_mat(2,2,48)
 REAL(DP) :: diff, diff1

 CALL set_group_desc(group_desc,nsym,code_group_ext)
 DO isym=1,nsym
    CALL set_sym_matrix(group_desc(isym), group_mat(1,1,isym), sinv(isym))
!    WRITE(6,*) 'isym', isym, group_desc(isym)
 ENDDO

 epos=1
 prd=0
 DO isym=1, nsym
    a_mat(:,:) = group_mat(:,:,isym)
    sinv1=sinv(isym)
    DO jsym=1, nsym
       b_mat(:,:) = group_mat(:,:,jsym)
       sinv2=sinv(jsym)
       sinvp = sinv1*sinv2
       c_mat(1,1)=a_mat(1,1) * b_mat(1,1) + a_mat(1,2) * b_mat(2,1)
       c_mat(1,2)=a_mat(1,1) * b_mat(1,2) + a_mat(1,2) * b_mat(2,2)
!       WRITE(6,*) 'isym,jsym a', isym, jsym,c_mat(1,1)
!       WRITE(6,*) 'isym,jsym b', isym, jsym,c_mat(1,2)
       pcount=0
       DO ksym = 1, nsym
          IF (sinv(ksym)==sinvp) THEN
             diff=ABS(c_mat(1,1)-group_mat(1,1,ksym))+ &
                  ABS(c_mat(1,2)-group_mat(1,2,ksym)) 
             IF (diff < 1.D-6) THEN
                prd(isym,jsym)=ksym
                pcount=pcount+1
             ELSE
                diff1=ABS(c_mat(1,1)+group_mat(1,1,ksym))+ &
                      ABS(c_mat(1,2)+group_mat(1,2,ksym)) 
                IF (diff1 < 1.D-6) THEN
                   prd(isym,jsym)=ksym
                   epos(isym,jsym)=-1
                   pcount=pcount+1
                ENDIF
             ENDIF
          END IF
       END DO
       IF (pcount/=1) &
          CALL errore('find_double_product_table','problem with matrices',1)
    END DO
 END DO

 RETURN
 END SUBROUTINE find_double_product_table
!

 FUNCTION nsym_group(code_group)
!
!  this function receives the code of the group and gives as output
!  the number of symmetry operations of the group
!
!   1  "C_1 " 1    11 "D_6 " 12   21 "D_3h" 12    31 "O   " 24
!   2  "C_i " 2    12 "C_2v" 4    22 "D_4h" 16    32 "O_h " 48
!   3  "C_s " 2    13 "C_3v" 6    23 "D_6h" 24 
!   4  "C_2 " 2    14 "C_4v" 8    24 "D_2d" 8
!   5  "C_3 " 3    15 "C_6v" 12   25 "D_3d" 12
!   6  "C_4 " 4    16 "C_2h" 4    26 "S_4 " 4
!   7  "C_6 " 6    17 "C_3h" 6    27 "S_6 " 6
!   8  "D_2 " 4    18 "C_4h" 8    28 "T   " 12
!   9  "D_3 " 6    19 "C_6h" 12   29 "T_h " 24
!   10 "D_4 " 8    20 "D_2h" 8    30 "T_d " 24

  IMPLICIT NONE
  INTEGER :: nsym_group
  INTEGER, INTENT(IN) :: code_group

  INTEGER :: nelem(32)
  DATA nelem / 1, 2,  2, 2,  3,  4,  6, 4,  6, 8, 12,  4,  6,  8, 12, 4, &
               6, 8, 12, 8, 12, 16, 24, 8, 12, 4,  6, 12, 24, 24, 24, 48 /

  nsym_group= nelem(code_group)

  RETURN
  END FUNCTION nsym_group

  END MODULE point_group
