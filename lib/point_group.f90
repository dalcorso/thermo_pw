!
! Copyright (C) 2014-2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE point_group
!
!  This module contains variables and routines to deal with the 
!  crystallographic point group symmetry. It complements the routines 
!  in find_group.f90, divide_class.f90 and divide_class_so.f90 in the 
!  PW/src directory of the QE package.
!  The conventions, such as the code group, the symmetry operation types,
!  the irreducible representations etc. are the same.
!
!  Presently it has routines to perform the following tasks:
!
!  Given a point group it can find which group it is, not only its name,
!  but also the sequence of symmetry operations.
!  Given a factor system there are routines that find to which p-equivalence
!  class the representation belong. There is a routine that finds the
!  gauge transformation that converts a given factor system in the
!  standard ones and a routine that sets the characters of the
!  irreducible (also projective) representations.
!
!  Given two point groups, the second a subgroup of the first,
!  and a list of representations of the first point group, there is a routine
!  that transforms it in a list of representations of the second group.
!  Given two point groups, the second a subgroup of the first, 
!  find which type it is. The different cases are discussed in the point-group
!  manual in the thermo_pw/Doc directory. The routine find_aux_ind_two_groups,
!  receives the rotation matrices of the two groups and gives an index that
!  correspond to the case.
!  Double groups are supported, however in this case the distinction between
!  different subgroup cases is irrelevant and not used.
!
!  There is also a routine able to transform the projective representations
!  into the appropriate representations of the subgroup. A routine sets
!  the standard factor system of the group and finds the factor system
!  of the subgroup. The user specifies in input both the ptype of the 
!  representation (of the point group, of the double point group or
!  projective) and must also specify in which type of representations
!  it expects that the representation can be decomposed. The routine
!  checks that this is actually the case and gives the decomposition.
!  See the routine find_projection_type for the definition of ptype.
!
!  The module provides also auxiliary routines to find the multiplication
!  table of the point group and of the double point group, to 
!  print on output the character tables or a factor system or a 
!  multiplication table, to set the o3 or the su2 rotations matrices
!  that correspond to each operation.
!
!  Among the variables that are offered by the modulus there is a 
!  list of colors for each irreducible representation and a short
!  name of the symmetry operation.
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
  DATA sym_label / 'E',   '2z', '2y',  '2x',   '2xy', '2x-y', '4-z', '4z',    &
              '2xz',   '2x-z', '4y', '4-y',   '2yz', '2y-z', '4-x', '4x',     &
              '3-x-y-z', '3-xyz',  '3xy-z', '3x-yz', '3xyz', '3-xy-z',        &
              '3x-y-z', '3-x-yz',  '6z', '6-z', '3z', '3-z', '21-10', '2210', & 
              '2010', '2110', &
              'i',   'i2z',   'i2y', 'i2x',  'i2xy', 'i2x-y', 'i4-z', 'i4z',  &
              'i2xz', 'i2x-z', 'i4y', 'i4-y', 'i2yz', 'i2y-z', 'i4-x', 'i4x', &
              'i3-x-y-z', 'i3-xyz', 'i3xy-z', 'i3x-yz', 'i3xyz', 'i3-xy-z',   &
              'i3x-y-z', 'i3-x-yz', 'i6z', 'i6-z', 'i3z', 'i3-z', 'i21-10',   &
              'i2210', 'i2010', 'i2110' /

  CHARACTER(LEN=40) :: sym_jones(64)
!
!  The jones symbol contains the crystal coordinates of the point  
!  x,y,z transformed by the symmetry. For cubic symetries the coordinates
!  coincide with the cartesian coordinates. For hexagonal symmetries 
!  they depend on the definition of the direct lattice vectors. 
!  We give first the  transformation using the definition of hexagonal axis 
!  of QE, than those using the ITA definition. For 2x and 2y symmetries 
!  we give three sets of transformed coordinates: cubic, hex QE and hex ITA.
!  NB: For QE  a_1 = a (1,0,0),          a_2 = a (-0.5,0.866,0)
!      for ITA a_1 = a (-0.5,-0.866,0)   a_2 = a (1,0,0)
!
   DATA sym_jones / ' x,y,z,         x,y,z,         x,y,z',      &
                    '-x,-y,z        -x,-y,z        -x,-y,z',     &
                    '-x,y,-z,        x-y,y,-z,     -x,y-x,-z',   &
                    ' x,-y,-z,       y-x,y,-z,      x,x-y,-z',   &
                    ' y,x,-z', '-y,-x,-z', ' y,-x,z',   &
                    '-y,x,z', ' z,-y,x', '-z,-y,-x',    &
                    ' z,y,-x', '-z,y,x', '-x,z,y',      &
                    '-x,-z,-y',' x,z,-y',' x,-z,y',      &
                    ' y,z,x','-y,z,-x',' y,-z,-x',       &
                    '-y,-z,x',' z,x,y',' z,-x,-y',       &
                    '-z,x,y','-z,x,-y',                &
                    '                x-y,x,z        x-y,x,z',            &
                    '                y,y-x,z        y,y-x,z',            &
                    '               -y,x-y,z       -y,x-y,z',            &
                    '                y-x,-x,z       y-x,-x,z',           &
                    '               -y,-x,-z,       y-x,y,-z',   &
                    '                x,x-y,-z,     -y,-x,-z',    &
                    '               -x,y-x,-z       y,x,-z',     &
                    '                y,x,-z,        x-y,-y,-z',  &  
                    '-x,-y,-z       -x,-y,-z,      -x,-y,-z',    &
                    ' x,y,-z         x,y,-z,        x,y,-z',     &
                    ' x,-y,z,        y-x,-y,z,      x,x-y,z',    &
                    '-x,y,z,         x-y,-y,z,     -x,y-x,z',    &
                    '-y,-x,z', ' y,x,z', '-y,x,-z',      &
                    ' y,-x,-z', '-z,y,-x', ' z,y,x',     &
                    '-z,-y,x', ' z,-y,-x', ' x,-z,-y',   &
                    ' x,z,y','-x,-z,y','-x,z,-y',       &
                    '-y,-z,-x',' y,-z,x','-y,z,x',      &
                    ' y,z,-x','-z,-x,-y','-z,x,y',      &
                    ' z,-x,-y',' z,-x,y',                &
                    '                y-x,-x,-z      y-x,-x,-z',  &
                    '               -y,x-y,-z      -y,x-y,-z',   &
                    '                y,y-x,-z       y,y-x,-z',   &
                    '                x-y,x,-z       x-y,x,-z',   &
                    '                y,x,z,         x-y,-y,z', &
                    '               -x,y-x,z,       y,x,z',    &
                    '                x,x-y,z       -y,-x,z',   &
                    '               -y,-x,z,        y-x,y,z' /

  LOGICAL :: cub_op(64)

  DATA cub_op / .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,   .TRUE.,   .TRUE.,  &     
                .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,   .TRUE.,   .TRUE.,  &
                .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,   .TRUE.,   .TRUE.,  &
                .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,   .TRUE.,   .TRUE.,  &
               .FALSE., .FALSE., .FALSE., .FALSE.,  .FALSE.,  .FALSE.,  &
               .FALSE., .FALSE.,  .TRUE.,  .TRUE.,   .TRUE.,   .TRUE.,  &
                .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,   .TRUE.,   .TRUE.,  &
                .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,   .TRUE.,   .TRUE.,  &
                .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,   .TRUE.,   .TRUE.,  &
                .TRUE.,  .TRUE., .FALSE., .FALSE.,  .FALSE.,  .FALSE.,  &
               .FALSE., .FALSE., .FALSE.,  .FALSE. /

  LOGICAL :: hex_op(64)

  DATA hex_op / .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .FALSE.,  .FALSE., &     
               .FALSE., .FALSE., .FALSE., .FALSE.,  .FALSE.,  .FALSE., &
               .FALSE., .FALSE., .FALSE., .FALSE.,  .FALSE.,  .FALSE., &
               .FALSE., .FALSE., .FALSE., .FALSE.,  .FALSE.,  .FALSE., &
                .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,   .TRUE.,   .TRUE., &
                .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,   .TRUE.,   .TRUE., &
               .FALSE., .FALSE., .FALSE., .FALSE.,  .FALSE.,  .FALSE., &
               .FALSE., .FALSE., .FALSE., .FALSE.,  .FALSE.,  .FALSE., &
               .FALSE., .FALSE., .FALSE., .FALSE.,  .FALSE.,  .FALSE., &
               .FALSE., .FALSE.,  .TRUE.,  .TRUE.,   .TRUE.,   .TRUE., &
                .TRUE.,  .TRUE.,  .TRUE.,  .TRUE. /

  PUBLIC convert_rap, find_aux_ind_two_groups, has_sigma_h, is_right_oriented,&
         color_rap, find_group_info_ext, find_irr_proj, sym_label,    &
         group_generators,  print_element_list, find_projection_type, &
         convert_rap_new, nsym_group, group_index_from_ext, is_subgroup, &
         find_double_product_table, set_group_desc, print_character_table, &
         print_compatibility_table, find_group_ext, &
         find_double_product_table_from_sym, set_sym_o3, set_sym_su2, &
         product_sym_su2, compute_classes, compute_classes_double, &
         print_kronecker_table, write_group_table,  &
         print_ptype_info, find_factor_system, sym_jones, transform_group, &
         hex_op, cub_op, point_group_bravais, find_group_tags,  &
         group_name_schoenflies, group_name_international

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
 
  SUBROUTINE convert_rap_new(n, list_in, list_out, group_ext_in,  &
                    group_ext_out, aux_ind, ptype_in, ptype_out, &
                    gauge_in, gauge_out)

  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: group_ext_in, group_ext_out, n, aux_ind
  INTEGER, INTENT(IN) :: list_in(n), ptype_in(3), ptype_out(3)
  REAL(DP) :: gauge_in(48), gauge_out(48)
  INTEGER, INTENT(OUT) :: list_out(n)

  INTEGER :: i, j, ndeg, group_in, group_out
  LOGICAL :: done(n)
  INTEGER :: rap_list(4)

  IF (SUM(ABS(gauge_in))< 1.D-8.AND.SUM(ABS(gauge_out))<1.D-8) THEN
!
!  standard case, no projective representations
!
     group_in = group_index_from_ext(group_ext_in)
     group_out = group_index_from_ext(group_ext_out)
     done=.FALSE.
     DO i=1,n
        IF (done(i)) CYCLE
        IF (list_in(i)<=0) THEN
            list_out(i)=list_in(i)
            CYCLE
        ENDIF
 
        IF (ptype_in(1)==-1) THEN
           CALL convert_one_rap_so(list_in(i), ndeg, rap_list, &
                                                     group_in, group_out )
        ELSE
           CALL convert_one_rap(list_in(i), ndeg, rap_list, group_in, &
                                group_out, aux_ind)
        ENDIF
        list_out(i) = rap_list(1)
        DO j=2, ndeg
           IF (list_in(i+j-1) /= list_in(i)) &
              CALL errore('convert_rap_new','Problem with degeneracy',1)
           list_out(i+j-1) = rap_list(j)
           done(i+j-1) = .TRUE.
        END DO
     END DO
  ELSE
!
!  projective representations 
!
     CALL convert_rap_proj(n, list_in, list_out, group_ext_in, &
                           group_ext_out, ptype_in, ptype_out, gauge_in, &
                           gauge_out)
  ENDIF

  RETURN
  END SUBROUTINE convert_rap_new

  SUBROUTINE convert_rap_proj(n, list_in, list_out, group_ext_in,    &
                                 group_ext_out, ptype_in, ptype_out, &
                                 gauge_in, gauge_out)
!
!   This routine converts a list of irreducible representations of 
!   group_ext_in into a list of irreducible representations of group_ext_out.
!   Note that in the input list each representation must appear as many
!   times as its dimension.
!   The routine supports also projective representations, as specified
!   by ptype_in. Note that only the decomposition in the representations
!   that belong to ptype_out is attempted and the two set of indeces must
!   be consistent.
!

  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE constants, ONLY : pi
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: list_in(n), group_ext_in, group_ext_out, &
                         ptype_in(3), ptype_out(3)
  INTEGER, INTENT(OUT) :: list_out(n)
  REAL(DP), INTENT(IN) :: gauge_in(48), gauge_out(48)

  INTEGER :: group_desc_in(48), group_desc_out(48), b_in_a(48)
  INTEGER :: nsym_in, nsym_out
  INTEGER :: group_in, group_out
  COMPLEX(DP) :: char_mat_proj_in(48,48), char_mat_proj_out(48,48), &
                 char_mat_sub(48,48), asum, phase(48), pha
  REAL(DP) :: arguments_in(48,48), arguments_out(48,48), gauge(48), arg
  CHARACTER(LEN=45) :: name_rap_in(48), name_rap_out(48)
  INTEGER :: counter, jrap, irap, isym, jsym, ideg, ndeg, rap, &
             nrap_proj_in, nrap_proj_out, sub_table(12,7), ptype(3), &
             prd(48,48), epos(48,48)
  INTEGER :: start, last, itables, ntables, irot
  LOGICAL :: done(n), switched_group

  group_in = group_index_from_ext(group_ext_in)
  group_out = group_index_from_ext(group_ext_out)

  CALL set_group_desc(group_desc_in, nsym_in, group_ext_in)
  CALL set_group_desc(group_desc_out, nsym_out, group_ext_out)

  b_in_a=0
  DO isym=1,nsym_out
     DO jsym=1,nsym_in
        IF (group_desc_in(jsym)==group_desc_out(isym)) b_in_a(isym)=jsym
     ENDDO
     IF (b_in_a(isym)==0) &
        CALL errore('convert_rap_proj','group_out not subgroup of group_in',1)
  ENDDO
!
!  Check that the input and output ptype are compatible
!
  CALL set_factors(group_in, ptype_in, arguments_in)
  IF (ptype_in(1)==-1) THEN
     CALL find_double_product_table(prd, epos, group_ext_in)
  ELSE
     CALL find_product_table(prd, group_ext_in)
     epos=1
  ENDIF
!
! pass to the gauge_in factors and take the subgroup
!   
  DO isym=1,nsym_out
     DO jsym=1,nsym_out
        arguments_out(isym,jsym) = arguments_in(b_in_a(isym),b_in_a(jsym)) &
                   + gauge_in(prd(b_in_a(isym), b_in_a(jsym))) &
                   - gauge_in(b_in_a(isym)) &
                   - gauge_in(b_in_a(jsym)) 
         IF (epos(b_in_a(isym),b_in_a(jsym))==-1) &
            arguments_out(isym,jsym)=arguments_out(isym,jsym) + pi
      END DO
   END DO

   CALL find_projection_type(group_out, group_ext_out, arguments_out, &
                              ptype, gauge, .FALSE.)
  !
  ! special case groups 1 {E} and 28 {E, I} double group cannot be
  ! recognized by the factor system of the point group.
  ! The representations are the same, but we use different names
  ! for the point group and double group so we set here ptype in this cases.
  !
  IF (((group_ext_in==1.OR.group_ext_in==28).OR.   &
                    (group_ext_out==1.OR.group_ext_out==28))) THEN
     IF (ptype_in(1)==-1) ptype(1)=-1
     IF (ptype_out(1)==-1) ptype(1)=-1
!
!  check that the expected ptype coincides with that found here
!
  ELSEIF (ptype(1) /= ptype_out(1) .OR. ptype(2) /= ptype_out(2) .OR. &
      ptype(3) /= ptype_out(3) ) THEN
      WRITE(stdout,'(5x,"Group in ",4i5)') group_ext_in, ptype_in(1:3)
      WRITE(stdout,'(5x,"ptype_out",5x,3i5)')  ptype(1:3)
      WRITE(stdout,'(5x,"Group out",4i5)') group_ext_out, ptype_out(1:3)
      CALL errore('convert_rap_proj','Decomposition not possible',1)
  ENDIF
!
!  If the code arrives here the two representations are compatible.
!  We can take the representations and bring them in the common factor
!  system
!

  CALL set_stand_irr_proj(group_ext_out, ptype_out, char_mat_proj_out, &
                      name_rap_out, nrap_proj_out, nsym_out)
  !
  !   apply the possible phase
  !
  DO isym=1,nsym_out
     arg=gauge_out(isym)
     pha=CMPLX(COS(arg),SIN(arg))
     DO irap=1,nrap_proj_out
        char_mat_proj_out(irap,isym) = char_mat_proj_out(irap,isym) / pha 
     END DO 
  END DO

  CALL set_stand_irr_proj(group_ext_in, ptype_in, char_mat_proj_in, &
                      name_rap_in, nrap_proj_in, nsym_in)
  !
  !   apply the possible phase
  !
  DO isym=1,nsym_in
     arg = gauge_in(isym)
     pha = CMPLX(COS(arg),SIN(arg))
     DO irap=1,nrap_proj_in
        char_mat_proj_in(irap,isym) = char_mat_proj_in(irap,isym) / pha
     END DO 
  END DO
!
!  take the subgroup
!
  DO irap=1, nrap_proj_in
     DO isym=1,nsym_out
        char_mat_sub(irap,isym)=char_mat_proj_in(irap, b_in_a(isym))
     ENDDO
  ENDDO
!
!  representation to decompose written in output
!
!  CALL write_group_char_mat(group_desc_out, nsym_out, char_mat_sub, &
!                                            name_rap_in, nrap_proj_in)
!
!  representations of group_out written in output
!
!  CALL write_group_char_mat(group_desc_out, nsym_out, char_mat_proj_out, &
!                                            name_rap_out, nrap_proj_out)
!
!  Do the actual decomposition
!
  DO irap=1, nrap_proj_in
     asum=0.0_DP
     DO isym=1,nsym_out
        asum = asum + char_mat_sub(irap,isym) * CONJG(char_mat_sub(irap,isym))
     ENDDO
     ndeg = NINT(DBLE(asum) / DBLE(nsym_out))  
!     write(6,*) ndeg,  asum / DBLE(nsym_out), ndeg - asum / DBLE(nsym_out)
     IF (ABS(ndeg - asum / DBLE(nsym_out)) > 1.D-6) &
        CALL errore('convert_rap_proj','problem with rap',irap)

     sub_table(irap,1) = NINT(DBLE(char_mat_sub(irap,1)))

     counter=1
     DO jrap=1, nrap_proj_out

        asum=(0.0_DP,0.0_DP)
        DO isym=1,nsym_out
           asum = asum + char_mat_sub(irap,isym) * &
                           CONJG(char_mat_proj_out(jrap,isym))
        ENDDO
        ndeg= NINT(DBLE(asum) / DBLE(nsym_out))
!        WRITE(6,*) irap, jrap, asum, ndeg
        IF (ABS(ndeg - asum / DBLE(nsym_out)) > 1.D-6) &
            CALL errore('convert_rap_proj','problem with ndeg',1)

        DO ideg =1, ndeg * NINT(DBLE(char_mat_proj_out(jrap,1)))
           counter=counter+1 
           sub_table(irap,counter) = jrap 
        END DO
     END DO
     IF (counter /= sub_table(irap,1)+1) &
        CALL errore('convert_rap_proj','problem with counter',1)
  END DO

  done = .FALSE.
  DO rap=1,n 
     IF (done(rap)) CYCLE
     irap=list_in(rap)
     IF (irap<=0) THEN
        list_out(rap)=list_in(rap)
     ELSE
       ndeg = sub_table(irap, 1)
       IF (ndeg==0) THEN
          WRITE(stdout,'(5x, "case not programmed group in",i5," group_out", i5, &
        " ptype in",3i5," ptype out", 3i5," rap",i5)') group_ext_in, &
                        group_ext_out, ptype_in, ptype_out, rap
          CALL errore('convert_one_rap_proj','problem representation not found',1)
       END IF
       DO ideg=1,ndeg
          IF (list_in(rap+ideg-1) /= irap) &
             CALL errore('convert_rap_proj','problem with representations',1)
          done(rap+ideg-1) = .TRUE.
          list_out(rap+ideg-1) = sub_table(irap, ideg+1)
       ENDDO
     ENDIF
  ENDDO

  RETURN
  END SUBROUTINE convert_rap_proj
 

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

        sub_table(12,2,1,1)=1
        sub_table(12,2,1,2)=1
        sub_table(12,2,2,1)=1
        sub_table(12,2,2,2)=2
        sub_table(12,2,3,1)=1
        sub_table(12,2,3,2)=4
        sub_table(12,2,4,1)=1
        sub_table(12,2,4,2)=3
        sub_table(12,2,5,1)=2
        sub_table(12,2,5,2)=3
        sub_table(12,2,5,3)=4
        sub_table(12,2,6,1)=2
        sub_table(12,2,6,2)=1
        sub_table(12,2,6,3)=2
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

        sub_table(12,4,1,1)=1
        sub_table(12,4,1,2)=1
        sub_table(12,4,2,1)=1
        sub_table(12,4,2,2)=4
        sub_table(12,4,3,1)=2
        sub_table(12,4,3,2)=1
        sub_table(12,4,3,3)=4
        sub_table(12,4,4,1)=3
        sub_table(12,4,4,2)=2
        sub_table(12,4,4,3)=3
        sub_table(12,4,4,4)=4
        sub_table(12,4,5,1)=3
        sub_table(12,4,5,2)=1
        sub_table(12,4,5,3)=2
        sub_table(12,4,5,4)=3
        sub_table(12,4,6,1)=1
        sub_table(12,4,6,2)=2
        sub_table(12,4,7,1)=1
        sub_table(12,4,7,2)=3
        sub_table(12,4,8,1)=2
        sub_table(12,4,8,2)=2
        sub_table(12,4,8,3)=3
        sub_table(12,4,9,1)=3
        sub_table(12,4,9,2)=1
        sub_table(12,4,9,3)=3
        sub_table(12,4,9,4)=4
        sub_table(12,4,10,1)=3
        sub_table(12,4,10,2)=1
        sub_table(12,4,10,3)=2
        sub_table(12,4,10,4)=4
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
  INTEGER :: ga, gb, group_a_ext, group_b_ext
  INTEGER :: group_desc_a(48), which_elem_a(48)
  INTEGER :: group_desc_b(48), which_elem_b(48)
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

  ga=group_a
  CALL find_group_info_ext(nsym_a, sr_a, ga, group_a_ext, &
                                                  which_elem_a, group_desc_a)
  gb=group_b
  CALL find_group_info_ext(nsym_b, sr_b, gb, group_b_ext, &
                                                  which_elem_b, group_desc_b)
  
  IF (.NOT. is_subgroup(group_a_ext, group_b_ext)) THEN
     aux_ind=-1
     RETURN
  ENDIF

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
              !
              aux_ind=0
              IF (group_desc_a(3)==group_desc_b(2)) aux_ind=1
              IF (group_desc_a(4)==group_desc_b(2)) aux_ind=2
              IF (aux_ind==0) CALL errore('find_aux_ind_two_groups',&
                               'C_2v and C_s have no common mirror', 1)
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

        aux_ind=0
        IF (group_desc_a(5)==group_desc_b(2).OR.  &
            group_desc_a(7)==group_desc_b(2)) aux_ind=1
        IF (group_desc_a(6)==group_desc_b(2) .OR. &
            group_desc_a(8)==group_desc_b(2)) aux_ind=2
        IF (aux_ind==0) CALL errore('find_aux_ind_two_groups',&
                               'C_4v and C_s have no common mirror', 1)

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
      IF (group_desc_a(5)==group_desc_b(3).OR.  &
          group_desc_a(7)==group_desc_b(3)) aux_ind=1
      IF (group_desc_a(6)==group_desc_b(3) .OR. &
          group_desc_a(8)==group_desc_b(3)) aux_ind=2
    
!              DO isym=1,nsym_b
!
!   find one of the mirrors
!
!                 IF (tipo_sym(sr_b(1,1,isym))==5) imirror=isym
!              ENDDO
!              CALL mirror_axis(sr_b(1,1,imirror),bx)
!              IF (is_axis(bx,1).OR.is_axis(bx,2).OR.is_axis(bx,3)) THEN
!                 aux_ind=1
!              ELSE
!                 aux_ind=2
!              ENDIF
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
                 aux_ind=2
              ELSEIF ( ibx==1 .OR. ibx==10 .OR. ibx==11 ) THEN
                 aux_ind=1
              ELSE
                 CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
              ENDIF
           CASE(4,5,7)
      !
      !    C_2, C_3, C_6
      !
              aux_ind=1
           CASE(12)
      !
      !    C_2v
      !
              CALL mirror_axis(sr_b(1,1,3),bx)
              CALL which_c2(bx, ibx)
              IF ( ibx==2 .OR. ibx==12 .OR. ibx==13 ) THEN
                 aux_ind=2
              ELSEIF ( ibx==1 .OR. ibx==10 .OR. ibx==11 ) THEN
                 aux_ind=1
              ELSE
                 CALL errore('find_aux_ind_two_groups',&
                                      ' Problem with axis direction ',1)
              ENDIF
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
                 aux_ind=2
              ELSEIF ( ibx==1 .OR. ibx==10 .OR. ibx==11 ) THEN
                 aux_ind=1
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
                   CALL versor(sr_b(1,1,isym),bx)
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
      !     C_2v  to choose between 1 2 3 or 4
      !
               IF (group_desc_b(2)==2.OR.group_desc_b(2)==3.OR.&
                   group_desc_b(2)==4) THEN

                  IF (group_desc_b(3)==34.OR.group_desc_b(3)==35.OR. &
                      group_desc_b(3)==36) THEN
                     aux_ind=1
                  ELSE
                     aux_ind=2
                  ENDIF
               ELSE
                  IF (group_desc_b(3)==34.OR.group_desc_b(3)==35.OR. &
                      group_desc_b(3)==36) THEN
                     aux_ind=3
                  ELSE
                     aux_ind=4
                  ENDIF
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
!  it corresponds to the standard ordering (that must match that used
!  for the projective representations).
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
! 38)   D_2  E  2z  2x  2y      1  2   4   3   = 30
! 39)   D_2  E  2z  2x-y 2xy    1  2   6   5   = 66
! 40)   D_2  E  2y  2x-z 2xz    1  3   10  9   = 191
! 41)   D_2  E  2x  2yz 2y-z    1  4   13  14  = 382
! 42)   D_2  E  2z  21-10 2110 1  2   29  32  = 1870
! 43)   D_2  E  2z  2010  2210 1  2   31  30  = 1866
! 44)   D_3  E  3z  3-z 2x 2110 2010   1  27 28 4 32 31 = 3515
! 45)   D_3  E  3z  3-z 2210  2y 21-10  1  27 28 30 3 29 = 3264
! 46)   D_3  E  3xyz   3-x-y-z 2y-z 2x-y 2x-z  1  21 17 14 6  10 = 1063
! 47)   D_3  E  3x-yz  3-xy-z 2yz 2x-z 2xy 1  20 22 13 10 5 = 1179 
! 48)   D_3  E  3-xyz  3x-y-z 2xz 2xy 2y-z   1  18 23 9 5  14 = 1156
! 49)   D_3  E  3-x-yz 3xy-z  2xz 2yz 2x-y 1  24 19 9  13 6 = 1224
! 50)   D_4  E  4z 2z  4-z 2x 2xy 2y 2x-y 1 8 2 7 4 5 3 6 = 204
! 51)   D_4  E  4y 2y  4-y 2z 2xz 2x 2x-z 1 11 3 12 2 9 4 10 = 476
! 52)   D_4  E  4x 2x  4-x 2yz 2z 2y-z 2y 1 15 4 16 13 2 14 3 =876
! 53)   D_6  E  6z 3z 2z 3-z 6-z 2x 2210 2110 2y 2010 21-10 
!                                   1 25 27 2 28 26 4 30 32 3 31 29 = 6570
! 54)   C_2v E  2x   i2y    i2z     1   4   35  34  = 2398
! 55)   C_2v E  2x   i2yz   i2y-z   1   4   45  46  = 4158
! 56)   C_2v E  2y   i2z    i2x     1   3   34  36  = 2462
! 57)   C_2v E  2y   i2x-z  i2xz    1   3   42  41  = 3455
! 58)   C_2v E  2z   i2x    i2y     1   2   36  35  = 2526
! 59)   C_2v E  2z   i2x-y  i2xy    1   2   38  37  = 2818
! 60)   C_2v E  2xy  i2z    i2x-y   1   5   34  38  = 2626
! 61)   C_2v E  2x-y i2xy   i2z     1   6   37  34  = 2562
! 62)   C_2v E  2xz  i2y    i2x-z   1   9   35  42  = 3071
! 63)   C_2v E  2x-z i2xz   i2y     1   10  41  35  = 3007
! 64)   C_2v E  2yz  i2y-z  i2x     1   13  46  36  = 3582
! 65)   C_2v E  2y-z i2x    i2yz    1   14  36  45  = 3518
! 66)   C_2v E  2z   i2010  i2210   1   2   63  62  = 7818
! 67)   C_2v E  2z   i21-10 i2110   1   2   61  64  = 7822
! 68)   C_2v E  2210  i2z   i2010   1   30  34  63  = 6026
! 69)   C_2v E  21-10 i2110 i2z     1   29  64  34  = 6094
! 70)   C_2v E  2110  i2z   i21-10  1   32  34  61  = 5902
! 71)   C_2v E  2010  i2210 i2z     1   31  62  34  = 5962
! 72)   C_3v E  3z 3-z i2x i2110 i2010 1 27 28 36 64 63 = 10875
! 73)   C_3v E  3z 3-z i2210 i2y i21-10  1 27 28 62 35 61 = 10304
! 74)   C_3v E  3xyz 3-x-y-z i2y-z i2x-y i2x-z 1 21 17 46 38 42 = 6055
! 75)   C_3v E  3x-yz 3-xy-z i2yz  i2x-z i2xy 1 20 22 45 42 37 = 6043
! 76)   C_3v E  3-xyz 3x-y-z i2xz i2xy  i2y-z  1 18 23 41 37 46 = 6020
! 77)   C_3v E  3-x-yz 3xy-z i2xz  i2yz i2x-y  1 24 19 41 45 38 = 6088
! 78)   C_4v E  4z 2z 4-z  i2x  i2xy  i2y i2x-y  1 8 2 7 36 37 35 38 = 5452
! 79)   C_4v E  4y 2y 4-y  i2z  i2xz  i2x i2x-z  1 11 3 12 34 41 36 42 = 6172
! 80)   C_4v E  4x 2x 4-x  i2yz  i2z i2y-z i2y   1 16 4 15 45 34 46 35 = 7020
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
! 96)   C_4h E  4z 2z 4-z i i4z i2z i4-z 1 8 2 7  33 40 34 39 = 5484
! 97)   C_4h E  4y 2y 4-y i i4y i2y i4-y 1 11 3 12 33 43 35 44 = 6374
! 98)   C_4h E  4x 2x 4-x i i4x i2x i4-x 1 16 4 15 33 48 36 47 = 7396
! 99)   C_6h E  6z 3z 2z 3-z 6-z i i6z i3z i2z i3-z i6-z 
!               1  25 27 2 28 26 33 57 59 34 60 58 = 18758 
! 100)  D_2h E  2z 2x   2y    i i2z i2x   i2y  1 2 4 3 33 34 36 35 = 4796
! 101)  D_2h E  2x 2yz  2y-z  i i2x i2yz  i2y-z  1 4 13 14 33 36 45 46 = 6908
! 102)  D_2h E  2y 2x-z 2xz   i i2y i2x-z i2xz  1 3 10 9 33 35 42 41 = 5950
! 103)  D_2h E  2z 2x-y 2xy   i i2z i2x-y i2xy  1 2 6 5 33 34 38 37 = 5124
! 104)  D_2h E  2z 21-10 2110 i i2z i21-10 i2110 1 2 29 32 33 34 61 64 = 11932
! 105)  D_2h E  2z 2010 2210  i i2z i2010 i2210  1 2 31 30 33 34 63 62 = 11924
! 106)  D_3h E  i6z 3z i2z 3-z i6-z 2x i2210 2110 i2y 2010 i21-10   
!                      1  57  27  34  28  58  4  62  32  35 31  61 = 20074
! 107)  D_3h E  i6z 3z i2z 3-z i6-z i2x 2210 i2119 2y i2010  21-10    
!                      1  57  27  34  28  58  36 30 64 3  63  29   = 20394 
! 108)  D_4h E  4z 2z 4-z 2x 2xy 2y 2x-y i i4z i2z i4-z i2x 
!            s_xy s_y s_x-y 1 8 2 7 4  5  3  6  33 40 34 39 36 37 35 38 = 10904
! 109)  D_4h E  4y 2y 4-y 2z 2xz 2x 2x-z  i i4y i2y i4-y i2z i2xz i2x i2x-z 
!                           1 11 3 12  2  9 4 10 33 43 35 44 34 41 36 42 =12472
! 110)  D_4h E 4x 2x 4-x 2yz 2z 2y-z 2y i i4x i2x i4-x i2yz i2z i2y-z i2y
!                          1 16 4 15  13 2 14 3 33 48 36 47 45 34 46 35 =  14296
! 111)  D_6h E 6z 3z 2z 3-z 6-z 2x 2210 2110 2y 2010 21-10
!            i i6z i3z i2z i3-z i6-z i2x  i2210  i2110  i2y i2010 i21-10
!            1 25 27 2 28 26 4 30 32 3 31 29 33 57 59 34 60 58 36 62 64 35
!            63 61 = 40660
! 112)  D_2d E i4z 2z i4-z 2x i2xy 2y i2x-y 1 40 2 39 4 37 3 38 = 5964
! 113)  D_2d E i4z 2z i4-z i2x 2xy i2y 2x-y 1 40 2 39 36 5 35 6  = 5708
! 114)  D_2d E i4y 2y i4-y 2z  i2xz  2x  i2x-z  1 43 3 44 2 41 4 42   = 7260
! 115)  D_2d E i4y 2y i4-y i2z 2xz i2x  2x-z 1 43 3 44 34 9 36 10 = 6428
! 116)  D_2d E i4x 2x i4-x i2yz 2z  i2y-z 2y 1 48 4 47 45 2 46 3   = 8684
! 117)  D_2d E i4x 2x i4-x 2yz i2z  2y-z i2y  1 48 4 47 13 34 14 35 = 7276
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
!            22 17 18 19 20 33 34 35 36 53 55 56 54 49 50 51 52 = 30292
! 134)  T_d  E 2x 2y 2z 3xyz 3x-y-z 3-xy-z 3-x-yz 3-x-yz 3-x-y-z 3-xyz
!            3x-yz 3xy-z i4-x i4x i4-y i4y i4-z i4z i2y-z i2yz i2x-z i2xz 
!            i2x-y i2xy
!            1 4 3 2 21 23 22 24 17 18 20 19  47 48 44 43 39 40 46 45 42 
!            41 38 37  = 25252
! 135)  O    E 2x 2y 2z 3xyz 3x-y-z 3-xy-z 3-x-yz 3-x-yz 3-x-y-z 3-xyz
!            3x-yz 3xy-z 4-x 4x 4-y 4y 4-z 4z 2y-z 2yz 2x-z 2xz 2x-y 2xy
!            1 4 3 2 21 23 22 24 17 18 20 19 15 16 12 11 7 8 14 13 10  
!            9 6 5  = 4900
! 136)  O_h  E 2x 2y 2z 3xyz 3x-y-z 3-xy-z 3-x-yz 3-x-yz 3-x-y-z 3-xyz
!            3x-yz 3xy-z 4-x 4x 4-y 4y 4-z 4z 2y-z 2yz 2x-z 2xz 2x-y 2xy
!            i i2x i2y i2z i3xyz i3x-y-z i3-xy-z i3-x-yz i3-x-yz i3-x-y-z i3-xyz
!            i3x-yz i3xy-z i4-x i4x i4-y i4y i4-z i4z i2y-z i2yz i2x-z 
!            i2xz i2x-y i2xy
!            1 4 3 2 21 23 22 24 17 18 20 19 15 16 12 11 7 8 14 13 10  
!            9 6 5 32 36 35 53 54 56 49 50 52 51 47 48 44 43 39 40 46 45 42
!            41 38 37 = 53576

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nsym
  REAL(DP), INTENT(IN) :: smat(3,3,nsym)
  INTEGER, INTENT(OUT) :: code_group, code_group_ext 
  INTEGER, INTENT(OUT) :: which_elem(nsym)
  INTEGER, INTENT(OUT) :: group_desc(48)


  INTEGER :: group_tags(48), group_tag

  INTEGER :: isym, jsym, i, nsym_

  CALL find_group_tags(nsym, smat, group_tags)

  CALL find_group_ext(group_tags, nsym, code_group_ext)

  code_group=group_index_from_ext(code_group_ext)  
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

  SUBROUTINE find_group_ext(group_tags, nsym, code_group_ext)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nsym
  INTEGER, INTENT(IN) :: group_tags(48)
  INTEGER, INTENT(OUT) :: code_group_ext

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

  INTEGER :: isym, i, group_tag

  group_tag=0
  DO isym=1,nsym
     group_tag = group_tag + group_tags(isym)**2
  ENDDO

  code_group_ext=0
  DO i=1,npg
     IF (group_tag==hash_tab(i)) code_group_ext=i
  ENDDO
!
!  groups 76 and 92 have the same hash value.  The routine here finds always 92.
!  We distinguish them from the number of operations.
!
  IF (code_group_ext==92.AND.nsym==6) code_group_ext=76

  IF (code_group_ext==0) &
     CALL errore('find_group_info_ext','input group unknown',1)

  RETURN
  END SUBROUTINE find_group_ext

  FUNCTION group_index_from_ext(code_group_ext)
!
!   This function converts the extended group code (1-136) into the standard
!   group code (1-32)
!
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
!
!   This routine contains the list of symmetry operations for each
!   extended group code (1-136) and sets them into group_desc
!
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
       group_desc(3)=4
       group_desc(4)=3
     CASE (39)
       group_desc(2)=2
       group_desc(3)=6
       group_desc(4)=5
     CASE (40)
       group_desc(2)=3
       group_desc(3)=10
       group_desc(4)=9
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
       group_desc(3)=31
       group_desc(4)=30
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
       group_desc(5)=13
       group_desc(6)=2
       group_desc(7)=14
       group_desc(8)=3
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
       group_desc(3)=42
       group_desc(4)=41
     CASE (58)
       group_desc(2)=2
       group_desc(3)=36
       group_desc(4)=35
     CASE (59)
       group_desc(2)=2
       group_desc(3)=38
       group_desc(4)=37
     CASE (60)
       group_desc(2)=5
       group_desc(3)=34
       group_desc(4)=38
     CASE (61)
       group_desc(2)=6
       group_desc(3)=37
       group_desc(4)=34
     CASE (62)
       group_desc(2)=9
       group_desc(3)=35
       group_desc(4)=42
     CASE (63)
       group_desc(2)=10
       group_desc(3)=41
       group_desc(4)=35
     CASE (64)
       group_desc(2)=13
       group_desc(3)=46
       group_desc(4)=36
     CASE (65)
       group_desc(2)=14
       group_desc(3)=36
       group_desc(4)=45
     CASE (66)
       group_desc(2)=2
       group_desc(3)=63
       group_desc(4)=62
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
       group_desc(3)=64
       group_desc(4)=34
     CASE (70)
       group_desc(2)=32
       group_desc(3)=34
       group_desc(4)=61
     CASE (71)
       group_desc(2)=31
       group_desc(3)=62
       group_desc(4)=34
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
       group_desc(5)=45
       group_desc(6)=34
       group_desc(7)=46
       group_desc(8)=35
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
       group_desc(3)=4
       group_desc(4)=3
       group_desc(5)=33
       group_desc(6)=34
       group_desc(7)=36
       group_desc(8)=35
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
       group_desc(3)=10
       group_desc(4)=9
       group_desc(5)=33
       group_desc(6)=35
       group_desc(7)=42
       group_desc(8)=41
     CASE (103)
       group_desc(2)=2
       group_desc(3)=6
       group_desc(4)=5
       group_desc(5)=33
       group_desc(6)=34
       group_desc(7)=38
       group_desc(8)=37
     CASE (104)
       group_desc(2)=2
       group_desc(3)=29
       group_desc(4)=32
       group_desc(5)=33
       group_desc(6)=34
       group_desc(7)=61
       group_desc(8)=64
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
       group_desc(2)=57
       group_desc(3)=27
       group_desc(4)=34
       group_desc(5)=28
       group_desc(6)=58
       group_desc(7)=36
       group_desc(8)=30
       group_desc(9)=64
       group_desc(10)=3
       group_desc(11)=63
       group_desc(12)=29
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
       group_desc(5)=13
       group_desc(6)=2
       group_desc(7)=14
       group_desc(8)=3
       group_desc(9)=33
       group_desc(10)=48
       group_desc(11)=36
       group_desc(12)=47
       group_desc(13)=45
       group_desc(14)=34
       group_desc(15)=46
       group_desc(16)=35
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
       group_desc(2)=40
       group_desc(3)=2
       group_desc(4)=39
       group_desc(5)=4
       group_desc(6)=37
       group_desc(7)=3
       group_desc(8)=38
     CASE (113)
       group_desc(2)=40
       group_desc(3)=2
       group_desc(4)=39
       group_desc(5)=36
       group_desc(6)=5
       group_desc(7)=35
       group_desc(8)=6
     CASE (114)
       group_desc(2)=43
       group_desc(3)=3
       group_desc(4)=44
       group_desc(5)=2
       group_desc(6)=41
       group_desc(7)=4
       group_desc(8)=42
     CASE (115)
       group_desc(2)=43
       group_desc(3)=3
       group_desc(4)=44
       group_desc(5)=34
       group_desc(6)=9
       group_desc(7)=36
       group_desc(8)=10
     CASE (116)
       group_desc(2)=48
       group_desc(3)=4
       group_desc(4)=47
       group_desc(5)=45
       group_desc(6)=2
       group_desc(7)=46
       group_desc(8)=3
     CASE (117)
       group_desc(2)=48
       group_desc(3)=4
       group_desc(4)=47
       group_desc(5)=13
       group_desc(6)=34
       group_desc(7)=14
       group_desc(8)=35
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
       group_desc(2)=4
       group_desc(3)=3
       group_desc(4)=2
       group_desc(5)=21
       group_desc(6)=23
       group_desc(7)=22
       group_desc(8)=24
       group_desc(9)=17
       group_desc(10)=18
       group_desc(11)=20
       group_desc(12)=19
       group_desc(13)=47
       group_desc(14)=48
       group_desc(15)=44
       group_desc(16)=43
       group_desc(17)=39
       group_desc(18)=40
       group_desc(19)=46
       group_desc(20)=45
       group_desc(21)=42
       group_desc(22)=41
       group_desc(23)=38
       group_desc(24)=37

     CASE (135,136)
!
!   O, O_h
!
       group_desc(2)=4
       group_desc(3)=3
       group_desc(4)=2
       group_desc(5)=21
       group_desc(6)=23
       group_desc(7)=22
       group_desc(8)=24
       group_desc(9)=17
       group_desc(10)=18
       group_desc(11)=20
       group_desc(12)=19
       group_desc(13)=15
       group_desc(14)=16
       group_desc(15)=12
       group_desc(16)=11
       group_desc(17)=7
       group_desc(18)=8
       group_desc(19)=14
       group_desc(20)=13
       group_desc(21)=10
       group_desc(22)=9
       group_desc(23)=6
       group_desc(24)=5
       IF (code_group_ext==136) THEN
          DO isym=1,24
             group_desc(isym+24) = group_desc(isym) + 32
          ENDDO
       ENDIF
     CASE DEFAULT
       CALL errore('find_group_info_ext','group index not found',1)
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

SUBROUTINE group_generators(code_group, row, column, n, linv, columni)
!
!  Since each point group has two generators A and B and is defined by
!  A^n1=E' and B^n2=E' and (AB)^2=E',
!  this routine gives, for each point group, the following information:
!  In column(i) the index in the list of the group operations (ordered
!  as specified in find_group_info_ext)
!      i         operations
!      1           A
!      2           B
!      3           AB
!  For cyclic groups B and AB are not set. 
!
!  In row(i) we put the index in the list of the group operations of
!  the element
!      i        operations
!      1          A^(-1)
!      2          B^(-1)
!      3          AB
!
!  In n(1) we put n1, in n(2) we put n2.
!
!  If the group contains inversion, linv is set to .true. and
!  columni(i) will contain the index in the list of the group operations
!      i       operation
!      1           I
!      2           IA
!      3           IB
!
!  This should be sufficient to determine the factor system and the type
!  of projective representation.
!  In order to apply possible phase factors the group elements of uniaxial
!  groups are all ordered in the form 
!  E A A^2 A^3 ... A^(n1-1) B  AB A^2B .... A^{n1-1}B ..
!  I IA IA^2 ...
!
!  Note that for improper point groups isomorphous to a proper point
!  group, A or B might be improper operations.
!
!  For T the elements order keeps elements in the same class close.
!  The order is 
!  
!  E  AB  BA   ABBA  A  B^2  BA^2  A^2B  A^2  B  AB^2  B^2A
!                     
!  For O the elements order keeps elements in the same class close. Note that
!  A does not coincide with that of T. This order for the cubic group
!  is taken from the book Kim (see below the complete reference).
!
!    E BA^2B-1  ABA-1B A^2 B A-1A-1B B-1ABA-1 AB-1A B-1 B-1A^2 A^2B-1 A-1BA-1
!    B-1A  A-1B  AB-1 BA-1 A-1 A A^2B-1A  AB A^2B-1AB  BA BA^2B-1A  BA-1B
!
!  T_d is isomorphic to O and ordered in the same way.
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: code_group
INTEGER, INTENT(OUT) :: row(3), column(3), columni(3), n(3)
LOGICAL, INTENT(OUT) :: linv


row=1
column=1
columni=1
n=1
n(3)=2
linv=.FALSE.

SELECT CASE (code_group)

   CASE(1)
      n(3)=1
   CASE(2)
      linv=.TRUE.
      columni(1) = 2 
      columni(2) = 2 
      columni(3) = 2 
   CASE(3,4)
!
!  C_s, C_2
!
      row(1)=2
      column(1)=2
      n(1)=2

   CASE(5)
!
!  C_3
!
      row(1)=3
      column(1)=2
      n(1)=3

   CASE(6)
!
!  C_4
!
      row(1)=4
      column(1)=2
      n(1)=4

   CASE(7)
!
!  C_6
!
      row(1)=6
      column(1)=2
      n(1)=6

   CASE(8)
!
!  D_2
!
      row(1)=2
      row(2)=3
      row(3)=4
      column(1)=2
      column(2)=3
      column(3)=4
      n(1)=2
      n(2)=2

   CASE(9)
!
!  D_3
!
      row(1)=3
      row(2)=4
      row(3)=5
      column(1)=2
      column(2)=4
      column(3)=5
      n(1)=3
      n(2)=2

   CASE(10)
!
!  D_4
!
      row(1)=4
      row(2)=5
      row(3)=6
      column(1)=2
      column(2)=5
      column(3)=6
      n(1)=4
      n(2)=2

   CASE(11)
!
!  D_6
!
      row(1)=6
      row(2)=7
      row(3)=8
      column(1)=2
      column(2)=7
      column(3)=8
      n(1)=6
      n(2)=2

   CASE(12)
!
!  C_2v
!
      row(1)=2
      row(2)=3
      row(3)=4
      column(1)=2
      column(2)=3
      column(3)=4
      n(1)=2
      n(2)=2

   CASE(13)
!
!  C_3v
!
      row(1)=3
      row(2)=4
      row(3)=5
      column(1)=2
      column(2)=4
      column(3)=5
      n(1)=3
      n(2)=2

   CASE(14)
!
!  C_4v
!
      row(1)=4
      row(2)=5
      row(3)=6
      column(1)=2
      column(2)=5
      column(3)=6
      n(1)=4
      n(2)=2

   CASE(15)
!
!  C_6v
!
      row(1)=6
      row(2)=7
      row(3)=8
      column(1)=2
      column(2)=7
      column(3)=8
      n(1)=6
      n(2)=2

   CASE(16)
!
!  C_2h
!
      row(1)=2
      column(1)=2
      n(1)=2

      linv=.TRUE. 
      columni(1)=3
      columni(2)=4
      columni(3)=3

   CASE(17)
!
!  C_3h
!
      row(1)=6
      column(1)=2
      n(1)=6

   CASE(18)
!
!  C_4h
!
      row(1)=4
      column(1)=2
      n(1)=4

      linv=.TRUE. 
      columni(1)=5
      columni(2)=6
      columni(3)=5


   CASE(19)
!
!  C_6h
!
      row(1)=6
      column(1)=2
      n(1)=6
!
      linv=.TRUE. 
      columni(1)=7
      columni(2)=8
      columni(3)=7
   CASE(20)
!
!  D_2h
!
      row(1)=2
      row(2)=3
      row(3)=4
      column(1)=2
      column(2)=3
      column(3)=4
      n(1)=2
      n(2)=2

      linv=.TRUE. 
      columni(1)=5
      columni(2)=6
      columni(3)=7

   CASE(21)
!
!  D_3h
!
      row(1)=6
      row(2)=7
      row(3)=8
      column(1)=2
      column(2)=7
      column(3)=8
      n(1)=6
      n(2)=2

   CASE(22)
!
!  D_4h
!
      row(1)=4
      row(2)=5
      row(3)=6
      column(1)=2
      column(2)=5
      column(3)=6

      n(1)=4
      n(2)=2

      linv=.TRUE. 
      columni(1)=9
      columni(2)=10
      columni(3)=13

   CASE(23)
!
!  D_6h
!
      row(1)=6
      row(2)=7
      row(3)=8
      column(1)=2
      column(2)=7
      column(3)=8

      n(1)=6
      n(2)=2

      linv=.TRUE. 
      columni(1)=13
      columni(2)=14
      columni(3)=19

   CASE(24)
!
!  D_2d
!
      row(1)=4
      row(2)=5
      row(3)=6
      column(1)=2
      column(2)=5
      column(3)=6
      n(1)=4
      n(2)=2

   CASE(25)
!
!  D_3d
!
      row(1)=3
      row(2)=4
      row(3)=5
      column(1)=2
      column(2)=4
      column(3)=5
      n(1)=3
      n(2)=2

      linv=.TRUE. 
      columni(1)=7
      columni(2)=8
      columni(3)=10

   CASE(26)
!
!  S_4
!
      row(1)=4
      column(1)=2
      n(1)=4

   CASE(27)
!
!  S_6
!
      row(1)=3
      column(1)=2
      n(1)=3

      linv=.TRUE. 
      columni(1)=4
      columni(2)=5
      columni(3)=4

   CASE(28)
!
!  T
!
      row(1)=9
      row(2)=6
      row(3)=2
      column(1)=5
      column(2)=10
      column(3)=2
      n(1)=3
      n(2)=3

   CASE(29)
!
!  T_h
!
      row(1)=9
      row(2)=6
      row(3)=2
      column(1)=5
      column(2)=10
      column(3)=2
      n(1)=3
      n(2)=3

      linv=.TRUE. 
      columni(1)=13
      columni(2)=17
      columni(3)=22

   CASE(30)
!
!  T_d
!
      row(1)=17
      row(2)=9
      row(3)=20
      column(1)=18
      column(2)=5
      column(3)=20
      n(1)=4
      n(2)=3

   CASE(31)
!
!  O
!
      row(1)=17
      row(2)=9
      row(3)=20
      column(1)=18
      column(2)=5
      column(3)=20
      n(1)=4
      n(2)=3

   CASE(32)
!
!  O_h
!
      row(1)=17
      row(2)=9
      row(3)=20
      column(1)=18
      column(2)=5
      column(3)=20
      n(1)=4
      n(2)=3

      linv=.TRUE. 
      columni(1)=25
      columni(2)=42
      columni(3)=29

CASE DEFAULT
      CALL errore('group_generators','point group not available',1)

END SELECT

RETURN
END SUBROUTINE group_generators

SUBROUTINE set_stand_irr_proj(cge, ptype, char_mat_proj, name_rap, nrap_proj, &
                                                                   nsym_proj)
!
!  This subroutine receives in input the extended group code and
!  sets the standard projective irreducible representations for that group.
!  One out of three types of representations are set depending on the input
!  ptype. The representations of the point group, the representations of the
!  double point group, or the projective representations corresponding to the
!  beta and gamma given by ptype.
!  The routine does not use the class and set the character for all elements
!  of the group. The order of the elements of the group is assumed to be that
!  written in the comments of the routine find_group_info_ext.
!
!  The character tables given in this routine have been derived 
!  following the ideas reported in the book:
!
!  Shoon K. Kim, Group Theoretical Methods and applications to molecules and
!             crystals, Cambridge University Press (1999).
!
!  The name of the representations follow the standard conventions, as
!  reported in the point group manual of thermo_pw, or the names introduced 
!  in the book.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: cge    ! code group extended
INTEGER, INTENT(IN) :: ptype(3)  ! type of the representation requested

INTEGER, INTENT(OUT) :: nrap_proj ! the number of representation
INTEGER, INTENT(OUT) :: nsym_proj ! the number of symmetry elements

CHARACTER(LEN=45), INTENT(OUT) :: name_rap(48)  
                                   ! Output: name of the representations

COMPLEX(DP) :: char_mat_proj(48,48) ! Output: character matrix of projective 
                                    !         rap (class not used)

REAL(DP), PARAMETER :: sqrt3=SQRT(3.0_DP), sqrt2=SQRT(2.0_DP)
COMPLEX(DP) :: w, w1
INTEGER :: group_code

group_code=group_index_from_ext(cge)
char_mat_proj=(0.0_DP,0.0_DP)
nrap_proj=0
nsym_proj=0

SELECT CASE (cge)
   CASE (1)
      nsym_proj=1     
      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
         nrap_proj=1
         name_rap(1)='A'
         char_mat_proj(1,1)=(1.0_DP,0.0_DP)
      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
         nrap_proj=1
         name_rap(1)='G_2'
         char_mat_proj(1,1)=(1.0_DP,0.0_DP)
      ENDIF
   CASE (2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,&
                                                          27)

      nsym_proj=2
      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=2
         name_rap(1)='A'''
         IF (cge<15) name_rap(1)='A'
         IF (cge==28) name_rap(1)='A_g'
         char_mat_proj(1,1:2)=(1.0_DP,0.0_DP)

         name_rap(2)='A'''''
         IF (cge<15) name_rap(2)='B'
         IF (cge==28) name_rap(2)='A_u'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)=(-1.0_DP,0.0_DP)
      ELSEIF(ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations
!
         nrap_proj=2

         name_rap(1)='G_3     M-12'
         IF (cge==28) name_rap(1)='G_2+    M-12'
         char_mat_proj(1,1)=(1.0_DP,0.0_DP)
         char_mat_proj(1,2)=(0.0_DP,1.0_DP)

         name_rap(2)='G_4     M12'
         IF (cge==28) name_rap(2)='G_2-    M12'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)=(0.0_DP,-1.0_DP)

      ENDIF
   CASE (28)
      nsym_proj=2
      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=2
         name_rap(1)='A_g'
         char_mat_proj(1,1:2)=(1.0_DP,0.0_DP)

         name_rap(2)='A_u'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)=(-1.0_DP,0.0_DP)

      ELSEIF(ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations, identical to those of the
!   point group, but we keep the different name
!
         nrap_proj=2

         name_rap(1)='G_2+'
         char_mat_proj(1,1)=(1.0_DP,0.0_DP)
         char_mat_proj(1,2)=(1.0_DP,0.0_DP)

         name_rap(2)='G_2-'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)=(-1.0_DP,0.0_DP)

      ENDIF

   CASE (29,30,31,32,33)
!
!   C_3
!

      nsym_proj=3
      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=3

         w=CMPLX(-0.5_DP,-sqrt3*0.5_DP)
         name_rap(1)='A'
         char_mat_proj(1,1:3)=(1.0_DP,0.0_DP)

         name_rap(2)='E       M2'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)= w **2
         char_mat_proj(2,3)= w 
     
         name_rap(3)='E*      M1'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)= w 
         char_mat_proj(3,3)= w **2
      ELSEIF(ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations
!
         w=CMPLX(0.5_DP,sqrt3*0.5_DP)
         nrap_proj=3

         name_rap(1)='G_4      M-12'
         char_mat_proj(1,1)=(1.0_DP,0.0_DP)
         char_mat_proj(1,2)= w 
         char_mat_proj(1,3)= CONJG(w)

         name_rap(2)='G_5      M12'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)= CONJG(w)
         char_mat_proj(2,3)= w

         name_rap(3)='G_6      M32'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)=(-1.0_DP,0.0_DP) 
         char_mat_proj(3,3)=(-1.0_DP,0.0_DP)
      ENDIF

   CASE (34,35,36,124,125,126)
!
!  C_4 or S_4
!

      w=CMPLX(sqrt2*0.5_DP,-sqrt2*0.5_DP)

      nsym_proj=4
      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=4

         name_rap(1)='A'
         char_mat_proj(1,1:4)=(1.0_DP,0.0_DP)

         name_rap(2)='B       M_2'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,3)=(1.0_DP,0.0_DP)
         char_mat_proj(2,4)=(-1.0_DP,0.0_DP)

         name_rap(3)='E       M-1'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)=(0.0_DP,1.0_DP)
         char_mat_proj(3,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,4)=(0.0_DP,-1.0_DP)

         name_rap(4)='E*      M1'
         char_mat_proj(4,1)=(1.0_DP,0.0_DP)
         char_mat_proj(4,2)=(0.0_DP,-1.0_DP)
         char_mat_proj(4,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,4)=(0.0_DP,1.0_DP)
     
      ELSEIF(ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations
!
         nrap_proj=4

         name_rap(1)='G_5     M-12'
         char_mat_proj(1,1)=(1.0_DP,0.0_DP)
         char_mat_proj(1,2)= CONJG(w)
         char_mat_proj(1,3)=(0.0_DP,1.0_DP)
         char_mat_proj(1,4)= w

         name_rap(2)='G_6     M_12'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)= w
         char_mat_proj(2,3)= w**2
         char_mat_proj(2,4)= CONJG(w)

         name_rap(3)='G_7     M_32'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)= w**3
         char_mat_proj(3,3)=(0.0_DP,1.0_DP)
         char_mat_proj(3,4)= CONJG(w**3)

         name_rap(4)='G_8     M-32'
         char_mat_proj(4,1)=(1.0_DP,0.0_DP)
         char_mat_proj(4,2)= -w
         char_mat_proj(4,3)=(0.0_DP,-1.0_DP)
         char_mat_proj(4,4)= -CONJG(w)

      ENDIF

   CASE (37,95)
!
!  C_6 or C_3h 
!

      nsym_proj=6
      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         w=CMPLX(0.5_DP,-sqrt3*0.5_DP)

         nrap_proj=6

         name_rap(1)='A'
         char_mat_proj(1,1:6)=(1.0_DP,0.0_DP)

         name_rap(2)='B       M_3'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,3)=(1.0_DP,0.0_DP)
         char_mat_proj(2,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,5)=(1.0_DP,0.0_DP)
         char_mat_proj(2,6)=(-1.0_DP,0.0_DP)

         w1=CONJG(w)
         name_rap(3)='E_1     M-1'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)= w1
         char_mat_proj(3,3)= w1**2
         char_mat_proj(3,4)= w1**3
         char_mat_proj(3,5)= CONJG(w1**2)
         char_mat_proj(3,6)= CONJG(w1)

         name_rap(4)='E_1*    M_1'
         char_mat_proj(4,1)=(1.0_DP,0.0_DP)
         char_mat_proj(4,2)= w
         char_mat_proj(4,3)= w**2
         char_mat_proj(4,4)= w**3
         char_mat_proj(4,5)= CONJG(w**2)
         char_mat_proj(4,6)= CONJG(w)

         w1=CONJG(w**2)
         name_rap(5)='E_2     M-2'
         char_mat_proj(5,1)=(1.0_DP,0.0_DP)
         char_mat_proj(5,2)= w1
         char_mat_proj(5,3)= w1**2
         char_mat_proj(5,4)= w1**3
         char_mat_proj(5,5)= CONJG(w1**2)
         char_mat_proj(5,6)= CONJG(w1)

         w1=w**2
         name_rap(6)='E_2*    M_2'
         char_mat_proj(6,1)=(1.0_DP,0.0_DP)
         char_mat_proj(6,2)= w1
         char_mat_proj(6,3)= w1**2
         char_mat_proj(6,4)= w1**3
         char_mat_proj(6,5)= CONJG(w1**2)
         char_mat_proj(6,6)= CONJG(w1)

      ELSEIF(ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!  double point group representations
!
         w=CMPLX(sqrt3*0.5_DP,-0.5_DP)

         nrap_proj=6

         w1=CONJG(w)
         name_rap(1)='G_7     M-12'
         char_mat_proj(1,1) = (1.0_DP,0.0_DP)
         char_mat_proj(1,2) = w1
         char_mat_proj(1,3) = w1**2
         char_mat_proj(1,4) = w1**3
         char_mat_proj(1,5) = CONJG(w1**2)
         char_mat_proj(1,6) = CONJG(w1)

         name_rap(2)='G_8     M_12'
         char_mat_proj(2,1) = (1.0_DP,0.0_DP)
         char_mat_proj(2,2) = w
         char_mat_proj(2,3) = w**2
         char_mat_proj(2,4) = w**3
         char_mat_proj(2,5) = CONJG(w**2)
         char_mat_proj(2,6) = CONJG(w)

         w1=w**5
         name_rap(3)='G_9     M_52'
         char_mat_proj(3,1) = (1.0_DP,0.0_DP)
         char_mat_proj(3,2) = w1
         char_mat_proj(3,3) = w1**2
         char_mat_proj(3,4) = w1**3
         char_mat_proj(3,5) = CONJG(w1**2)
         char_mat_proj(3,6) = CONJG(w1)

         w1=CONJG(w**5)
         name_rap(4)='G_10     M-52'
         char_mat_proj(4,1) = (1.0_DP,0.0_DP)
         char_mat_proj(4,2) = w1
         char_mat_proj(4,3) = w1**2
         char_mat_proj(4,4) = w1**3
         char_mat_proj(4,5) = CONJG(w1**2)
         char_mat_proj(4,6) = CONJG(w1)


         w1=CONJG(w**3)
         name_rap(5)='G_11    M-32'
         char_mat_proj(5,1) = (1.0_DP,0.0_DP)
         char_mat_proj(5,2) = w1
         char_mat_proj(5,3) = w1**2
         char_mat_proj(5,4) = w1**3
         char_mat_proj(5,5) = CONJG(w1**2)
         char_mat_proj(5,6) = CONJG(w1)
     
         w1=w**3
         name_rap(6)='G_12     M_32'
         char_mat_proj(6,1) = (1.0_DP,0.0_DP)
         char_mat_proj(6,2) = w1
         char_mat_proj(6,3) = w1**2
         char_mat_proj(6,4) = w1**3
         char_mat_proj(6,5) = CONJG(w1**2)
         char_mat_proj(6,6) = CONJG(w1)

      ENDIF
     
   CASE (38,39,40,41,42,43,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71)
!
!  D_2  or C_2v
!
      nsym_proj=4

      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=4

         name_rap(1)='A       A_1'
         IF (group_code==12) name_rap(1)='A_1'
         char_mat_proj(1,1:4)=(1.0_DP,0.0_DP)

         name_rap(2)='B_1     A_2'
         IF (group_code==12) name_rap(2)='A_2'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)=(1.0_DP,0.0_DP)
         char_mat_proj(2,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,4)=(-1.0_DP,0.0_DP)

         name_rap(3)='B_2     B_1'
         IF (group_code==12) name_rap(3)='B_1'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,3)=(1.0_DP,0.0_DP)
         char_mat_proj(3,4)=(-1.0_DP,0.0_DP)

         name_rap(4)='B_3     B_2'
         IF (group_code==12) name_rap(4)='B_2'
         char_mat_proj(4,1)=(1.0_DP,0.0_DP)
         char_mat_proj(4,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,4)=(1.0_DP,0.0_DP)

     ELSEIF(ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations
!
         nrap_proj=1

         name_rap(1)='G_5     E_12'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
      ENDIF

   CASE (44,45,46,47,48,49,72,73,74,75,76,77)
!
!  D_3 or C_3v 
!
      nsym_proj=6

      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=3
         name_rap(1)='A_1'
         char_mat_proj(1,1:6)=(1.0_DP,0.0_DP)

         name_rap(2)='A_2'
         char_mat_proj(2,1:3)=(1.0_DP,0.0_DP)
         char_mat_proj(2,4:6)=(-1.0_DP,0.0_DP)

         name_rap(3)='E       E_1'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,2:3)=(-1.0_DP,0.0_DP)

      ELSEIF(ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations
!
         nrap_proj=3

         name_rap(1)='G_4     E_12'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,2:3)=(1.0_DP,0.0_DP)

         name_rap(2)='G_5     B2'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2:3)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,4)=(0.0_DP,-1.0_DP)
         char_mat_proj(2,5:6)=(0.0_DP,1.0_DP)

         name_rap(3)='G_6     B1'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2:3)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,4)=(0.0_DP,1.0_DP)
         char_mat_proj(3,5:6)=(0.0_DP,-1.0_DP)

      ENDIF

   CASE (50,51,52,78,79,80,112,113,114,115,116,117)
!
!  D_4 or C_4v or D_2d
!
      nsym_proj=8

      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=5
         name_rap(1)='A_1'
         char_mat_proj(1,1:8)=(1.0_DP,0.0_DP)

         name_rap(2)='A_2'
         char_mat_proj(2,1:4)=(1.0_DP,0.0_DP)
         char_mat_proj(2,5:8)=(-1.0_DP,0.0_DP)

!
!  for the point group 52, 80, 116 and 117 the two representations are
!  reversed. Actually we stick with the standard definition that \sigma_v
!  are the mirror perpendicular to x, y, z and \sigma_d those perpendicular
!  to x=y, x=-y etc. and similar definition for C_2' and C_2''.
!
         IF (cge==52.OR.cge==80.OR.cge==116.OR.cge==117) THEN
            name_rap(3)='B_1'
            char_mat_proj(3,1)=(1.0_DP,0.0_DP)
            char_mat_proj(3,2)=(-1.0_DP,0.0_DP)
            char_mat_proj(3,3)=(1.0_DP,0.0_DP)
            char_mat_proj(3,4)=(-1.0_DP,0.0_DP)
            char_mat_proj(3,5)=(-1.0_DP,0.0_DP)
            char_mat_proj(3,6)=(1.0_DP,0.0_DP)
            char_mat_proj(3,7)=(-1.0_DP,0.0_DP)
            char_mat_proj(3,8)=(1.0_DP,0.0_DP)

            name_rap(4)='B_2'
            char_mat_proj(4,1)=(1.0_DP,0.0_DP)
            char_mat_proj(4,2)=(-1.0_DP,0.0_DP)
            char_mat_proj(4,3)=(1.0_DP,0.0_DP)
            char_mat_proj(4,4)=(-1.0_DP,0.0_DP)
            char_mat_proj(4,5)=(1.0_DP,0.0_DP)
            char_mat_proj(4,6)=(-1.0_DP,0.0_DP)
            char_mat_proj(4,7)=(1.0_DP,0.0_DP)
            char_mat_proj(4,8)=(-1.0_DP,0.0_DP)
         ELSE
            name_rap(3)='B_1'
            char_mat_proj(3,1)=(1.0_DP,0.0_DP)
            char_mat_proj(3,2)=(-1.0_DP,0.0_DP)
            char_mat_proj(3,3)=(1.0_DP,0.0_DP)
            char_mat_proj(3,4)=(-1.0_DP,0.0_DP)
            char_mat_proj(3,5)=(1.0_DP,0.0_DP)
            char_mat_proj(3,6)=(-1.0_DP,0.0_DP)
            char_mat_proj(3,7)=(1.0_DP,0.0_DP)
            char_mat_proj(3,8)=(-1.0_DP,0.0_DP)

            name_rap(4)='B_2'
            char_mat_proj(4,1)=(1.0_DP,0.0_DP)
            char_mat_proj(4,2)=(-1.0_DP,0.0_DP)
            char_mat_proj(4,3)=(1.0_DP,0.0_DP)
            char_mat_proj(4,4)=(-1.0_DP,0.0_DP)
            char_mat_proj(4,5)=(-1.0_DP,0.0_DP)
            char_mat_proj(4,6)=(1.0_DP,0.0_DP)
            char_mat_proj(4,7)=(-1.0_DP,0.0_DP)
            char_mat_proj(4,8)=(1.0_DP,0.0_DP)
         ENDIF

         name_rap(5)='E'
         char_mat_proj(5,1)=(2.0_DP,0.0_DP)
         char_mat_proj(5,3)=(-2.0_DP,0.0_DP)
      ELSEIF(ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations
!
         nrap_proj=2

         name_rap(1)='G_6     E_12'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,2)=CMPLX(sqrt2,0.0_DP)
         char_mat_proj(1,4)=CMPLX(sqrt2,0.0_DP)

         name_rap(2)='G_7     E_32'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,2)=CMPLX(-sqrt2,0.0_DP)
         char_mat_proj(2,4)=CMPLX(-sqrt2,0.0_DP)
      ENDIF

   CASE (53,81)
!
!   D_6 or D_3h or C_6v
!
      nsym_proj=12

      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=6
         name_rap(1)='A_1'
         char_mat_proj(1,1:12)=(1.0_DP,0.0_DP)

         name_rap(2)='A_2'
         char_mat_proj(2,1:6)=(1.0_DP,0.0_DP)
         char_mat_proj(2,7:12)=(-1.0_DP,0.0_DP)

         name_rap(3)='B_1'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,3)=(1.0_DP,0.0_DP)
         char_mat_proj(3,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,5)=(1.0_DP,0.0_DP)
         char_mat_proj(3,6)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,7)=(1.0_DP,0.0_DP)
         char_mat_proj(3,8)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,9)=(1.0_DP,0.0_DP)
         char_mat_proj(3,10)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,11)=(1.0_DP,0.0_DP)
         char_mat_proj(3,12)=(-1.0_DP,0.0_DP)

         name_rap(4)='B_2'
         char_mat_proj(4,1)=(1.0_DP,0.0_DP)
         char_mat_proj(4,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,3)=(1.0_DP,0.0_DP)
         char_mat_proj(4,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,5)=(1.0_DP,0.0_DP)
         char_mat_proj(4,6)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,7)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,8)=(1.0_DP,0.0_DP)
         char_mat_proj(4,9)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,10)=(1.0_DP,0.0_DP)
         char_mat_proj(4,11)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,12)=(1.0_DP,0.0_DP)


         name_rap(5)='E_1'
         char_mat_proj(5,1)=(2.0_DP,0.0_DP)
         char_mat_proj(5,2)=(1.0_DP,0.0_DP)
         char_mat_proj(5,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,4)=(-2.0_DP,0.0_DP)
         char_mat_proj(5,5)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,6)=(1.0_DP,0.0_DP)

         name_rap(6)='E_2'
         char_mat_proj(6,1)=(2.0_DP,0.0_DP)
         char_mat_proj(6,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,4)=(2.0_DP,0.0_DP)
         char_mat_proj(6,5)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,6)=(-1.0_DP,0.0_DP)
      ELSEIF(ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations
!
         nrap_proj=3

         name_rap(1)='G_7     E_12'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,2)=(sqrt3,0.0_DP)
         char_mat_proj(1,3)=(1.0_DP,0.0_DP)
         char_mat_proj(1,4)=(0.0_DP,0.0_DP)
         char_mat_proj(1,5)=(1.0_DP,0.0_DP)
         char_mat_proj(1,6)=(sqrt3,0.0_DP)

         name_rap(2)='G_8     E_52'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,2)=CMPLX(-sqrt3,0.0_DP)
         char_mat_proj(2,3)=(1.0_DP,0.0_DP)
         char_mat_proj(2,4)=(0.0_DP,0.0_DP)
         char_mat_proj(2,5)=(1.0_DP,0.0_DP)
         char_mat_proj(2,6)=CMPLX(-sqrt3,0.0_DP)

         name_rap(3)='G_9     E_32'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,2)=(0.0_DP,0.0_DP)
         char_mat_proj(3,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,4)=(0.0_DP,0.0_DP)
         char_mat_proj(3,5)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,6)=(0.0_DP,0.0_DP)

      ENDIF

   CASE (106,107)
!
!   D_3h  temporaly divided from D_6 for compatibility with
!   previous character tables
!
      nsym_proj=12

      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=6

         name_rap(1)='A_1'''
         char_mat_proj(1,1:12)=(1.0_DP,0.0_DP)

         name_rap(2)='A_2'''
         char_mat_proj(2,1:6)=(1.0_DP,0.0_DP)
         char_mat_proj(2,7:12)=(-1.0_DP,0.0_DP)

         name_rap(3)='E'''
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,4)=(2.0_DP,0.0_DP)
         char_mat_proj(3,5)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,6)=(-1.0_DP,0.0_DP)

!
         name_rap(4)='A_1'''''
         char_mat_proj(4,1)=(1.0_DP,0.0_DP)
         char_mat_proj(4,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,3)=(1.0_DP,0.0_DP)
         char_mat_proj(4,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,5)=(1.0_DP,0.0_DP)
         char_mat_proj(4,6)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,7)=(1.0_DP,0.0_DP)
         char_mat_proj(4,8)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,9)=(1.0_DP,0.0_DP)
         char_mat_proj(4,10)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,11)=(1.0_DP,0.0_DP)
         char_mat_proj(4,12)=(-1.0_DP,0.0_DP)
!
         name_rap(5)='A_2'''''
         char_mat_proj(5,1)=(1.0_DP,0.0_DP)
         char_mat_proj(5,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,3)=(1.0_DP,0.0_DP)
         char_mat_proj(5,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,5)=(1.0_DP,0.0_DP)
         char_mat_proj(5,6)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,7)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,8)=(1.0_DP,0.0_DP)
         char_mat_proj(5,9)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,10)=(1.0_DP,0.0_DP)
         char_mat_proj(5,11)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,12)=(1.0_DP,0.0_DP)

!
         name_rap(6)='E'''''
         char_mat_proj(6,1)=(2.0_DP,0.0_DP)
         char_mat_proj(6,2)=(1.0_DP,0.0_DP)
         char_mat_proj(6,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,4)=(-2.0_DP,0.0_DP)
         char_mat_proj(6,5)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,6)=(1.0_DP,0.0_DP)
      ELSEIF(ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN

!   double point group representations
!
         nrap_proj=3

         name_rap(1)='G_7     E_12'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,2)=(sqrt3,0.0_DP)
         char_mat_proj(1,3)=(1.0_DP,0.0_DP)
         char_mat_proj(1,4)=(0.0_DP,0.0_DP)
         char_mat_proj(1,5)=(1.0_DP,0.0_DP)
         char_mat_proj(1,6)=(sqrt3,0.0_DP)

         name_rap(2)='G_8     E_52'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,2)=CMPLX(-sqrt3,0.0_DP)
         char_mat_proj(2,3)=(1.0_DP,0.0_DP)
         char_mat_proj(2,4)=(0.0_DP,0.0_DP)
         char_mat_proj(2,5)=(1.0_DP,0.0_DP)
         char_mat_proj(2,6)=CMPLX(-sqrt3,0.0_DP)

         name_rap(3)='G_9     E_32'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,2)=(0.0_DP,0.0_DP)
         char_mat_proj(3,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,4)=(0.0_DP,0.0_DP)
         char_mat_proj(3,5)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,6)=(0.0_DP,0.0_DP)

      ENDIF

   CASE (82,83,84,85,86,87,88,89,90,91,92,93,94)
!
!   C_2h
!
      nsym_proj=4

      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=4

         name_rap(1)='A_g'
         char_mat_proj(1,1:2)=(1.0_DP,0.0_DP)
         char_mat_proj(1,3:4)=char_mat_proj(1,1:2)

         name_rap(2)='B_g'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,3:4)=char_mat_proj(2,1:2)

         name_rap(3)='A_u'
         char_mat_proj(3,1:2)=(1.0_DP,0.0_DP)
         char_mat_proj(3,3:4)=-char_mat_proj(3,1:2)

         name_rap(4)='B_u'
         char_mat_proj(4,1)=(1.0_DP,0.0_DP)
         char_mat_proj(4,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,3:4)=-char_mat_proj(4,1:2)


      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!  double point group representations
!
         nrap_proj=4

         name_rap(1)='G_3+    M-12g'
         char_mat_proj(1,1)=(1.0_DP,0.0_DP)
         char_mat_proj(1,2)=(0.0_DP,1.0_DP)
         char_mat_proj(1,3:4)=char_mat_proj(1,1:2)

         name_rap(2)='G_4+    M12g'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)=(0.0_DP,-1.0_DP)
         char_mat_proj(2,3:4)=char_mat_proj(2,1:2)

         name_rap(3)='G_3-    M-12u'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)=(0.0_DP,1.0_DP)
         char_mat_proj(3,3:4)=-char_mat_proj(3,1:2)


         name_rap(4)='G_4-    M12u'
         char_mat_proj(4,1)=(1.0_DP,0.0_DP)
         char_mat_proj(4,2)=(0.0_DP,-1.0_DP)
         char_mat_proj(4,3:4)=-char_mat_proj(4,1:2)

      ELSEIF (ptype(1)==1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN
!
!  projective representations (beta=-1)
!
         nrap_proj=1

         name_rap(1)='BA      M_1M_0'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN
!
!  projective representations (beta=-1)
!
         nrap_proj=1

         name_rap(1)='G_4G_3  M_12M-12'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)

      ENDIF
   CASE (127,128,129,130,131)
!
!  S_6
!     
      w=CMPLX(-0.5_DP,-sqrt3*0.5_DP)

      nsym_proj=6
      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=6

         name_rap(1)='A_g'
         char_mat_proj(1,1:3)=(1.0_DP,0.0_DP)
         char_mat_proj(1,4:6)=char_mat_proj(1,1:3)

         name_rap(2)='E_g     M2'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)= w **2
         char_mat_proj(2,3)= w 
         char_mat_proj(2,4:6)=char_mat_proj(2,1:3)
     
         name_rap(3)='E*_g    M1'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)= w 
         char_mat_proj(3,3)= w **2
         char_mat_proj(3,4:6)=char_mat_proj(3,1:3)

         name_rap(4)='A_u'
         char_mat_proj(4,1:3)=(1.0_DP,0.0_DP)
         char_mat_proj(4,4:6)=-char_mat_proj(4,1:3)

         name_rap(5)='E_u     M2'
         char_mat_proj(5,1)=(1.0_DP,0.0_DP)
         char_mat_proj(5,2)= w **2
         char_mat_proj(5,3)= w 
         char_mat_proj(5,4:6)=-char_mat_proj(5,1:3)
     
         name_rap(6)='E*_u    M1'
         char_mat_proj(6,1)=(1.0_DP,0.0_DP)
         char_mat_proj(6,2)= w 
         char_mat_proj(6,3)= w **2
         char_mat_proj(6,4:6)=-char_mat_proj(6,1:3)

      ELSEIF(ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations
!
         w=CMPLX(0.5_DP,sqrt3*0.5_DP)
         nrap_proj=6

         name_rap(1)='G_4+     M-12g'
         char_mat_proj(1,1)=(1.0_DP,0.0_DP)
         char_mat_proj(1,2)=  w 
         char_mat_proj(1,3)= CONJG(w) 
         char_mat_proj(1,4:6)=char_mat_proj(1,1:3)

         name_rap(2)='G_5+     M12g'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)= CONJG(w) 
         char_mat_proj(2,3)= w 
         char_mat_proj(2,4:6)=char_mat_proj(2,1:3)

         name_rap(3)='G_6+     M32g'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)=(-1.0_DP,0.0_DP) 
         char_mat_proj(3,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,4:6)=char_mat_proj(3,1:3)

         name_rap(4)='G_4-     M-12u'
         char_mat_proj(4,1)=(1.0_DP,0.0_DP)
         char_mat_proj(4,2)=  w 
         char_mat_proj(4,3)= CONJG(w) 
         char_mat_proj(4,4:6)=-char_mat_proj(4,1:3)

         name_rap(5)='G_5-     M12u'
         char_mat_proj(5,1)=(1.0_DP,0.0_DP)
         char_mat_proj(5,2)= CONJG(w) 
         char_mat_proj(5,3)= w 
         char_mat_proj(5,4:6)=-char_mat_proj(5,1:3)

         name_rap(6)='G_6-     M32u'
         char_mat_proj(6,1)=(1.0_DP,0.0_DP)
         char_mat_proj(6,2)=(-1.0_DP,0.0_DP) 
         char_mat_proj(6,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,4:6)=-char_mat_proj(6,1:3)

      ENDIF
   CASE (96,97,98)
!
!  C_4h
!
      nsym_proj=8
      w=CMPLX(sqrt2/2.0_DP,-sqrt2/2.0_DP)

      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=8

         name_rap(1)='A_g'
         char_mat_proj(1,1:4)=(1.0_DP,0.0_DP)
         char_mat_proj(1,5:8)=char_mat_proj(1,1:4)

         name_rap(2)='B_g     M_2g'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,3)=(1.0_DP,0.0_DP)
         char_mat_proj(2,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,5:8)=char_mat_proj(2,1:4)

         name_rap(3)='E_g     M-1g'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)=(0.0_DP,1.0_DP)
         char_mat_proj(3,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,4)=(0.0_DP,-1.0_DP)
         char_mat_proj(3,5:8)=char_mat_proj(3,1:4)
     
         name_rap(4)='E_g*    M_1g'
         char_mat_proj(4,1)=(1.0_DP,0.0_DP)
         char_mat_proj(4,2)=(0.0_DP,-1.0_DP)
         char_mat_proj(4,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,4)=(0.0_DP,1.0_DP)
         char_mat_proj(4,5:8)=char_mat_proj(4,1:4)

         name_rap(5)='A_u'
         char_mat_proj(5,1:4)=(1.0_DP,0.0_DP)
         char_mat_proj(5,5:8)=-char_mat_proj(5,1:4)

         name_rap(6)='B_u     M_2u'
         char_mat_proj(6,1)=(1.0_DP,0.0_DP)
         char_mat_proj(6,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,3)=(1.0_DP,0.0_DP)
         char_mat_proj(6,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,5:8)=-char_mat_proj(6,1:4)

         name_rap(7)='E_u     M-1u'
         char_mat_proj(7,1)=(1.0_DP,0.0_DP)
         char_mat_proj(7,2)=(0.0_DP,1.0_DP)
         char_mat_proj(7,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(7,4)=(0.0_DP,-1.0_DP)
         char_mat_proj(7,5:8)=-char_mat_proj(7,1:4)
     
         name_rap(8)='E_u*    M_1u'
         char_mat_proj(8,1)=(1.0_DP,0.0_DP)
         char_mat_proj(8,2)=(0.0_DP,-1.0_DP)
         char_mat_proj(8,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(8,4)=(0.0_DP,1.0_DP)
         char_mat_proj(8,5:8)=-char_mat_proj(8,1:4)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations
!
         nrap_proj=8

         name_rap(1)='G_5+    M-12g'
         char_mat_proj(1,1)=(1.0_DP,0.0_DP)
         char_mat_proj(1,2)= -w**3
         char_mat_proj(1,3)=(0.0_DP,1.0_DP)
         char_mat_proj(1,4)= -CONJG(w**3)
         char_mat_proj(1,5:8)=char_mat_proj(1,1:4)

         name_rap(2)='G_6+    M12g'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)= w
         char_mat_proj(2,3)= w**2
         char_mat_proj(2,4)= CONJG(w)
         char_mat_proj(2,5:8)=char_mat_proj(2,1:4)

         name_rap(3)='G_7+    M32g'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)= w**3
         char_mat_proj(3,3)=(0.0_DP,1.0_DP)
         char_mat_proj(3,4)= CONJG(w**3)
         char_mat_proj(3,5:8)=char_mat_proj(3,1:4)

         name_rap(4)='G_8+    M-32g'
         char_mat_proj(4,1)=(1.0_DP,0.0_DP)
         char_mat_proj(4,2)= -w
         char_mat_proj(4,3)=(0.0_DP,-1.0_DP)
         char_mat_proj(4,4)= -CONJG(w)
         char_mat_proj(4,5:8)=char_mat_proj(4,1:4)

         name_rap(5)='G_5-    M-12u'
         char_mat_proj(5,1)=(1.0_DP,0.0_DP)
         char_mat_proj(5,2)= -w**3
         char_mat_proj(5,3)=(0.0_DP,1.0_DP)
         char_mat_proj(5,4)= -CONJG(w**3)
         char_mat_proj(5,5:8)=-char_mat_proj(5,1:4)

         name_rap(6)='G_6-    M12u'
         char_mat_proj(6,1)=(1.0_DP,0.0_DP)
         char_mat_proj(6,2)= w
         char_mat_proj(6,3)= w**2
         char_mat_proj(6,4)= CONJG(w)
         char_mat_proj(6,5:8)=-char_mat_proj(6,1:4)

         name_rap(7)='G_7-    M32u'
         char_mat_proj(7,1)=(1.0_DP,0.0_DP)
         char_mat_proj(7,2)= w**3
         char_mat_proj(7,3)=(0.0_DP,1.0_DP)
         char_mat_proj(7,4)= CONJG(w**3)
         char_mat_proj(7,5:8)=-char_mat_proj(7,1:4)

         name_rap(8)='G_8-    M-32u'
         char_mat_proj(8,1)=(1.0_DP,0.0_DP)
         char_mat_proj(8,2)= -w
         char_mat_proj(8,3)=(0.0_DP,-1.0_DP)
         char_mat_proj(8,4)= -CONJG(w)
         char_mat_proj(8,5:8)=-char_mat_proj(8,1:4)

      ELSEIF (ptype(1)==1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN
!
!  projective representations (beta=-1) point group
!
         nrap_proj=2

         name_rap(1)='E*E     M_1M-1'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,3)=(-2.0_DP,0.0_DP)

         name_rap(2)='AB      M_2M_0'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,3)=(2.0_DP,0.0_DP)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN
!
!  projective representations (beta=-1) double point group
!
         nrap_proj=2

         name_rap(1)='G_6G_8  M12M-32'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,3)=(0.0_DP,-2.0_DP)

         name_rap(2)='G_7G_5  M_32M-12'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,3)=(0.0_DP,2.0_DP)
      ENDIF

   CASE (99)
!
!  C_6h
!
      nsym_proj=12

      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         w=CMPLX(0.5_DP,-sqrt3*0.5_DP)

         nrap_proj=12

         name_rap(1)='A_g'
         char_mat_proj(1,1:6)=(1.0_DP,0.0_DP)
         char_mat_proj(1,7:12)=char_mat_proj(1,1:6)

         name_rap(2)='B_g     M3g'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,3)=(1.0_DP,0.0_DP)
         char_mat_proj(2,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,5)=(1.0_DP,0.0_DP)
         char_mat_proj(2,6)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,7:12)=char_mat_proj(2,1:6)

         w1=CONJG(w)
         name_rap(3)='E_1g    M-1g'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)= w1
         char_mat_proj(3,3)= w1**2
         char_mat_proj(3,4)= w1**3
         char_mat_proj(3,5)= w1**4
         char_mat_proj(3,6)= w1**5
         char_mat_proj(3,7:12)=char_mat_proj(3,1:6)

         name_rap(4)='E_1*g   M_1g'
         char_mat_proj(4,1)=(1.0_DP,0.0_DP)
         char_mat_proj(4,2)= w
         char_mat_proj(4,3)= w**2
         char_mat_proj(4,4)= w**3
         char_mat_proj(4,5)= w**4
         char_mat_proj(4,6)= w**5
         char_mat_proj(4,7:12)=char_mat_proj(4,1:6)

         w1=CONJG(w**2)
         name_rap(5)='E_2g    M-2g'
         char_mat_proj(5,1)=(1.0_DP,0.0_DP)
         char_mat_proj(5,2)= w1
         char_mat_proj(5,3)= w1**2
         char_mat_proj(5,4)= w1**3
         char_mat_proj(5,5)= w1**4
         char_mat_proj(5,6)= w1**5
         char_mat_proj(5,7:12)=char_mat_proj(5,1:6)

         w1=w**2
         name_rap(6)='E_2*g   M_2g'
         char_mat_proj(6,1)=(1.0_DP,0.0_DP)
         char_mat_proj(6,2)= w1
         char_mat_proj(6,3)= w1**2
         char_mat_proj(6,4)= w1**3
         char_mat_proj(6,5)= w1**4
         char_mat_proj(6,6)= w1**5
         char_mat_proj(6,7:12)=char_mat_proj(6,1:6)

         name_rap(7)='A_u'
         char_mat_proj(7,1:6)=(1.0_DP,0.0_DP)
         char_mat_proj(7,7:12)=-char_mat_proj(7,1:6)

         name_rap(8)='B_u'
         char_mat_proj(8,1)=(1.0_DP,0.0_DP)
         char_mat_proj(8,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(8,3)=(1.0_DP,0.0_DP)
         char_mat_proj(8,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(8,5)=(1.0_DP,0.0_DP)
         char_mat_proj(8,6)=(-1.0_DP,0.0_DP)
         char_mat_proj(8,7:12)=-char_mat_proj(8,1:6)

         w1=CONJG(w)
         name_rap(9)='E_1u    M-1u'
         char_mat_proj(9,1)=(1.0_DP,0.0_DP)
         char_mat_proj(9,2)= w1
         char_mat_proj(9,3)= w1**2
         char_mat_proj(9,4)= w1**3
         char_mat_proj(9,5)= w1**4
         char_mat_proj(9,6)= w1**5
         char_mat_proj(9,7:12)=-char_mat_proj(9,1:6)

         name_rap(10)='E_1*u   M_1u'
         char_mat_proj(10,1)=(1.0_DP,0.0_DP)
         char_mat_proj(10,2)= w
         char_mat_proj(10,3)= w**2
         char_mat_proj(10,4)= w**3
         char_mat_proj(10,5)= w**4
         char_mat_proj(10,6)= w**5
         char_mat_proj(10,7:12)=-char_mat_proj(10,1:6)

         w1=CONJG(w**2)
         name_rap(11)='E_2u    M-2u'
         char_mat_proj(11,1)=(1.0_DP,0.0_DP)
         char_mat_proj(11,2)= w1
         char_mat_proj(11,3)= w1**2
         char_mat_proj(11,4)= w1**3
         char_mat_proj(11,5)= w1**4
         char_mat_proj(11,6)= w1**5
         char_mat_proj(11,7:12)=-char_mat_proj(11,1:6)

         w1=w**2
         name_rap(12)='E_2*u   M_2u'
         char_mat_proj(12,1)=(1.0_DP,0.0_DP)
         char_mat_proj(12,2)= w1
         char_mat_proj(12,3)= w1**2
         char_mat_proj(12,4)= w1**3
         char_mat_proj(12,5)= w1**4
         char_mat_proj(12,6)= w1**5
         char_mat_proj(12,7:12)=-char_mat_proj(12,1:6)


      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations
!
         nrap_proj=12
         w=CMPLX(sqrt3*0.5_DP,-0.5_DP)

         w1=CONJG(w)
         name_rap(1)='G_7+    M-12g'
         char_mat_proj(1,1)=(1.0_DP,0.0_DP)
         char_mat_proj(1,2)= w1
         char_mat_proj(1,3)= w1**2
         char_mat_proj(1,4)= w1**3
         char_mat_proj(1,5)= CONJG(w1**2)
         char_mat_proj(1,6)= CONJG(w1)
         char_mat_proj(1,7:12)=char_mat_proj(1,1:6)

         name_rap(2)='G_8+    M_12g'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)= w
         char_mat_proj(2,3)= w**2
         char_mat_proj(2,4)= w**3
         char_mat_proj(2,5)= CONJG(w**2)
         char_mat_proj(2,6)= CONJG(w)
         char_mat_proj(2,7:12)=char_mat_proj(2,1:6)

         w1=w**5
         name_rap(3)='G_9+    M_52g'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)= w1
         char_mat_proj(3,3)= w1**2
         char_mat_proj(3,4)= w1**3
         char_mat_proj(3,5)= CONJG(w1**2)
         char_mat_proj(3,6)= CONJG(w1)
         char_mat_proj(3,7:12)=char_mat_proj(3,1:6)

         w1=CONJG(w**5)
         name_rap(4)='G_10+   M-52g'
         char_mat_proj(4,1)=(1.0_DP,0.0_DP)
         char_mat_proj(4,2)= w1
         char_mat_proj(4,3)= w1**2
         char_mat_proj(4,4)= w1**3
         char_mat_proj(4,5)= CONJG(w1**2)
         char_mat_proj(4,6)= CONJG(w1)
         char_mat_proj(4,7:12)=char_mat_proj(4,1:6)

         w1=CONJG(w**3)
         name_rap(5)='G_11+   M-32g'
         char_mat_proj(5,1)=(1.0_DP,0.0_DP)
         char_mat_proj(5,2)= w1
         char_mat_proj(5,3)= w1**2
         char_mat_proj(5,4)= w1**3
         char_mat_proj(5,5)= CONJG(w1**2)
         char_mat_proj(5,6)= CONJG(w1)
         char_mat_proj(5,7:12)=char_mat_proj(5,1:6)
     
         w1=w**3
         name_rap(6)='G_12+   M_32g'
         char_mat_proj(6,1)=(1.0_DP,0.0_DP)
         char_mat_proj(6,2)= w1
         char_mat_proj(6,3)= w1**2
         char_mat_proj(6,4)= w1**3
         char_mat_proj(6,5)= CONJG(w1**2)
         char_mat_proj(6,6)= CONJG(w1)
         char_mat_proj(6,7:12)=char_mat_proj(6,1:6)
     

         w1=CONJG(w)
         name_rap(7)='G_7-    M-12u'
         char_mat_proj(7,1)=(1.0_DP,0.0_DP)
         char_mat_proj(7,2)= w1
         char_mat_proj(7,3)= w1**2
         char_mat_proj(7,4)= w1**3
         char_mat_proj(7,5)= CONJG(w1**2)
         char_mat_proj(7,6)= CONJG(w1)
         char_mat_proj(7,7:12)=-char_mat_proj(7,1:6)

         name_rap(8)='G_8-    M_12u'
         char_mat_proj(8,1)=(1.0_DP,0.0_DP)
         char_mat_proj(8,2)= w
         char_mat_proj(8,3)= w**2
         char_mat_proj(8,4)= w**3
         char_mat_proj(8,5)= CONJG(w**2)
         char_mat_proj(8,6)= CONJG(w)
         char_mat_proj(8,7:12)=-char_mat_proj(8,1:6)

         w1=w**5
         name_rap(9)='G_9-    M_52u'
         char_mat_proj(9,1)=(1.0_DP,0.0_DP)
         char_mat_proj(9,2)= w1
         char_mat_proj(9,3)= w1**2
         char_mat_proj(9,4)= w1**3
         char_mat_proj(9,5)= CONJG(w1**2)
         char_mat_proj(9,6)= CONJG(w1)
         char_mat_proj(9,7:12)=-char_mat_proj(9,1:6)

         w1=CONJG(w**5)
         name_rap(10)='G_10-   M-52u'
         char_mat_proj(10,1)=(1.0_DP,0.0_DP)
         char_mat_proj(10,2)= w1
         char_mat_proj(10,3)= w1**2
         char_mat_proj(10,4)= w1**3
         char_mat_proj(10,5)= CONJG(w1**2)
         char_mat_proj(10,6)= CONJG(w1)
         char_mat_proj(10,7:12)=-char_mat_proj(10,1:6)

         w1=CONJG(w**3)
         name_rap(11)='G_11-   M-32u'
         char_mat_proj(11,1)=(1.0_DP,0.0_DP)
         char_mat_proj(11,2)= w1
         char_mat_proj(11,3)= w1**2
         char_mat_proj(11,4)= w1**3
         char_mat_proj(11,5)= CONJG(w1**2)
         char_mat_proj(11,6)= CONJG(w1)
         char_mat_proj(11,7:12)=-char_mat_proj(11,1:6)

         w1=w**3
         name_rap(12)='G_12-   M_32u'
         char_mat_proj(12,1)=(1.0_DP,0.0_DP)
         char_mat_proj(12,2)= w1
         char_mat_proj(12,3)= w1**2
         char_mat_proj(12,4)= w1**3
         char_mat_proj(12,5)= CONJG(w1**2)
         char_mat_proj(12,6)= CONJG(w1)
         char_mat_proj(12,7:12)=-char_mat_proj(12,1:6)
     

      ELSEIF (ptype(1)==1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN
!
!   projective representations (beta=-1) point group
!
         nrap_proj=3

         name_rap(1)='E_1*E_2 M_1M_-2'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,3)=CMPLX(-1.0_DP,-sqrt3)
         char_mat_proj(1,5)=CMPLX(-1.0_DP,sqrt3) 

         name_rap(2)='E_2*E_1 M_2M-1'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,3)=CMPLX(-1.0_DP,sqrt3)
         char_mat_proj(2,5)=CMPLX(-1.0_DP,-sqrt3) 

         name_rap(3)='BA      M_3M_0'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,3)=(2.0_DP,0.0_DP)
         char_mat_proj(3,5)=(2.0_DP,0.0_DP) 

      ELSEIF (ptype(1)==-1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN
!
!   projective representations (beta=-1) double point group
!
         nrap_proj=3

         name_rap(1)='G_8G_10 M_12M-52'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,3)=CMPLX(1.0_DP,-sqrt3)
         char_mat_proj(1,5)=CMPLX(1.0_DP, sqrt3) 

         name_rap(2)='G12G11  M_32M-32'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(2,5)=(-2.0_DP,0.0_DP) 

         name_rap(3)='G_9G_7  M_52M-12'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,3)=CMPLX(1.0_DP,sqrt3)
         char_mat_proj(3,5)=CMPLX(1.0_DP,-sqrt3) 

      ENDIF

   CASE (100,101,102,103,104,105)
!
!  D_2h
!
      nsym_proj=8
      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=8

         name_rap(1)='A_1g'
         char_mat_proj(1,1:4)=(1.0_DP,0.0_DP)
         char_mat_proj(1,5:8)=char_mat_proj(1,1:4)

         name_rap(2)='B_1g    A_2g'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2)=(1.0_DP,0.0_DP)
         char_mat_proj(2,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,5:8)=char_mat_proj(2,1:4)

         name_rap(3)='B_2g    B_1g'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,3)=(1.0_DP,0.0_DP)
         char_mat_proj(3,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,5:8)=char_mat_proj(3,1:4)

         name_rap(4)='B_3g    B_2g'
         char_mat_proj(4,1)=(1.0_DP,0.0_DP)
         char_mat_proj(4,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,4)=(1.0_DP,0.0_DP)
         char_mat_proj(4,5:8)=char_mat_proj(4,1:4)

         name_rap(5)='A_1u'
         char_mat_proj(5,1:4)=(1.0_DP,0.0_DP)
         char_mat_proj(5,5:8)=-char_mat_proj(5,1:4)

         name_rap(6)='B_1u    A_2u'
         char_mat_proj(6,1)=(1.0_DP,0.0_DP)
         char_mat_proj(6,2)=(1.0_DP,0.0_DP)
         char_mat_proj(6,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,5:8)=-char_mat_proj(6,1:4)

         name_rap(7)='B_2u    B_1u'
         char_mat_proj(7,1)=(1.0_DP,0.0_DP)
         char_mat_proj(7,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(7,3)=(1.0_DP,0.0_DP)
         char_mat_proj(7,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(7,5:8)=-char_mat_proj(7,1:4)

         name_rap(8)='B_3u    B_2u'
         char_mat_proj(8,1)=(1.0_DP,0.0_DP)
         char_mat_proj(8,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(8,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(8,4)=(1.0_DP,0.0_DP)
         char_mat_proj(8,5:8)=-char_mat_proj(8,1:4)


      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations
!
         nrap_proj=2
         name_rap(1)='G_5+    E_12g'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,5:8)=char_mat_proj(1,1:4)


         name_rap(2)='G_5-    E_12u'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,5:8)=-char_mat_proj(2,1:4)

      ELSEIF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==-1) THEN
!
!   projective representations (gamma=-1)
!
         nrap_proj=2

         name_rap(1)='AB_1    A_1A_2'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,2)=(2.0_DP,0.0_DP)

         name_rap(2)='B_2B_3  B_1B_2'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,2)=(-2.0_DP,0.0_DP)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==-1) THEN

         nrap_proj=2

         name_rap(1)='G_5y    E_12y'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,6)=(0.0_DP,-2.0_DP)

         name_rap(2)='G_5-y   E_12-y'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,6)=(0.0_DP,2.0_DP)

      ELSEIF (ptype(1)==1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN
!
!   projective representations (beta=-1)
!

         nrap_proj=2

         name_rap(1)='AB_2    A_1B_1'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,3)=(2.0_DP,0.0_DP)

         name_rap(2)='B_1B_3  A_2B_2'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,3)=(-2.0_DP,0.0_DP)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN

         nrap_proj=2

         name_rap(1)='G_5z    E_12z'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,7)=(0.0_DP,-2.0_DP)

         name_rap(2)='G_5-z   E_12-z'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,7)=(0.0_DP,2.0_DP)

      ELSEIF (ptype(1)==1.AND.ptype(2)==-1.AND.ptype(3)==-1) THEN
!
!   projective representations (beta=-1, gamma=-1)
!

         nrap_proj=2

         name_rap(1)='AB_3    A_1B_2'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,4)=(2.0_DP,0.0_DP)

         name_rap(2)='B_1B_2  A_2B_1'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,4)=(-2.0_DP,0.0_DP)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==-1.AND.ptype(3)==-1) THEN

         nrap_proj=2

         name_rap(1)='G_5x    E_12x'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,8)=(0.0_DP,-2.0_DP)

         name_rap(2)='G_5-x   E_12-x'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,8)=(0.0_DP,2.0_DP)

      END IF

   CASE (118,119,120,121,122,123)
!
!  D_3d 
!
      nsym_proj=12
      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=6

         name_rap(1)='A_1g'
         char_mat_proj(1,1:12)=(1.0_DP,0.0_DP)

         name_rap(2)='A_2g'
         char_mat_proj(2,1:3)=(1.0_DP,0.0_DP)
         char_mat_proj(2,4:6)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,7:12)=char_mat_proj(2,1:6)

         name_rap(3)='E_g     E_1g'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,2:3)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,7:12)=char_mat_proj(3,1:6)

         name_rap(4)='A_1u'
         char_mat_proj(4,1:6)=(1.0_DP,0.0_DP)
         char_mat_proj(4,7:12)=-char_mat_proj(4,1:6)

         name_rap(5)='A_2u'
         char_mat_proj(5,1:3)=(1.0_DP,0.0_DP)
         char_mat_proj(5,4:6)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,7:12)=-char_mat_proj(5,1:6)

         name_rap(6)='E_u     E_1u'
         char_mat_proj(6,1)=(2.0_DP,0.0_DP)
         char_mat_proj(6,2:3)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,7:12)=-char_mat_proj(6,1:6)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations
!
         nrap_proj=6

         name_rap(1)='G4+     E12g'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,2:3)=(1.0_DP,0.0_DP)
         char_mat_proj(1,7:12)=char_mat_proj(1,1:6)

         name_rap(2)='G5+     B2g'
         char_mat_proj(2,1)=(1.0_DP,0.0_DP)
         char_mat_proj(2,2:3)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,4)=(0.0_DP,1.0_DP)
         char_mat_proj(2,5:6)=(0.0_DP,-1.0_DP)
         char_mat_proj(2,7:12)=char_mat_proj(2,1:6)

         name_rap(3)='G6+     B1g'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2:3)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,4)=(0.0_DP,-1.0_DP)
         char_mat_proj(3,5:6)=(0.0_DP,1.0_DP)
         char_mat_proj(3,7:12)=char_mat_proj(3,1:6)


         name_rap(4)='G4-     E12u'
         char_mat_proj(4,1)=(2.0_DP,0.0_DP)
         char_mat_proj(4,2:3)=(1.0_DP,0.0_DP)
         char_mat_proj(4,7:12)=-char_mat_proj(4,1:6)

         name_rap(5)='G5-     B2u'
         char_mat_proj(5,1)=(1.0_DP,0.0_DP)
         char_mat_proj(5,2:3)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,4)=(0.0_DP,1.0_DP)
         char_mat_proj(5,5:6)=(0.0_DP,-1.0_DP)
         char_mat_proj(5,7:12)=-char_mat_proj(5,1:6)

         name_rap(6)='G6-     B1u'
         char_mat_proj(6,1)=(1.0_DP,0.0_DP)
         char_mat_proj(6,2:3)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,4)=(0.0_DP,-1.0_DP)
         char_mat_proj(6,5:6)=(0.0_DP,1.0_DP)
         char_mat_proj(6,7:12)=-char_mat_proj(6,1:6)


      ELSEIF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==-1) THEN
!
!   projective representations (gamma=-1)  point group
!
         nrap_proj=3

         name_rap(1)='A_1A_2'
         char_mat_proj(1,1:3)=(2.0_DP,0.0_DP)

         name_rap(2)='Ey      E_1y'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,2:3)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,8)=CMPLX(0.0_DP,-sqrt3)
         char_mat_proj(2,9)=CMPLX(0.0_DP,sqrt3)

         name_rap(3)='E-y     E_1-y'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,2:3)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,8)=CMPLX(0.0_DP,sqrt3)
         char_mat_proj(3,9)=CMPLX(0.0_DP,-sqrt3)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==-1) THEN
!
!   projective representations (gamma=-1) double point group
!
         nrap_proj=3

         name_rap(1)='G_6G_5  B_1B_2'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,2:3)=CMPLX(-2.0_DP,0.0_DP)

         name_rap(2)='G_4y    E_12y'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,2:3)=(1.0_DP,0.0_DP)
         char_mat_proj(2,8)=CMPLX(0.0_DP,-sqrt3)
         char_mat_proj(2,9)=CMPLX(0.0_DP,sqrt3)

         name_rap(3)='G_4-y   E_12-y'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,2:3)=(1.0_DP,0.0_DP)
         char_mat_proj(3,8)=CMPLX(0.0_DP,sqrt3)
         char_mat_proj(3,9)=CMPLX(0.0_DP,-sqrt3)

      ENDIF

   CASE (108,109,110)
!
!  D_4h
!
      nsym_proj=16
      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=10

         name_rap(1)='A_1g'
         char_mat_proj(1,1:8)=(1.0_DP,0.0_DP)
         char_mat_proj(1,9:16)=char_mat_proj(1,1:8)

         name_rap(2)='A_2g'
         char_mat_proj(2,1:4)=(1.0_DP,0.0_DP)
         char_mat_proj(2,5:8)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,9:16)=char_mat_proj(2,1:8)
!
!  See the explanation for the D_4 group
!
         IF (cge==110) THEN
            name_rap(3)='B_1g'
            char_mat_proj(3,1)=(1.0_DP,0.0_DP)
            char_mat_proj(3,2)=(-1.0_DP,0.0_DP)
            char_mat_proj(3,3)=(1.0_DP,0.0_DP)
            char_mat_proj(3,4)=(-1.0_DP,0.0_DP)
            char_mat_proj(3,5)=(-1.0_DP,0.0_DP)
            char_mat_proj(3,6)=(1.0_DP,0.0_DP)
            char_mat_proj(3,7)=(-1.0_DP,0.0_DP)
            char_mat_proj(3,8)=(1.0_DP,0.0_DP)
            char_mat_proj(3,9:16)=char_mat_proj(3,1:8)

            name_rap(4)='B_2g'
            char_mat_proj(4,1)=(1.0_DP,0.0_DP)
            char_mat_proj(4,2)=(-1.0_DP,0.0_DP)
            char_mat_proj(4,3)=(1.0_DP,0.0_DP)
            char_mat_proj(4,4)=(-1.0_DP,0.0_DP)
            char_mat_proj(4,5)=(1.0_DP,0.0_DP)
            char_mat_proj(4,6)=(-1.0_DP,0.0_DP)
            char_mat_proj(4,7)=(1.0_DP,0.0_DP)
            char_mat_proj(4,8)=(-1.0_DP,0.0_DP)
            char_mat_proj(4,9:16)=char_mat_proj(4,1:8)
         ELSE
            name_rap(3)='B_1g'
            char_mat_proj(3,1)=(1.0_DP,0.0_DP)
            char_mat_proj(3,2)=(-1.0_DP,0.0_DP)
            char_mat_proj(3,3)=(1.0_DP,0.0_DP)
            char_mat_proj(3,4)=(-1.0_DP,0.0_DP)
            char_mat_proj(3,5)=(1.0_DP,0.0_DP)
            char_mat_proj(3,6)=(-1.0_DP,0.0_DP)
            char_mat_proj(3,7)=(1.0_DP,0.0_DP)
            char_mat_proj(3,8)=(-1.0_DP,0.0_DP)
            char_mat_proj(3,9:16)=char_mat_proj(3,1:8)

            name_rap(4)='B_2g'
            char_mat_proj(4,1)=(1.0_DP,0.0_DP)
            char_mat_proj(4,2)=(-1.0_DP,0.0_DP)
            char_mat_proj(4,3)=(1.0_DP,0.0_DP)
            char_mat_proj(4,4)=(-1.0_DP,0.0_DP)
            char_mat_proj(4,5)=(-1.0_DP,0.0_DP)
            char_mat_proj(4,6)=(1.0_DP,0.0_DP)
            char_mat_proj(4,7)=(-1.0_DP,0.0_DP)
            char_mat_proj(4,8)=(1.0_DP,0.0_DP)
            char_mat_proj(4,9:16)=char_mat_proj(4,1:8)
         ENDIF
         name_rap(5)='E_2g'
         char_mat_proj(5,1)=(2.0_DP,0.0_DP)
         char_mat_proj(5,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(5,9:16)=char_mat_proj(5,1:8)

         name_rap(6)='A_1u'
         char_mat_proj(6,1:8)=(1.0_DP,0.0_DP)
         char_mat_proj(6,9:16)=-char_mat_proj(6,1:8)

         name_rap(7)='A_2u'
         char_mat_proj(7,1:4)=(1.0_DP,0.0_DP)
         char_mat_proj(7,5:8)=(-1.0_DP,0.0_DP)
         char_mat_proj(7,9:16)=-char_mat_proj(7,1:8)

         IF (cge==110) THEN
            name_rap(8)='B_1u'
            char_mat_proj(8,1)=(1.0_DP,0.0_DP)
            char_mat_proj(8,2)=(-1.0_DP,0.0_DP)
            char_mat_proj(8,3)=(1.0_DP,0.0_DP)
            char_mat_proj(8,4)=(-1.0_DP,0.0_DP)
            char_mat_proj(8,5)=(-1.0_DP,0.0_DP)
            char_mat_proj(8,6)=(1.0_DP,0.0_DP)
            char_mat_proj(8,7)=(-1.0_DP,0.0_DP)
            char_mat_proj(8,8)=(1.0_DP,0.0_DP)
            char_mat_proj(8,9:16)=-char_mat_proj(8,1:8)

            name_rap(9)='B_2u'
            char_mat_proj(9,1)=(1.0_DP,0.0_DP)
            char_mat_proj(9,2)=(-1.0_DP,0.0_DP)
            char_mat_proj(9,3)=(1.0_DP,0.0_DP)
            char_mat_proj(9,4)=(-1.0_DP,0.0_DP)
            char_mat_proj(9,5)=(1.0_DP,0.0_DP)
            char_mat_proj(9,6)=(-1.0_DP,0.0_DP)
            char_mat_proj(9,7)=(1.0_DP,0.0_DP)
            char_mat_proj(9,8)=(-1.0_DP,0.0_DP)
            char_mat_proj(9,9:16)=-char_mat_proj(9,1:8)
         ELSE
            name_rap(8)='B_1u'
            char_mat_proj(8,1)=(1.0_DP,0.0_DP)
            char_mat_proj(8,2)=(-1.0_DP,0.0_DP)
            char_mat_proj(8,3)=(1.0_DP,0.0_DP)
            char_mat_proj(8,4)=(-1.0_DP,0.0_DP)
            char_mat_proj(8,5)=(1.0_DP,0.0_DP)
            char_mat_proj(8,6)=(-1.0_DP,0.0_DP)
            char_mat_proj(8,7)=(1.0_DP,0.0_DP)
            char_mat_proj(8,8)=(-1.0_DP,0.0_DP)
            char_mat_proj(8,9:16)=-char_mat_proj(8,1:8)

            name_rap(9)='B_2u'
            char_mat_proj(9,1)=(1.0_DP,0.0_DP)
            char_mat_proj(9,2)=(-1.0_DP,0.0_DP)
            char_mat_proj(9,3)=(1.0_DP,0.0_DP)
            char_mat_proj(9,4)=(-1.0_DP,0.0_DP)
            char_mat_proj(9,5)=(-1.0_DP,0.0_DP)
            char_mat_proj(9,6)=(1.0_DP,0.0_DP)
            char_mat_proj(9,7)=(-1.0_DP,0.0_DP)
            char_mat_proj(9,8)=(1.0_DP,0.0_DP)
            char_mat_proj(9,9:16)=-char_mat_proj(9,1:8)
         ENDIF
         name_rap(10)='E_2u'
         char_mat_proj(10,1)=(2.0_DP,0.0_DP)
         char_mat_proj(10,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(10,9:16)=-char_mat_proj(10,1:8)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations
!
         nrap_proj=4

         name_rap(1)='G_6+    E_12g'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,2)=CMPLX(sqrt2,0.0_DP)
         char_mat_proj(1,4)=CMPLX(sqrt2,0.0_DP)
         char_mat_proj(1,9:16)=char_mat_proj(1,1:8)

         name_rap(2)='G_7+    E_32g'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,2)=CMPLX(-sqrt2,0.0_DP)
         char_mat_proj(2,4)=CMPLX(-sqrt2,0.0_DP)
         char_mat_proj(2,9:16)=char_mat_proj(2,1:8)

         name_rap(3)='G_6-    E_12u'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,2)=CMPLX(sqrt2,0.0_DP)
         char_mat_proj(3,4)=CMPLX(sqrt2,0.0_DP)
         char_mat_proj(3,9:16)=-char_mat_proj(3,1:8)

         name_rap(4)='G_7-    E_32u'
         char_mat_proj(4,1)=(2.0_DP,0.0_DP)
         char_mat_proj(4,2)=CMPLX(-sqrt2,0.0_DP)
         char_mat_proj(4,4)=CMPLX(-sqrt2,0.0_DP)
         char_mat_proj(4,9:16)=-char_mat_proj(4,1:8)

      ELSEIF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==-1) THEN
!
!   projective representations (gamma=-1)
!
         nrap_proj=4

         name_rap(1)='A_1A_2'
         char_mat_proj(1,1:4)=(2.0_DP,0.0_DP)

         name_rap(2)='B_1B_2'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,2)=(-2.0_DP,0.0_DP)
         char_mat_proj(2,3)=(2.0_DP,0.0_DP)
         char_mat_proj(2,4)=(-2.0_DP,0.0_DP)

         name_rap(3)='E_y     E_1y'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,10)=(0.0_DP,-2.0_DP)
         char_mat_proj(3,12)=(0.0_DP,2.0_DP)

         name_rap(4)='E-y     E_1-y'
         char_mat_proj(4,1)=(2.0_DP,0.0_DP)
         char_mat_proj(4,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(4,10)=(0.0_DP,2.0_DP)
         char_mat_proj(4,12)=(0.0_DP,-2.0_DP)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==-1) THEN

         nrap_proj=4

         name_rap(1)='G_6y    E_12y'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,2)=CMPLX(sqrt2,0.0_DP)
         char_mat_proj(1,4)=CMPLX(sqrt2,0.0_DP)
         char_mat_proj(1,10)=CMPLX(0.0_DP,-sqrt2)
         char_mat_proj(1,11)=(0.0_DP,-2.0_DP)
         char_mat_proj(1,12)=CMPLX(0.0_DP,sqrt2)

         name_rap(2)='G_6-y   E_12-y'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,2)=CMPLX(sqrt2,0.0_DP)
         char_mat_proj(2,4)=CMPLX(sqrt2,0.0_DP)
         char_mat_proj(2,10)=CMPLX(0.0_DP,sqrt2)
         char_mat_proj(2,11)=(0.0_DP,2.0_DP)
         char_mat_proj(2,12)=CMPLX(0.0_DP,-sqrt2)


         name_rap(3)='G_7y     E_32y'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,2)=CMPLX(-sqrt2,0.0_DP)
         char_mat_proj(3,4)=CMPLX(-sqrt2,0.0_DP)
         char_mat_proj(3,10)=CMPLX(0.0_DP,-sqrt2)
         char_mat_proj(3,11)=(0.0_DP,2.0_DP)
         char_mat_proj(3,12)=CMPLX(0.0_DP,sqrt2)

         name_rap(4)='G_7-y    E_32-y'
         char_mat_proj(4,1)=(2.0_DP,0.0_DP)
         char_mat_proj(4,2)=CMPLX(-sqrt2,0.0_DP)
         char_mat_proj(4,4)=CMPLX(-sqrt2,0.0_DP)
         char_mat_proj(4,10)=CMPLX(0.0_DP,sqrt2)
         char_mat_proj(4,11)=(0.0_DP,-2.0_DP)
         char_mat_proj(4,12)=CMPLX(0.0_DP,-sqrt2)

      ELSEIF (ptype(1)==1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN
!
!   projective representations (beta=-1)
!
         nrap_proj=4

         IF (cge==110) THEN
            name_rap(1)='A_1B_2'
         ELSE
            name_rap(1)='A_1B_1'
         ENDIF
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,3)=(2.0_DP,0.0_DP)
         char_mat_proj(1,5)=(2.0_DP,0.0_DP)
         char_mat_proj(1,7)=(2.0_DP,0.0_DP)

         IF (cge==110) THEN
            name_rap(2)='A_2B_1'
         ELSE
            name_rap(2)='A_2B_2'
         ENDIF
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,3)=(2.0_DP,0.0_DP)
         char_mat_proj(2,5)=(-2.0_DP,0.0_DP)
         char_mat_proj(2,7)=(-2.0_DP,0.0_DP)

         name_rap(3)='E_z      E_1z'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,13)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,15)=(2.0_DP,0.0_DP)

         name_rap(4)='E-z      E_1-z'
         char_mat_proj(4,1)=(2.0_DP,0.0_DP)
         char_mat_proj(4,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(4,13)=(2.0_DP,0.0_DP)
         char_mat_proj(4,15)=(-2.0_DP,0.0_DP)


      ELSEIF (ptype(1)==-1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN

         nrap_proj=1

         name_rap(1)='G_6G_7  E_12E_32'
         char_mat_proj(1,1)=(4.0_DP,0.0_DP)

      ELSEIF (ptype(1)==1.AND.ptype(2)==-1.AND.ptype(3)==-1) THEN
!
!   projective representations (beta=-1, gamma=-1)
!
         nrap_proj=4

         IF (cge==110) THEN
            name_rap(1)='A_1B_1'
         ELSE
            name_rap(1)='A_1B_2'
         ENDIF
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,3)=(2.0_DP,0.0_DP)
         char_mat_proj(1,6)=(2.0_DP,0.0_DP)
         char_mat_proj(1,8)=(2.0_DP,0.0_DP)

         IF (cge==110) THEN
            name_rap(2)='A_2B_2'
         ELSE
            name_rap(2)='A_2B_1'
         ENDIF
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,3)=(2.0_DP,0.0_DP)
         char_mat_proj(2,6)=(-2.0_DP,0.0_DP)
         char_mat_proj(2,8)=(-2.0_DP,0.0_DP)

         name_rap(3)='E_x     E_1x'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,14)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,16)=(2.0_DP,0.0_DP)

         name_rap(4)='E-x     E_1-x'
         char_mat_proj(4,1)=(2.0_DP,0.0_DP)
         char_mat_proj(4,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(4,14)=(2.0_DP,0.0_DP)
         char_mat_proj(4,16)=(-2.0_DP,0.0_DP)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==-1.AND.ptype(3)==-1) THEN

         nrap_proj=1
         name_rap(1)='G_6G_7   E_12E_32'
         char_mat_proj(1,1)=(4.0_DP,0.0_DP)

      ENDIF

   CASE (111)
!
!  D_6h
!
      nsym_proj=24
      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations
!
         nrap_proj=12

         name_rap(1)='A_1g'
         char_mat_proj(1,1:12)=(1.0_DP,0.0_DP)
         char_mat_proj(1,13:24)=char_mat_proj(1,1:12)

         name_rap(2)='A_2g'
         char_mat_proj(2,1:6)=(1.0_DP,0.0_DP)
         char_mat_proj(2,7:12)=(-1.0_DP,0.0_DP)
         char_mat_proj(2,13:24)=char_mat_proj(2,1:12)

         name_rap(3)='B_1g'
         char_mat_proj(3,1)=(1.0_DP,0.0_DP)
         char_mat_proj(3,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,3)=(1.0_DP,0.0_DP)
         char_mat_proj(3,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,5)=(1.0_DP,0.0_DP)
         char_mat_proj(3,6)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,7)=(1.0_DP,0.0_DP)
         char_mat_proj(3,8)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,9)=(1.0_DP,0.0_DP)
         char_mat_proj(3,10)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,11)=(1.0_DP,0.0_DP)
         char_mat_proj(3,12)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,13:24)=char_mat_proj(3,1:12)


         name_rap(4)='B_2g'
         char_mat_proj(4,1)=(1.0_DP,0.0_DP)
         char_mat_proj(4,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,3)=(1.0_DP,0.0_DP)
         char_mat_proj(4,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,5)=(1.0_DP,0.0_DP)
         char_mat_proj(4,6)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,7)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,8)=(1.0_DP,0.0_DP)
         char_mat_proj(4,9)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,10)=(1.0_DP,0.0_DP)
         char_mat_proj(4,11)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,12)=(1.0_DP,0.0_DP)
         char_mat_proj(4,13:24)=char_mat_proj(4,1:12)

         name_rap(5)='E_1g'
         char_mat_proj(5,1)=(2.0_DP,0.0_DP)
         char_mat_proj(5,2)=(1.0_DP,0.0_DP)
         char_mat_proj(5,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,4)=(-2.0_DP,0.0_DP)
         char_mat_proj(5,5)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,6)=(1.0_DP,0.0_DP)
         char_mat_proj(5,13:24)=char_mat_proj(5,1:12)

         name_rap(6)='E_2g'
         char_mat_proj(6,1)=(2.0_DP,0.0_DP)
         char_mat_proj(6,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,4)=(2.0_DP,0.0_DP)
         char_mat_proj(6,5)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,6)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,13:24)=char_mat_proj(6,1:12)

         name_rap(7)='A_1u'
         char_mat_proj(7,1:12)=(1.0_DP,0.0_DP)
         char_mat_proj(7,13:24)=-char_mat_proj(7,1:12)

         name_rap(8)='A_2u'
         char_mat_proj(8,1:6)=(1.0_DP,0.0_DP)
         char_mat_proj(8,7:12)=(-1.0_DP,0.0_DP)
         char_mat_proj(8,13:24)=-char_mat_proj(8,1:12)

         name_rap(9)='B_1u'
         char_mat_proj(9,1)=(1.0_DP,0.0_DP)
         char_mat_proj(9,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(9,3)=(1.0_DP,0.0_DP)
         char_mat_proj(9,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(9,5)=(1.0_DP,0.0_DP)
         char_mat_proj(9,6)=(-1.0_DP,0.0_DP)
         char_mat_proj(9,7)=(1.0_DP,0.0_DP)
         char_mat_proj(9,8)=(-1.0_DP,0.0_DP)
         char_mat_proj(9,9)=(1.0_DP,0.0_DP)
         char_mat_proj(9,10)=(-1.0_DP,0.0_DP)
         char_mat_proj(9,11)=(1.0_DP,0.0_DP)
         char_mat_proj(9,12)=(-1.0_DP,0.0_DP)
         char_mat_proj(9,13:24)=-char_mat_proj(9,1:12)

         name_rap(10)='B_2u'
         char_mat_proj(10,1)=(1.0_DP,0.0_DP)
         char_mat_proj(10,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(10,3)=(1.0_DP,0.0_DP)
         char_mat_proj(10,4)=(-1.0_DP,0.0_DP)
         char_mat_proj(10,5)=(1.0_DP,0.0_DP)
         char_mat_proj(10,6)=(-1.0_DP,0.0_DP)
         char_mat_proj(10,7)=(-1.0_DP,0.0_DP)
         char_mat_proj(10,8)=(1.0_DP,0.0_DP)
         char_mat_proj(10,9)=(-1.0_DP,0.0_DP)
         char_mat_proj(10,10)=(1.0_DP,0.0_DP)
         char_mat_proj(10,11)=(-1.0_DP,0.0_DP)
         char_mat_proj(10,12)=(1.0_DP,0.0_DP)
         char_mat_proj(10,13:24)=-char_mat_proj(10,1:12)

         name_rap(11)='E_1u'
         char_mat_proj(11,1)=(2.0_DP,0.0_DP)
         char_mat_proj(11,2)=(1.0_DP,0.0_DP)
         char_mat_proj(11,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(11,4)=(-2.0_DP,0.0_DP)
         char_mat_proj(11,5)=(-1.0_DP,0.0_DP)
         char_mat_proj(11,6)=(1.0_DP,0.0_DP)
         char_mat_proj(11,13:24)=-char_mat_proj(11,1:12)

         name_rap(12)='E_2u'
         char_mat_proj(12,1)=(2.0_DP,0.0_DP)
         char_mat_proj(12,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(12,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(12,4)=(2.0_DP,0.0_DP)
         char_mat_proj(12,5)=(-1.0_DP,0.0_DP)
         char_mat_proj(12,6)=(-1.0_DP,0.0_DP)
         char_mat_proj(12,13:24)=-char_mat_proj(12,1:12)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations
!
         nrap_proj=6


         name_rap(1)='G_7+    E_12g'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,2)=CMPLX(sqrt3,0.0_DP)
         char_mat_proj(1,3)=(1.0_DP,0.0_DP)
         char_mat_proj(1,4)=(0.0_DP,0.0_DP)
         char_mat_proj(1,5)=(1.0_DP,0.0_DP)
         char_mat_proj(1,6)=CMPLX(sqrt3,0.0_DP)
         char_mat_proj(1,13:24)=char_mat_proj(1,1:12)

         name_rap(2)='G_8+    E_52g'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,2)=CMPLX(-sqrt3,0.0_DP)
         char_mat_proj(2,3)=(1.0_DP,0.0_DP)
         char_mat_proj(2,4)=(0.0_DP,0.0_DP)
         char_mat_proj(2,5)=(1.0_DP,0.0_DP)
         char_mat_proj(2,6)=CMPLX(-sqrt3,0.0_DP)
         char_mat_proj(2,13:24)=char_mat_proj(2,1:12)

         name_rap(3)='G_9+    E_32g'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,2)=(0.0_DP,0.0_DP)
         char_mat_proj(3,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,4)=(0.0_DP,0.0_DP)
         char_mat_proj(3,5)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,6)=(0.0_DP,0.0_DP)
         char_mat_proj(3,13:24)=char_mat_proj(3,1:12)

         name_rap(4)='G_7-    E_12u'
         char_mat_proj(4,1)=(2.0_DP,0.0_DP)
         char_mat_proj(4,2)=CMPLX(sqrt3,0.0_DP)
         char_mat_proj(4,3)=(1.0_DP,0.0_DP)
         char_mat_proj(4,4)=(0.0_DP,0.0_DP)
         char_mat_proj(4,5)=(1.0_DP,0.0_DP)
         char_mat_proj(4,6)=CMPLX(sqrt3,0.0_DP)
         char_mat_proj(4,13:24)=-char_mat_proj(4,1:12)

         name_rap(5)='G_8-    E_52u'
         char_mat_proj(5,1)=(2.0_DP,0.0_DP)
         char_mat_proj(5,2)=CMPLX(-sqrt3,0.0_DP)
         char_mat_proj(5,3)=(1.0_DP,0.0_DP)
         char_mat_proj(5,4)=(0.0_DP,0.0_DP)
         char_mat_proj(5,5)=(1.0_DP,0.0_DP)
         char_mat_proj(5,6)=CMPLX(-sqrt3,0.0_DP)
         char_mat_proj(5,13:24)=-char_mat_proj(5,1:12)

         name_rap(6)='G_9-    E_32u'
         char_mat_proj(6,1)=(2.0_DP,0.0_DP)
         char_mat_proj(6,2)=(0.0_DP,0.0_DP)
         char_mat_proj(6,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(6,4)=(0.0_DP,0.0_DP)
         char_mat_proj(6,5)=(-2.0_DP,0.0_DP)
         char_mat_proj(6,6)=(0.0_DP,0.0_DP)
         char_mat_proj(6,13:24)=-char_mat_proj(6,1:12)


      ELSEIF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==-1) THEN
!
!   projective representations (gamma=-1) point group
!
         nrap_proj=6

         name_rap(1)='A_1A_2'
         char_mat_proj(1,1:6)=(2.0_DP,0.0_DP)

         name_rap(2)='B_1B_2'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,2)=(-2.0_DP,0.0_DP)
         char_mat_proj(2,3)=(2.0_DP,0.0_DP)
         char_mat_proj(2,4)=(-2.0_DP,0.0_DP)
         char_mat_proj(2,5)=(2.0_DP,0.0_DP)
         char_mat_proj(2,6)=(-2.0_DP,0.0_DP)

         name_rap(3)='E_1y'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,2)=(1.0_DP,0.0_DP)
         char_mat_proj(3,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,4)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,5)=(-1.0_DP,0.0_DP)
         char_mat_proj(3,6)=(1.0_DP,0.0_DP)
         char_mat_proj(3,14)=CMPLX(0.0_DP,-sqrt3)
         char_mat_proj(3,15)=CMPLX(0.0_DP,-sqrt3)
         char_mat_proj(3,17)=CMPLX(0.0_DP,sqrt3)
         char_mat_proj(3,18)=CMPLX(0.0_DP,sqrt3)
   
         name_rap(4)='E_1-y'
         char_mat_proj(4,1)=(2.0_DP,0.0_DP)
         char_mat_proj(4,2)=(1.0_DP,0.0_DP)
         char_mat_proj(4,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,4)=(-2.0_DP,0.0_DP)
         char_mat_proj(4,5)=(-1.0_DP,0.0_DP)
         char_mat_proj(4,6)=(1.0_DP,0.0_DP)
         char_mat_proj(4,14)=CMPLX(0.0_DP,sqrt3)
         char_mat_proj(4,15)=CMPLX(0.0_DP,sqrt3)
         char_mat_proj(4,17)=CMPLX(0.0_DP,-sqrt3)
         char_mat_proj(4,18)=CMPLX(0.0_DP,-sqrt3)

         name_rap(5)='E_2y'
         char_mat_proj(5,1)=(2.0_DP,0.0_DP)
         char_mat_proj(5,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,4)=(2.0_DP,0.0_DP)
         char_mat_proj(5,5)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,6)=(-1.0_DP,0.0_DP)
         char_mat_proj(5,14)=CMPLX(0.0_DP,-sqrt3)
         char_mat_proj(5,15)=CMPLX(0.0_DP,sqrt3)
         char_mat_proj(5,17)=CMPLX(0.0_DP,-sqrt3)
         char_mat_proj(5,18)=CMPLX(0.0_DP,sqrt3)

         name_rap(6)='E_2-y'
         char_mat_proj(6,1)=(2.0_DP,0.0_DP)
         char_mat_proj(6,2)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,3)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,4)=(2.0_DP,0.0_DP)
         char_mat_proj(6,5)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,6)=(-1.0_DP,0.0_DP)
         char_mat_proj(6,14)=CMPLX(0.0_DP,sqrt3)
         char_mat_proj(6,15)=CMPLX(0.0_DP,-sqrt3)
         char_mat_proj(6,17)=CMPLX(0.0_DP,sqrt3)
         char_mat_proj(6,18)=CMPLX(0.0_DP,-sqrt3)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==-1) THEN
!
!   projective representations (gamma=-1) double point group
!
         nrap_proj=6

         name_rap(1)='G_7y    E_12y'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,2)=CMPLX(sqrt3,0.0_DP)
         char_mat_proj(1,3)=(1.0_DP,0.0_DP)
         char_mat_proj(1,5)=(1.0_DP,0.0_DP)
         char_mat_proj(1,6)=CMPLX(sqrt3,0.0_DP)
         char_mat_proj(1,14)=(0.0_DP,-1.0_DP)
         char_mat_proj(1,15)=CMPLX(0.0_DP,-sqrt3)
         char_mat_proj(1,16)=(0.0_DP,-2.0_DP)
         char_mat_proj(1,17)=CMPLX(0.0_DP,sqrt3)
         char_mat_proj(1,18)=(0.0_DP,1.0_DP)

         name_rap(2)='G_7-y   E_12-y'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,2)=CMPLX(sqrt3,0.0_DP)
         char_mat_proj(2,3)=(1.0_DP,0.0_DP)
         char_mat_proj(2,5)=(1.0_DP,0.0_DP)
         char_mat_proj(2,6)=CMPLX(sqrt3,0.0_DP)
         char_mat_proj(2,14)=(0.0_DP,1.0_DP)
         char_mat_proj(2,15)=CMPLX(0.0_DP,sqrt3)
         char_mat_proj(2,16)=(0.0_DP,2.0_DP)
         char_mat_proj(2,17)=CMPLX(0.0_DP,-sqrt3)
         char_mat_proj(2,18)=(0.0_DP,-1.0_DP)
   

         name_rap(3)='G_9y    E_32y'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,5)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,14)=(0.0_DP,-2.0_DP)
         char_mat_proj(3,16)=(0.0_DP,2.0_DP)
         char_mat_proj(3,18)=(0.0_DP,2.0_DP)

         name_rap(4)='G_9-y   E_32-y'
         char_mat_proj(4,1)=(2.0_DP,0.0_DP)
         char_mat_proj(4,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(4,5)=(-2.0_DP,0.0_DP)
         char_mat_proj(4,14)=(0.0_DP,2.0_DP)
         char_mat_proj(4,16)=(0.0_DP,-2.0_DP)
         char_mat_proj(4,18)=(0.0_DP,-2.0_DP)
   
   
         name_rap(5)='G_8y    E_52y'
         char_mat_proj(5,1)=(2.0_DP,0.0_DP)
         char_mat_proj(5,2)=CMPLX(-sqrt3,0.0_DP)
         char_mat_proj(5,3)=(1.0_DP,0.0_DP)
         char_mat_proj(5,5)=(1.0_DP,0.0_DP)
         char_mat_proj(5,6)=CMPLX(-sqrt3,0.0_DP)
         char_mat_proj(5,14)=(0.0_DP,-1.0_DP)
         char_mat_proj(5,15)=CMPLX(0.0_DP,sqrt3)
         char_mat_proj(5,16)=(0.0_DP,-2.0_DP)
         char_mat_proj(5,17)=CMPLX(0.0_DP,-sqrt3)
         char_mat_proj(5,18)=(0.0_DP,1.0_DP)

         name_rap(6)='G_8-y   E_52-y'
         char_mat_proj(6,1)=(2.0_DP,0.0_DP)
         char_mat_proj(6,2)=CMPLX(-sqrt3,0.0_DP)
         char_mat_proj(6,3)=(1.0_DP,0.0_DP)
         char_mat_proj(6,5)=(1.0_DP,0.0_DP)
         char_mat_proj(6,6)=CMPLX(-sqrt3,0.0_DP)
         char_mat_proj(6,14)=(0.0_DP,1.0_DP)
         char_mat_proj(6,15)=CMPLX(0.0_DP,-sqrt3)
         char_mat_proj(6,16)=(0.0_DP,2.0_DP)
         char_mat_proj(6,17)=CMPLX(0.0_DP,sqrt3)
         char_mat_proj(6,18)=(0.0_DP,-1.0_DP)
   
      ELSEIF (ptype(1)==1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN
!
!   projective representations (beta=-1) point group
!
         nrap_proj=3

         name_rap(1)='A_1B_1'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,3)=(2.0_DP,0.0_DP)
         char_mat_proj(1,5)=(2.0_DP,0.0_DP)
         char_mat_proj(1,7)=(2.0_DP,0.0_DP)
         char_mat_proj(1,9)=(2.0_DP,0.0_DP)
         char_mat_proj(1,11)=(2.0_DP,0.0_DP)

         name_rap(2)='A_2B_2'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,3)=(2.0_DP,0.0_DP)
         char_mat_proj(2,5)=(2.0_DP,0.0_DP)
         char_mat_proj(2,7)=(-2.0_DP,0.0_DP)
         char_mat_proj(2,9)=(-2.0_DP,0.0_DP)
         char_mat_proj(2,11)=(-2.0_DP,0.0_DP)

         name_rap(3)='E_1E_2'
         char_mat_proj(3,1)=(4.0_DP,0.0_DP)
         char_mat_proj(3,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,5)=(-2.0_DP,0.0_DP)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN
!
!   projective representations (beta=-1) double point group
!
         nrap_proj=3

         name_rap(1)='G_9z    E_32z'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(1,5)=(-2.0_DP,0.0_DP)

         char_mat_proj(1,19)=(0.0_DP,2.0_DP)
         char_mat_proj(1,21)=(0.0_DP,-2.0_DP)
         char_mat_proj(1,23)=(0.0_DP,-2.0_DP)

         name_rap(2)='G_9-z   E_32-z'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(2,5)=(-2.0_DP,0.0_DP)
         char_mat_proj(2,19)=(0.0_DP,-2.0_DP)
         char_mat_proj(2,21)=(0.0_DP,2.0_DP)
         char_mat_proj(2,23)=(0.0_DP,2.0_DP)

         name_rap(3)='G_7G_8  E_12E_52'
         char_mat_proj(3,1)=(4.0_DP,0.0_DP)
         char_mat_proj(3,3)=(2.0_DP,0.0_DP)
         char_mat_proj(3,5)=(2.0_DP,0.0_DP)
  
      ELSEIF (ptype(1)==1.AND.ptype(2)==-1.AND.ptype(3)==-1) THEN
!
!   projective representations (beta=-1, gamma=-1) point group
!
         nrap_proj=3

         name_rap(1)='A_1B_2'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,3)=(2.0_DP,0.0_DP)
         char_mat_proj(1,5)=(2.0_DP,0.0_DP)
         char_mat_proj(1,8)=(2.0_DP,0.0_DP)
         char_mat_proj(1,10)=(2.0_DP,0.0_DP)
         char_mat_proj(1,12)=(2.0_DP,0.0_DP)

         name_rap(2)='A_2B_1'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,3)=(2.0_DP,0.0_DP)
         char_mat_proj(2,5)=(2.0_DP,0.0_DP)
         char_mat_proj(2,8)=(-2.0_DP,0.0_DP)
         char_mat_proj(2,10)=(-2.0_DP,0.0_DP)
         char_mat_proj(2,12)=(-2.0_DP,0.0_DP)

         name_rap(3)='E_1E_2'
         char_mat_proj(3,1)=(4.0_DP,0.0_DP)
         char_mat_proj(3,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(3,5)=(-2.0_DP,0.0_DP)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==-1.AND.ptype(3)==-1) THEN
!
!   projective representations (beta=-1, gamma=-1) double point group
!
         nrap_proj=3

         name_rap(1)='G_9x    E_32x'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(1,5)=(-2.0_DP,0.0_DP)
         char_mat_proj(1,20)=(0.0_DP,2.0_DP)
         char_mat_proj(1,22)=(0.0_DP,-2.0_DP)
         char_mat_proj(1,24)=(0.0_DP,-2.0_DP)

         name_rap(2)='G_9-x   E_32-x'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,3)=(-2.0_DP,0.0_DP)
         char_mat_proj(2,5)=(-2.0_DP,0.0_DP)
         char_mat_proj(2,20)=(0.0_DP,-2.0_DP)
         char_mat_proj(2,22)=(0.0_DP,2.0_DP)
         char_mat_proj(2,24)=(0.0_DP,2.0_DP)

         name_rap(3)='G_7G_8  E_12E_52'
         char_mat_proj(3,1)=(4.0_DP,0.0_DP)
         char_mat_proj(3,3)=(2.0_DP,0.0_DP)
         char_mat_proj(3,5)=(2.0_DP,0.0_DP)

      END IF

   CASE (132)
!
!  T
!
      w=CMPLX(-0.5_DP,-sqrt3/2.0_DP)

      nsym_proj=12
      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations 
!
         nrap_proj=4

         name_rap(1)='A_1'
         char_mat_proj(1,1:12)=(1.0_DP,0.0_DP)

         name_rap(2)='A'''
         char_mat_proj(2,1:4)=(1.0_DP,0.0_DP)
         char_mat_proj(2,5:8)= w
         char_mat_proj(2,9:12)= w ** 2

         name_rap(3)='A'''''
         char_mat_proj(3,1:4)=(1.0_DP,0.0_DP)
         char_mat_proj(3,5:8)= w**2
         char_mat_proj(3,9:12)= w 

         name_rap(4)='T'
         char_mat_proj(4,1)=(3.0_DP,0.0_DP)
         char_mat_proj(4,2:4)=(-1.0_DP,0.0_DP) 

      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations 
!
         nrap_proj=3

         name_rap(1)='G_5     E12'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,5:8)=(1.0_DP,0.0_DP) 
         char_mat_proj(1,9:12)=(1.0_DP,0.0_DP) 

         name_rap(2)='G_6     E''''12'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,5:8)= w ** 2
         char_mat_proj(2,9:12)= CONJG(w**2)

         name_rap(3)='G_7     E''12'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,5:8)= w
         char_mat_proj(3,9:12)= CONJG(w)

      ENDIF

   CASE (133)
!
!  T_h
!
      w=CMPLX(-0.5_DP,-sqrt3/2.0_DP)

      nsym_proj=24
      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations 
!
         nrap_proj=8

         name_rap(1)='A_1g'
         char_mat_proj(1,1:12)=(1.0_DP,0.0_DP)
         char_mat_proj(1,13:24)=char_mat_proj(1,1:12)

         name_rap(2)='A''g'
         char_mat_proj(2,1:4)=(1.0_DP,0.0_DP)
         char_mat_proj(2,5:8)= w
         char_mat_proj(2,9:12)= w ** 2
         char_mat_proj(2,13:24)=char_mat_proj(2,1:12)

         name_rap(3)='A''''g'
         char_mat_proj(3,1:4)=(1.0_DP,0.0_DP)
         char_mat_proj(3,5:8)= w**2
         char_mat_proj(3,9:12)= w 
         char_mat_proj(3,13:24)=char_mat_proj(3,1:12)

         name_rap(4)='Tg'
         char_mat_proj(4,1)=(3.0_DP,0.0_DP)
         char_mat_proj(4,2:4)=(-1.0_DP,0.0_DP) 
         char_mat_proj(4,13:24)=char_mat_proj(4,1:12)

         name_rap(5)='A_1u'
         char_mat_proj(5,1:12)=(1.0_DP,0.0_DP)
         char_mat_proj(5,13:24)=-char_mat_proj(5,1:12)

         name_rap(6)='A''u'
         char_mat_proj(6,1:4)=(1.0_DP,0.0_DP)
         char_mat_proj(6,5:8)= w
         char_mat_proj(6,9:12)= w ** 2
         char_mat_proj(6,13:24)=-char_mat_proj(6,1:12)

         name_rap(7)='A''''u'
         char_mat_proj(7,1:4)=(1.0_DP,0.0_DP)
         char_mat_proj(7,5:8)= w**2
         char_mat_proj(7,9:12)= w 
         char_mat_proj(7,13:24)=-char_mat_proj(7,1:12)

         name_rap(8)='Tu'
         char_mat_proj(8,1)=(3.0_DP,0.0_DP)
         char_mat_proj(8,2:4)=(-1.0_DP,0.0_DP) 
         char_mat_proj(8,13:24)=-char_mat_proj(8,1:12)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations 
!
         nrap_proj=6

         name_rap(1)='G_5+    E12g'
         char_mat_proj(1,1)=(2.0_DP,0.0_DP)
         char_mat_proj(1,5:8)=(1.0_DP,0.0_DP) 
         char_mat_proj(1,9:12)=(1.0_DP,0.0_DP) 
         char_mat_proj(1,13:24)=char_mat_proj(1,1:12)

         name_rap(2)='G_6+    E''''g'
         char_mat_proj(2,1)=(2.0_DP,0.0_DP)
         char_mat_proj(2,5:8)= w ** 2
         char_mat_proj(2,9:12)= CONJG(w**2)
         char_mat_proj(2,13:24)=char_mat_proj(2,1:12)

         name_rap(3)='G_7+    E''g'
         char_mat_proj(3,1)=(2.0_DP,0.0_DP)
         char_mat_proj(3,5:8)= w
         char_mat_proj(3,9:12)= CONJG(w)
         char_mat_proj(3,13:24)=char_mat_proj(3,1:12)

         name_rap(4)='G_5-    E12u'
         char_mat_proj(4,1)=(2.0_DP,0.0_DP)
         char_mat_proj(4,5:8)=(1.0_DP,0.0_DP) 
         char_mat_proj(4,9:12)=(1.0_DP,0.0_DP) 
         char_mat_proj(4,13:24)=-char_mat_proj(4,1:12)

         name_rap(5)='G_6-    E''''u'
         char_mat_proj(5,1)=(2.0_DP,0.0_DP)
         char_mat_proj(5,5:8)= w ** 2
         char_mat_proj(5,9:12)= CONJG(w**2)
         char_mat_proj(5,13:24)=-char_mat_proj(5,1:12)

         name_rap(6)='G_7-    E''u'
         char_mat_proj(6,1)=(2.0_DP,0.0_DP)
         char_mat_proj(6,5:8)= w
         char_mat_proj(6,9:12)= CONJG(w)
         char_mat_proj(6,13:24)=-char_mat_proj(6,1:12)

      ENDIF

   CASE (134,135)
!
!  T_d or O
!
      nsym_proj=24

      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations 
!
         nrap_proj=5

         name_rap(1)='A_1'
         char_mat_proj(1,1:24)=(1.0_DP,0.0_DP)

         name_rap(2)='A_2'
         char_mat_proj(2,1) = (1.0_DP,0.0_DP)
         char_mat_proj(2,2:4) = (1.0_DP,0.0_DP)
         char_mat_proj(2,5:12) = (1.0_DP,0.0_DP)
         char_mat_proj(2,13:24) = (-1.0_DP,0.0_DP)

         name_rap(3)='E'
         char_mat_proj(3,1) = (2.0_DP,0.0_DP)
         char_mat_proj(3,2:4) = (2.0_DP,0.0_DP)
         char_mat_proj(3,5:12) = (-1.0_DP,0.0_DP)

         name_rap(4)='T_1'
         char_mat_proj(4,1) = (3.0_DP,0.0_DP)
         char_mat_proj(4,2:4) = (-1.0_DP,0.0_DP)
         char_mat_proj(4,13:18) = (1.0_DP,0.0_DP)
         char_mat_proj(4,19:24) = (-1.0_DP,0.0_DP)

         name_rap(5)='T_2'
         char_mat_proj(5,1) = (3.0_DP,0.0_DP)
         char_mat_proj(5,2:4) = (-1.0_DP,0.0_DP)
         char_mat_proj(5,13:18) = (-1.0_DP,0.0_DP)
         char_mat_proj(5,19:24) = (1.0_DP,0.0_DP)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations 
!
         nrap_proj=3

         name_rap(1)='G_6     E_12'
         char_mat_proj(1,1) = (2.0_DP,0.0_DP)
         char_mat_proj(1,5:12) = (1.0_DP,0.0_DP)
         char_mat_proj(1,13:18) = CMPLX(sqrt2,0.0_DP)

         name_rap(2)='G_7     E''12'
         char_mat_proj(2,1) = (2.0_DP,0.0_DP)
         char_mat_proj(2,5:12) = (1.0_DP,0.0_DP)
         char_mat_proj(2,13:18) = CMPLX(-sqrt2,0.0_DP)

         name_rap(3)='G_8     Q'
         char_mat_proj(3,1) = (4.0_DP,0.0_DP)
         char_mat_proj(3,5:12) = (-1.0_DP,0.0_DP)

      ENDIF

   CASE (136)
!
!   O_h
!
      nsym_proj=48

      IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   point group representations 
!
         nrap_proj=10

         name_rap(1)='A_1g'
         char_mat_proj(1,1:24)=(1.0_DP,0.0_DP)
         char_mat_proj(1,25:48)=char_mat_proj(1,1:24)

         name_rap(2)='A_2g'
         char_mat_proj(2,1) = (1.0_DP,0.0_DP)
         char_mat_proj(2,2:4) = (1.0_DP,0.0_DP)
         char_mat_proj(2,5:12) = (1.0_DP,0.0_DP)
         char_mat_proj(2,13:24) = (-1.0_DP,0.0_DP)
         char_mat_proj(2,25:48)=char_mat_proj(2,1:24)

         name_rap(3)='Eg'
         char_mat_proj(3,1) = (2.0_DP,0.0_DP)
         char_mat_proj(3,2:4) = (2.0_DP,0.0_DP)
         char_mat_proj(3,5:12) = (-1.0_DP,0.0_DP)
         char_mat_proj(3,25:48)=char_mat_proj(3,1:24)

         name_rap(4)='T_1g'
         char_mat_proj(4,1) = (3.0_DP,0.0_DP)
         char_mat_proj(4,2:4) = (-1.0_DP,0.0_DP)
         char_mat_proj(4,13:18) = (1.0_DP,0.0_DP)
         char_mat_proj(4,19:24) = (-1.0_DP,0.0_DP)
         char_mat_proj(4,25:48)=char_mat_proj(4,1:24)

         name_rap(5)='T_2g'
         char_mat_proj(5,1) = (3.0_DP,0.0_DP)
         char_mat_proj(5,2:4) = (-1.0_DP,0.0_DP)
         char_mat_proj(5,13:18) = (-1.0_DP,0.0_DP)
         char_mat_proj(5,19:24) = (1.0_DP,0.0_DP)
         char_mat_proj(5,25:48)=char_mat_proj(5,1:24)

         name_rap(6)='A_1u'
         char_mat_proj(6,1:24)=(1.0_DP,0.0_DP)
         char_mat_proj(6,25:48)=-char_mat_proj(6,1:24)

         name_rap(7)='A_2u'
         char_mat_proj(7,1) = (1.0_DP,0.0_DP)
         char_mat_proj(7,2:4) = (1.0_DP,0.0_DP)
         char_mat_proj(7,5:12) = (1.0_DP,0.0_DP)
         char_mat_proj(7,13:24) = (-1.0_DP,0.0_DP)
         char_mat_proj(7,25:48)=-char_mat_proj(7,1:24)

         name_rap(8)='Eu'
         char_mat_proj(8,1) = (2.0_DP,0.0_DP)
         char_mat_proj(8,2:4) = (2.0_DP,0.0_DP)
         char_mat_proj(8,5:12) = (-1.0_DP,0.0_DP)
         char_mat_proj(8,25:48) = -char_mat_proj(8,1:24)

         name_rap(9)='T_1u'
         char_mat_proj(9,1) = (3.0_DP,0.0_DP)
         char_mat_proj(9,2:4) = (-1.0_DP,0.0_DP)
         char_mat_proj(9,13:18) = (1.0_DP,0.0_DP)
         char_mat_proj(9,19:24) = (-1.0_DP,0.0_DP)
         char_mat_proj(9,25:48) = -char_mat_proj(9,1:24)

         name_rap(10)='T_2u'
         char_mat_proj(10,1) = (3.0_DP,0.0_DP)
         char_mat_proj(10,2:4) = (-1.0_DP,0.0_DP)
         char_mat_proj(10,13:18) = (-1.0_DP,0.0_DP)
         char_mat_proj(10,19:24) = (1.0_DP,0.0_DP)
         char_mat_proj(10,25:48) = -char_mat_proj(10,1:24)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
!
!   double point group representations 
!
         nrap_proj=6

         name_rap(1)='G_6+    E_12g'
         char_mat_proj(1,1) = (2.0_DP,0.0_DP)
         char_mat_proj(1,5:12) = (1.0_DP,0.0_DP)
         char_mat_proj(1,13:18) = CMPLX(sqrt2,0.0_DP)
         char_mat_proj(1,25:48)=char_mat_proj(1,1:24)

         name_rap(2)='G_7+    E''12g'
         char_mat_proj(2,1) = (2.0_DP,0.0_DP)
         char_mat_proj(2,5:12) = (1.0_DP,0.0_DP)
         char_mat_proj(2,13:18) = CMPLX(-sqrt2,0.0_DP)
         char_mat_proj(2,25:48)=char_mat_proj(2,1:24)

         name_rap(3)='G_8+    Qg'
         char_mat_proj(3,1) = (4.0_DP,0.0_DP)
         char_mat_proj(3,5:12) = (-1.0_DP,0.0_DP)
         char_mat_proj(3,25:48)=char_mat_proj(3,1:24)

         name_rap(4)='G_6-    E_12u'
         char_mat_proj(4,1) = (2.0_DP,0.0_DP)
         char_mat_proj(4,5:12) = (1.0_DP,0.0_DP)
         char_mat_proj(4,13:18) = CMPLX(sqrt2,0.0_DP)
         char_mat_proj(4,25:48) = -char_mat_proj(4,1:24)
 
         name_rap(5)='G_7-    E''12u'
         char_mat_proj(5,1) = (2.0_DP,0.0_DP)
         char_mat_proj(5,5:12) = (1.0_DP,0.0_DP)
         char_mat_proj(5,13:18) = CMPLX(-sqrt2,0.0_DP)
         char_mat_proj(5,25:48) = -char_mat_proj(5,1:24)

         name_rap(6)='G_8-     Qu'
         char_mat_proj(6,1) = (4.0_DP,0.0_DP)
         char_mat_proj(6,5:12) = (-1.0_DP,0.0_DP)
         char_mat_proj(6,25:48) = -char_mat_proj(6,1:24)

      ELSEIF (ptype(1)==1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN
!
!   projective representations (beta=-1) point group
!
         nrap_proj=4

         name_rap(1)='A_1A_2'
         char_mat_proj(1,1) = (2.0_DP,0.0_DP)
         char_mat_proj(1,2:4) = (2.0_DP,0.0_DP)
         char_mat_proj(1,5:12) = (2.0_DP,0.0_DP)

         name_rap(2)='Ey'
         char_mat_proj(2,1) = (2.0_DP,0.0_DP)
         char_mat_proj(2,2:4) = (2.0_DP,0.0_DP)
         char_mat_proj(2,5:12) = (-1.0_DP,0.0_DP)
         char_mat_proj(2,29:32) = CMPLX(0.0_DP,-sqrt3)
         char_mat_proj(2,33:36) = CMPLX(0.0_DP,sqrt3)

         name_rap(3)='E-y'
         char_mat_proj(3,1) = (2.0_DP,0.0_DP)
         char_mat_proj(3,2:4) = (2.0_DP,0.0_DP)
         char_mat_proj(3,5:12) = (-1.0_DP,0.0_DP)
         char_mat_proj(3,29:32) = CMPLX(0.0_DP,sqrt3)
         char_mat_proj(3,33:36) = CMPLX(0.0_DP,-sqrt3)

         name_rap(4)='T_1T_2'
         char_mat_proj(4,1) = (6.0_DP,0.0_DP)
         char_mat_proj(4,2:4) = (-2.0_DP,0.0_DP)

      ELSEIF (ptype(1)==-1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN
!
!   projective representations (beta=-1) double point group
!
         nrap_proj=3

         name_rap(1)='G_6G_7  E12E12'
         char_mat_proj(1,1) = (4.0_DP,0.0_DP)
         char_mat_proj(1,5:12) = (2.0_DP,0.0_DP)

         name_rap(2)='G_8y    Qy'
         char_mat_proj(2,1) = (4.0_DP,0.0_DP)
         char_mat_proj(2,5:12) = (-1.0_DP,0.0_DP)
         char_mat_proj(2,29:32) = CMPLX(0.0_DP,sqrt3)
         char_mat_proj(2,33:36) = CMPLX(0.0_DP,-sqrt3)

         name_rap(3)='G_8-y   Q-y'
         char_mat_proj(3,1) = (4.0_DP,0.0_DP)
         char_mat_proj(3,5:12) = (-1.0_DP,0.0_DP)
         char_mat_proj(3,29:32) = CMPLX(0.0_DP,-sqrt3)
         char_mat_proj(3,33:36) = CMPLX(0.0_DP,sqrt3)

      ENDIF

   CASE DEFAULT
      CALL errore('set_stand_irr_proj','group representations not implemented',1)
END SELECT

RETURN
END SUBROUTINE set_stand_irr_proj

SUBROUTINE  transform_s_to_cart( sk, sr, nsym, at, bg )
  !----------------------------------------------------------------------
  !
  !     This routine transforms symmetry matrices expressed in the
  !     basis of the crystal axis into rotations in cartesian axis
  !     A similar routine is already in symm_base, but has not arguments.
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

SUBROUTINE find_projection_type(code_group, code_group_ext, argument,  &
                                     ptype, gauge, verbose)
!
!  The inputs of the routine are the codes of the point group, standard and
!  extended, and a set of arguments for the phases of the factor system.
!  The outputs are the ptype array (see below) and the gauge factors 
!  (arguments of the phases) to be applied to each group element 
!  to bring the input phase factors into the standard one through a 
!  p-equivalence (dividing the standard characters by these gauges,
!  they become the characters of the matrices that correspond to the input 
!  factor system). The arguments of the factor systems must be ordered
!  with the standard order of the point groups.
!
!  ptype is an array of three integers +-1:
!  -E beta gamma where the first number says if the single (1) or
!   double (-1) group representations are used. 
!   beta and gamma correspond to the factors of
!   IAI = beta A and IBI = gamma B, where I is the inversion and A and B
!   are the group generators.
!
!  The tables of the irreducible representations provided by 
!  set_irr_proj are written for the standard factor system.
!
!  Set verbose to .TRUE. if you want to print detailed information of
!  the calculations of this routine. 

  USE kinds, ONLY : DP
  USE constants, ONLY : tpi, pi
  USE io_global, ONLY : stdout
  IMPLICIT NONE

  INTEGER,  INTENT(IN) :: code_group, code_group_ext
  REAL(DP), INTENT(IN) :: argument(48,48)
  LOGICAL,  INTENT(IN) :: verbose  ! if true write a lot of messages on output
  INTEGER,  INTENT(OUT) :: ptype(3)
  REAL(DP), INTENT(INOUT) :: gauge(48)

  INTEGER :: isym, jsym, nsym, tau, kappa, i
  INTEGER :: group_desc(48), epos(48,48), prd(48,48), tab(48,48)
  INTEGER :: row(3), column(3), n(3), columni(3)
  COMPLEX(DP) :: pha
  REAL(DP) :: gauge_a, gauge_b, gauge_ap, gauge_bp, gauge_as, gauge_bs, &
              gauge_i, gauge_e
  REAL(DP) :: arg_alpha, arg_beta, arg_gmma, arg_gmma1, arg_eta,  &
              arg_alphai, arg_betai, arge1, arg_i
  REAL(DP) :: arg
  COMPLEX(DP) :: factor(48,48), phase(48)
  CHARACTER(LEN=11) :: group_name
  LOGICAL :: linv
!
!  set up the group elements and optionaly write the input factor system
!
  CALL set_group_desc(group_desc,nsym,code_group_ext)

  IF (verbose) THEN
     WRITE(stdout,'(/,5x,"Find projection type, group  ", a)') &
                                                   TRIM(group_name(code_group))
     DO isym=1,nsym
        DO jsym=1,nsym
           arg=argument(isym,jsym)
           factor(isym,jsym)=CMPLX( COS(arg), SIN(arg) )   
        ENDDO
     ENDDO

     WRITE(stdout,'(/,5x, "The input factor system of the point group:")')
     CALL write_group_table(group_desc, nsym, factor)

     WRITE(stdout,'(/,5x, "A^n    = l1 E''")') 
     WRITE(stdout,'(5x,   "B^m    = l2 E''")') 
     WRITE(stdout,'(5x,   "(AB)^2 = l3 E''")') 
     WRITE(stdout,'(5x,   "(E'')^2 = eta E")') 
  ENDIF

  CALL group_generators(code_group, row, column, n, linv, columni)

  IF (code_group < 28) tau=n(1)
  IF (code_group == 28 .OR. code_group == 29 ) tau = 6
  IF (code_group == 30 .OR. code_group == 31 .OR. code_group == 32 ) tau = 12
!
!  From the arguments find l1, l2, l3 and find tentative phases for the gauge
!
  arg_beta=0.0_DP
  arg_gmma=0.0_DP
  arg_alpha=0.0_DP
  IF (n(1)==2.AND.n(2)==1) THEN
     arg_alpha=argument(2,2) 
  ELSEIF (n(1)==3.AND.n(2)==1) THEN
     arg_alpha=argument(2,2) + argument(3,2)
  ELSEIF (n(1)==4.AND.n(2)==1) THEN
     arg_alpha=argument(2,2) * 2.0_DP + argument(3,3)
  ELSEIF (n(1)==6.AND.n(2)==1) THEN
     arg_alpha=argument(2,2) * 3.0_DP + argument(3,3) + argument(5,3)
  ELSEIF (n(1)==2.AND.n(2)==2) THEN
     arg_alpha=argument(2,2) 
     arg_beta=argument(row(2),column(2)) 
     arg_gmma= argument(column(1),column(2))*2.0_DP + argument(row(3),column(3))
  ELSEIF (n(1)==3.AND.n(2)==2) THEN
     arg_alpha=argument(2,2) + argument(3,2)
     arg_beta=argument(row(2),column(2)) 
     arg_gmma= argument(column(1),column(2))*2.0_DP + argument(row(3),column(3))
  ELSEIF (n(1)==4.AND.n(2)==2) THEN
     arg_alpha=argument(2,2) * 2.0_DP + argument(3,3)
     arg_beta=argument(row(2),column(2)) 
     arg_gmma= argument(column(1),column(2))*2.0_DP + argument(row(3),column(3))
  ELSEIF (n(1)==6.AND.n(2)==2) THEN
     arg_alpha=argument(2,2) * 3.0_DP + argument(3,3) + argument(5,3)
     arg_beta=argument(row(2),column(2)) 
     arg_gmma= argument(column(1),column(2))*2.0_DP + argument(row(3),column(3))
  ELSEIF (n(1)==4.AND.n(2)==3) THEN
     arg_alpha=argument(18,18) * 2.0_DP + argument(4,4)
     arg_beta= argument(5,5) + argument(9,5)
     arg_gmma= argument(column(1),column(2))*2.0_DP + argument(row(3),column(3))
  ELSEIF (n(1)==3.AND.n(2)==3) THEN
     arg_alpha=argument(5,5) + argument(5,9)
     arg_beta=argument(10,10) + argument(10,6)
     arg_gmma= argument(column(1),column(2))*2.0_DP + argument(row(3),column(3))
  ENDIF   

  CALL zero_tpi(arg_alpha)
  CALL zero_tpi(arg_beta)
  CALL zero_tpi(arg_gmma)

  IF (verbose) THEN
     WRITE(stdout,'(/,5x,"Input l1    phi= ",f11.6, 2x, " l1 =",2f14.6)')  &
               arg_alpha * 180._DP / pi, COS(arg_alpha), SIN(arg_alpha)
     WRITE(stdout,'(5x,"Input l2    phi= ",f11.6, 2x, " l2 =",2f14.6)')  &
               arg_beta * 180._DP / pi, COS(arg_beta), SIN(arg_beta)
     WRITE(stdout,'(5x,"Input l3    phi= ",f11.6, 2x, " l3 =",2f14.6)')  &
               arg_gmma * 180._DP / pi, COS(arg_gmma), SIN(arg_gmma)
  ENDIF

  ptype=1
  arg=arg_alpha / DBLE ( n(1) )
  gauge_ap=-arg
  
  arg=arg_beta/ DBLE( n(2) )
  gauge_bp=-arg

  arg_gmma1=arg_gmma-arg_alpha*2.0_DP/DBLE(n(1))-arg_beta*2.0_DP/ DBLE( n(2) )

  IF (verbose) THEN
     WRITE(stdout,'(/,5x,"Tentative phase_a", 4x, f12.6, 3x, 2f16.6)') &
                       gauge_ap * 180.0_DP/ pi, COS(gauge_ap), SIN(gauge_ap)
     WRITE(stdout,'(5x,"Tentative phase_b", 4x, f12.6, 3x, 2f16.6)') &
                       gauge_bp * 180.0_DP/ pi, COS(gauge_bp), SIN(gauge_bp)
     WRITE(stdout,'(/,5x,"new l3      ",f12.6, 3x, 2f16.6)') &
                      arg_gmma1 * 180._DP/ pi, COS(arg_gmma1), SIN(arg_gmma1)
  ENDIF
  !
  ! Determine if a switch between point group and double point group
  ! is necessary and adjust the phases
  !
  kappa = - NINT(arg_gmma1 * DBLE(tau) / tpi)
  kappa = MOD (kappa, tau)
  IF (kappa < 0 ) kappa=kappa+tau

  gauge_e=0.0_DP
  arge1 = pi * kappa 
 
  IF (verbose) THEN
     WRITE(stdout,'(5x,"kappa= ",i5)') kappa
     WRITE(stdout,'(5x,"D''(E'') = eta D(E'') where eta ", f11.6, 3x, &
                             &2f13.6)') arge1 * 180._DP/ pi, COS(arge1), &
                                                             SIN(arge1)
  END IF

  IF (ABS(SIN(arge1)) > 1.D-6 ) &
              CALL errore('find_projection_type','problem with E''',1)

  IF (ABS(COS(arge1)+1.0_DP) < 1.D-6) THEN
!     IF (n(1)==3.AND.n(2)==2) THEN
!
!   D_3 or C_3v can always be represented with the point group operations
!
!        WRITE(stdout,'(5x,"Odd order group, double group exchange is not &
!                                             &necessary")') 
!        gauge_ap = gauge_ap + pi
!        gauge_bp = gauge_bp + pi / 2.0_DP
!     ELSE
        ptype(1)=-ptype(1)
        gauge_e=0.0_DP
        IF (ptype(1)==-1) gauge_e=pi
!     ENDIF
  ELSEIF (ABS(COS(arge1)-1.0_DP) > 1.D-6) THEN
     IF (ABS(SIN(arge1)) > 1.D-6 ) &
              CALL errore('find_projection_type','problem 1 with E''',1)
  END IF

  arg= arge1 / DBLE(n(1))
  gauge_as = arg
  gauge_a = gauge_as + gauge_ap

  arg=arge1 / DBLE(n(2))
  gauge_bs = arg
  gauge_b = gauge_bp + gauge_bs

  IF (verbose) THEN
     WRITE(stdout,'(/,5x,"Second part phase_a", 2x, f12.6, 3x, 2f16.6)') &
                          gauge_as * 180.0_DP/ pi, COS(gauge_as), SIN(gauge_as)
     WRITE(stdout,'(5x,"Final phase_a", 23x, 2f16.6)') COS(gauge_a), &
                                                       SIN(gauge_a)
     WRITE(stdout,'(/,5x,"Second part phase_b", 2x, f12.6, 3x, 2f16.6)') &
                          gauge_bs * 180.0_DP/ pi, COS(gauge_bs), SIN(gauge_bs)
     WRITE(stdout,'(5x,"Final phase_b", 23x, 2f16.6,/)') COS(gauge_b), &
                                                         SIN(gauge_b)
  ENDIF

  IF (linv) THEN
!
!   the group has inversion. Determine here the projective type of the
!   input factor system and possibly find the gauge for inversion.
!
     IF (verbose) THEN
        WRITE(stdout,'(/,5x,"The group has inversion")')

        WRITE(stdout,'(/,5x, "I^2    = aii I")') 
        WRITE(stdout,'(5x,   "IAI    = betai A")') 
        WRITE(stdout,'(5x,   "IBI    = gammai B")') 
     ENDIF

     arg_i=argument(columni(1),columni(1))
     arg_alphai=argument(columni(1),column(1))+argument(columni(2),columni(1))  
     IF (n(2)>1) THEN
        arg_betai=argument(columni(1),column(2))+ &
                          argument(columni(3),columni(1))    
     ELSE
        arg_betai=0.0_DP
     ENDIF
     
     CALL zero_tpi(arg_i)
     CALL zero_tpi(arg_alphai)
     CALL zero_tpi(arg_betai)
 
     IF (verbose) THEN
        WRITE(stdout,'(/,5x,"Input i        phi= ",f11.6, 2x, " aii =   ",&
                         2f14.6)')  arg_i*180._DP/pi, COS(arg_i), SIN(arg_i)
        WRITE(stdout,'(5x,"Input betai    phi= ",f11.6, 2x, " betai  =",&
                         2f14.6)')  arg_alphai*180._DP/pi, COS(arg_alphai), &
                                                           SIN(arg_alphai)
        WRITE(stdout,'(5x,"Input gammai   phi= ",f11.6, 2x, " gammai =",&
                         2f14.6)')  arg_betai*180._DP/pi, COS(arg_betai), &
                                                          SIN(arg_betai)
     END IF

     arg=arg_i / 2.0_DP
     gauge_i = -arg

     arg_alphai = arg_alphai - arg_i 
     IF (n(2)>1) THEN
        arg_betai = arg_betai - arg_i 
     ELSE
        arg_betai = 0.0_DP
     ENDIF

     IF (verbose) THEN
        WRITE(stdout,'(/,5x,"Final phase_i", 29x, 2f14.6)') COS(gauge_i), &
                                                            SIN(gauge_i)
        WRITE(stdout,'(5x,"Output betai   phi= ",f11.6, 2x, " betai  =",&
                         2f14.6)')  arg_alphai*180._DP/pi, COS(arg_alphai), &
                                                           SIN(arg_alphai)
        WRITE(stdout,'(5x,"Output gammai  phi= ",f11.6, 2x, " gammai =",&
                         2f14.6)')  arg_betai*180._DP/pi, COS(arg_betai), &
                                                           SIN(arg_betai)
     END IF
     IF (ABS(COS(arg_alphai)+1.0_DP)<1.D-6) ptype(2)=-1
     IF (ABS(COS(arg_betai)+1.0_DP)<1.D-6)  ptype(3)=-1
     !
     !   Check that inversion is ok
     !
     IF ( ABS(SIN(arg_alphai))>1.D-6 .OR. ABS(SIN(arg_betai)) > 1.D-6) &
         CALL errore('find_projection_type','problem with inversion',1)
     IF ( ABS(COS(arg_alphai)-1.0_DP)>1.D-6 .AND. & 
          ABS(COS(arg_alphai)+1.0_DP)>1.D-6 .AND. & 
          ABS(COS(arg_betai)-1.0_DP)>1.D-6 .AND. & 
          ABS(COS(arg_betai)+1.0_DP)>1.D-6 ) &
             CALL errore('find_projection_type','problem with inversion1',1)

  END IF
!
!  Now compute the actual gauge factors using the group relationships
!
  gauge(1)=0.0_DP
  SELECT CASE (code_group)
    CASE (1)
    CASE (2)
!
! C_i
!
        gauge(2)=gauge_i

     CASE(3,4)
!
!  C_2, C_s
!
        gauge(2)=gauge_a

     CASE(5)
!
!  C_3
!
        gauge(2)=gauge_a
        gauge(3) = -argument(2,3) - gauge_a 

     CASE(6,26)
!
!  C_4 or S_4
!
        gauge(2) = gauge_a
        gauge(3) = argument(2,2) + 2.0_DP * gauge_a 
        gauge(4) = -argument(2,4) - gauge_a 

     CASE(7,17)
!
!  C_6 or C_3h 
!
        gauge(2) = gauge_a
        gauge(3) = argument(2,2) + 2.0_DP * gauge_a
        gauge(4) = argument(2,3) + gauge(3) + gauge_a
        gauge(5) = -argument(3,5) - gauge(3)
        gauge(6) = -argument(2,6) - gauge_a

     CASE(8,12)
!
!   D_2 or C_2v
!
        gauge(2) = gauge_a
        gauge(3) = gauge_b
        gauge(4) = argument(2,3) + gauge_a + gauge_b

     CASE(9,13)
!
!   D_3 or C_3v
!
        gauge(2) = gauge_a 
        gauge(3) = -argument(2,3) - gauge_a 
        gauge(4) = gauge_b 
        gauge(5) = argument(2,4) + gauge(2) + gauge_b 
        gauge(6) = argument(3,4) + gauge(3) + gauge_b 

     CASE(10,14,24)
!
!   D_4 or C_4v or D_2d
!
        gauge(2) = gauge_a
        gauge(3) = argument(2,2) + 2.0_DP * gauge_a 
        gauge(4) = - argument(2,4) - gauge_a
        gauge(5) = gauge_b
        gauge(6) = argument(2,5) + gauge_a + gauge_b
        gauge(7) = argument(3,5) + gauge(3) + gauge_b
        gauge(8) = argument(4,5) + gauge(4) + gauge_b

     CASE(11,15,21)
!
!   D_6 or C_6v or D_3h
!
        gauge(2) = gauge_a
        gauge(3) = argument(2,2) + 2.0_DP * gauge_a 
        gauge(4) = argument(2,3) + gauge(3) + gauge_a
        gauge(5) = -argument(3,5) - gauge(3)
        gauge(6) = -argument(2,6) - gauge_a
        gauge(7) = gauge_b
        gauge(8) = argument(2,7) + gauge_a + gauge_b
        gauge(9) = argument(3,7) + gauge(3) + gauge_b
        gauge(10) = argument(4,7) + gauge(4) + gauge_b
        gauge(11) = argument(5,7) + gauge(5) + gauge_b
        gauge(12) = argument(6,7) + gauge(6) + gauge_b

     CASE(16)
!
!   C_2h
!
        gauge(2) = gauge_a
        gauge(3) = gauge_i
        gauge(4) = argument(3,2) + gauge_a + gauge_i

     CASE(18)
!
!   C_4h
!
        gauge(2) = gauge_a
        gauge(3) = argument(2,2) + 2.0_DP * gauge_a 
        gauge(4) = -argument(2,4) - gauge_a
        gauge(5) = gauge_i
        gauge(6) = argument(5,2) + gauge_a + gauge_i
        gauge(7) = argument(5,3) + gauge(3) + gauge_i
        gauge(8) = argument(5,4) + gauge(4) + gauge_i


     CASE(19)
!
!   C_6h
!
        gauge(2)= gauge_a
        gauge(3) = argument(2,2) + 2.0_DP * gauge_a 
        gauge(4) = argument(2,3) + gauge(3) + gauge_a
        gauge(5) = - argument(3,5) - gauge(3)
        gauge(6) = - argument(2,6) - gauge_a
        gauge(7) = gauge_i
        gauge(8) = argument(7,2) + gauge_a + gauge_i
        gauge(9) = argument(7,3) + gauge(3) + gauge_i
        gauge(10) = argument(7,4) + gauge(4) + gauge_i
        gauge(11) = argument(7,5) + gauge(5) + gauge_i
        gauge(12) = argument(7,6) + gauge(6) + gauge_i

     CASE(20)
!
!    D_2h
!
        gauge(2) = gauge_a
        gauge(3) = gauge_b
        gauge(4) = argument(2,3) + gauge_a + gauge_b
        gauge(5) = gauge_i
        gauge(6) = argument(5,2) + gauge_a + gauge_i
        gauge(7) = argument(5,3) + gauge_b + gauge_i
        gauge(8) = argument(5,4) + gauge(4) + gauge_i 

     CASE(22)
!
!   D_4h
!
        gauge(2) = gauge_a
        gauge(3) = argument(2,2) + 2.0_DP * gauge_a 
        gauge(4) = - argument(2,4) - gauge_a 
        gauge(5) = gauge_b
        gauge(6) = argument(2,5) + gauge_a + gauge_b
        gauge(7) = argument(3,5) + gauge(3) + gauge_b
        gauge(8) = argument(4,5) + gauge(4) + gauge_b
        gauge(9) = gauge_i
        gauge(10) = argument(9,2) + gauge_a + gauge_i
        gauge(11) = argument(9,3) + gauge(3) + gauge_i
        gauge(12) = argument(9,4) + gauge(4) + gauge_i 
        gauge(13) = argument(9,5) + gauge_b + gauge_i 
        gauge(14) = argument(9,6) + gauge(6) + gauge_i 
        gauge(15) = argument(9,7) + gauge(7) + gauge_i 
        gauge(16) = argument(9,8) + gauge(8) + gauge_i 

     CASE(23)
!
!   D_6h
!
        gauge(2) = gauge_a
        gauge(3) = argument(2,2) + 2.0_DP* gauge_a 
        gauge(4) = argument(2,3) + gauge_a + gauge(3)
        gauge(5) = -argument(3,5) - gauge(3) 
        gauge(6) = -argument(2,6) - gauge_a
        gauge(7) = gauge_b
        gauge(8) = argument(2,7) + gauge_a + gauge_b
        gauge(9) = argument(3,7) + gauge(3) + gauge_b
        gauge(10) = argument(4,7) + gauge(4) + gauge_b 
        gauge(11) = argument(5,7) + gauge(5) + gauge_b
        gauge(12) = argument(6,7) + gauge(6) + gauge_b 
        gauge(13) = gauge_i
        gauge(14) = argument(13,2) + gauge_a + gauge_i
        gauge(15) = argument(13,3) + gauge(3) + gauge_i
        gauge(16) = argument(13,4) + gauge(4) + gauge_i
        gauge(17) = argument(13,5) + gauge(5) + gauge_i
        gauge(18) = argument(13,6) + gauge(6) + gauge_i
        gauge(19) = argument(13,7) + gauge_b + gauge_i
        gauge(20) = argument(13,8) + gauge(8) + gauge_i
        gauge(21) = argument(13,9) + gauge(9) + gauge_i
        gauge(22) = argument(13,10) + gauge(10) + gauge_i
        gauge(23) = argument(13,11) + gauge(11) + gauge_i
        gauge(24) = argument(13,12) + gauge(12) + gauge_i

     CASE(25)
!
!   D_3d
!
        gauge(2) = gauge_a
        gauge(3) = - argument(2,3) - gauge_a 
        gauge(4) = gauge_b 
        gauge(5) = argument(2,4) + gauge_a + gauge_b
        gauge(6) = argument(3,4) + gauge(3) + gauge_b
        gauge(7) = gauge_i
        gauge(8) = argument(7,2) + gauge_a + gauge_i
        gauge(9) = argument(7,3) + gauge(3) + gauge_i
        gauge(10) = argument(7,4) + gauge(4) + gauge_i 
        gauge(11) = argument(7,5) + gauge(5) + gauge_i 
        gauge(12) = argument(7,6) + gauge(6) + gauge_i 

     CASE(27)
!
!    S_6
!
        gauge(2) = gauge_a
        gauge(3) = -argument(2,3) - gauge_a
        gauge(4) = gauge_i
        gauge(5) = argument(4,2) + gauge_a + gauge_i
        gauge(6) = argument(4,3) + gauge(3) + gauge_i

     CASE(28,29)
!
!   T and T_h
!
        gauge(5) = gauge_a
        gauge(10) = gauge_b
        gauge(6) = -argument(10,6) - gauge_b 
        gauge(9) = -argument(5,9) - gauge_a 
        gauge(4) = argument(5,10) + gauge_a + gauge_b 
        gauge(3) = argument(10,5) + gauge_a + gauge_b 
        gauge(2) = argument(4,3) + gauge(4) + gauge(3)
        gauge(7) = argument(10,9) + gauge(9) + gauge(10)
        gauge(8) = argument(9,10) + gauge(9) + gauge(10)
        gauge(11) = argument(5,6) + gauge(5) + gauge(6)
        gauge(12) = argument(6,5) + gauge(5) + gauge(6)

        IF (code_group==29) THEN

           gauge(13) = gauge_i
           DO i=2,12
              gauge(12+i) = argument(13,i) + gauge(i) + gauge_i
           END DO

        END IF

     CASE(30,31,32)
!
!   T_d, O, and O_h
!
        gauge(18) = gauge_a
        gauge(5)  = gauge_b
        gauge(20) = argument(18,5) + gauge_a + gauge_b
        gauge(4)  = argument(18,18) + 2.0_DP*gauge_a 
        gauge(17) = -argument(18,17) - gauge_a 
        gauge(22) = argument(5,18) + gauge_b + gauge_a
        gauge(9)  = -argument(5,9) - gauge_b 
        gauge(15) = argument(18,9) + gauge_a + gauge(9) 
        gauge(14) = argument(17,5) + gauge_b + gauge(17) 
        gauge(13) = argument(9,18) + gauge(9) + gauge_a
        gauge(16) = argument(5,17) + gauge(17) + gauge_b
        gauge(11) = argument(4,9) + gauge(4) + gauge(9)
        gauge(10) = argument(9,4) + gauge(4) + gauge(9)
        gauge(8) = argument(15,18) + gauge(15) + gauge_a
        gauge(24) = argument(16,5) + gauge(5) + gauge(16)
        gauge(19) = argument(11,18) + gauge(11) + gauge(18)
        gauge(21) = argument(11,20) + gauge(11) + gauge(20)
        gauge(3) = argument(18,24) + gauge(18) + gauge(24)
        gauge(2) = argument(5,11) + gauge(5) + gauge(11)
        gauge(23) = argument(5,19) + gauge(5) + gauge(19)
        gauge(7) = argument(13,16) + gauge(13) + gauge(16)
        gauge(6) = argument(17,14) + gauge(17) + gauge(14)
        gauge(12) = argument(17,16) + gauge(17) + gauge(16)

        IF (code_group==32) THEN

           gauge(25) = gauge_i

           DO i=2,24
              gauge(24+i) = argument(25,i) + gauge(i) + gauge_i
           END DO

        END IF

   CASE DEFAULT
     CALL errore('find_projection_type','point group not available',1)
END SELECT
!
!   Print on output several checks, the group tables, the list of gauge
!   factors and the factor system after the gauge transformation
!
IF (verbose) THEN
   IF (ptype(1)==-1) THEN
      CALL find_double_product_table(prd, epos, code_group_ext)
   ELSE
      CALL find_product_table(prd, code_group_ext)
      epos=1
   ENDIF

!   WRITE(stdout,'(/,5x, "The product table:")')
!   tab(:,:)=prd(:,:)*epos(:,:)
!   CALL write_group_table_integer(group_desc, nsym, tab)

   WRITE(stdout,'(/,5x, "The following phases applied to the input factors")')
   WRITE(stdout,'(5x,"make them p-equivalent to the standard ones:",/)')

   DO isym=1,nsym
      arg=gauge(isym)
      phase(isym)=CMPLX(COS(arg),SIN(arg))
      WRITE(stdout,'(5x,i3,3x,a8,2f14.5)') isym, &
                                   sym_label(group_desc(isym)), phase(isym)
   END DO

   DO isym=1,nsym
      DO jsym=1,nsym
          arg=argument(isym,jsym)
          pha = CMPLX(COS(arg),SIN(arg))
          factor(isym,jsym)= pha * phase(isym) * phase(jsym) /  &
                             phase(prd(isym,jsym)) * DBLE(epos(isym,jsym))
      END DO
   END DO

   WRITE(stdout,'(/,5x, "The factor system after the gauge transformation:")')
   CALL write_group_table(group_desc, nsym, factor)

END IF

RETURN
END SUBROUTINE find_projection_type

SUBROUTINE find_irr_proj(code_group_ext,char_mat_proj,name_rap,nrap_proj,&
                                              nsym,ptype,gauge,verbosity)
!
!  This routine receives the arguments of the phases (gauge) of a given 
!  gauge. It sets the character tables that correspond to the standard factor
!  system and transform them into the characters of the symmetry operations 
!  written in the new gauge. 
!
!  ptype selects which representations are set:
!
!  those of the standard point group,
!  those of the double point group,
!  those of the projective set with beta and gamma given in ptype.
!
!  If verbosity=.TRUE. the character table is written on output
!
USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE 
INTEGER,  INTENT(IN) :: code_group_ext, nsym
INTEGER,  INTENT(IN) :: ptype(3)
REAL(DP), INTENT(IN) :: gauge(nsym)
LOGICAL,  INTENT(IN) :: verbosity

INTEGER, INTENT(OUT) :: nrap_proj
COMPLEX(DP), INTENT(INOUT) :: char_mat_proj(48,48)
CHARACTER(LEN=45), INTENT(INOUT) :: name_rap(48)

INTEGER :: group_desc(48)
INTEGER :: ntables, itables, start, last, irap, irot, nsym_
REAL(DP) :: arg
COMPLEX(DP) :: pha
!
! start by setting the standard representation
!
CALL set_stand_irr_proj(code_group_ext, ptype, char_mat_proj, name_rap, &
                               nrap_proj, nsym_)
!
! and apply the gauge that corresponds to the input phases
!
DO irot=1, nsym
   arg=gauge(irot)
   pha=CMPLX(COS(arg),SIN(arg))
   DO irap=1,nrap_proj
      char_mat_proj(irap,irot) = char_mat_proj(irap,irot) / pha
   END DO
END DO
!
!  write the representation on output if requested
!
IF (verbosity) THEN

   CALL set_group_desc(group_desc,nsym_,code_group_ext)

   CALL write_gauge(gauge, group_desc, nsym)

   CALL print_ptype_info(ptype, code_group_ext)

   WRITE(stdout,'(5x, "after the gauge transformation:")')

   CALL write_group_char_mat(group_desc, nsym, char_mat_proj, name_rap, &
                                                                nrap_proj)
ENDIF
RETURN
END SUBROUTINE find_irr_proj

SUBROUTINE print_element_list(code_group_ext)
!
!  This routine print the list of elements of a point group, given its
!  extended code
!
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: code_group_ext
INTEGER :: i, group_desc(48), nsym

CALL set_group_desc(group_desc,nsym,code_group_ext)

WRITE(stdout,'((5x,4(i2," - ",i2,1x,a8)))') (i, group_desc(i), &
                                     sym_label(group_desc(i)), i=1,nsym)

RETURN
END SUBROUTINE print_element_list


SUBROUTINE compute_classes(cge, nclasses, nelem, elem)
!
!  This subroutine finds the classes of a point group
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: cge  ! extended code group index
INTEGER, INTENT(OUT) :: nclasses, nelem(24), elem(18,24)

INTEGER :: group_desc(48), nsym, prd(48,48), invs(48), done(48)

INTEGER :: isym, jsym, ind1, ind2

CALL set_group_desc(group_desc,nsym,cge)
CALL find_product_table(prd, cge)

DO isym=1,nsym
   DO jsym=1,nsym
      IF (prd(isym,jsym)==1) invs(isym)=jsym
   ENDDO
ENDDO

done=0
nelem=0
nclasses=0
DO isym=1,nsym
   IF (done(isym)/=0) CYCLE
   nclasses=nclasses+1
   nelem(nclasses)=1
   elem(nelem(nclasses),nclasses)=isym
   done(isym)=1
   DO jsym=1,nsym
      ind1=prd(jsym,isym)
      ind2=prd(ind1,invs(jsym))
      IF (done(ind2)==0) THEN
         nelem(nclasses)=nelem(nclasses)+1
         elem(nelem(nclasses),nclasses)=ind2
         done(ind2)=1
      ENDIF
   ENDDO
ENDDO  

RETURN
END SUBROUTINE compute_classes

SUBROUTINE compute_classes_double(cge, nclasses, nelem, elem, has_e)
!
!  This subroutine finds the classes of a double point group
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: cge  ! extended code group index
INTEGER, INTENT(OUT) :: nclasses, nelem(24), elem(18,24), has_e(18,24)

INTEGER :: group_desc(48), nsym, prd(48,48), epos(48,48), invs(48), done(96)

INTEGER :: isym, jsym, ind1, ind2, has_e1, has_e2, has_einv(48)

CALL set_group_desc(group_desc,nsym,cge)
CALL find_double_product_table(prd, epos, cge)

DO isym=1,nsym
   DO jsym=1,nsym
      IF (prd(isym,jsym)==1) THEN
         invs(isym)=jsym
         has_einv(isym)=epos(isym,jsym)
      ENDIF
   ENDDO
ENDDO

done=0
nelem=0
nclasses=0
DO isym=1,nsym
   IF (done(isym)/=0) CYCLE
   nclasses=nclasses+1
   IF (nclasses > 18) CALL errore('compute_classes_double', &
                                    'Too many classes',1)
   nelem(nclasses)=1
   elem(nelem(nclasses),nclasses)=isym
   has_e(nelem(nclasses),nclasses)=1
   done(isym)=1
!
!  make the conjugation with all the operation of the point group
!
   DO jsym=1,nsym
      ind1=prd(jsym,isym)
      has_e1 = epos(jsym,isym)
      ind2=prd(ind1,invs(jsym))
      has_e2 = epos(ind1,invs(jsym)) * has_e1 * has_einv(jsym)
      IF (has_e2==1) THEN
         IF (done(ind2)==0) THEN
            nelem(nclasses)=nelem(nclasses)+1
            elem(nelem(nclasses),nclasses)=ind2
            has_e(nelem(nclasses),nclasses)=1
            done(ind2)=1
         ENDIF
      ELSE
         IF (done(ind2+nsym)==0) THEN
            nelem(nclasses)=nelem(nclasses)+1
            elem(nelem(nclasses),nclasses)=ind2
            has_e(nelem(nclasses),nclasses)=-1
            done(ind2+nsym)=1
         ENDIF
      ENDIF
   ENDDO
!
!  make the conjugation with all the operations of the point group multiplied
!  by -E
!
   DO jsym=1,nsym
      ind1=prd(jsym,isym)
      has_e1 = -epos(jsym,isym) 
      ind2=prd(ind1,invs(jsym))
      has_e2 = -epos(ind1,invs(jsym)) * has_e1 * has_einv(jsym)
      IF (has_e2==1) THEN
         IF (done(ind2)==0) THEN
            nelem(nclasses)=nelem(nclasses)+1
            elem(nelem(nclasses),nclasses)=ind2
            has_e(nelem(nclasses),nclasses)=1
            done(ind2)=1
         ENDIF
      ELSE
         IF (done(ind2+nsym)==0) THEN
            nelem(nclasses)=nelem(nclasses)+1
            elem(nelem(nclasses),nclasses)=ind2
            has_e(nelem(nclasses),nclasses)=-1
            done(ind2+nsym)=1
         ENDIF
      ENDIF
   ENDDO
ENDDO  
!
!  here classes starting from the elements multiplied by -E
!
DO isym=1,nsym
   IF (done(isym+nsym)/=0) CYCLE
   nclasses=nclasses+1
   IF (nclasses > 18) CALL errore('compute_classes_double', &
                                 'Too many classes 1',1)
   nelem(nclasses)=1
   elem(nelem(nclasses),nclasses)=isym
   has_e(nelem(nclasses),nclasses)=-1
   done(isym+nsym)=1
!
!  make the conjugation with all the operation of the point group
!
   DO jsym=1,nsym
      ind1=prd(jsym,isym)
!
!   the minus is the -E of isym
!
      has_e1 = -epos(jsym,isym)
      ind2=prd(ind1,invs(jsym))
      has_e2 = epos(ind1,invs(jsym)) * has_e1 * has_einv(jsym)
      IF (has_e2==1) THEN
         IF (done(ind2)==0) THEN
            nelem(nclasses)=nelem(nclasses)+1
            elem(nelem(nclasses),nclasses)=ind2
            has_e(nelem(nclasses),nclasses)=1
            done(ind2)=1
         ENDIF
      ELSE
         IF (done(ind2+nsym)==0) THEN
            nelem(nclasses)=nelem(nclasses)+1
            elem(nelem(nclasses),nclasses)=ind2
            has_e(nelem(nclasses),nclasses)=-1
            done(ind2+nsym)=1
         ENDIF
      ENDIF
   ENDDO
!
!  make the conjugation with all the operations of the point group multiplied
!  by -E
!
   DO jsym=1,nsym
      ind1=prd(jsym,isym)
!
!   the plus sign is due to the fact that both isym and jsym have -E 
!
      has_e1 = epos(jsym,isym) 
      ind2=prd(ind1,invs(jsym))
      has_e2 = -epos(ind1,invs(jsym)) * has_e1 * has_einv(jsym)
      IF (has_e2==1) THEN
         IF (done(ind2)==0) THEN
            nelem(nclasses)=nelem(nclasses)+1
            elem(nelem(nclasses),nclasses)=ind2
            has_e(nelem(nclasses),nclasses)=1
            done(ind2)=1
         ENDIF
      ELSE
         IF (done(ind2+nsym)==0) THEN
            nelem(nclasses)=nelem(nclasses)+1
            elem(nelem(nclasses),nclasses)=ind2
            has_e(nelem(nclasses),nclasses)=-1
            done(ind2+nsym)=1
         ENDIF
      ENDIF
   ENDDO
ENDDO  

RETURN
END SUBROUTINE compute_classes_double

FUNCTION is_subgroup(group_in_ext, group_out_ext)
!
!  The function returns .TRUE. if group_out is a subgroup of group_in
!  The two extended codes of the groups are given in input
!
IMPLICIT NONE
LOGICAL :: is_subgroup
INTEGER, INTENT(IN) :: group_in_ext, group_out_ext

INTEGER :: nsym_in, nsym_out
INTEGER :: group_desc_in(48), group_desc_out(48)
INTEGER :: isym, jsym, irot

IF (group_in_ext<1.OR.group_in_ext>136.OR.group_out_ext<1.OR.group_out_ext>136)&
   THEN
   is_subgroup=.FALSE.
   RETURN
ENDIF

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

SUBROUTINE set_sym_su2(sym_num, smat, sinv)
!
!  This routine uses the Cayley-Klein parameters to set the su2 rotation
!  matrices for the 32 proper rotations defined in the module. 
!  See sym_label for the name of each rotation. The integer sinv is set 
!  to 1 if the operation has not inversion, to -1 if the operation is to 
!  be multiplied by inversion. Note that inversion does not act on spin, 
!  so the su2 matrices remain the same. To multiply by E' just change
!  the sign of the matrix. The choices of the su2 matrix corresponding
!  to each rotation is done following S.K. Kim, Group theoretical methods.
!  except for D_3 and C_3v for which the same matrices used for the
!  other groups are used. The character tables are changed accordingly.
!
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: sym_num
  INTEGER, INTENT(OUT) :: sinv
  COMPLEX(DP), INTENT(OUT) :: smat(2,2)

  INTEGER :: snum
  COMPLEX(DP) :: a(32), b(32)
  REAL(DP) :: sqrt2=SQRT(2.0_DP), sqrt3=SQRT(3.0_DP)

  IF (sym_num < 1 .OR. sym_num > 64) CALL errore('set_sym_su2', &
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
  END SUBROUTINE set_sym_su2

  SUBROUTINE find_product_table(prd, code_group_ext)

  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: prd(48,48)
  INTEGER, INTENT(IN) :: code_group_ext

  REAL(DP) :: group_mat(3,3,48), a_mat(3,3), b_mat(3,3), &
              c_mat(3,3) 
  INTEGER :: group_desc(48)
  INTEGER :: isym, jsym, ksym, nsym
  !
  ! Find the list of operations of the group and sets their 3x3 matrices
  !
  CALL set_group_desc(group_desc,nsym,code_group_ext)

  DO isym=1, nsym
     CALL set_sym_o3(group_mat(:,:,isym),group_desc(isym))
  ENDDO
  !
  !  make the product and compare with the other matrices
  !
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

  SUBROUTINE set_sym_o3(mat,isym)
!
!  Set the symmetry matrix isym as a 3x3 orthogonal matrix in cartesian i
!  coordinates
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: isym
  REAL(DP), INTENT(OUT) :: mat(3,3)

  real(DP), PARAMETER :: sin3 = 0.866025403784438597d0, cos3 = 0.5d0, &
                        msin3 =-0.866025403784438597d0, mcos3 = -0.5d0

  REAL(DP) :: s0(3,3,32)

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

  IF (isym < 33 ) THEN
     mat(:,:)=s0(:,:,isym)
  ELSE
     mat(:,:)=-s0(:,:,isym-32)
  ENDIF
 
  END SUBROUTINE set_sym_o3

 SUBROUTINE find_double_product_table(prd, epos, code_group_ext)

 IMPLICIT NONE
 INTEGER, INTENT(IN)  :: code_group_ext
 INTEGER, INTENT(INOUT) :: prd(48,48), epos(48,48)
 
 INTEGER :: group_desc(48), nsym

 CALL set_group_desc(group_desc,nsym,code_group_ext)

 CALL find_double_product_table_from_sym(prd, epos, group_desc, nsym)

 RETURN
 END SUBROUTINE find_double_product_table

 SUBROUTINE find_double_product_table_from_sym(prd, epos, group_desc, nsym)
!
!  This routine provides the product table prd of a given double group.
!  Each entry prd(i,j) is the index of the operation that results from
!  the product of S_i S_j (in this order).
!  The symmetry operations are only those without E'. If necessary
!  the table can be extended to the entire group using EE=E, EE'=E'E=E', 
!  E'E'=E.
!
!  The routine provides also the factor system epos of the point group
!  considered as a projective representation of the double group.
!
 USE kinds, ONLY : DP
 IMPLICIT NONE
 INTEGER, INTENT(OUT) :: prd(48,48), epos(48,48)

 INTEGER :: group_desc(48), sinv(48), isym, jsym, ksym, lsym, nsym
 COMPLEX(DP) :: group_mat(2,2,48)
  
 epos=1
 prd=0
 DO isym=1, nsym
    DO jsym=1, nsym
       CALL product_sym_su2(group_desc(isym), group_desc(jsym), &
                            ksym, epos(isym,jsym))
       DO lsym=1,nsym
          IF (group_desc(lsym)==ksym) prd(isym,jsym)=lsym
       ENDDO
    END DO
 END DO

 RETURN
 END SUBROUTINE find_double_product_table_from_sym

 SUBROUTINE product_sym_su2(isym, jsym, prd, epos)
 !
 !  This routine recives the indeces of two symmetry operations, the
 !  list of symmetry operations in su2 form and gives the index of the
 !  product inside the group list and the possible -E operation
 !
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: isym, jsym
 INTEGER, INTENT(OUT) :: prd, epos

 COMPLEX(DP) :: a_mat(2,2), b_mat(2,2), c_mat(2,2), group_mat(2,2,64)
 INTEGER :: sinv(64)
 REAL(DP) :: diff, diff1
 INTEGER :: sinv1, sinv2, sinvp
 INTEGER :: pcount, ksym

 DO ksym=1,64
    CALL set_sym_su2(ksym, group_mat(1,1,ksym), sinv(ksym))
!    WRITE(6,*) 'ksym', ksym, group_desc(ksym)
 ENDDO

 a_mat(:,:)=group_mat(:,:,isym)
 sinv1=sinv(isym)
 b_mat(:,:)=group_mat(:,:,jsym)
 sinv2=sinv(jsym)

 sinvp = sinv1*sinv2
 c_mat(1,1)=a_mat(1,1) * b_mat(1,1) + a_mat(1,2) * b_mat(2,1)
 c_mat(1,2)=a_mat(1,1) * b_mat(1,2) + a_mat(1,2) * b_mat(2,2)

 epos=1
 prd=0
 pcount=0
 DO ksym = 1, 64
    IF (sinv(ksym)==sinvp) THEN
       diff=ABS(c_mat(1,1)-group_mat(1,1,ksym))+ &
            ABS(c_mat(1,2)-group_mat(1,2,ksym))
       IF (diff < 1.D-6) THEN
          prd=ksym
          pcount=pcount+1
       ELSE
          diff1=ABS(c_mat(1,1)+group_mat(1,1,ksym))+ &
                ABS(c_mat(1,2)+group_mat(1,2,ksym))
          IF (diff1 < 1.D-6) THEN
             prd=ksym
             epos=-1
             pcount=pcount+1
          ENDIF
       ENDIF
    END IF
 END DO
 IF (pcount/=1) &
    CALL errore('product_sym_su2','The product of these matrices is &
                                   &not among the allowed symmetries',1)

 
 RETURN
 END SUBROUTINE product_sym_su2
!
 SUBROUTINE  set_factors(group_code, ptype, argument)
!
!  This routine set the standard factor system of a point group that
!  correspond to a given projection type. Note that only the factors
!  corresponding to beta and gamma for groups that have inversion 
!  are set here. The double group factors are not set by this routine and
!  can be added with the epos factors calculated by find_double_product_table.
!
 USE constants, ONLY : pi
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: group_code
 INTEGER, INTENT(IN) :: ptype(3)
 REAL(DP), INTENT(OUT) :: argument(48,48)

 INTEGER :: factors(48,48), isym, jsym, nsym

 nsym=0 
 factors=1
 SELECT CASE (group_code)

  CASE (16,18,19,28)
!
!  C_2h, C_4h, C_6h, S_6
!
      IF (group_code==16) nsym=4
      IF (group_code==18) nsym=8
      IF (group_code==19) nsym=12
      IF (group_code==28) nsym=6

      DO isym=2, nsym/2
         DO jsym=nsym/2+1,nsym
            factors(isym, jsym) = ptype(2)**(isym-1)
            factors(isym+nsym/2, jsym) = ptype(2)**(isym-1)
         ENDDO
      ENDDO

  CASE (20,25,22,23)
!
!  D_2h, D_3d, D_4h, D_6h
!
      IF (group_code==20) nsym=8
      IF (group_code==25) nsym=12
      IF (group_code==22) nsym=16
      IF (group_code==23) nsym=24

      DO isym=1, nsym/4
         DO jsym=nsym/2+1,nsym
            factors(isym, jsym) = ptype(2)**(isym-1)
            factors(isym+nsym/4, jsym) = ptype(2)**(isym-1) * ptype(3) 
            factors(isym+nsym/2, jsym) = ptype(2)**(isym-1)
            factors(isym+3*nsym/4, jsym) = ptype(2)**(isym-1) * ptype(3) 
         ENDDO
      ENDDO
  CASE (32)
!
!  O_h
!
    nsym=48
    DO isym=13,24
       DO jsym=nsym/2+1,nsym
          factors(isym, jsym) = ptype(2)
          factors(isym+nsym/2, jsym) = ptype(2)
       END DO
    END DO
  CASE DEFAULT

 END SELECT

 argument=0.0_DP
 DO isym=1,nsym
    DO jsym=1,nsym
       IF (factors(isym,jsym)==-1) argument(isym,jsym)=pi
    END DO
 END DO
 RETURN
 END SUBROUTINE set_factors

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

SUBROUTINE zero_tpi(arg)
!
!  This subroutine receives the argument of a phase and brings it in the
!  [0,tpi] range.
!
USE kinds, ONLY : DP
USE constants, ONLY : tpi
IMPLICIT NONE

REAL(DP), INTENT(INOUT) :: arg

arg = MOD(arg,tpi)
IF (arg < 0.0_DP) arg=arg+tpi

RETURN
END SUBROUTINE zero_tpi

SUBROUTINE write_group_table(group_desc, nsym, factor)

USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE

INTEGER, INTENT(IN) :: nsym, group_desc(48)
INTEGER :: isym, jsym, itables, ntables, start, last
COMPLEX(DP) :: factor(48, 48)
LOGICAL :: ctable

ctable=SUM(ABS(AIMAG(factor))) > 1.D-4
ntables= nsym / 8
IF (MOD(nsym,8) /= 0 ) ntables=ntables+1
DO itables=1,ntables
   start=(itables-1)*8+1
   last=MIN(itables*8,nsym)
   WRITE(stdout,'(/,5x,"real part")')
   WRITE(stdout,'(8x,8(1x,a8))') (sym_label(group_desc(jsym)), &
                                  jsym=start,last)
   DO isym=1,nsym
      WRITE(stdout,'(a8,8(f7.2,2x))') sym_label(group_desc(isym)), &
       (REAL(factor(isym, jsym)),jsym=start,last)
   ENDDO
   IF (ctable) THEN
      WRITE(stdout,'(5x,"imaginary part")')
      WRITE(stdout,'(8x,8(1x,a8))') (sym_label(group_desc(jsym)), &
                                         jsym=start,last)
      DO isym=1,nsym
         WRITE(stdout,'(a8,8(f7.2,2x))') sym_label(group_desc(isym)), &
           (AIMAG(factor(isym, jsym)),jsym=start,last)
      ENDDO
   END IF
END DO

RETURN
END SUBROUTINE write_group_table

SUBROUTINE write_group_table_integer(group_desc, nsym, factor)

USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE

INTEGER, INTENT(IN) :: nsym, group_desc(48)
INTEGER :: isym, jsym, itables, ntables, start, last
INTEGER :: factor(48, 48)

ntables= nsym / 8
IF (MOD(nsym,8) /= 0 ) ntables=ntables+1
DO itables=1,ntables
   start=(itables-1)*8+1
   last=MIN(itables*8,nsym)
   WRITE(stdout,'(8x,8(1x,a8))') (sym_label(group_desc(jsym)), &
                                  jsym=start,last)
   DO isym=1,nsym
      WRITE(stdout,'(a8,8(i7,2x))') sym_label(group_desc(isym)), &
       (factor(isym, jsym),jsym=start,last)
   ENDDO
ENDDO

RETURN
END SUBROUTINE write_group_table_integer

SUBROUTINE write_group_char_mat(group_desc, nsym, char_mat, name_rap, nrap)
USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: nsym, nrap, group_desc(48)
COMPLEX(DP), INTENT(IN) :: char_mat(48, 48)
CHARACTER(LEN=45), INTENT(IN) :: name_rap(48)

INTEGER :: irap, irot, itables, ntables, start, last
LOGICAL :: ctable

ctable=SUM(ABS(AIMAG(char_mat))) > 1.D-4
ntables= nsym / 8
IF (MOD(nsym,8) /= 0 ) ntables=ntables+1
DO itables=1,ntables
   start=(itables-1)*8+1
   last=MIN(itables*8,nsym)
   WRITE(stdout,'(/,5x,"real part")')
   WRITE(stdout,'(8x,8(1x,a8))') (sym_label(group_desc(irot)), irot=start,last)
   DO irap=1,nrap
      WRITE(stdout,'(a8,9f8.2)') name_rap(irap), &
      (REAL(char_mat(irap,irot)),irot=start,last)
   ENDDO

   IF (ctable) THEN
      WRITE(stdout,'(5x,"imaginary part")')
      WRITE(stdout,'(8x,8(1x,a8))') (sym_label(group_desc(irot)), &
                                                          irot=start,last)
      DO irap=1,nrap
         WRITE(stdout,'(a8,9f8.2)') name_rap(irap), &
         (AIMAG(char_mat(irap,irot)),irot=start,last)
      ENDDO
   END IF
ENDDO
WRITE(stdout,*)

RETURN
END SUBROUTINE write_group_char_mat

SUBROUTINE print_character_table(cge,ptype)
  
USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER,INTENT(IN) :: cge, ptype(3)

INTEGER :: irot, irap, nsym, nrap, itables, ntables, start, last
COMPLEX(DP) :: char_mat_proj(48,48)
CHARACTER(LEN=45) :: name_rap(48)
INTEGER :: group_desc(48)
LOGICAL :: ctable

CALL set_group_desc(group_desc,nsym,cge)
CALL set_stand_irr_proj(cge, ptype, char_mat_proj, name_rap, nrap, nsym)

IF (nrap==0) RETURN

CALL print_ptype_info(ptype, cge)

ctable=SUM(ABS(AIMAG(char_mat_proj))) > 1.D-4
ntables= nsym / 8
IF (MOD(nsym,8) /= 0 ) ntables=ntables+1
DO itables=1,ntables
   start=(itables-1)*8+1
   last=MIN(itables*8,nsym)
   WRITE(stdout,'(/,5x,"real part")')
   WRITE(stdout,'(13x,8(a7,1x))') (sym_label(group_desc(irot)),&
                                                           irot=start,last)
   DO irap=1,nrap
      WRITE(stdout,'(a8,9f8.2)') TRIM(name_rap(irap)), &
      (REAL(char_mat_proj(irap,irot)),irot=start,last)
   ENDDO
   IF (ctable) THEN
      WRITE(stdout,'(5x,"imaginary part")')
      WRITE(stdout,'(13x,8(a7,1x))') (sym_label(group_desc(irot)),&
                                                             irot=start,last)
      DO irap=1,nrap
         WRITE(stdout,'(a8,9f8.2)') TRIM(name_rap(irap)), &
            (AIMAG(char_mat_proj(irap,irot)),irot=start,last)
      ENDDO
   ENDIF
ENDDO

RETURN
END SUBROUTINE print_character_table

SUBROUTINE print_compatibility_table(cge_in,ptype_in,cge_out)
!
!  This routine writes the compatibility table between the irreducible 
!  representations of the group cge_in with the projection type ptype_in, 
!  and the irreducible representations of one of its subgroups cge_out.
!  It determines also which representations type of cge_out must be
!  used to decompose those of cge_in. Projective representations are
!  allowed.
!
USE kinds, ONLY : DP
USE constants, ONLY : pi
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER,INTENT(IN) :: cge_in, ptype_in(3), cge_out

COMPLEX(DP) :: char_mat_proj_in(48,48), char_mat_proj_out(48,48)
REAL(DP) :: arguments_in(48,48), arguments_out(48,48), gauge_in(48), &
                                           gauge_out(48)
CHARACTER(LEN=45) :: name_rap_in(48), name_rap_out(48)

INTEGER :: group_desc_in(48), group_desc_out(48), b_in_a(48), list_in(100), &
           list_out(100), rap(12), prd(48,48), epos(48,48), ptype_out(3)

INTEGER :: group_in, group_out, nsym_in, nsym_out, nrap_in, nrap_out
INTEGER :: n, i, ndeg, isym, jsym, irap, jrap

CHARACTER(LEN=256) :: rap_name
CHARACTER(LEN=6) :: int_to_char

IF (.NOT.(is_subgroup(cge_in,cge_out))) RETURN

CALL set_stand_irr_proj(cge_in, ptype_in, char_mat_proj_in, &
                                name_rap_in, nrap_in, nsym_in)

IF (nrap_in==0) RETURN

WRITE(stdout,*)
CALL print_ptype_info(ptype_in, cge_in)

group_in = group_index_from_ext(cge_in)
group_out = group_index_from_ext(cge_out)
CALL set_group_desc(group_desc_in,nsym_in,cge_in)
CALL set_group_desc(group_desc_out,nsym_out,cge_out)

arguments_in=0.0_DP
CALL set_factors(group_in, ptype_in, arguments_in)

IF (ptype_in(1)==-1) THEN
   CALL find_double_product_table(prd, epos, cge_in)
ELSE
   CALL find_product_table(prd, cge_in)
   epos=1
ENDIF

b_in_a=0
DO isym=1,nsym_out
   DO jsym=1,nsym_in
      IF (group_desc_in(jsym)==group_desc_out(isym)) b_in_a(isym)=jsym
   ENDDO
   IF (b_in_a(isym)==0) &
      CALL errore('convert_rap_proj','group_out not subgroup of group_in',1)
ENDDO
   
DO isym=1,nsym_out
   DO jsym=1,nsym_out
      arguments_out(isym,jsym) = arguments_in(b_in_a(isym),b_in_a(jsym)) 
      IF (epos(b_in_a(isym),b_in_a(jsym))==-1) &
            arguments_out(isym,jsym)=arguments_out(isym,jsym) + pi
   END DO
END DO

CALL find_projection_type(group_out, cge_out, arguments_out, &
                          ptype_out, gauge_out, .FALSE.)
!
! special case groups 1 {E} and 28 {E, I} the factor system does not
! allow to determine if the point group or the double point group
! must be used (actually they have the same rappresentation), but we keep two
! different names and force here the name to use).
!
IF (((cge_in==1.OR.cge_in==28).OR.(cge_out==1.OR.cge_out==28))&
                                    .AND.ptype_in(1)==-1) ptype_out(1)=-1
n=0
DO irap=1,nrap_in
   ndeg=NINT(DBLE(char_mat_proj_in(irap,1)))
   DO i=1,ndeg
      n=n+1
      list_in(n)=irap
   ENDDO
ENDDO

CALL convert_rap_proj(n, list_in, list_out, cge_in, cge_out, ptype_in, &
                     ptype_out, gauge_in, gauge_out )

CALL set_stand_irr_proj(cge_out, ptype_out, char_mat_proj_out, &
                                name_rap_out, nrap_out, nsym_out)

WRITE(stdout,'(5x,"decomposed into ")')
CALL print_ptype_info(ptype_out, cge_out)
WRITE(stdout,*)

n=0
DO irap=1,nrap_in
   ndeg=NINT(DBLE(char_mat_proj_in(irap,1)))
   rap=0
   DO i=1,ndeg
      n=n+1
      rap(list_out(n))=rap(list_out(n))+1
   ENDDO
   rap_name=''
   DO jrap=1,nrap_out
      ndeg=NINT(DBLE(char_mat_proj_out(jrap,1)))
      rap(jrap)=rap(jrap)/ndeg
      CALL add_rap_name(rap_name, rap(jrap), name_rap_out(jrap)(1:8))
   ENDDO
   WRITE(6,'(5x, a8," = ",a)') name_rap_in(irap), TRIM(rap_name)
ENDDO

RETURN
END SUBROUTINE print_compatibility_table

SUBROUTINE print_kronecker_table(cge_in,ptype1_in,ptype2_in,cge_out,lcomp)
!
!  This routine writes the compatibility table between the product of
!  two irreducible representations of the group cge_in with the projection 
!  type ptype1_in and ptype2_in and the irreducible representations of one 
!  of its subgroups cge_out.
!  It determines which representations type of cge_out must be
!  used to decompose those of cge_in. Projective representations are
!  allowed.
!
!  When lcomp=.TRUE. it decomposes chi1^*(S) chi2(S)
!  When lcomp=.FALSE. it decomposes chi1(S) chi2(S)
!
USE kinds, ONLY : DP
USE constants, ONLY : pi
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER,INTENT(IN) :: cge_in, ptype1_in(3), ptype2_in(3), cge_out
LOGICAL, INTENT(IN) :: lcomp

COMPLEX(DP) :: char_mat_proj1_in(48,48), char_mat_proj2_in(48,48), &
               char_mat_proj_out(48,48), &
               char_mat_prod(48), char_mat_sub(48), pha, asum
REAL(DP) :: arguments1_in(48,48), arguments2_in(48,48), arguments_out(48,48), &
            arguments1_out(48,48), arguments2_out(48,48), gauge_in(48), &
            gauge_out(48), arg
CHARACTER(LEN=45) :: name_rap1_in(48), name_rap2_in(48), name_rap_out(48)

INTEGER :: group_desc_in(48), group_desc_out(48), b_in_a(48), &
           prd(48,48), epos1(48,48), epos2(48,48), ptype_out(3)

INTEGER :: group_in, group_out, nsym_in, nsym_out, nrap1_in, nrap2_in,nrap_out
INTEGER :: n, i, l, m,  ndeg, isym, jsym, irap, jrap, krap

CHARACTER(LEN=256) :: rap_name
CHARACTER(LEN=6) :: int_to_char

IF (.NOT.(is_subgroup(cge_in,cge_out))) RETURN

CALL set_stand_irr_proj(cge_in, ptype1_in, char_mat_proj1_in, &
                                name_rap1_in, nrap1_in, nsym_in)
CALL set_stand_irr_proj(cge_in, ptype2_in, char_mat_proj2_in, &
                                name_rap2_in, nrap2_in, nsym_in)

IF (nrap1_in==0.OR.nrap2_in==0) RETURN

WRITE(stdout,'(/,5x,"The first representation from")')
CALL print_ptype_info(ptype1_in, cge_in)
WRITE(stdout,'(5x,"the second representation from")')
CALL print_ptype_info(ptype2_in, cge_in)

group_in = group_index_from_ext(cge_in)
group_out = group_index_from_ext(cge_out)
CALL set_group_desc(group_desc_in,nsym_in,cge_in)
CALL set_group_desc(group_desc_out,nsym_out,cge_out)

arguments1_in=0.0_DP
CALL set_factors(group_in, ptype1_in, arguments1_in)
arguments2_in=0.0_DP
CALL set_factors(group_in, ptype2_in, arguments2_in)

IF (ptype1_in(1)==-1) THEN
   CALL find_double_product_table(prd, epos1, cge_in)
ELSE
   CALL find_product_table(prd, cge_in)
   epos1=1
ENDIF

IF (ptype2_in(1)==-1) THEN
   CALL find_double_product_table(prd, epos2, cge_in)
ELSE
   CALL find_product_table(prd, cge_in)
   epos2=1
ENDIF

b_in_a=0
DO isym=1,nsym_out
   DO jsym=1,nsym_in
      IF (group_desc_in(jsym)==group_desc_out(isym)) b_in_a(isym)=jsym
   ENDDO
   IF (b_in_a(isym)==0) &
      CALL errore('convert_rap_proj','group_out not subgroup of group_in',1)
ENDDO
   
DO isym=1,nsym_out
   DO jsym=1,nsym_out
      arguments1_out(isym,jsym) = arguments1_in(b_in_a(isym),b_in_a(jsym)) 
      IF (epos1(b_in_a(isym),b_in_a(jsym))==-1) &
            arguments1_out(isym,jsym)=arguments1_out(isym,jsym) + pi
      arguments2_out(isym,jsym) = arguments2_in(b_in_a(isym),b_in_a(jsym)) 
      IF (epos2(b_in_a(isym),b_in_a(jsym))==-1) &
            arguments2_out(isym,jsym)=arguments2_out(isym,jsym) + pi
      arguments_out(isym,jsym) = arguments1_out(isym,jsym) + &
                                 arguments2_out(isym,jsym)
   END DO
END DO

CALL find_projection_type(group_out, cge_out, arguments_out, &
                          ptype_out, gauge_out, .FALSE.)
!
! special case groups 1 {E} and 28 {E, I} the factor system does not
! allow to determine if the point group or the double point group
! must be used (actually they have the same rappresentation), but we keep two
! different names and force here the name to use).
!
IF (((cge_in==1.OR.cge_in==28).OR.(cge_out==1.OR.cge_out==28))&
                                    .AND.ptype1_in(1)==-1) ptype_out(1)=-1

CALL set_stand_irr_proj(cge_out, ptype_out, char_mat_proj_out, &
                                name_rap_out, nrap_out, nsym_out)

WRITE(stdout,'(5x,"decomposed into ")')
CALL print_ptype_info(ptype_out, cge_out)
WRITE(stdout,*)

DO isym=1,nsym_out
   arg=gauge_out(isym)
   pha=CMPLX(COS(arg),SIN(arg))
   DO irap=1,nrap_out
      char_mat_proj_out(irap,isym) = char_mat_proj_out(irap,isym) / pha
   END DO
END DO

IF (lcomp) THEN
   WRITE(stdout,'(5x,"Decomposing the product chi^*(S) x chi(S)")')
ELSE
   WRITE(stdout,'(5x,"Decomposing the product chi(S) x chi(S)")')
ENDIF

DO irap=1,nrap1_in
   DO krap=1,nrap2_in
!
!   make the product of the two representations
!
      IF (lcomp) THEN
         char_mat_prod(1:nsym_in) = CONJG(char_mat_proj1_in(irap,1:nsym_in)) * &
                                    char_mat_proj2_in(krap,1:nsym_in)
      ELSE
         char_mat_prod(1:nsym_in) = char_mat_proj1_in(irap,1:nsym_in) * &
                                    char_mat_proj2_in(krap,1:nsym_in)
      ENDIF
!
!  project on the subgroup
!
      DO isym=1,nsym_out
         char_mat_sub(isym)=char_mat_prod(b_in_a(isym))
      ENDDO
!
!  Do the projection on the subgroup representations
!
      rap_name=''
      DO jrap=1, nrap_out

         asum=(0.0_DP,0.0_DP)
         DO isym=1,nsym_out
            asum = asum + char_mat_sub(isym) * &
                           CONJG(char_mat_proj_out(jrap,isym))
         ENDDO
         ndeg= NINT(DBLE(asum) / DBLE(nsym_out))
!         WRITE(6,*) irap, jrap, asum, ndeg
         IF (ABS(ndeg - asum / DBLE(nsym_out)) > 1.D-6) &
             CALL errore('compute_kronecker_table','problem with ndeg',1)

         CALL add_rap_name(rap_name, ndeg, name_rap_out(jrap)(1:8))
      ENDDO
!
!  and write the output
!
      WRITE(6,'(5x, a8," x ",a8," =",3x,a)') TRIM(name_rap1_in(irap)), &
                                TRIM(name_rap2_in(krap)), TRIM(rap_name)
   ENDDO
ENDDO

RETURN
END SUBROUTINE print_kronecker_table

SUBROUTINE add_rap_name(rap_name, ndeg, add_name)
!
! This routines adds the string add_name to rap_name. When rap_name
! is not empty its adds a + sign. Moreover if nged > 1 it puts it 
! as a coefficient to add_name
!
IMPLICIT NONE
CHARACTER(LEN=256), INTENT(INOUT) :: rap_name
INTEGER, INTENT(IN) :: ndeg
CHARACTER(LEN=8), INTENT(IN) :: add_name
CHARACTER(LEN=6) :: int_to_char

IF (rap_name=='') THEN
   IF (ndeg>1) &
      rap_name=TRIM(int_to_char(ndeg))//' '//TRIM(add_name)
   IF (ndeg==1) &
         rap_name=TRIM(add_name)
ELSE
   IF (ndeg>1) &
      rap_name=TRIM(rap_name)//' + '//TRIM(int_to_char(ndeg))&
                                  //' '//TRIM(add_name)
   IF (ndeg==1) &
      rap_name=TRIM(rap_name)//' + '//TRIM(add_name)
ENDIF

RETURN
END SUBROUTINE add_rap_name

SUBROUTINE write_gauge(gauge, group_desc, nsym)
!
!  This routine prints the gauge factors. In input gauge are the
!  arguments of the phases.
!
USE kinds, ONLY : DP
USE constants, ONLY : pi
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: nsym
INTEGER, INTENT(IN) :: group_desc(48)
REAL(DP), INTENT(IN) :: gauge(nsym)
REAL(DP) :: arg
COMPLEX(DP) :: pha

INTEGER :: isym

WRITE(stdout,'(/,5x,"Gauge factors:   argument (deg)         phase",/)')
DO isym=1, nsym
   arg=gauge(isym)
   pha=CMPLX(COS(arg), SIN(arg))
   WRITE(stdout,'(5x,a8,6x,f10.2,6x,2f10.4)') sym_label(group_desc(isym)), &
                                              gauge(isym)*180._DP/pi, pha
END DO
WRITE(stdout,*)

RETURN
END SUBROUTINE write_gauge

SUBROUTINE print_ptype_info(ptype, cge)
!
!  This routine can be used to print information on ptype on output.
!
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: ptype(3), cge

INTEGER :: code_group
CHARACTER(LEN=11) :: group_name, gname

code_group=group_index_from_ext(cge)
gname= group_name(code_group)

IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
   WRITE(stdout,'(5x,"the irreducible reprepresentations of the point group")')
   WRITE(stdout,'(5x,"number ",i3,2x,a)') cge, TRIM(gname)
ENDIF
IF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==1) THEN
   WRITE(stdout,'(5x,"the irreducible representations of the double point ")')
   WRITE(stdout,'(5x,"group number ",i3,2x,a)') cge, TRIM(gname)
ENDIF
IF (ptype(1)==1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN
   WRITE(stdout,'(5x,"the projective irreducible representations with beta=-1")')
   WRITE(stdout,'(5x,"point group number ", i3,2x,a)') cge, TRIM(gname)
END IF
IF (ptype(1)==1.AND.ptype(2)==1.AND.ptype(3)==-1) THEN
   WRITE(stdout,'(5x,"the projective irreducible representations with gamma=-1 ")')
   WRITE(stdout,'(5x,"point group number ", i3,2x,a)') cge, TRIM(gname)
END IF
IF (ptype(1)==1.AND.ptype(2)==-1.AND.ptype(3)==-1) THEN
   WRITE(stdout,'(5x,"the projective irreducible representations with &
             beta=-1 and gamma=-1")')
   WRITE(stdout,'(5x,"point group number ", i3,2x,a)') cge, TRIM(gname)
END IF
IF (ptype(1)==-1.AND.ptype(2)==-1.AND.ptype(3)==1) THEN
   WRITE(stdout,'(5x,"the projective irreducible representations with beta=-1")')
   WRITE(stdout,'(5x,"double group number ", i3,2x,a)') cge, TRIM(gname)
ENDIF
IF (ptype(1)==-1.AND.ptype(2)==1.AND.ptype(3)==-1) THEN
   WRITE(stdout,'(5x,"the projective irreducible representations with gamma=-1")')
   WRITE(stdout,'(5x,"double group number ", i3,2x,a)') cge, TRIM(gname)
ENDIF
IF (ptype(1)==-1.AND.ptype(2)==-1.AND.ptype(3)==-1) THEN
   WRITE(stdout,'(5x,"the projective irreducible representations with beta=-1 and &
          &gamma=-1")') 
   WRITE(stdout,'(5x,"double group number ", i3,2x,a)') cge, TRIM(gname)
ENDIF
RETURN
END SUBROUTINE print_ptype_info

SUBROUTINE find_factor_system(sym_mat, dim_rap, nsym, cge, phase, verbosity)
!
!  This routine receives a set of nsym routines of dimensions dim_rap that
!  are supposed to be a representation (possibly projective) of the point 
!  group given by the extended code group and gives as output the factor 
!  system of the representation. If verbosity=.TRUE. it writes on 
!  output both the input matrices and the factor system. The order of the
!  matrices must be the standard order as in find_info_group_ext.
!
USE kinds, ONLY : DP
USE io_global, ONLY : stdout

IMPLICIT NONE
INTEGER, INTENT(IN) :: dim_rap, nsym, cge
COMPLEX(DP), INTENT(IN) :: sym_mat(dim_rap,dim_rap,nsym)
COMPLEX(DP), INTENT(OUT) :: phase(48,48)
LOGICAL, INTENT(IN) :: verbosity

INTEGER :: group_desc(48), prd(48,48), epos(48,48)
INTEGER :: isym, jsym, ksym, lsym, nsym_, i, j
COMPLEX(DP) :: c_mat(dim_rap,dim_rap)

CALL set_group_desc(group_desc, nsym_, cge)
CALL find_double_product_table(prd, epos, cge)

IF (verbosity) THEN
   WRITE(stdout,*)
   DO isym=1,nsym
      WRITE(stdout,'(a8)') sym_label(group_desc(isym))
      DO i=1,dim_rap
         WRITE(stdout,'(6f13.5)') (sym_mat(i,j,isym),j=1,dim_rap)
      END DO
   END DO
END IF

DO isym=1,nsym
   DO jsym=1,nsym
      c_mat = MATMUL(sym_mat(:,:,isym), sym_mat(:,:,jsym))
      DO i=1, dim_rap
         IF (ABS(sym_mat(1, i, prd(isym,jsym)))>1.D-8) THEN
            phase(isym,jsym) = c_mat(1,i) / sym_mat(1,i,prd(isym,jsym)) 
            EXIT
         END IF
         IF (i==dim_rap) CALL errore('find_factor_system',&
                                     'one row of the matrix is zero',1)
      END DO
      !
      !   check that all the rest of the matrix adjusted with the phase
      !   coincides with the product matrix
      !
      DO ksym=1,dim_rap
         DO lsym=1,dim_rap
            IF (ABS(phase(isym,jsym)*sym_mat(ksym,lsym,prd(isym,jsym))-   &
                        c_mat(ksym,lsym)) > 1.D-4)  THEN
               WRITE(stdout,'(5x,"isym, jsym",4i5)') isym, jsym, ksym, lsym
               WRITE(stdout,'(5x,"Phase",2f15.5)') phase(isym,jsym)
               WRITE(stdout,'(5x,"phase*mat",2f15.5," product",2f15.5)')  &
                     phase(isym,jsym)*sym_mat(ksym,lsym,prd(isym,jsym)),  &
                     c_mat(ksym,lsym) 
               CALL errore('find_factor_system','wrong phase',1)
            END IF
         END DO
      END DO
   END DO
END DO

IF (verbosity) CALL write_group_table(group_desc, nsym, phase)

RETURN
END SUBROUTINE find_factor_system

SUBROUTINE point_group_bravais(code_group, bl_code_group)
!
!   This routine receives the code of a point group and gives the
!   code_group of the Bravais lattice compatible with the input point group
!
!   For each row it gives the code of the last point group of the row
!
!   C_1  C_i
!   C_s  C_2  C_2h
!   D_2  C_2v D_2h
!   S_4  D_2d C_4 C_4h C_4v  D_4  D_4h
!   C_3  S_6  C_3v D_3  D_3d
!   C_3h D_3h C_6 C_6h C_6v D_6 D_6h
!   T  T_h T_d O O_h
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: code_group
INTEGER, INTENT(OUT) :: bl_code_group

INTEGER :: bl_cg(32)

DATA bl_cg  / 2,         2,         16,          16,   &
             25,        22,         23,          20,   &
             25,        22,         23,          20,   &
             25,        22,         23,          16,   &
             23,        22,         23,          20,   &
             23,        22,         23,          22,   &
             25,        22,         25,          32,   &
             32,        32,         32,          32    /

bl_code_group = bl_cg(code_group)

RETURN
END SUBROUTINE point_group_bravais

SUBROUTINE transform_group(code_group_ext_in, op_code, code_group_ext_out)
!
!  This routine find the group conjugate of the input group according to
!  the operation op_code
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: code_group_ext_in, op_code
INTEGER, INTENT(OUT) :: code_group_ext_out

INTEGER :: group_desc_in(48), nsym_in, group_desc_out(48), which_elem(48), &
           isym, code_group
REAL(DP) :: sr(3,3,48), sr_out(3,3,48), s_op(3,3)

code_group_ext_out=0

CALL set_group_desc(group_desc_in, nsym_in, code_group_ext_in)

CALL set_sym_o3(s_op, op_code)
DO isym=1,nsym_in
   CALL set_sym_o3(sr(:,:,isym),group_desc_in(isym))
   sr_out(:,:,isym)=MATMUL(TRANSPOSE(s_op), MATMUL(sr(:,:,isym),s_op)) 
ENDDO

CALL find_group_info_ext(nsym_in, sr_out, code_group, code_group_ext_out, &
                                               which_elem, group_desc_out)
RETURN
END SUBROUTINE transform_group

SUBROUTINE group_name_schoenflies(code, group_name)
!
!   This routine receives the code of a point group and writes its
!   Schoenflies name in a form that can be used to write a gnuplot
!   scripts.
!
INTEGER :: code
CHARACTER(LEN=12) :: group_name

CHARACTER(LEN=12) :: gname(32)

data gname  / "C_1         ", "C_i         ", "C_s         ", "C_2         ", &
              "C_3         ", "C_4         ", "C_6         ", "D_2         ", &
              "D_3         ", "D_4         ", "D_6         ", "C_{2v}      ", &
              "C_{3v}      ", "C_{4v}      ", "C_{6v}      ", "C_{2h}      ", &
              "C_{3h}      ", "C_{4h}      ", "C_{6h}      ", "D_{2h}      ", &
              "D_{3h}      ", "D_{4h}      ", "D_{6h}      ", "D_{2d}      ", &
              "D_{3d}      ", "S_4         ", "S_6         ", "T           ", &
              "T_h         ", "T_d         ", "O           ", "O_h         "  /

IF (code < 1 .OR. code > 32 ) CALL errore('group_name_schoenflies', &
                                               'code is out of range',1)
group_name=gname(code)

RETURN
END SUBROUTINE group_name_schoenflies

SUBROUTINE group_name_international(code, group_name)
!
!   This routine receives the code of a point group and writes its
!   international name in a form that can be used to write a gnuplot
!   scripts.
!
INTEGER :: code
CHARACTER(LEN=12) :: group_name

CHARACTER(LEN=12) :: gname(32)

data gname  /"1           ", "@^{/=24-}1  ", "m           ", "2           ", &
             "3           ", "4           ", "6           ", "222         ", &
             "32          ", "422         ", "622         ", "mm2         ", &
             "3m          ", "4mm         ", "6mm         ", "2/m         ", &
             "@^{/=24-}6  ", "4/m         ", "6/m         ", "mmm         ", &
             "@^{/=24-}62m", "4/mmm       ", "6/mmm       ", "@^{/=24-}42m", &
             "@^{/=24-}3m ", "@^{/=24-}4  ", "@^{/=24-}3  ", "23          ", &
             "m@^{/=24-}3 ", "@^{/=24-}43m", "432         ", "m@^{/=24-}3m"  /

IF (code < 1 .OR. code > 32 ) CALL errore('group_name_international', &
                                               'code is out of range',1)
group_name=gname(code)

RETURN
END SUBROUTINE group_name_international

END MODULE point_group

