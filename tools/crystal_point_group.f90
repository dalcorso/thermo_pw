!
! Copyright (C) 2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
PROGRAM crystal_point_group
!--------------------------------------------------------------------
!
!
USE kinds, ONLY : DP
USE mp_global,        ONLY : mp_startup, mp_global_end
USE environment,      ONLY : environment_start, environment_end
USE point_group,      ONLY : print_element_list, group_index_from_ext, &
                             is_subgroup, find_double_product_table_from_sym, &
                             sym_label, set_group_desc, print_character_table, &
                             print_compatibility_table, set_sym_o3, &
                             set_sym_su2, product_sym_su2, compute_classes, &
                             compute_classes_double, print_kronecker_table, &
                             sym_jones, transform_group, hex_op, cub_op, &
                             group_intersection, find_irreducible, &
                             set_rotations_character, set_irr_times_d
USE io_global,        ONLY : stdout

IMPLICIT NONE

CHARACTER(LEN=11) :: group_name
CHARACTER(LEN=9) :: code='CPG', label(48)
INTEGER :: work_choice, igroup, group_index_ext, group_index_in_ext, &
           group_index_out_ext
INTEGER :: isym, jsym, nsym, itables, ntables, start, last
INTEGER :: group_desc(48), epos(48,48), prd(48,48), ptype(3), ptype1(3), &
           ptype2(3), linvs
INTEGER :: i, j, k, l, m, n, iepos, ksym, iclass, ielem, giin_ext, giout_ext
INTEGER :: nclasses, nelem(24), elem(18,24), has_e(18,24)
INTEGER :: group_index_a_ext, group_index_b_ext, group_index_c_ext, nrap, &
           nrap_gp, irap, ndim, ncount, subwork_choice
INTEGER :: stdin
INTEGER, ALLOCATABLE :: rap_list(:)
CHARACTER(LEN=45) :: name_rap(48)
CHARACTER(LEN=45), ALLOCATABLE :: name_rap_list(:)
CHARACTER(LEN=1) :: chparity
REAL(DP):: mat(3,3), jang
LOGICAL :: do_parity
COMPLEX(DP) :: cmat(2,2), character_rap(48), char_mat_proj(48,48)

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

stdin=5
WRITE(stdout,'(/,5x,"Choose what to write")')
WRITE(stdout,'(5x,"1) List the point groups (short)")')
WRITE(stdout,'(5x,"2) List the symmetry operations")')
WRITE(stdout,'(5x,"3) Product of two rotations")')
WRITE(stdout,'(5x,"4) List the point groups (long)")')
WRITE(stdout,'(5x,"5) List the elements of a point group")')
WRITE(stdout,'(5x,"6) Write the matrices of a point group elements")')
WRITE(stdout,'(5x,"7) Write the product table")')
WRITE(stdout,'(5x,"8) List the point group classes")')
WRITE(stdout,'(5x,"9) List the double group classes")')
WRITE(stdout,'(5x,"10) Write the character table")')
WRITE(stdout,'(5x,"11) List subgroups")')
WRITE(stdout,'(5x,"12) List supergroups")')
WRITE(stdout,'(5x,"13) Write one compatibility table")')
WRITE(stdout,'(5x,"14) Write all compatibility tables for one group")')
WRITE(stdout,'(5x,"15) Write all compatibility tables")')
WRITE(stdout,'(5x,"16) Decompose representations")')
WRITE(stdout,'(5x,"17) List conjugate groups")')
WRITE(stdout,'(5x,"18) Find intersection between two groups")')


READ(stdin,*) work_choice

IF (work_choice == 1) THEN
   WRITE(stdout,*)
   WRITE(stdout,'(4(i5,2x,a12))') (igroup, TRIM(group_name(igroup)), &
                                           igroup=1,32)
ELSEIF (work_choice == 2) THEN
      WRITE(stdout,'(/,5x,"Available symmetry operations:")')
      WRITE(stdout,'(5x,"Number  Label  Jones symbol cubic   hex QE        &
                                                          &hex ITA")')
      WRITE(stdout,'(5x,i3,5x,a8,2x,a41)') (isym, sym_label(isym), &
                                               sym_jones(isym), isym=1,64)
      WRITE(stdout,'(/,5x,"The number indicates the angle (2-180, 3-120, &
                                             4-90, 6-60)")')
      WRITE(stdout,'(5x,"The letters indicate the rotation axis")')
      WRITE(stdout,'(5x,"i means multiplication by inversion")')
      WRITE(stdout,'(5x,"QE indicates the hexagonal axes of Quantum ESPRESSO")')
      WRITE(stdout,'(5x,"ITA indicates those of &
                         &the International Tables for Crystallography A")')
ELSEIF (work_choice == 3) THEN
   WRITE(stdout,'(5x,"Give the number of the two rotations &
                                     &(point 2 gives the list)")')
   READ(stdin,*) isym, jsym
   CALL product_sym_su2(isym,jsym,ksym,iepos)
   CALL set_sym_o3(mat, ksym)
   CALL set_sym_su2(ksym,cmat,linvs)
   IF (iepos > 0) THEN
      WRITE(stdout,'(5x,"The product of ",a8," and ",a8," is ", a8)') &
              sym_label(isym), sym_label(jsym), sym_label(ksym)
   ELSE
      WRITE(stdout,'(5x,"The product of ",a8," and ",a8," is -", a8)') &
              sym_label(isym), sym_label(jsym), sym_label(ksym)
   ENDIF
   WRITE(stdout,'(i5,3x,a8)') ksym, sym_label(ksym)
   WRITE(stdout,'(5x,3f8.2,5x,"a=",2f6.2,3x," b=",2f6.2, i5 )') &
                                (mat(1,j),j=1,3), iepos * cmat(1,1), &
                                      iepos * cmat(1,2), linvs
   DO i=2,3
      WRITE(stdout,'(5x,3f8.2)') (mat(i,j),j=1,3)
   ENDDO
   WRITE(stdout,'(/,5x,"A minus in front of the operation name means &
                                multiplication by -E")')
   WRITE(stdout,'(5x,"a and b are the Cayley-Klein parameters")')
   WRITE(stdout,*)
ELSEIF (work_choice == 4) THEN
   DO igroup=1,136
      WRITE(stdout,'(5x,"Group number",i5,3x,a11)') igroup, &
                          group_name(group_index_from_ext(igroup)) 
      CALL print_element_list(igroup)
      WRITE(stdout,*) 
   END DO
ELSEIF (work_choice == 5) THEN
   CALL read_group_index(group_desc,nsym,group_index_ext)
   CALL print_element_list(group_index_ext)
ELSEIF (work_choice == 6) THEN
   CALL read_group_index(group_desc,nsym,group_index_ext)
   DO isym=1,nsym
      CALL set_sym_o3(mat, group_desc(isym))
      CALL set_sym_su2(group_desc(isym),cmat,linvs)
      WRITE(stdout,'(i5," - ",i2,3x,a8)') isym, group_desc(isym), &
                                        sym_label(group_desc(isym))
      WRITE(stdout,'(5x,3f8.2,5x,"a=",2f6.2,3x," b=",2f6.2, i5 )') &
                                (mat(1,j),j=1,3), cmat(1,1), cmat(1,2), linvs
      DO i=2,3
         WRITE(stdout,'(5x,3f8.2)') (mat(i,j),j=1,3)
      ENDDO
      WRITE(stdout,*)
   ENDDO
   WRITE(stdout,'(5x,"a and b are the Cayley-Klein parameters")')
ELSEIF (work_choice == 7) THEN
   CALL read_group_index(group_desc, nsym, group_index_ext)

   CALL find_double_product_table_from_sym(prd, epos, group_desc, nsym)

   WRITE(stdout,'(/,5x, "The double group product table",/)')
   ntables= nsym / 8
   IF (MOD(nsym,8) /= 0 ) ntables=ntables+1
   DO itables=1,ntables
      start=(itables-1)*8+1
      last=MIN(itables*8,nsym)
      WRITE(stdout,'(8x,8(1x,a8))') (sym_label(group_desc(jsym)),jsym=start,last)
      DO isym=1,nsym
         WRITE(stdout,'(a8,i3,2x,7(i7,2x))') sym_label(group_desc(isym)), &
        (prd(isym, jsym)*epos(isym,jsym),jsym=start,last)
      ENDDO
      WRITE(stdout,*)
   ENDDO
   WRITE(stdout,'(/,5x,"Row x column multiplication table")')
   WRITE(stdout,'(5x,"- means multiplication by -E:")')
   WRITE(stdout,'(5x, "Point group product table: neglect the minus signs")')
ELSEIF (work_choice == 8) THEN
   CALL read_group_index(group_desc, nsym, group_index_ext)
   CALL compute_classes(group_index_ext, nclasses, nelem, elem)
   WRITE(stdout,'(/,5x,"There are",i4," classes")') nclasses
   DO iclass=1,nclasses
      WRITE(stdout,'(/,5x,"Class number ",i4," has ",i4," element(s)")') &
                                          iclass, nelem(iclass)
      WRITE(stdout,'(4(i5," - ",i2,3x,a8))') (ielem, &
            group_desc(elem(ielem,iclass)), sym_label(group_desc(&
                         elem(ielem,iclass))), ielem=1, nelem(iclass))
   END DO
ELSEIF (work_choice == 9) THEN
   CALL read_group_index(group_desc, nsym, group_index_ext)
   CALL compute_classes_double(group_index_ext, nclasses, nelem, elem, has_e)
   WRITE(stdout,'(/,5x,"There are",i4," classes")') nclasses
   DO iclass=1,nclasses
      WRITE(stdout,'(/,5x,"Class number ",i4," has ",i4," element(s)")') &
                                          iclass, nelem(iclass)
      DO ielem=1, nelem(iclass)
         IF (has_e(ielem,iclass)==-1) THEN
            label(ielem)='-'//sym_label(group_desc(elem(ielem,iclass)))
         ELSE
            label(ielem)=sym_label(group_desc(elem(ielem,iclass)))
         ENDIF
      ENDDO
      WRITE(stdout,'(4(i3," - ",i2,3x,a9))') (ielem, &
            group_desc(elem(ielem,iclass)), label(ielem), ielem=1, &
                                                          nelem(iclass))
   END DO
ELSEIF (work_choice == 10) THEN
   CALL read_group_index(group_desc, nsym, group_index_ext)
   DO k=1,-1,-2
      DO j=1,-1,-2
         DO i=1,-1,-2
            ptype(1)=i
            ptype(2)=j
            ptype(3)=k
            CALL print_character_table(group_index_ext,ptype)
         END DO
      END DO
   END DO
ELSEIF (work_choice == 11) THEN
   CALL read_group_index(group_desc, nsym, group_index_ext)
   CALL print_element_list(group_index_ext)

   WRITE(stdout,'(/,5x,"Subgroup list:")') 
   DO igroup=1,136
      IF (is_subgroup(group_index_ext,igroup)) THEN
         WRITE(stdout,'(5x,"Group number",i5,3x,a11)') igroup, &
                          group_name(group_index_from_ext(igroup)) 
         CALL print_element_list(igroup)
         WRITE(stdout,*) 
      ENDIF
   ENDDO   
ELSEIF (work_choice == 12) THEN
   CALL read_group_index(group_desc, nsym, group_index_ext)
   CALL print_element_list(group_index_ext)

   WRITE(stdout,'(/,5x,"Supergroups list:")') 
   DO igroup=1,136
      IF (is_subgroup(igroup,group_index_ext)) THEN
         WRITE(stdout,'(5x,"Group number",i5,3x,a11)') igroup, &
                          group_name(group_index_from_ext(igroup)) 
         CALL print_element_list(igroup)
         WRITE(stdout,*) 

      ENDIF
   ENDDO   
ELSEIF (work_choice == 13 .OR. work_choice == 14 .OR. work_choice==15) THEN
   giin_ext=0
   giout_ext=0
   IF (work_choice /= 15) &
      CALL read_group_index(group_desc,nsym,giin_ext)
   IF (work_choice ==13) &
      CALL read_group_index(group_desc,nsym,giout_ext)
   IF (work_choice==13.AND..NOT.is_subgroup(giin_ext, giout_ext)) THEN
      WRITE(stdout,'(/,5x,"Group ",i3,&
                " is not a subgroup of group ",i3)') giout_ext, giin_ext
      WRITE(stdout,'(/,5x,"Elements list of group ",i3,":")') giin_ext
      CALL print_element_list(giin_ext)
      WRITE(stdout,'(/,5x,"Elements list of group ",i3,":")') giout_ext
      CALL print_element_list(giout_ext)
   ENDIF
   DO group_index_in_ext=1,136
      DO group_index_out_ext=1,136
         IF ( group_index_in_ext /= giin_ext .AND. giin_ext /= 0) CYCLE
         IF ( group_index_out_ext /= giout_ext .AND. giout_ext /= 0) CYCLE
         IF (.NOT.is_subgroup(group_index_in_ext, group_index_out_ext)) CYCLE

         WRITE(stdout,'(/,5x,70("-"))')
         WRITE(stdout,'(/,5x,"Group number   ", i4, 2x, a11, &
                               &" : Subgroup number", i4, 2x, a11)')  &
                 group_index_in_ext,   &
                 group_name(group_index_from_ext(group_index_in_ext)), &
                 group_index_out_ext,    &
                 group_name(group_index_from_ext(group_index_out_ext))

         DO k=1,-1,-2
            DO j=1,-1,-2
               DO i=1,-1,-2
                  ptype(1)=i
                  ptype(2)=j
                  ptype(3)=k
                  CALL print_compatibility_table(group_index_in_ext,ptype,&
                                                          group_index_out_ext)
               END DO
            END DO
         END DO
      END DO
   END DO
ELSEIF (work_choice == 16) THEN

   WRITE(stdout,'(5x,"1) Decompose Kronecker products table (chi x chi)")')
   WRITE(stdout,'(5x,"2) Decompose Kronecker products table (chi^* x chi)")')
   WRITE(stdout,'(5x,"3) Find the representations of translations")')
   WRITE(stdout,'(5x,"4) Find the representations of rotations")')
   WRITE(stdout,'(5x,"5) Decompose the representations of the full &
                                           &rotation group D(j)(+-)")')
   WRITE(stdout,'(5x,"6) Decompose the product of the group representations &
                                               &chi and D(1/2)(+)")')
   WRITE(stdout,'(5x,"7) Decompose the product of the group representations &
                                               &chi and D(j)(+-)")')
   READ(stdin,*) subwork_choice 
   IF (subwork_choice ==1 .OR. subwork_choice==2) THEN
      CALL read_group_index(group_desc,nsym,group_index_in_ext)
      CALL read_group_index(group_desc,nsym,group_index_out_ext)

!     DO group_index_in_ext=1,136
!        DO group_index_out_ext=1,136
            IF (is_subgroup(group_index_in_ext, group_index_out_ext)) THEN
               WRITE(stdout,'(/,5x,"Group number   ", i4, 2x, a11, &
                                  &" : Subgroup number", i4, 2x, a11)')  &
                       group_index_in_ext,   &
                       group_name(group_index_from_ext(group_index_in_ext)), &
                       group_index_out_ext,    &
                       group_name(group_index_from_ext(group_index_out_ext))

               DO k=1,-1,-2
                  DO j=1,-1,-2
                     DO i=1,-1,-2
                        DO l=1,-1,-2
                           DO m=1,-1,-2
                              DO n=1,-1,-2
                                 ptype1(1)=n
                                 ptype1(2)=l
                                 ptype1(3)=j
                                 ptype2(1)=m
                                 ptype2(2)=i
                                 ptype2(3)=k
                                 IF (m>n) CYCLE
                                 IF (i>l) CYCLE
                                 IF (k>j) CYCLE
                                 IF (work_choice==14) THEN
                                    CALL print_kronecker_table(group_index_in_ext,&
                                   ptype1,ptype2,group_index_out_ext,.FALSE.)
                                 ELSE
                                    CALL print_kronecker_table(group_index_in_ext,&
                                      ptype1,ptype2,group_index_out_ext,.TRUE.)
                                 END IF
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            ENDIF
!        ENDDO
!     ENDDO
   ELSEIF (subwork_choice==3.OR.subwork_choice==4) THEN
      CALL read_group_index(group_desc,nsym,group_index_ext)

!DO group_index_ext=1,136
!   WRITE(stdout,'(/,5x,"Group number   ", i4, 2x, a11)') &
!                    group_index_ext,   &
!                    group_name(group_index_from_ext(group_index_ext))
!
!
!   translations are vectors, so they have the same character table of the rotations 
!   with l=1 and odd parity, rotations are axial vectors so the parity is even
!
      IF (subwork_choice==3) THEN
         CALL set_rotations_character(group_index_ext, character_rap,1.0_DP,.TRUE.)
      ELSE
         CALL set_rotations_character(group_index_ext, character_rap,1.0_DP,.FALSE.)
      ENDIF

      ndim=3
      ALLOCATE(rap_list(ndim))
      ALLOCATE(name_rap_list(ndim))
      ptype=1
      CALL find_irreducible(group_index_ext, ptype, character_rap, &
                                rap_list, name_rap_list, nrap, ndim)

      IF (subwork_choice==3) THEN
         WRITE(stdout,'(/,5x,"Translations are:")') 
      ELSE
         WRITE(stdout,'(/,5x,"Rotations are:")') 
      ENDIF
      CALL print_rap(nrap, name_rap_list)
      DEALLOCATE(name_rap_list)
      DEALLOCATE(rap_list)
!ENDDO
   ELSEIF (subwork_choice==5) THEN
      CALL read_group_index(group_desc,nsym,group_index_ext)

      WRITE(stdout,'(5x,"Angular momentum j (0, 1, 2, 3, ... for s,p,d,f... &
             or 0.5, 1.5, 2.5, ...) ?")')
      READ(stdin,*) jang

      WRITE(stdout,'(/,5x,"Parity of D(j): + or - ?")') 
      READ(stdin,*) chparity
      do_parity=.FALSE.
      IF (chparity=='-') do_parity=.TRUE.

!DO group_index_ext=1,136
!   WRITE(stdout,'(/,5x,"Group number   ", i4, 2x, a11)') &
!                    group_index_ext,   &
!                    group_name(group_index_from_ext(group_index_ext))

      CALL set_rotations_character(group_index_ext, character_rap, jang, do_parity)

      ptype=1
      ndim=NINT(2.0_DP*jang+1.0_DP)
      ALLOCATE(rap_list(ndim))
      ALLOCATE(name_rap_list(ndim))
!
!  Use the double group to decompose the half integer representations
!
      IF ((NINT(jang)-jang)> 1.D-8) ptype(1)=-1
      CALL find_irreducible(group_index_ext, ptype, character_rap, &
                                    rap_list, name_rap_list, nrap, ndim)

      IF (do_parity) THEN
         WRITE(stdout,'(/,5x,"D(",f6.2,"(-) decomposes into:")') jang
      ELSE
         WRITE(stdout,'(/,5x,"D(",f6.2,"(+) decomposes into:")') jang
      ENDIF
   
      CALL print_rap(nrap, name_rap_list)

      DEALLOCATE(rap_list)
      DEALLOCATE(name_rap_list)
!ENDDO
   ELSEIF (subwork_choice==6.OR.subwork_choice==7) THEN
      CALL read_group_index(group_desc,nsym,group_index_ext)
!
!   Spin does not change for parity so use the (+) representation
!
      ptype1=1
      IF (subwork_choice==7) THEN
         WRITE(stdout,'(/,5x,"Point group (1) or double point group (-1) &
                                                      &representations?")')
         READ(stdin,*) ptype1(1)
         IF (ptype1(1)/=1.AND.ptype1(1)/=-1) ptype1(1)=1

         WRITE(stdout,'(5x,"Angular momentum j (0, 1, 2, 3, ... for s,p,d,f... &
             or 0.5, 1.5, 2.5, ...) ?")')
         READ(stdin,*) jang

         WRITE(stdout,'(5x,"Parity of D(j): + or - ?")')
         READ(stdin,*) chparity
         do_parity=.FALSE.
         IF (chparity=='-') do_parity=.TRUE.
      ELSE
         jang=0.5_DP
         do_parity=.FALSE.
      ENDIF
!DO group_index_ext=1,136
!   WRITE(stdout,'(/,5x,"Group number   ", i4, 2x, a11)') &
!                    group_index_ext,   &
!                    group_name(group_index_from_ext(group_index_ext))

      CALL set_irr_times_d(group_index_ext, ptype1, char_mat_proj, name_rap, &
                                            nrap_gp, nsym, jang, do_parity)
      ptype=ptype1
      IF (MOD(NINT(2.0_DP*jang),2)==1) ptype(1)=-1*ptype1(1)
      DO irap=1,nrap_gp
         IF (subwork_choice==6) THEN
            WRITE(stdout,'(/,5x,a," times D(1/2)(+) decomposes into:")') &
                              TRIM(name_rap(irap)(1:8))
         ELSE
            IF (do_parity) THEN
               WRITE(stdout,'(/,5x,a," times D(",f6.1,")(-) &
                         &decomposes into:")') TRIM(name_rap(irap)(1:8)), jang
            ELSE
               WRITE(stdout,'(/,5x,a," times D(",f6.1,")(+) &
                             &decomposes into:")')  TRIM(name_rap(irap)(1:8)),& 
                             jang
            ENDIF 
         ENDIF
         character_rap(1:nsym) = char_mat_proj(irap,1:nsym)
         ndim=NINT(REAL(char_mat_proj(irap,1)))
         ALLOCATE(rap_list(ndim))
         ALLOCATE(name_rap_list(ndim))
         CALL find_irreducible(group_index_ext, ptype, character_rap, &
                                    rap_list, name_rap_list, nrap, ndim)
         CALL print_rap(nrap, name_rap_list)

         DEALLOCATE(rap_list)
         DEALLOCATE(name_rap_list)
      ENDDO
!    ENDDO
   ENDIF
ELSEIF (work_choice==17) THEN
   CALL read_group_index(group_desc,nsym,group_index_in_ext)

   WRITE(stdout,*)
   DO ielem=1,64
      IF (hex_op(ielem).AND.is_subgroup(111,group_index_in_ext) &
        .OR. (cub_op(ielem).AND.is_subgroup(136,group_index_in_ext))) THEN
         CALL transform_group(group_index_in_ext, ielem, group_index_out_ext)
         IF (group_index_out_ext > 0) THEN
            WRITE(stdout,'(5x,"Symmetry ",i4,2x,a8," => ", i4,2x,a11)') &
            ielem, TRIM(sym_label(ielem)), group_index_out_ext, &
                  group_name(group_index_from_ext(group_index_out_ext)) 
         END IF
      END IF
   END DO
ELSEIF (work_choice==18) THEN
   CALL read_group_index(group_desc,nsym,group_index_a_ext)
   CALL read_group_index(group_desc,nsym,group_index_b_ext)

   CALL group_intersection(group_index_a_ext, group_index_b_ext, &
                                              group_index_c_ext)
   WRITE(stdout,'(/,5x,"The intersection between group number",i5,3x,a11)') &
                         group_index_a_ext,  &
                         group_name(group_index_from_ext(group_index_a_ext))
   WRITE(stdout,'(5x,"with elements")') 
   CALL print_element_list(group_index_a_ext)
   WRITE(stdout,'(/,5x,"and group number",i5,3x,a11)') group_index_b_ext, &
                         group_name(group_index_from_ext(group_index_b_ext))
   WRITE(stdout,'(5x,"with elements")') 
   CALL print_element_list(group_index_b_ext)
   WRITE(stdout,'(/,5x,"is group number",i5,3x,a11)') group_index_c_ext, &
                         group_name(group_index_from_ext(group_index_c_ext))
   WRITE(stdout,'(/,5x,"with elements")') 
   CALL print_element_list(group_index_c_ext)

ENDIF

CALL environment_end( code )
CALL mp_global_end ()

END PROGRAM crystal_point_group

!--------------------------------------------------------------------
SUBROUTINE read_group_index(group_desc, nsym, group_index_ext)
!--------------------------------------------------------------------
!
USE point_group, ONLY : group_index_from_ext
USE io_global, ONLY : stdout
USE point_group, ONLY : find_group_ext, set_group_desc

IMPLICIT NONE
INTEGER :: group_index_ext, nsym
INTEGER :: group_desc(48)
INTEGER :: isym 
INTEGER :: stdin
CHARACTER(LEN=11) :: group_name

stdin=5
WRITE(stdout,'(/,5x,"Extended point group code (1-136) (see the list using 4)?")')

READ(stdin,*) group_index_ext
CALL set_group_desc(group_desc,nsym,group_index_ext)

WRITE(stdout,'(/,5x,"Group number",i5,3x,a11)') group_index_ext, &
                          group_name(group_index_from_ext(group_index_ext)) 

RETURN
END SUBROUTINE read_group_index

!-----------------------------------------------------------------------
SUBROUTINE print_rap(nrap, name_rap_list)
!-----------------------------------------------------------------------

USE io_global, ONLY : stdout
IMPLICIT NONE

INTEGER, INTENT(IN) :: nrap
CHARACTER(LEN=45) :: name_rap_list(nrap)

INTEGER :: ncount, irap
CHARACTER(LEN=45) :: current_name

IF (nrap>0) THEN
   ncount=1
   current_name=name_rap_list(1)(1:8)
ENDIF
IF (nrap==1) WRITE(stdout,'(9x,a)') TRIM(current_name)
DO irap=2, nrap
   IF (name_rap_list(irap)(1:8) /= current_name) THEN
      IF (ncount==1) THEN
         WRITE(stdout,'(9x,a)') TRIM(current_name)
      ELSE
         WRITE(stdout,'(5x,i3,1x,a)') ncount, TRIM(current_name)
      ENDIF
      ncount=1
      current_name=name_rap_list(irap)(1:8)
      IF (irap==nrap) WRITE(stdout,'(9x,a)') TRIM(current_name)
   ELSEIF (irap==nrap) THEN
      ncount=ncount+1
      WRITE(stdout,'(5x,i3,1x,a)') ncount, TRIM(current_name)
   ELSE
      ncount=ncount+1
   ENDIF
ENDDO
RETURN
END SUBROUTINE print_rap
