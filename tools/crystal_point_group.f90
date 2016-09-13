!
! Copyright (C) 2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM crystal_point_group
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
                             compute_classes_double, print_kroneker_table, &
                             kovalev_cubic, kovalev_hexagonal
USE io_global,        ONLY : stdout

IMPLICIT NONE

CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=11) :: group_name
CHARACTER(LEN=9) :: code='CPG', label(48)
INTEGER :: work_choice, igroup, group_index_ext, group_index_in_ext, &
           group_index_out_ext
INTEGER :: isym, jsym, nsym, itables, ntables, start, last
INTEGER :: group_desc(48), epos(48,48), prd(48,48), ptype(3), ptype1(3), &
           ptype2(3), linvs
INTEGER :: i, j, k, l, m, n, iepos, ksym, iclass, ielem, giin_ext, giout_ext
INTEGER :: nclasses, nelem(24), elem(18,24), has_e(18,24)
REAL(DP):: mat(3,3)
COMPLEX(DP) :: cmat(2,2)

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

WRITE(stdout,'(/,5x,"Choose what to write")')
WRITE(stdout,'(5x,"1) List point groups (short)")')
WRITE(stdout,'(5x,"2) List symmetry operations")')
WRITE(stdout,'(5x,"3) Product of two rotations")')
WRITE(stdout,'(5x,"4) List point groups (long)")')
WRITE(stdout,'(5x,"5) List the elements of a point group")')
WRITE(stdout,'(5x,"6) Write the matrices of a point group elements")')
WRITE(stdout,'(5x,"7) Write product table")')
WRITE(stdout,'(5x,"8) List the point group classes")')
WRITE(stdout,'(5x,"9) List the double group classes")')
WRITE(stdout,'(5x,"10) Write character table")')
WRITE(stdout,'(5x,"11) List subgroups")')
WRITE(stdout,'(5x,"12) List supergroups")')
WRITE(stdout,'(5x,"13) Write one compatibility table")')
WRITE(stdout,'(5x,"14) Write all compatibility tables for one group")')
WRITE(stdout,'(5x,"15) Write all compatibility tables")')
WRITE(stdout,'(5x,"16) Decompose Kroneker products table (chi x chi)")')
WRITE(stdout,'(5x,"17) Decompose Kroneker products table (chi^* x chi)")')
WRITE(stdout,'(5x,"18) Write Kovalev symmetry operations list")')
WRITE(stdout,'(5x,"19) Write Kovalev cubic product table")')
WRITE(stdout,'(5x,"20) Write Kovalev hexagonal product table")')

READ(5,*) work_choice

IF (work_choice == 1) THEN
   WRITE(stdout,*)
   WRITE(stdout,'(4(i5,2x,a12))') (igroup, TRIM(group_name(igroup)), &
                                           igroup=1,32)
ELSEIF (work_choice == 2) THEN
      WRITE(stdout,'(/,5x,"Available symmetry operations:",/)')
      WRITE(stdout,'(6(i3,2x,a8))') (isym, sym_label(isym), isym=1,64)
      WRITE(stdout,'(/,5x,"The number indicates the angle (2-180, 3-120, &
                                             4-90, 6-60)")')
      WRITE(stdout,'(5x,"The letters indicate the rotation axis")')
      WRITE(stdout,'(5x,"i means multiplication by inversion")')
ELSEIF (work_choice == 3) THEN
   WRITE(stdout,'(5x,"Give the number of the two rotations &
                                     &(point 2 gives the list)")')
   READ(5,*) isym, jsym
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
ELSEIF (work_choice == 16 .OR. work_choice==17) THEN
   CALL read_group_index(group_desc,nsym,group_index_in_ext)
   CALL read_group_index(group_desc,nsym,group_index_out_ext)

!   DO group_index_in_ext=1,136
!      DO group_index_out_ext=1,136
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
                                 CALL print_kroneker_table(group_index_in_ext,&
                                   ptype1,ptype2,group_index_out_ext,.FALSE.)
                              ELSE
                                 CALL print_kroneker_table(group_index_in_ext,&
                                   ptype1,ptype2,group_index_out_ext,.TRUE.)
                              END IF
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         ENDIF
!      ENDDO
!   ENDDO
ELSEIF (work_choice == 18 ) THEN
WRITE(stdout,'(/,5x, "Kovalev symmetries - cubic groups",/)')
DO isym=1,48
   WRITE(stdout,'(5x,i5," - ",i3,3x,a8)') isym, kovalev_cubic(isym), &
                                        sym_label(kovalev_cubic(isym))
ENDDO
WRITE(stdout,'(/,5x, "Kovalev symmetries - hexagonal groups",/)')
DO isym=1,24
   WRITE(stdout,'(5x,i5," - ",i3,3x,a8)') isym, kovalev_hexagonal(isym), &
                                        sym_label(kovalev_hexagonal(isym))
ENDDO

ELSEIF (work_choice == 19 .OR. work_choice==20 ) THEN
IF (work_choice==19) THEN
   nsym=48
   group_index_ext=136
   group_desc=kovalev_cubic
ELSE
   nsym=24
   group_index_ext=111
   group_desc(1:24)=kovalev_hexagonal
ENDIF

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

END IF

CALL environment_end( code )
CALL mp_global_end ()

END PROGRAM crystal_point_group

SUBROUTINE read_group_index(group_desc, nsym, group_index_ext)
USE point_group, ONLY : group_index_from_ext
USE io_global, ONLY : stdout
USE point_group, ONLY : find_group_ext, set_group_desc

IMPLICIT NONE
INTEGER :: group_index_ext, nsym
INTEGER :: group_desc(48)
INTEGER :: isym 
CHARACTER(LEN=11) :: group_name
WRITE(stdout,'(/,5x,"Extended point group code (see the list using 4)?")')
WRITE(stdout,'(5x,"Set 0 to give the number of symmetry operations and &
                                                    &their list.")')
READ(5,*) group_index_ext

IF (group_index_ext==0) THEN
   READ(5,*) nsym
   READ(5,*) (group_desc(isym), isym=1,nsym)
   CALL find_group_ext(group_desc,nsym,group_index_ext)
ELSE
   CALL set_group_desc(group_desc,nsym,group_index_ext)
ENDIF

WRITE(stdout,'(/,5x,"Group number",i5,3x,a11)') group_index_ext, &
                          group_name(group_index_from_ext(group_index_ext)) 

RETURN
END SUBROUTINE read_group_index
