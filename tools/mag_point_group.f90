!
! Copyright (C) 2017 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------
PROGRAM mag_point_group
!-------------------------------------------------------------------
!
!
USE kinds, ONLY : DP
USE mp_global,        ONLY : mp_startup, mp_global_end
USE environment,      ONLY : environment_start, environment_end
USE point_group,      ONLY : print_element_list, is_subgroup,    &
                             group_index_from_ext
USE magnetic_point_group,  ONLY : mag_group_name, magnetic_type, &
                             find_mag_code_group, is_mag_group,  &
                             find_group_subgroup_ext, a_birss_code_group, &
                             b_birss_code_group
USE io_global,        ONLY : stdout

IMPLICIT NONE

CHARACTER(LEN=9) :: code='MPG'
CHARACTER(LEN=11) :: group_name
INTEGER :: work_choice, igroup, code_group, code_subgroup, mag_group_count, &
           group_code_ext, subgroup_code_ext, mag_group_code, &
           group_desc(48), subgroup_desc(48), mag_group_code_ext, nsym, &
           mag_code
INTEGER :: stdin

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

stdin=5
WRITE(stdout,'(/,5x,"Choose what to write")')
WRITE(stdout,'(5x,"1) List the magnetic point groups (short)")')
WRITE(stdout,'(5x,"2) List the magnetic point groups (long)")')
WRITE(stdout,'(5x,"3) List the elements of a magnetic point group ")')
WRITE(stdout,'(5x,"4) List the A and B groups of magnetic point groups")')
WRITE(stdout,'(5x,"5) Find the A and B Birss groups of a magnetic group ")')
WRITE(stdout,'(5x,"6) For each point group list magnetic groups that have it as A group ")')
WRITE(stdout,'(5x,"7) For each point group list magnetic groups that have it as B group ")')

READ(stdin,*) work_choice

IF (work_choice == 1) THEN
   WRITE(stdout,*)
   WRITE(stdout,'(2(i5,3x,a21,2x,a2,4x))') (igroup, &
            TRIM(mag_group_name(igroup)), magnetic_type(igroup), igroup=1,122)
   WRITE(stdout,'(/,5x, "DP = Diamagnetic, Paramagnetic")')
   WRITE(stdout,'(5x, "F = Ferromagnetic")')
   WRITE(stdout,'(5x, "AF = AntiFerromagnetic")')
ELSEIF (work_choice == 2) THEN
   mag_group_count=0
   DO group_code_ext=1,136
      DO subgroup_code_ext=1,136
         IF (is_subgroup(group_code_ext, subgroup_code_ext)) THEN
            code_group=group_index_from_ext(group_code_ext)
            code_subgroup=group_index_from_ext(subgroup_code_ext)
            IF (is_mag_group(code_group,code_subgroup)) THEN
               CALL find_mag_code_group(code_group,code_subgroup,mag_group_code)
               mag_group_count=mag_group_count+1
               WRITE(stdout,'(2x,i4,2x,i4,3x,a21,2x,&
                             &a2,2(i5,3x,a11))') &
                       mag_group_count, mag_group_count+272,&
                      TRIM(mag_group_name(mag_group_code)),&
                       magnetic_type(mag_group_code), &
                       group_code_ext, &
                       group_name(group_index_from_ext(group_code_ext)), &
                       subgroup_code_ext, &
                       group_name(group_index_from_ext(subgroup_code_ext))

            ENDIF
         ENDIF
      ENDDO
   ENDDO
   WRITE(stdout,'(/,5x, "F = Ferromagnetic")')
   WRITE(stdout,'(5x, "AF = AntiFerromagnetic")')
   WRITE(stdout,'(5x, "Codes from 1-136: gray groups")')
   WRITE(stdout,'(5x, "Codes from 137-272: groups with no time reversal")')
   WRITE(stdout,'(5x, "Use crystal_point_group to list them")')
ELSEIF(work_choice==3) THEN
   CALL read_group_code_mag(group_desc,subgroup_desc,nsym,mag_group_code_ext)
   CALL find_group_subgroup_ext(group_code_ext,subgroup_code_ext, &
                                    mag_group_code_ext)
   WRITE(stdout,'(/,5x,"Group number ", i5)') group_code_ext
   CALL print_element_list(group_code_ext)
   IF (mag_group_code_ext>272) THEN
      WRITE(stdout,'(/,5x,"Proper subgroup of operations without time &
                     &reversal", i5)') subgroup_code_ext
      CALL print_element_list(subgroup_code_ext)
   ELSEIF(mag_group_code_ext>136) THEN
      WRITE(stdout,'(/,5x,"No operation requires time reversal")') 
   ELSEIF(mag_group_code_ext>1) THEN
      WRITE(stdout,'(/,5x,"Gray group, time reversal is a symmetry operation")') 
   ENDIF
ELSEIF(work_choice==4) THEN
   WRITE(stdout,'(/,17x,"Magnetic group", 12x, " A Group", 10x, " B Group")')
   DO mag_code=1,122
      WRITE(stdout,'(5x,i5,3x,a20,7x,a12,7x,a12)') mag_code, &
                    TRIM(mag_group_name(mag_code)), &
                    TRIM(group_name(a_birss_code_group(mag_code))), &
                    TRIM(group_name(b_birss_code_group(mag_code)))
   ENDDO
ELSEIF(work_choice==5) THEN
   WRITE(stdout,'(5x," Code of the magnetic point group (see 1) ")')
   READ(stdin,*) mag_code

   IF (mag_code<1.OR.mag_code>122) CALL errore('mag_point_group',&
                              'Wrong magnetic group code',1)
   IF (mag_code<33) THEN
      WRITE(stdout,'(/,5x,"Magnetic group (Group I): ",a)') &
                                              TRIM(mag_group_name(mag_code))
   ELSEIF (mag_code<65) THEN
      WRITE(stdout,'(/,5x,"Nonmagnetic group (Group II): ",a)') &
                                              TRIM(mag_group_name(mag_code))
   ELSE
      WRITE(stdout,'(/,5x,"Magnetic group (Group III): ",a)') &
                                              TRIM(mag_group_name(mag_code))
   ENDIF
   WRITE(stdout,'(/,5x,"A Birss group: ",a)') TRIM(group_name(&
                                              a_birss_code_group(mag_code)))
   WRITE(stdout,'(5x,"B Birss group: ",a)') TRIM(group_name(&
                                              b_birss_code_group(mag_code)))
ELSEIF(work_choice==6) THEN
   DO igroup=1,32
      WRITE(stdout,'(/,5x,a,":")') TRIM(group_name(igroup))
      DO mag_code=1,122
         IF (mag_code>32.AND.mag_code<65) CYCLE
         IF (a_birss_code_group(mag_code)==igroup) &
            WRITE(stdout,'(5x,a)') TRIM(mag_group_name(mag_code))
      ENDDO
   ENDDO
ELSEIF(work_choice==7) THEN
   DO igroup=1,32
      WRITE(stdout,'(/,5x,a,":")') TRIM(group_name(igroup))
      DO mag_code=1,122
         IF (mag_code>32.AND.mag_code<65) CYCLE
         IF (b_birss_code_group(mag_code)==igroup) &
            WRITE(stdout,'(5x,a)') TRIM(mag_group_name(mag_code))
      ENDDO
   ENDDO
ENDIF

CALL environment_end( code )
CALL mp_global_end ()

END PROGRAM mag_point_group

!--------------------------------------------------------------------------
SUBROUTINE read_group_code_mag(group_desc, subgroup_desc, nsym, &
                                                         mag_group_code_ext)
!--------------------------------------------------------------------------
USE magnetic_point_group, ONLY : find_group_subgroup_ext, mag_group_name, &
                                 mag_group_index_from_ext
USE io_global,   ONLY : stdout
USE point_group, ONLY : find_group_ext, set_group_desc, group_index_from_ext

IMPLICIT NONE
INTEGER :: mag_group_code_ext, nsym
INTEGER :: group_desc(48), subgroup_desc(48)
INTEGER :: isym, subnsym, group_code_ext, subgroup_code_ext
INTEGER :: stdin
CHARACTER(LEN=11) :: group_name

WRITE(stdout,'(/,5x,"Extended magnetic point group code &
                                              &(see the list using 2)?")')
stdin=5
READ(stdin,*) mag_group_code_ext

IF (mag_group_code_ext>589) &
   CALL errore('read_group_code_mag','input magnetic code out of range',1)

CALL find_group_subgroup_ext(group_code_ext, subgroup_code_ext, &
                                             mag_group_code_ext)

CALL set_group_desc(group_desc,nsym,group_code_ext)
subnsym=nsym/2
IF (subgroup_code_ext > 0) &
   CALL set_group_desc(subgroup_desc,subnsym,subgroup_code_ext)

WRITE(stdout,'(/,5x,"Magnetic Group number",i5,3x,a21)') mag_group_code_ext, &
            mag_group_name(mag_group_index_from_ext(mag_group_code_ext))

RETURN
END SUBROUTINE read_group_code_mag
