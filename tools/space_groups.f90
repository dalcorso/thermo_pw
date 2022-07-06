!
! Copyright (C) 2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM crystal_space_groups
!
!  This is a simple code that reads three direct lattice vectors and
!  finds ibrav and celldm of that lattice. It can be used to transform
!  any triplet of primitive vectors in the input needed by thermo_pw or
!  QE.
!  Alternatively it checks the routine find_ibrav of the 
!  thermo_pw module by generating all the Bravais lattice vectors 
!  and checking that the module identify correctly the bravais lattice, then
!  it applies an arbitrary rotation and again checks that 
!  the routine find_ibrav identify correctly the lattice.
!
USE kinds, ONLY : DP
USE mp_global,        ONLY : mp_startup, mp_global_end
USE environment,      ONLY : environment_start, environment_end
USE point_group,      ONLY : set_group_desc, sym_label
USE space_groups,     ONLY : set_point_group_code, set_standard_sg,    &
                             find_space_group_names, find_space_group, &
                             set_add_info, find_space_group_number,    &
                             find_sg_name => sg_name, symmorphic_sg
USE io_global,        ONLY : stdout

IMPLICIT NONE

CHARACTER(LEN=9) :: code='Space'
CHARACTER(LEN=11) :: gname
REAL(DP), PARAMETER :: eps1=1.D-8, eps2=1.D-5
INTEGER :: group_desc(48), group_code, group_code_ext, work_choice, ibrav, &
           sgc, sgc_, ssgc, nsym, isym, ivec, counter, aux_sg, ftau(3,48)
INTEGER :: stdin
REAL(DP) :: celldm(6), sr(3,3,48), ft(3,48), s01(3), s02(3), at(3,3), bg(3,3), &
            omega
CHARACTER(LEN=11) :: group_name, sg_name
CHARACTER(LEN=12) :: sg_name1

CALL mp_startup (start_images=.TRUE.)
CALL environment_start(code)

stdin=5
WRITE(stdout,'(/,5x,"Choose what to do")')
WRITE(stdout,'(5x,"1) Write the space group names")')
WRITE(stdout,'(5x,"2) Write the space group number")')
WRITE(stdout,'(5x,"3) Write the symmetry elements of a space group")')
WRITE(stdout,'(5x,"4) Write symmorphic space groups")')

READ(stdin,*) work_choice

celldm(1)=1.0_DP
celldm(2)=1.2_DP
celldm(3)=1.4_DP
celldm(4)=-0.4_DP
celldm(5)=-0.3_DP
celldm(6)=-0.2_DP


IF (work_choice==1) THEN
   WRITE(stdout,'(5x,"Space group number?")')
   READ(stdin,*) sgc

   CALL set_add_info()

   IF (sgc > 0.AND. sgc <= 230 ) THEN
      WRITE(stdout,'(5x,"The name(s) of the space group ",i4, " is (are):",/)')&
                                                                  sgc
      CALL find_space_group_names(sgc)
      WRITE(stdout,'(/,5x,"The first name is in the ITA tables. ")')
      WRITE(stdout,'(5x,"s short name, S Schoenflies symbol. ")') 
      WRITE(stdout,'(5x,"* (#) name used in 1935 (1952) edition of IT.")')
      IF (sgc > 15 .AND. sgc < 75) THEN
         WRITE(stdout,'(5x,"Orthorombic group.")')
         WRITE(stdout,'(5x,"Apply the second transformation to the &
                                         &coordinates and")')
         WRITE(stdout,'(5x,"to celldm to pass to the ITA group, and the first &
                    &transformation")')
         WRITE(stdout,'(5x,"to return to the given group.")')
      ENDIF
   ELSE
      WRITE(stdout,'(5x,"Group number must be between 1 and 230" )')
   ENDIF
ELSEIF (work_choice==2) THEN
   WRITE(stdout,'(5x,"Space group name?")')
   READ(stdin,'(a)') sg_name

   CALL set_add_info()
   CALL find_space_group_number(sg_name,sgc)
   IF (sgc > 0 ) THEN
      WRITE(stdout,'(5x,"The space group number of ",a," is ",i5)') &
                                             TRIM(sg_name), sgc
      WRITE(stdout,'(5x,"Other names of the same group:",/)')
      CALL find_space_group_names(sgc)
      WRITE(stdout,'(/,5x,"The first name is in the ITA tables. ")')
      WRITE(stdout,'(5x,"s short name, S Schoenflies symbol. ")')
      WRITE(stdout,'(5x,"* (#) name used in 1935 (1952) edition of IT.")')
      IF (sgc > 15 .AND. sgc < 75) THEN
        WRITE(stdout,'(5x,"Orthorombic group.")')
        WRITE(stdout,'(5x,"Apply the second transformation to the coordinates &
                                &and")')
         WRITE(stdout,'(5x,"to celldm to pass to the ITA group, and the first &
                    &transformation")')
         WRITE(stdout,'(5x,"to return to the given group.")')
      ENDIF
   ELSE
      WRITE(stdout,'(5x,"This group name is unknown" )')
   ENDIF

ELSEIF (work_choice==3) THEN
   WRITE(stdout,'(5x,"Input space_group_number? ")')
   READ(stdin,*) sgc

   CALL set_add_info()

!   DO sgc=1,230

   WRITE(stdout,'(/,5x,"This is space group number",i5,/)') sgc
   CALL set_point_group_code(sgc, group_code, group_code_ext)
   CALL set_group_desc(group_desc, nsym, group_code_ext)
   CALL set_standard_sg(sgc,ibrav,celldm,nsym,sr,ft)

   CALL latgen(ibrav,celldm,at(1,1), at(1,2), at(1,3), omega)
   CALL recips(at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))

   CALL find_space_group_names(sgc)

   CALL find_space_group(ibrav, nsym, sr, ft, at, bg, sgc_, aux_sg, &
                                                            s01, s02,.TRUE.)

   CALL find_sg_name(sgc_, aux_sg, sg_name1)
   WRITE(stdout,'(/,5x,"input/output space group",i5," / ",i4, 2x, a)') sgc, &
                                                      sgc_, TRIM(sg_name1)

   IF (sgc /= sgc_) CALL errore('space_group','problem with space',1)

!   ENDDO

ELSEIF (work_choice==4) THEN
   WRITE(stdout,*)
   counter=0
   DO sgc=1,230
      IF (symmorphic_sg(sgc,ssgc)) THEN
         counter=counter+1
         CALL find_sg_name(sgc, 1, sg_name1)
         WRITE(stdout,'(5x,i5,"  space group number", i5, 2x, a)') counter, &
                                                  sgc, TRIM(sg_name1)
      END IF
   END DO
ENDIF

CALL environment_end(code)
CALL mp_global_end ()

END PROGRAM crystal_space_groups

