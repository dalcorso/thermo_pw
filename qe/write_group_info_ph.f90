!
! Copyright (C) 2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE write_group_info_ph(flag)
!----------------------------------------------------------------------------
!
! This routine writes on output the main information on the point group
! If flag is .false. writes only the character table. If flag is .true.
! writes also the elements of each class.
! This routine differs from the one contained in divide_class_so because
! it never uses the information on the double point group, since the phonon
! are classified with the representations of the point group
!
!
USE rap_point_group,      ONLY : code_group, nclass, nelem, elem, which_irr, &
                                 char_mat, name_rap, name_class, gname,      &
                                 elem_name
USE io_global,            ONLY : stdout

IMPLICIT NONE

INTEGER :: iclass, irot, i, idx, irap
LOGICAL :: is_complex, is_complex_so, flag

WRITE(stdout,'(/,5x,"point group ",a11)') gname
WRITE(stdout,'(5x, "there are", i3," classes")') nclass

WRITE(stdout,'(5x, "the character table:")')
WRITE(stdout,'(/,7x,12(a5,1x))') (name_class(irot),irot=1,nclass)
DO iclass=1,nclass
   WRITE(stdout,'(a5,12f6.2)') name_rap(iclass), &
      (REAL(char_mat(iclass,irot)),irot=1,nclass)
ENDDO
idx=code_group
IF (is_complex(idx)) THEN
   WRITE(stdout,'(5x,"imaginary part")')
   DO iclass=1,nclass
      WRITE(stdout,'(a5,12f6.2)') name_rap(iclass), &
           (AIMAG(char_mat(iclass,irot)),irot=1,nclass)
   ENDDO
ENDIF
IF (flag) THEN
   WRITE(stdout,'(/5x, "the symmetry operations in each class and &
                        &the name of the first element:",/)')
   DO irap = 1, nclass
      DO iclass=1,nclass
         IF (which_irr(iclass)/=irap) CYCLE
         WRITE(stdout,'(5x,a5,12i5)') name_class(which_irr(iclass)), &
            (elem(i,iclass), i=1,nelem(iclass))
!
!    The name of the first element of each class is written on output
!
         WRITE(stdout,'(10x,a)') elem_name(1,iclass)
      ENDDO
   ENDDO
END IF

RETURN
END SUBROUTINE write_group_info_ph

