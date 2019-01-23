!
! Copyright (C) 2010 Quantum ESPRESSO group
! Copyright (C) 2018 Andrea Dal Corso (extended to noncollinear magnetic
!                  systems)
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
  SUBROUTINE prepare_sym_analysis_tpw(nsym,sr,sname,t_rev,magnetic_sym)

  USE kinds,    ONLY : DP
  USE rap_point_group,  ONLY : code_group, nclass, nelem, elem, which_irr,  &
                               char_mat, name_rap, gname, name_class, ir_ram, &
                               elem_name
  USE rap_point_group_is, ONLY : code_group_is, gname_is

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nsym
  REAL(DP), INTENT(IN) :: sr(3,3,nsym)
  INTEGER, INTENT(IN) :: t_rev(nsym)
  LOGICAL, INTENT(IN) :: magnetic_sym
  CHARACTER(len=45), INTENT(IN) :: sname(48)

  INTEGER :: nsym_is, isym
  REAL(DP) :: sr_is(3,3,48)
  CHARACTER(len=45) :: sname_is(48)
!
!  Find the group name and sets its irreducible representation in the
!  rap_point_group module variables
!
  CALL find_group(nsym,sr,gname,code_group)
!
!  If some symmetry needs the time reversal check which group is formed
!  by the operations that do not need time reversal. Presently only this
!  group is used to classify the phonon modes.
!
  IF (magnetic_sym) THEN
     nsym_is=0
     DO isym=1,nsym
        IF (t_rev(isym)==0) THEN
           nsym_is=nsym_is+1
           sr_is(:,:,nsym_is) = sr(:,:,isym)
           sname_is(nsym_is) = sname(isym)
        ENDIF
     ENDDO
     CALL find_group(nsym_is,sr_is,gname_is,code_group_is)
     CALL set_irr_rap(code_group_is,nclass,char_mat,name_rap,name_class,ir_ram)
     CALL divide_class(code_group_is,nsym_is,sr_is,nclass,nelem,elem,which_irr)
     CALL set_class_el_name(nsym_is,sname_is,nclass,nelem,elem,elem_name)
  ELSE
     CALL set_irr_rap(code_group,nclass,char_mat,name_rap,name_class,ir_ram)
     CALL divide_class(code_group,nsym,sr,nclass,nelem,elem,which_irr)
     CALL set_class_el_name(nsym,sname,nclass,nelem,elem,elem_name)
  ENDIF

  RETURN
  END SUBROUTINE prepare_sym_analysis_tpw
