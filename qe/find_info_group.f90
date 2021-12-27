!
! Copyright (C) 2006-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE find_info_group_tpw(nsym,s,t_rev,ft,d_spink,gk,sname,  &
     s_is,d_spin_is,gk_is, invs_is,is_symmorphic,search_sym)
  !
  ! This routine receives as input a point group and sets the corresponding
  ! variables for the description of the classes and of the irreducible
  ! representations. It sets also the group name and code.
  ! In the magnetic case it selects the invariat subgroup.
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at, bg
  USE fft_base,             ONLY : dfftp
  USE noncollin_module,     ONLY : noncolin, domag
  USE rap_point_group,      ONLY : code_group, nclass, nelem, elem, which_irr, &
       char_mat, name_rap, name_class, gname, ir_ram
  USE rap_point_group_so,   ONLY : nrap, nelem_so, elem_so, has_e, &
       which_irr_so, char_mat_so, name_rap_so, &
       name_class_so, d_spin, name_class_so1
  USE rap_point_group_is,   ONLY : nsym_is, sr_is, ft_is, gname_is, &
       sname_is, code_group_is

  IMPLICIT NONE

  INTEGER, INTENT(in) :: nsym,        & ! dimension of the group
       s(3,3,48),   & ! rotation matrices
       t_rev(48),   & ! if time reversal is need
       gk(3,48)
  REAL(dp), INTENT(IN) :: ft(3,48)! fractionary translation

  INTEGER, INTENT(out) :: s_is(3,3,48),   & ! rotation matrices
       gk_is(3,48), invs_is(48)

  COMPLEX(DP),INTENT(out)   :: d_spink(2,2,48),  & ! rotation in spin space
       d_spin_is(2,2,48)   ! rotation in spin space

  LOGICAL, INTENT(out) :: is_symmorphic, &  ! true if the gruop is symmorphic
       search_sym        ! true if gk

  CHARACTER(len=45), INTENT(in) :: sname(48)

  REAL(DP) :: sr(3,3,48)
  INTEGER :: isym, jsym, ss(3,3)
  LOGICAL :: found


  is_symmorphic=.true.
  search_sym=.true.

  DO isym=1,nsym
     is_symmorphic=( is_symmorphic.and.(ft(1,isym)==0.d0).and.  &
          (ft(2,isym)==0.d0).and.  &
          (ft(3,isym)==0.d0) )
  ENDDO

  IF (.NOT.is_symmorphic) THEN

     DO isym=1,nsym
        DO jsym=1,nsym
           search_sym=search_sym.AND.(ABS(gk(1,isym)*ft(1,jsym)+ &
                     gk(2,isym)*ft(2,jsym)+gk(3,isym)*ft(3,jsym))<1.D-8) 
        ENDDO
     ENDDO
  ENDIF
  !
  !  Set the group name, divide it in classes and set the irreducible
  !  representations
  !
  nsym_is=0
  DO isym=1,nsym
     CALL s_axis_to_cart (s(1,1,isym), sr(1,1,isym), at, bg)
     IF (noncolin) THEN
        !
        !  In the noncollinear magnetic case finds the invariant subgroup of the point
        !  group of k. Presently we use only this subgroup to classify the levels.
        !
        IF (domag) THEN
           IF (t_rev(isym)==0) THEN
              nsym_is=nsym_is+1
              CALL s_axis_to_cart (s(1,1,isym), sr_is(1,1,nsym_is), at, bg)
              CALL find_u(sr_is(1,1,nsym_is),d_spin_is(1,1,nsym_is))
              s_is(:,:,nsym_is)=s(:,:,isym)
              gk_is(:,nsym_is)=gk(:,isym)
              ft_is(:,nsym_is)=ft(:,isym)
              sname_is(nsym_is)=sname(isym)
           ENDIF
        ELSE
           CALL find_u(sr(1,1,isym),d_spink(1,1,isym))
        ENDIF
     ENDIF
  ENDDO

  IF (noncolin.AND.domag) THEN
!
!   find the inverse of each element
!
     DO isym = 1, nsym_is
        found = .FALSE.
        DO jsym = 1, nsym_is
           !
           ss = MATMUL (s_is(:,:,jsym),s_is(:,:,isym))
           ! s(:,:,1) is the identity
           IF ( ALL ( s_is(:,:,1) == ss(:,:) ) ) THEN
              invs_is (isym) = jsym
              found = .TRUE.
           ENDIF
        END DO
        IF ( .NOT.found) CALL errore ('find_info_group', ' Not a group', 1)
     ENDDO
!
!   Recheck if we can compute the representations. The group is now smaller
!
     is_symmorphic=.TRUE.
     search_sym=.TRUE.

     DO isym=1,nsym_is
        is_symmorphic=( is_symmorphic.AND.(ft_is(1,isym)==0.d0).AND.  &
                                          (ft_is(2,isym)==0.d0).AND.  &
                                          (ft_is(3,isym)==0.d0) )
     ENDDO
     IF (.NOT.is_symmorphic) THEN
        DO isym=1,nsym_is
           DO jsym=1,nsym_is
              search_sym=search_sym.AND.(ABS(gk_is(1,isym)*ft_is(1,jsym)+ &
                                             gk_is(2,isym)*ft_is(2,jsym)+     &
                                             gk_is(3,isym)*ft_is(3,jsym))<1.D-8) 
           ENDDO
        ENDDO
     ENDIF
  END IF
  !
  !  Set the group name, divide it in classes and set the irreducible
  !
  CALL find_group(nsym,sr,gname,code_group)
  IF (noncolin) THEN
     IF (domag) THEN
        CALL find_group(nsym_is,sr_is,gname_is,code_group_is)
        CALL set_irr_rap_so(code_group_is,nclass,nrap,char_mat_so, &
             name_rap_so,name_class_so,name_class_so1)
        CALL divide_class_so_tpw(code_group_is,nsym_is,sr_is,d_spin_is,&
             has_e,nclass,nelem_so,elem_so,which_irr_so)
     ELSE
        CALL set_irr_rap_so(code_group,nclass,nrap,char_mat_so, &
             name_rap_so,name_class_so,name_class_so1)
        CALL divide_class_so_tpw(code_group,nsym,sr,d_spink, &
             has_e,nclass,nelem_so,elem_so,which_irr_so)
     ENDIF
  ELSE
     CALL set_irr_rap(code_group,nclass,char_mat,name_rap,name_class,ir_ram)
     CALL divide_class_tpw(code_group,nsym,sr,nclass,nelem,elem,which_irr)
  ENDIF

  RETURN
END SUBROUTINE find_info_group_tpw
