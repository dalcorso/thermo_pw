!
! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE ph_symmetry
  !
  USE kinds, ONLY : DP
  USE ions_base, ONLY : nat
  !
  SAVE
  PRIVATE

  INTEGER :: s_ph(3,3,48), ftau_ph(3,48), nsym_ph, invs_ph(48), t_rev_ph(48)
  INTEGER, ALLOCATABLE :: irt_ph(:,:)

  REAL(DP) :: sr_ph(3,3,48), gi_ph(3,48)
  REAL(DP), ALLOCATABLE :: rtau_ph(:,:,:)
  CHARACTER(len=45) :: sname_ph(48)

  INTEGER :: qcode_old=0

  PUBLIC  manage_ph_symmetry, initialize_gcode_old

CONTAINS

SUBROUTINE manage_ph_symmetry(dyn, w2, num_rap_mode, xq, search_sym, flag)
!
!  This routine is a driver that takes all the modes at a given point q
!  and writes on output the irreducible representations of each mode.
!  It works also for nonsymmorphic points group. The point group and the
!  symmetry information should be calculated before calling this 
!  routine and must be in the standard variables, s, sr, ftau, etc.
!  In the noncollinear magnetic system the routine removes the symmetries
!  that require time reversal and use the smaller group to identify the
!  symmetry of the mode.
!
USE kinds,         ONLY : DP
USE cell_base,     ONLY : at, bg
USE ions_base,     ONLY : nat, tau, ntyp => nsp, ityp, amass
USE control_ph,    ONLY : lgamma_gamma
USE fft_base,      ONLY : dfftp
USE control_lr,    ONLY : lgamma
USE noncollin_module, ONLY : nspin_mag
USE rap_point_group,  ONLY : code_group, nclass, nelem, elem, elem_name
USE proj_rap_point_group, ONLY : lqproj, qptype, which_elem, group_desc, &
                                 code_groupq_ext
USE lattices,      ONLY : zone_border
USE point_group,   ONLY : find_group_info_ext
USE io_global,     ONLY : stdout

IMPLICIT NONE

COMPLEX(DP) :: dyn(3*nat, 3*nat)
REAL(DP)    :: w2(3*nat), xq(3)
INTEGER     :: num_rap_mode(3*nat), flag
LOGICAL     :: search_sym

REAL(DP)    :: ft(3,48), wrk(3,48)
INTEGER     :: isym, ptype(3), ierr, gii(3,48)
LOGICAL     :: magnetic_sym
LOGICAL     :: symmorphic_or_nzb
!
IF (search_sym) THEN
   CALL set_symmetry_ph()
   IF (symmorphic_or_nzb()) THEN
      lqproj=0
      IF (zone_border(xq,at,bg,-1)) lqproj=2
      qptype=1
      magnetic_sym=(nspin_mag==4)
      CALL prepare_sym_analysis(nsym_ph,sr_ph,t_rev_ph,magnetic_sym)
      CALL find_group_info_ext(nsym_ph,sr_ph,code_group,code_groupq_ext, &
                                                    which_elem, group_desc)
      IF (code_group/=qcode_old) THEN
         CALL set_class_el_name(nsym_ph,sname_ph,nclass,nelem,elem,elem_name)
         CALL write_group_info_ph(.TRUE.)
      ENDIF
      CALL find_mode_sym_tpw (dyn, w2, tau, nat, nsym_ph, s_ph, sr_ph, irt_ph, &
            xq, rtau_ph, amass, ntyp, ityp, flag, lgamma_gamma, .FALSE., &
            num_rap_mode, ierr)
      CALL print_mode_sym(w2, num_rap_mode, lgamma)
   ELSE
      WRITE(stdout,'(/,5x,"Zone border point and nonsymmorphic &
                                                  &operations. Using")')
      lqproj=1
      DO isym = 1, nsym_ph
         ft(1,isym) = DBLE(ftau_ph(1,isym)) / DBLE(dfftp%nr1)
         ft(2,isym) = DBLE(ftau_ph(2,isym)) / DBLE(dfftp%nr2)
         ft(3,isym) = DBLE(ftau_ph(3,isym)) / DBLE(dfftp%nr3)
      END DO
      wrk(:,1:nsym_ph)=gi_ph(:,1:nsym_ph)
      CALL cryst_to_cart (nsym_ph, wrk, at, -1)
      gii(:,1:nsym_ph)=NINT(wrk(:,1:nsym_ph))

      CALL prepare_sym_analysis_proj(nsym_ph,s_ph,sr_ph,ft,gii,t_rev_ph,&
                                                                    ptype,-1)
      CALL find_mode_sym_proj (dyn, w2, tau, nat, nsym_ph, s_ph, sr_ph, ft, &
                            gii, invs_ph, irt_ph, xq, rtau_ph, amass, ntyp, &
                            ityp, flag, .FALSE., .FALSE., num_rap_mode, ierr)

      CALL print_mode_sym_proj(w2, num_rap_mode, ptype)
   ENDIF
   qcode_old=code_group
   CALL unset_symmetry_ph()
ENDIF

RETURN
END SUBROUTINE manage_ph_symmetry

SUBROUTINE set_symmetry_ph()
!
!  This routine sets the symmetry elements that must be used to
!  analyze the phonon modes. They are the small point group of q in the
!  standard case, or the small point group of q calculated only from the
!  symmetries that do not contain time reversal in the noncollinear
!  magnetic case.
!
USE ions_base,        ONLY : nat
USE symm_base,        ONLY : s, sr, irt, ftau, nsym, invs, t_rev, sname
USE lr_symm_base,     ONLY : gi, nsymq, rtau
USE noncollin_module, ONLY : nspin_mag

IMPLICIT NONE

LOGICAL :: sym(48)
INTEGER :: isym, jsym, ss(3,3)
LOGICAL :: found

ALLOCATE (irt_ph(48, nat))
ALLOCATE (rtau_ph(3, 48, nat))
IF (nspin_mag==4) THEN
!
!  In the noncollinear magnetic case use only the symmetries that do
!  not require time reversal
!
   nsym_ph=0
   DO isym=1,nsymq
      IF (t_rev(isym)==0) THEN
         nsym_ph=nsym_ph+1
         s_ph(:,:,nsym_ph) = s(:,:,isym)
         sr_ph(:,:,nsym_ph) = sr(:,:,isym)
         irt_ph(nsym_ph,1:nat) = irt(isym,1:nat)
         ftau_ph(:,nsym_ph) = ftau(:,isym)
         t_rev_ph(nsym_ph) = 0
         gi_ph(:,nsym_ph) = gi(:,isym)
         sname_ph(nsym_ph) = sname(isym)
         rtau_ph(:,nsym_ph,1:nat)=rtau(:,isym,1:nat)
      ENDIF
   ENDDO

   CALL find_inverse_s( nsym_ph, s_ph, invs_ph)
ELSE
   nsym_ph=nsymq
   s_ph=s
   sr_ph=sr
   irt_ph=irt
   ftau_ph=ftau
   invs_ph=invs
   t_rev_ph=t_rev
   gi_ph=gi
   sname_ph=sname
   rtau_ph=rtau
ENDIF

RETURN
END SUBROUTINE set_symmetry_ph

SUBROUTINE unset_symmetry_ph

IMPLICIT NONE

DEALLOCATE(irt_ph)
DEALLOCATE(rtau_ph)
RETURN
END SUBROUTINE unset_symmetry_ph

SUBROUTINE initialize_gcode_old(code)
IMPLICIT NONE
INTEGER :: code

qcode_old=code

RETURN
END SUBROUTINE initialize_gcode_old

END MODULE ph_symmetry
