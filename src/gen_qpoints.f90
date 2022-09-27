! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE gen_qpoints_tpw (ibrav, at_, bg_, nat, tau, ityp, nk1, nk2, nk3, &
     nqx, nq, q, wq)
  !-----------------------------------------------------------------------
  !
  !  This routine generates a q point mesh for doing integrals over
  !  the Brillouin zone or to compute the phonon density of states.
  !  The q points are symmetrized and only inequivalent q points
  !  with appropriate weights are given in output.
  !
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE symm_base,  ONLY : set_sym_bl, find_sym, s, irt, nsym, &
                         nrot, t_rev, time_reversal, allfrac, remove_sym
  USE ifc,        ONLY : m_loc
  USE noncollin_module, ONLY : nspin_mag
  USE initial_conf, ONLY : nr1_save, nr2_save, nr3_save
  USE mp_world,     ONLY : world_comm, mpime, nproc
  !
  IMPLICIT NONE
  ! input
  INTEGER :: ibrav, nat, nk1, nk2, nk3, iq, ityp(*)
  REAL(DP) :: at_(3,3), bg_(3,3), tau(3,nat)
  ! output
  INTEGER :: nqx, nq
  REAL(DP) :: q(3,nqx), wq(nqx)
  ! local
  REAL(DP) :: xqq(3)
  LOGICAL :: skip_equivalence=.FALSE., magnetic_sym, minus_q
  !
  time_reversal = .TRUE.
  magnetic_sym=(nspin_mag==4)
  minus_q=.TRUE.
  t_rev(:) = 0
  xqq (:) =0.d0
  at = at_
  bg = bg_
  CALL set_sym_bl ( )
  !
  !   Here the q points are reduced with the symmetry of the Bravais 
  !   lattice that has time_reversal=.TRUE. also in magnetic systems.
  !
  CALL kpoint_grid_tpw ( nrot, time_reversal, skip_equivalence, s, &
        t_rev, bg, nqx, 0,0,0, nk1,nk2,nk3, nq, q, wq, &
                                       world_comm, mpime, nproc)
  !
  !   The q points are then reopened with the point group of the
  !   crystal that accounts for magnetism if present in the
  !   dynamical matrices files
  !
  CALL find_sym ( nat, tau, ityp, magnetic_sym, m_loc )
  IF ( .NOT. allfrac ) CALL remove_sym ( nr1_save, nr2_save, nr3_save )
  !
  !   
  !  The points at q and -q have the same phonon frequencies also
  !  in the magnetic case and minus_q is always .TRUE.
  !  magnetic_sym has to be considered .FALSE. because .TRUE. should
  !  be used only for the electronic bands.
  !
  CALL irreducible_BZ (nrot, s, nsym, minus_q, .FALSE., &
                       at, bg, nqx, nq, q, wq, t_rev)
  !
  RETURN
END SUBROUTINE gen_qpoints_tpw
!
!-----------------------------------------------------------------------
SUBROUTINE gen_qpoints_serial_tpw (ibrav, at_, bg_, nat, tau, ityp, &
     nk1, nk2, nk3, nqx, nq, q, wq)
  !-----------------------------------------------------------------------
  !
  !  This routine generates a q point mesh for doing integrals over
  !  the Brillouin zone or to compute the phonon density of states.
  !  The q points are symmetrized and only inequivalent q points
  !  with appropriate weights are given in output.
  !
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE symm_base,  ONLY : set_sym_bl, find_sym, s, irt, nsym, &
                         nrot, t_rev, time_reversal, allfrac, remove_sym
  USE ifc,        ONLY : m_loc
  USE noncollin_module, ONLY : nspin_mag
  USE initial_conf, ONLY : nr1_save, nr2_save, nr3_save
  !
  IMPLICIT NONE
  ! input
  INTEGER :: ibrav, nat, nk1, nk2, nk3, iq, ityp(*)
  REAL(DP) :: at_(3,3), bg_(3,3), tau(3,nat)
  ! output
  INTEGER :: nqx, nq, n
  REAL(DP) :: q(3,nqx), wq(nqx)
  ! local
  REAL(DP) :: xqq(3)
  LOGICAL :: skip_equivalence=.FALSE., magnetic_sym, minus_q
  !
  time_reversal = .TRUE.
  magnetic_sym=(nspin_mag==4)
  minus_q=.TRUE.
  t_rev(:) = 0
  xqq (:) =0.d0
  at = at_
  bg = bg_
  CALL set_sym_bl ( )
  !
  !   Here the q points are reduced with the symmetry of the Bravais 
  !   lattice that has time_reversal=.TRUE. also in magnetic systems.
  !
  CALL kpoint_grid_serial_tpw( nrot, time_reversal, skip_equivalence, s, &
        t_rev, bg, nqx, 0, 0, 0, nk1,nk2,nk3, nq, q, wq )

!  CALL kpoint_grid( nrot, time_reversal, skip_equivalence, s, &
!        t_rev, bg, nqx, 0, 0, 0, nk1,nk2,nk3, nq, q, wq )

  ! WRITE(6,*) nq
  !
  !   The q points are then reopened with the point group of the
  !   crystal that accounts for magnetism if present in the
  !   dynamical matrices files
  !
  CALL find_sym ( nat, tau, ityp, magnetic_sym, m_loc )
  IF ( .NOT. allfrac ) CALL remove_sym ( nr1_save, nr2_save, nr3_save )
  !
  !   
  !  The points at q and -q have the same phonon frequencies also
  !  in the magnetic case and minus_q is always .TRUE.
  !  magnetic_sym has to be considered .FALSE. because .TRUE. should
  !  be used only for the electronic bands.
  !
  CALL irreducible_BZ (nrot, s, nsym, minus_q, .FALSE., &
                       at, bg, nqx, nq, q, wq, t_rev)
  !
  RETURN
END SUBROUTINE gen_qpoints_serial_tpw
