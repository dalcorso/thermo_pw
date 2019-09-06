! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE gen_qpoints_tpw (ibrav, at_, bg_, nat, tau, ityp, nk1, nk2, nk3, &
     ntetra, nqx, nq, q, wq)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE symm_base,  ONLY : set_sym_bl, find_sym, s, irt, nsym, &
                         nrot, t_rev, time_reversal,  sname, &
                         allfrac, remove_sym
  USE initial_conf, ONLY : nr1_save, nr2_save, nr3_save
  USE mp_world,     ONLY : world_comm, mpime, nproc
  !
  IMPLICIT NONE
  ! input
  INTEGER :: ibrav, nat, nk1, nk2, nk3, ntetra, iq, ityp(*)
  REAL(DP) :: at_(3,3), bg_(3,3), tau(3,nat)
  ! output
  INTEGER :: nqx, nq, tetra(4,ntetra)
  REAL(DP) :: q(3,nqx), wq(nqx)
  ! local
  REAL(DP) :: xqq(3), mdum(3,nat)
  LOGICAL :: magnetic_sym=.FALSE., skip_equivalence=.FALSE.
  !
  time_reversal = .true.
  t_rev(:) = 0
  xqq (:) =0.d0
  at = at_
  bg = bg_
  CALL set_sym_bl ( )
  !
  CALL kpoint_grid_tpw ( nrot, time_reversal, skip_equivalence, s, &
        t_rev, bg, nqx, 0,0,0, nk1,nk2,nk3, nq, q, wq, world_comm, mpime, nproc)
  !
  CALL find_sym ( nat, tau, ityp, .NOT.time_reversal, mdum )
  IF ( .NOT. allfrac ) CALL remove_sym ( nr1_save, nr2_save, nr3_save )
  !
  CALL irreducible_BZ (nrot, s, nsym, time_reversal, magnetic_sym, &
                       at, bg, nqx, nq, q, wq, t_rev)


  !
!  IF (ntetra /= 6 * nk1 * nk2 * nk3) &
!       CALL errore ('gen_qpoints','inconsistent ntetra',1)
  !
!  write(stdout,*) 'tetrahedra'
!  CALL tetrahedra (nsym, s, time_reversal, t_rev, at, bg, nqx, 0, 0, 0, &
!       nk1, nk2, nk3, nq, q, ntetra, tetra)
  !
  RETURN
END SUBROUTINE gen_qpoints_tpw
