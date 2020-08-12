!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE find_representations_mode_q ( nat, ntyp, xq, w2, u, tau, ityp, &
                  amass, num_rap_mode, nspin_mag, qcode_old )
!-----------------------------------------------------------------------
!
!  This routine find the symmetry of the phonon modes at a given q.
!
!
  USE kinds,        ONLY : DP
  USE cell_base,    ONLY : at, bg
  USE symm_base,    ONLY : s, irt, nsym, time_reversal, copy_sym, &
                           s_axis_to_cart, inverse_s
  USE lr_symm_base, ONLY : gi, nsymq, rtau
  USE control_ph,   ONLY : search_sym
  USE ph_symmetry,  ONLY : manage_ph_symmetry

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: nat, ntyp, nspin_mag, qcode_old
  REAL(DP), INTENT(IN) :: xq(3), amass(ntyp), tau(3,nat)
  REAL(DP), INTENT(IN) :: w2(3*nat)
  INTEGER, INTENT(IN)  :: ityp(nat)
  COMPLEX(DP), INTENT(IN) :: u(3*nat,3*nat)
  INTEGER, INTENT(OUT) :: num_rap_mode(3*nat)

  REAL(DP) :: gimq (3)
  INTEGER :: irotmq
  LOGICAL :: minus_q, sym(48), magnetic_sym
  LOGICAL :: symmorphic_or_nzb
!
!  This routine assumes that u contains the eigenvectors of the dynamical
!  matrix
!
!
!  find the small group of q and set the quantities needed to
!  symmetrize a phonon mode
!
  IF (.NOT.search_sym) RETURN
  time_reversal=(nspin_mag/=4)
  minus_q=.TRUE.
  IF (nspin_mag/=4) minus_q=.FALSE.

  sym(1:nsym)=.true.
  call smallg_q (xq, 0, at, bg, nsym, s, sym, minus_q)
  nsymq=copy_sym(nsym,sym)
  CALL s_axis_to_cart ()
  CALL set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
  CALL inverse_s()
  CALL sgam_lr (at, bg, nsymq, s, irt, tau, rtau, nat)
!
!  Search the symmetry of the phonon mode.
!
  CALL manage_ph_symmetry(u, w2, num_rap_mode, xq, .TRUE., -1)

  RETURN
END SUBROUTINE find_representations_mode_q
