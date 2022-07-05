!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE rotate_and_add_dyn_tpw (phi, phi2, nat, isym, s, invs, irt, &
     rtau, sxq)
  !-----------------------------------------------------------------------
  !! Rotates a dynamical matrix (\(\text{phi}\)) in crystal coordinates
  !! according to the specified symmetry operation and add the rotated
  !! matrix to \(\text{phi2}\). \(\text{phi}\) is left unmodified.
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : tpi
  USE symm_base, ONLY : t_rev
  IMPLICIT NONE
  ! input variables
  integer :: nat
  !! number of atoms in the unit cell
  integer :: isym
  !! index of the symm.op.
  integer :: s(3,3,48)
  !! the symmetry operations
  integer :: invs(48)
  !! index of the inverse operations
  integer :: irt(48,nat)
  !! index of the rotated atom
  complex(DP) :: phi(3,3,nat,nat)
  !! the input dyn.mat. in crystal coordinates
  complex(DP) :: phi2(3,3,nat,nat)
  !! the rotated dyn.mat. in crystal coordinates
  real(DP) :: rtau(3,48,nat)
  !! for eaxh atom and rotation gives the R vector involved
  real(DP) :: sxq(3)
  !! the rotated q involved in this sym.op.
  !
  ! ... local variables
  INTEGER :: na, nb, sna, snb, ism1, i, j, k, l
  ! counters on atoms
  ! indices of rotated atoms
  ! index of the inverse symm.op.
  ! generic counters
  REAL(DP) :: arg
  ! argument of the phase
  COMPLEX(DP) :: phase, work


  ism1 = invs (isym)
  DO na = 1, nat
     DO nb = 1, nat
        sna = irt (isym, na)
        snb = irt (isym, nb)
        arg = (sxq (1) * (rtau (1, isym, na) - rtau (1, isym, nb) ) &
             + sxq (2) * (rtau (2, isym, na) - rtau (2, isym, nb) ) &
             + sxq (3) * (rtau (3, isym, na) - rtau (3, isym, nb) ) ) * tpi
        phase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
        DO i = 1, 3
           DO j = 1, 3
              work = CMPLX(0.d0, 0.d0,kind=DP)
              DO k = 1, 3
                 DO l = 1, 3
                    IF (t_rev(isym)==1) THEN
                       work = work + s (i, k, ism1) * s (j, l, ism1) * &
                               CONJG(phi (k, l, na, nb)) * phase
                    ELSE
                       work = work + s (i, k, ism1) * s (j, l, ism1) * &
                               phi (k, l, na, nb) * phase
                    ENDIF
                 ENDDO
              ENDDO
              phi2 (i, j, sna, snb) = phi2 (i, j, sna, snb) + work
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  RETURN
END SUBROUTINE rotate_and_add_dyn_tpw
