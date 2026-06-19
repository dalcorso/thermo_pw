!
! Copyright (C) 2026 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE small_g_g (gvec, at, nrot, s, t_rev, sym )
  !-----------------------------------------------------------------------
  !
  ! This routine selects, among the symmetry matrices of the point group
  ! of a crystal, the symmetry operations which leave G unchanged.
  ! SG=G+G' is not allowed, only SG=G operations are kept.
  !
  USE KINDS, ONLY : DP
  
  IMPLICIT NONE

  REAL(DP), PARAMETER :: accep = 1.e-6_dp

  REAL(DP), INTENT(in) :: at (3, 3), gvec (3)
  ! input: the reciprocal lattice vectors
  ! input: the direct lattice vectors
  ! input: the G point 

  INTEGER, INTENT(IN) :: s (3, 3, 48), nrot, t_rev(48)
  ! input: the symmetry matrices
  ! input: number of symmetry operations

  LOGICAL, INTENT(INOUT) :: sym(48)
  ! input-output: .true. if symm. op. S G = G 
  !
  !  local variables
  !
  real(DP) :: ag (3), rag (3)
  ! q vector in crystal basis
  ! the rotated of the q vector
  ! the zero vector

  integer :: irot, ipol, jpol
  ! counter on symmetry op.
  ! counter on polarizations
  ! counter on polarizations

  !
  !   Transform G to the crystal basis
  !
  ag = gvec
  CALL cryst_to_cart (1, ag, at, - 1)
  !
  !   Test all symmetries of the crystal to see if this operation send SG in G 
  !
  DO irot = 1, nrot
     if (.NOT.sym (irot) ) CYCLE
     rag(:) = 0.d0
     DO ipol = 1, 3
        DO jpol = 1, 3
           rag(ipol) = rag(ipol) + DBLE( s(ipol,jpol,irot) ) * ag( jpol)
        ENDDO
     ENDDO
     IF (t_rev(irot)==1) rag=-rag
     sym (irot) = (SQRT (  (rag(1)-ag(1))**2 + (rag(2)-ag(2))**2 + &
                           (rag(3)-ag(3))**2 ) ) < accep
  ENDDO
  !
  RETURN
  !
END SUBROUTINE small_g_g

