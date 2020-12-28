!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE star_q_tpw (xq, at, bg, nsym, s, invs, t_rev, nq, sxq, isq, &
                                                         imq, verbosity )
  !-----------------------------------------------------------------------
  ! generate the star of q vectors that are equivalent to the input one
  ! NB: input s(:,:,1:nsym) must contain all crystal symmetries,
  ! i.e. not those of the small-qroup of q only
  !
  USE io_global,  ONLY : stdout
  USE kinds, ONLY : DP
  IMPLICIT NONE
  !
  REAL(DP), PARAMETER :: accep=1.e-5_dp

  INTEGER, INTENT(IN) :: nsym, s (3, 3, 48), invs(48), t_rev(48)
  ! nsym matrices of symmetry operations
  ! invs: list of inverse operation indices
  REAL(DP), INTENT(IN) :: xq (3), at (3, 3), bg (3, 3)
  ! xq: q vector
  ! at: direct lattice vectors
  ! bg: reciprocal lattice vectors
  !
  INTEGER, INTENT(OUT) :: nq, isq (48), imq
  ! nq  : degeneracy of the star of q
  ! isq : index of q in the star for a given sym
  ! imq : index of -q in the star (0 if not present)

  REAL(DP), INTENT(OUT) :: sxq (3, 48)
  ! list of vectors in the star of q
  LOGICAL, INTENT(IN) :: verbosity
  ! if true prints several messages.
  !
  INTEGER :: nsq (48), isym, ism1, iq, i
  ! number of symmetry ops. of bravais lattice
  ! counters on symmetry ops.
  ! index of inverse of isym
  ! counters
  REAL(DP) :: saq (3, 48), aq (3), raq (3), zero (3)
  ! auxiliary list of q (crystal coordinates)
  ! input q in crystal coordinates
  ! rotated q in crystal coordinates
  ! coordinates of fractionary translations
  ! a zero vector: used in eqvect

  LOGICAL, EXTERNAL :: eqvect
  ! function used to compare two vectors
  !
  zero(:) = 0.d0
  !
  ! go to  crystal coordinates
  !
  DO i = 1, 3
     aq(i) = xq(1) * at(1,i) + xq(2) * at(2,i) + xq(3) * at(3,i)
  ENDDO
  !
  ! create the list of rotated q
  !
  DO i = 1, 48
     nsq (i) = 0
     isq (i) = 0
  ENDDO
  nq = 0
  DO isym = 1, nsym
     ism1 = invs (isym)
     DO i = 1, 3
        raq (i) = s (i, 1, ism1) * aq (1) &
                + s (i, 2, ism1) * aq (2) &
                + s (i, 3, ism1) * aq (3)
     ENDDO
     IF (t_rev(isym)==1) raq=-raq
     DO i = 1, 3
        sxq (i, 48) = bg (i, 1) * raq (1) &
                    + bg (i, 2) * raq (2) &
                    + bg (i, 3) * raq (3)
     ENDDO
     DO iq = 1, nq
        IF (eqvect (raq, saq (1, iq), zero, accep) ) THEN
           isq (isym) = iq
           nsq (iq) = nsq (iq) + 1
        ENDIF
     ENDDO
     IF (isq (isym) == 0) THEN
        nq = nq + 1
        nsq (nq) = 1
        isq (isym) = nq
        saq(:,nq) = raq(:)
        DO i = 1, 3
           sxq (i, nq) = bg (i, 1) * saq (1, nq) &
                       + bg (i, 2) * saq (2, nq) &
                       + bg (i, 3) * saq (3, nq)
        ENDDO
     ENDIF
  ENDDO
  !
  ! set imq index if needed and check star degeneracy
  !
  raq (:) = - aq(:)
  imq = 0
  DO iq = 1, nq
     IF (eqvect (raq, saq (1, iq), zero, accep) ) imq = iq
     IF (nsq(iq)*nq /= nsym) CALL errore ('star_q', 'wrong degeneracy', iq)
  ENDDO
  !
  ! writes star of q
  !
  IF (verbosity) THEN
     WRITE( stdout, * )
     WRITE( stdout, '(5x,a,i4)') 'Number of q in the star = ', nq
     WRITE( stdout, '(5x,a)') 'List of q in the star:'
     WRITE( stdout, '(7x,i4,3f14.9)') (iq, (sxq(i,iq), i=1,3), iq=1,nq)
     IF (imq == 0) THEN
        WRITE( stdout, '(5x,a)') 'In addition there is the -q list: '
        WRITE( stdout, '(7x,i4,3f14.9)') (iq, (-sxq(i,iq), i=1,3), iq=1,nq)
     ENDIF
  ENDIF
  RETURN
END SUBROUTINE star_q_tpw
