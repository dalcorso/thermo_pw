!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE set_kplusq_tpw (xk, wk, xq, nks, npk, diago_bands, &
                                                    isym_bands, ik_origin)
  !-----------------------------------------------------------------------
  !! This routine sets the k and k+q points (with zero weight) used in
  !! the preparatory run for a linear response calculation.
  !
  !! * on input: xk and wk contain k-points and corresponding weights;
  !! * on output: the number of points is doubled and xk and wk in the
  !!              odd  positions are the original ones; those in the
  !!              even positions are the corresponding k+q values.
  !
  !! The gamma point is treated in a special way. No change is done
  !! to the k-points.
  !     This routine generalizes the one of Quantum espresso in two ways
  !     First because for each k+q it searches if among the k points there
  !     is one that can be used with a symmetry rotation of the point group
  !     to obtain k+q: k+q=S^-1 k + G, where possibly G is a reciprocal lattice
  !     vector. In that case it sets diago_bands to .FALSE. so that the
  !     wavefunctions at k+q can be obtained by rotation.
  !     Second when pools are used it rearranges the k and k+q couples so
  !     that symmetry related k and k+q points remain in the same pool
  !     as much as possible. The routine divide_et_impera then restore
  !     the diagonalization for those k+q for which the symmetry related
  !     k point is in another pool. 
  !
  USE kinds,   ONLY: DP
  USE cell_base, ONLY : at
  USE symm_base, ONLY : nsym, s, sr, t_rev
  USE start_k,   ONLY : nk1, nk2, nk3
  USE band_computation, ONLY : nks0
  USE mp_pools, ONLY : npool
  !
  IMPLICIT NONE
  !
  INTEGER :: npk
  !! inout: maximum allowed number of k
  INTEGER :: nks
  !! inout: starting and ending number of
  REAL(DP) :: xk(3,npk)
  !! inout: coordinates of k points
  REAL(DP) :: wk(npk)
  !! inout: weights of k points
  REAL(DP) :: xq(3)
  !! input: coordinates of a q-point

  INTEGER :: isym_bands(npk), ik_origin(npk)
  LOGICAL :: diago_bands(npk)
  ! 
  ! ... local variables
  !
  LOGICAL :: lgamma, samev
  ! true if xq is the gamma point
  INTEGER :: ik, jk, kpol, isym, jsym, j, nk
  ! counter on k
  ! counter
  REAL(DP), ALLOCATABLE :: xk_save(:,:), wk_new(:)
  REAL(DP) :: xkq(3), xks(3), eps
  INTEGER :: invs(3,3,48), table(48,48)
  LOGICAL :: diago_bands_save(nks)
  !
  eps = 1.d-12
  !
  ! ... shift the k points in the odd positions and fill the even ones with k+
  !
  diago_bands_save(1:nks)=diago_bands(1:nks)

  lgamma = ABS(xq(1))<eps .AND. ABS(xq(2))<eps .AND. ABS(xq(3))<eps

  IF (.NOT.lgamma) THEN

     IF (2 * nks > npk) CALL errore( 'set_kplusq', 'too many k points', nks )

     ALLOCATE(xk_save(3,nks))
     xk_save(:,1:nks)=xk(:,1:nks)
     CALL cryst_to_cart(nks, xk_save, at, -1)

     CALL multable (nsym, s, table)
     !
     !   And we set the matrices of the inverse
     !
     DO isym = 1, nsym
        DO jsym = 1, nsym
           IF (table (isym, jsym)==1) invs (:,:,isym) = s(:,:,jsym)
        ENDDO
     ENDDO

     DO ik = nks, 1, - 1
        DO j = 1, 3
           xk(j,2*ik-1) = xk(j,ik)
           xk(j,2*ik) = xk(j,ik) + xq(j)
        ENDDO
        wk(2*ik-1) = wk(ik)
        wk(2*ik) = 0.d0
        diago_bands( 2 * ik - 1 ) = diago_bands( ik )
        isym_bands( 2 * ik - 1 ) = isym_bands( ik )
        ik_origin( 2 * ik - 1) = 2 * ik_origin( ik )- 1
!
!   now check the k+q. Search if among the first nks0 k points
!   there is a k that rotated with a symmetry of the point group gives k+q+G
!
        xkq(:) = xk(:, 2*ik)
        CALL cryst_to_cart(1, xkq, at, -1)
        DO jk=1,nks
           IF (.NOT.diago_bands_save(jk)) CYCLE
           DO isym= 1, nsym
              DO kpol = 1, 3
                 xks (kpol) = invs(kpol, 1, isym)*xk_save (1,jk) + &
                              invs(kpol, 2, isym)*xk_save (2,jk) + &
                              invs(kpol, 3, isym)*xk_save (3,jk)
              ENDDO
              IF (t_rev(isym)==1)  xks (:)=-xks(:)
              samev = ABS( xks(1) - xkq(1) - &
                     NINT( xks(1) - xkq(1) ) ) < 1.0d-5 .AND. &
                      ABS( xks(2) - xkq(2) - &
                     NINT( xks(2) - xkq(2) ) ) < 1.0d-5 .AND. &
                      ABS( xks(3) - xkq(3) - &
                     NINT( xks(3) - xkq(3) ) ) < 1.0d-5
              IF (samev) THEN
                 diago_bands( 2 * ik ) = .FALSE.
                 isym_bands( 2 * ik ) = isym
                 ik_origin( 2 * ik ) = 2 * jk - 1
                 GOTO 100
              ENDIF
           ENDDO
        ENDDO            
        diago_bands( 2 * ik  ) = .TRUE.
        isym_bands( 2 * ik ) = 1
        ik_origin( 2 * ik ) = 2 * ik
100     CONTINUE
     ENDDO
     nks = 2 * nks

     DEALLOCATE(xk_save)
!
!  to skip the band symmetrization and diagonalize all as before
!  uncomment the following line
!
!     diago_bands=.TRUE.

  ENDIF
  RETURN
END SUBROUTINE set_kplusq_tpw
