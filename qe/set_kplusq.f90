!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine set_kplusq_tpw (xk, wk, xq, nks, npk, diago_bands, &
                                                    isym_bands, ik_origin)
  !-----------------------------------------------------------------------
  !     This routine sets the k and k+q points (with zero weight) used in
  !     the preparatory run for a linear response calculation.
  !
  !     on input: xk and wk contain k-points and corresponding weights
  !               xq the q point
  !
  !     on output: the number of points is doubled and xk and wk in the
  !                odd  positions are the original ones; those in the
  !                even positions are the corresponding k+q values.
  !     the gamma point is treated in a special way. No change is done
  !     to the k-points
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
  USE kinds, only : DP
  USE cell_base, ONLY : at
  USE symm_base, ONLY : nsym, s, sr, t_rev
  USE start_k,   ONLY : nk1, nk2, nk3
  USE band_computation, ONLY : nks0
  USE mp_pools, ONLY : npool
  implicit none
  !
  !    First the dummy variables
  !

  integer :: npk, nks
  ! input-output: maximum allowed number of k
  ! input-output: starting and ending number of
  INTEGER :: isym_bands(npk), ik_origin(npk)
  LOGICAL :: diago_bands(npk)
  real(DP) :: xk (3, npk), wk (npk), eps, xq (3)
  ! input-output: coordinates of k points
  ! input-output: weights of k points
  ! the smallest xq
  ! input: coordinates of a q-point
  !
  !    And then the local variables
  !
  logical :: lgamma, samev
  ! true if xq is the gamma point
  integer :: ik, jk, kpol, isym, jsym, j, nk, set_pos, ipool, resto
  ! counter on k
  ! counter
  REAL(DP), ALLOCATABLE :: xk_save(:,:), wk_new(:)
  INTEGER, ALLOCATABLE :: isym_bands_new(:), ik_origin_new(:), per(:), perm1(:)
  INTEGER, ALLOCATABLE :: start(:), nksp(:), nks0p(:), cur_pos(:)
  LOGICAL, ALLOCATABLE :: diago_bands_new(:), done(:)
  REAL(DP) :: xkq(3), xks(3)
  INTEGER :: invs(3,3,48), table(48,48)
  !
  eps = 1.d-12
  !
  ! shift the k points in the odd positions and fill the even ones with k+
  !

  lgamma = abs (xq (1) ) .lt.eps.and.abs (xq (2) ) .lt.eps.and.abs ( &
       xq (3) ) .lt.eps

  if (.not.lgamma) then

     if (2 * nks.gt.npk) call errore ('set_kplusq', 'too many k points', &
          & nks)

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

     do ik = nks, 1, - 1
        do j = 1, 3
           xk (j, 2 * ik - 1) = xk (j, ik)
           xk (j, 2 * ik) = xk (j, ik) + xq (j)
        enddo
        wk (2 * ik - 1) = wk (ik)
        wk (2 * ik) = 0.d0
        diago_bands( 2 * ik - 1 ) = diago_bands( ik )
        isym_bands( 2 * ik - 1 ) = isym_bands( ik )
        ik_origin( 2 * ik - 1) = 2 * ik_origin( ik )- 1
!
!   now check the k+q. Search if among the first nks0 k points
!   there is a k that rotated with a symmetry of the point group gives k+q+G
!
        xkq(:) = xk(:, 2*ik)
        CALL cryst_to_cart(1, xkq, at, -1)
        DO jk=1,nks0
           DO isym= 1, nsym
              DO kpol = 1, 3
                 xks (kpol) = invs(kpol, 1, isym)*xk_save (1,jk) + &
                              invs(kpol, 2, isym)*xk_save (2,jk) + &
                              invs(kpol, 3, isym)*xk_save (3,jk)
              ENDDO
              IF (t_rev(isym)==1)  xks (:)=-xks(:)
              samev = abs (xks (1) - xkq(1) - &
                     nint (xks (1) - xkq(1) ) ) < 1.0d-5 .AND. &
                      abs (xks (2) - xkq(2) - &
                     nint (xks (2) - xkq(2) ) ) < 1.0d-5 .AND. &
                      abs (xks (3) - xkq(3) - &
                     nint (xks (3) - xkq(3) ) ) < 1.0d-5
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
!   If there are pools, rearrange the k and k+q points so that the
!   largest number of symmetry related points goes in the same pool
!
    IF (npool > 1) THEN
!
!   First distribute the nks0 points. Count how many per pool, and how
!   many total k and k+q points per pools are computed
!       
       ALLOCATE(nks0p(npool))
       ALLOCATE(start(npool))
       ALLOCATE(nksp(npool))
       DO ipool=1, npool
          nks0p(ipool)= nks0 / npool
          resto = nks0 - nks0p(ipool) * npool
          IF (ipool <= resto) nks0p(ipool)=nks0p(ipool)+1
          nksp(ipool)= nks/2 / npool
          resto = nks/2 - nksp(ipool) * npool
          IF (ipool <= resto) nksp(ipool)=nksp(ipool)+1
          start(ipool) = nksp(ipool) * (ipool - 1)
          IF (ipool > resto) start(ipool) = start(ipool) + resto 
       ENDDO

       ALLOCATE(per(nks/2))
       ALLOCATE(perm1(nks/2))
       ALLOCATE(done(nks/2))
       ALLOCATE(cur_pos(npool))
!       DO ik=1,nks
!          WRITE(6,*) 'ik, ik_origin', ik, ik_origin(ik)
!       ENDDO
!
!   first distribute the k0
!
       done=.FALSE.
       cur_pos=0
       DO ipool=1, npool
          set_pos=0
          DO jk=1,nks0p(ipool)
             cur_pos(ipool)=cur_pos(ipool)+1
             IF (cur_pos(ipool) > nks0p(ipool)) EXIT
             IF (cur_pos(ipool) > set_pos) THEN
!
!   search the first free k0 and set it in the current pool
!
                DO ik=1, nks0
                   IF (.NOT.done(ik)) THEN
                      per(start(ipool)+cur_pos(ipool))=ik
                      perm1(ik)=start(ipool)+cur_pos(ipool)
                      done(ik)=.TRUE.
                      GOTO 200
                   ENDIF
                ENDDO
             ENDIF
200          CONTINUE                
!
!    now search new points until we reach nks0p(ipool)
!
             nk=0
             DO ik=1, nks0
                IF (cur_pos(ipool)+nk+1 > nks0p(ipool)) EXIT
                IF (ik_origin(2*per(start(ipool)+cur_pos(ipool)))==&
                                                                (2*ik-1)) THEN
                   IF (.NOT.done(ik)) THEN
                      nk=nk+1 
                      per(start(ipool)+cur_pos(ipool)+nk)=ik
                      perm1(ik)=start(ipool)+cur_pos(ipool)+nk
                      done(ik)=.TRUE.
                      set_pos=cur_pos(ipool)+nk
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!       DO ipool=1, npool
!          DO ik=1,nks0p(ipool)
!             WRITE(6,*) 'ipool, start+ik, per', ipool, start(ipool)+ik, &
!                                                    per(start(ipool)+ik)
!          ENDDO    
!       ENDDO
!
!   Then add all those that were not used but have the k0 in the pool
!   Start by seaching for couples satisfied if put in a given pool
!       
       DO ik=nks0+1, nks/2
couples: IF (.NOT.done(ik)) THEN
             DO ipool=1, npool
                DO jk=1, nks0p(ipool)
                   IF (cur_pos(ipool)>= nksp(ipool)) CYCLE
                   IF (ik_origin(2*ik-1)==2*per(start(ipool)+jk)-1 .AND. &
                              ik_origin(2*ik)==2*per(start(ipool)+jk)-1) THEN
                      cur_pos(ipool)=cur_pos(ipool)+1
                      per(start(ipool)+cur_pos(ipool))=ik
                      perm1(ik)=start(ipool)+cur_pos(ipool)
                      done(ik)=.TRUE. 
                      EXIT couples
                   ENDIF
                ENDDO
             ENDDO
          ENDIF couples
       ENDDO
!
!   And then sets were only one point needs to be diagonalized if in a
!   given pools.
!   We start from the last pools that are probably those with more points
!   to compute and therefore must be served first
!
       DO ik=nks0+1, nks/2
singles: IF (.NOT.done(ik)) THEN
             DO ipool=npool, 1, -1
                DO jk=1, nks0p(ipool)
                   IF (cur_pos(ipool)>= nksp(ipool)) CYCLE
                   IF (ik_origin(2*ik-1)==2*per(start(ipool)+jk)-1 .OR. &
                              ik_origin(2*ik)==2*per(start(ipool)+jk)-1) THEN
                      cur_pos(ipool)=cur_pos(ipool)+1
                      per(start(ipool)+cur_pos(ipool))=ik
                      perm1(ik)=start(ipool)+cur_pos(ipool)
                      done(ik)=.TRUE. 
                      EXIT singles
                   ENDIF
                ENDDO
             ENDDO
          ENDIF singles
       ENDDO
!
!   Finally add the lone pairs, those k and k+q sets whose k0 are in a pool 
!   which is already filled and must go in another pool where both
!   must be computed by diagonalization
!
       DO ik=nks0+1, nks/2
lone_pairs:IF (.NOT.done(ik)) THEN
             DO ipool=1, npool
                DO jk=1, nks0p(ipool)
                   IF (cur_pos(ipool)>= nksp(ipool)) CYCLE
                   cur_pos(ipool)=cur_pos(ipool)+1
                   per(start(ipool)+cur_pos(ipool))=ik
                   perm1(ik)=start(ipool)+cur_pos(ipool)
                   done(ik)=.TRUE. 
                   EXIT lone_pairs
                ENDDO
             ENDDO
          ENDIF lone_pairs
       ENDDO

       ALLOCATE(xk_save(3,nks))
       ALLOCATE(wk_new(nks))
       ALLOCATE(diago_bands_new(nks))
       ALLOCATE(isym_bands_new(nks))
       ALLOCATE(ik_origin_new(nks))

       xk_save(:,1:nks)=xk(:,1:nks)
       wk_new(1:nks)=wk(1:nks)
       diago_bands_new(1:nks)=diago_bands(1:nks)
       isym_bands_new(1:nks)=isym_bands(1:nks)
       ik_origin_new(1:nks)=ik_origin(1:nks)
       DO ik=1,nks/2
!          WRITE(6,*) ik, per(ik)
          xk(:,2*ik-1)=xk_save(:,2*per(ik)-1)
          xk(:,2*ik)=xk_save(:,2*per(ik))
          wk(2*ik-1)=wk_new(2*per(ik)-1)
          wk(2*ik)=wk_new(2*per(ik))
          diago_bands(2*ik-1)=diago_bands_new(2*per(ik)-1)
          diago_bands(2*ik)=diago_bands_new(2*per(ik))
          isym_bands(2*ik-1)=isym_bands_new(2*per(ik)-1)
          isym_bands(2*ik)=isym_bands_new(2*per(ik))
          ik_origin(2*ik-1)= 2*perm1((ik_origin_new(2*per(ik)-1)+1)/2)-1
          ik_origin(2*ik)=2*perm1((ik_origin_new(2*per(ik))+1)/2)-1
       ENDDO

       DEALLOCATE(wk_new)
       DEALLOCATE(xk_save)
       DEALLOCATE(diago_bands_new)
       DEALLOCATE(isym_bands_new)
       DEALLOCATE(ik_origin_new)
       DEALLOCATE(per)
       DEALLOCATE(perm1)
       DEALLOCATE(done)
       DEALLOCATE(nks0p)
       DEALLOCATE(nksp)
       DEALLOCATE(start)
       DEALLOCATE(cur_pos)

    ENDIF
!
!  to skip the band symmetrization and diagonalize all as before
!  uncomment the following line
!
!     diago_bands=.TRUE.

  endif
  return
end subroutine set_kplusq_tpw
