!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE divide_et_impera_tpw( nkstot, xk, wk, isk, nks, diago_bands, &
                                          isym_bands, ik_origin )
  !!
  !! This routine divides the k points across nodes, sets the variable
  !! nks equal to the local (on this processors) number of k-points
  !! (nkstot on input is the total number of k-points)
  !! The distributed has "granularity kunit", that is, kunit consecutive 
  !! points stay on the same processor. Usually kunit=1; kunit=2 is used 
  !! in phonon calculations, when one has interspersed k_i and k_i+q and
  !! it is needed that they stay on the same processor, kunit=4 when
  !! one wants in the same processor k_i and k_i+q of both spins
  !!
  USE kinds,     ONLY : DP
  USE mp_pools,  ONLY : my_pool_id, npool, kunit
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: nkstot
  !! total number of k-points
  INTEGER, INTENT(INOUT) :: isk(nkstot)
  !! spin index of each kpoint (used in LSDA calculations only)
  REAL (DP), INTENT(INOUT) :: xk(3,nkstot)
  !! k-points (on all processors)
  REAL (DP), INTENT(INOUT) :: wk(nkstot)
  !! k-point weights
  INTEGER, INTENT(OUT)  :: nks
  !! number of k-points per pool
  !
  INTEGER :: ik, nbase, rest

  INTEGER :: isym_bands(nkstot), ik_origin(nkstot)
  LOGICAL :: diago_bands(nkstot)
  !
  ! simple case: no pools
  !
  IF ( npool == 1 ) THEN
     nks = nkstot
     RETURN
  END IF
  !
  IF ( MOD( nkstot, kunit ) /= 0 ) &
     CALL errore( 'divide_et_impera', 'nkstot/kunit is not an integer', nkstot )
  !
  nks    = kunit * ( nkstot / kunit / npool )
  !
  IF (nks == 0) CALL errore('divide_et_impera','some nodes have no k-points', 1)
  !
  rest = ( nkstot - nks * npool ) / kunit
  !
  IF ( my_pool_id < rest ) nks = nks + kunit
  !
  ! ... calculates nbase = the position in the list of the first point
  ! ...                    that belongs to this pool, minus one
  !
  nbase = nks * my_pool_id
  IF ( my_pool_id >= rest ) nbase = nbase + rest * kunit
  !
  ! ... displace the nks points in the pool to the first positions of the list
  !
  IF ( nbase > 0 ) THEN
     !
     xk(:,1:nks) = xk(:,nbase+1:nbase+nks)
     wk (1:nks)  = wk(nbase+1:nbase+nks)
     isk(1:nks)  =isk(nbase+1:nbase+nks)
  ENDIF

  DO ik=1,nks
     IF ( ik_origin(nbase+ik)>nbase .AND. &
                          ik_origin(nbase+ik)<=nbase+nks ) THEN
        diago_bands(ik) = diago_bands(nbase+ik)
        isym_bands(ik) = isym_bands(nbase+ik)
        ik_origin(ik)=ik_origin(nbase+ik) - nbase
     ELSE
!
!   If the original k point is not in this pool, diagonalizes
!
        diago_bands(ik)=.TRUE.
     ENDIF
  ENDDO
  !
  RETURN
  !
END SUBROUTINE divide_et_impera_tpw
