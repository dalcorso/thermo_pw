!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE kpoint_grid_tpw ( nrot, time_reversal, skip_equivalence, s, t_rev, &
              bg, npk, k1,k2,k3, nk1,nk2,nk3, nks, xk, wk, comm, me_proc, nproc)
!-----------------------------------------------------------------------
!
!  Automatic generation of a uniform grid of k-points
!  Added the parallelization with all the processors
!
  USE kinds, ONLY: DP
  USE mp,    ONLY : mp_sum
  IMPLICIT NONE
  !
  INTEGER, INTENT(in):: nrot, npk, k1, k2, k3, nk1, nk2, nk3, &
                        t_rev(48), s(3,3,48), comm, nproc, me_proc
  LOGICAL, INTENT(in):: time_reversal, skip_equivalence
  real(DP), INTENT(in):: bg(3,3)
  !
  INTEGER, INTENT(out) :: nks
  real(DP), INTENT(out):: xk(3,npk)
  real(DP), INTENT(out):: wk(npk)
  ! LOCAL:
  real(DP), PARAMETER :: eps=1.0d-5
  real(DP) :: xkr(3), fact, xx, yy, zz
  real(DP), ALLOCATABLE:: xkg(:,:), wkk(:), xk_(:,:), wk_(:)
  INTEGER :: nkr, i,j,k, ns, n, nk, startk, lastk, pos, iproc
  INTEGER, ALLOCATABLE :: equiv(:), npos(:), nks_save(:)
  LOGICAL :: in_the_list
  !
  nkr=nk1*nk2*nk3

  CALL divide(comm, nkr, startk, lastk)

  ALLOCATE (xkg(3,startk:lastk)) 
  ALLOCATE (wkk(nkr))
  ALLOCATE (equiv( nkr))
  ALLOCATE (npos( nkr))
  !
  xkg=0.0_DP
  DO n=startk, lastk
     i= (n - 1) / nk2 / nk3 + 1
     j= (n - 1 - (i-1)*nk2*nk3) / nk3 + 1
     k= (n - 1 - (i-1)*nk2*nk3 - (j-1)*nk3) + 1
     xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
     xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
     xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
  ENDDO
     
!  DO i=1,nk1
!     DO j=1,nk2
!        DO k=1,nk3
           !  this is nothing but consecutive ordering
!           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
!           !  xkg are the components of the complete grid in crystal axis
!           xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
!           xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
!           xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
!        ENDDO
!     ENDDO
!  ENDDO

  !  equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
  !  equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)



  equiv=0
  npos=0
  wkk=0.0_DP
  DO nk=startk,lastk
     equiv(nk)=nk
     npos(nk)=nk
  ENDDO

  IF ( skip_equivalence ) THEN
    CALL infomsg('kpoint_grid', 'ATTENTION: skip check of k-points equivalence')
    wkk = 1.d0
  ELSE
    DO nk=startk,lastk
    !  check if this k-point has already been found equivalent to another
      IF (equiv(nk) == nk) THEN
        wkk(nk)   = 1.0d0
        !  check if there are equivalent k-point to this in the list
        !  (excepted those previously found to be equivalent to another)
        !  check both k and -k
        DO ns=1,nrot
           DO i=1,3
              xkr(i) = s(i,1,ns) * xkg(1,nk) &
                     + s(i,2,ns) * xkg(2,nk) &
                     + s(i,3,ns) * xkg(3,nk)
              xkr(i) = xkr(i) - nint( xkr(i) )
           ENDDO
           IF(t_rev(ns)==1) xkr = -xkr
           xx = xkr(1)*nk1 - 0.5d0*k1
           yy = xkr(2)*nk2 - 0.5d0*k2
           zz = xkr(3)*nk3 - 0.5d0*k3
           in_the_list = abs(xx-nint(xx))<=eps .and. &
                         abs(yy-nint(yy))<=eps .and. &
                         abs(zz-nint(zz))<=eps
           IF (in_the_list) THEN
              i = mod ( nint ( xkr(1)*nk1 - 0.5d0*k1 + 2*nk1), nk1 ) + 1
              j = mod ( nint ( xkr(2)*nk2 - 0.5d0*k2 + 2*nk2), nk2 ) + 1
              k = mod ( nint ( xkr(3)*nk3 - 0.5d0*k3 + 2*nk3), nk3 ) + 1
              n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
              IF (equiv(n)==0) THEN
!
!  this n is on another processor, set that nk must be copied on n
!  if this position is smaller than the actual one
!
                 IF ( n < npos(nk)) npos(nk)=n
              ELSE
                 IF (n>nk .AND. equiv(n)==n ) THEN
                    equiv(n) = nk
                    wkk(nk)=wkk(nk)+1.0d0
                 ELSE
                    IF (equiv(n)/=nk .or. n<nk) &
                     CALL errore('kpoint_grid', &
                        'something wrong in the checking algorithm',1)
                 ENDIF
              ENDIF
           ENDIF
           IF ( time_reversal ) THEN
              xx =-xkr(1)*nk1 - 0.5d0*k1
              yy =-xkr(2)*nk2 - 0.5d0*k2
              zz =-xkr(3)*nk3 - 0.5d0*k3
              in_the_list=abs(xx-nint(xx))<=eps.and.abs(yy-nint(yy))<=eps &
                                                 .and. abs(zz-nint(zz))<=eps
              IF (in_the_list) THEN
                 i = mod ( nint (-xkr(1)*nk1 - 0.5d0 * k1 + 2*nk1), nk1 ) + 1
                 j = mod ( nint (-xkr(2)*nk2 - 0.5d0 * k2 + 2*nk2), nk2 ) + 1
                 k = mod ( nint (-xkr(3)*nk3 - 0.5d0 * k3 + 2*nk3), nk3 ) + 1
                 n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                 IF (equiv(n)==0) THEN
                    IF (n < npos(nk)) npos(nk)=n
                 ELSE
                    IF (n>nk .and. equiv(n)==n) THEN
                       equiv(n) = nk
                       wkk(nk)=wkk(nk)+1.0d0
                    ELSE
                       IF ((equiv(n)/=nk.or.n<nk) .AND.equiv(n)/=0) &
                         CALL errore('kpoint_grid', &
                            'something wrong in the checking algorithm',2)
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
      ENDIF
    ENDDO
    CALL mp_sum(equiv, comm)
    CALL mp_sum(npos, comm)
    CALL mp_sum(wkk, comm)
!
!   We made only a partial collection in each processor, now collect
!   everything.
!
    DO nk=1,nkr
    !  check if this k-point has already been found equivalent to another
      IF (equiv(nk) == nk ) THEN
         IF (npos(nk)/=nk) THEN
            equiv(nk)=npos(nk)
            wkk(npos(nk))=wkk(npos(nk))+wkk(nk)
         ENDIF
      ENDIF
    ENDDO
  ENDIF

  !  count irreducible points and order them
  nks=0
  DO nk=startk,lastk
     IF (equiv(nk)==nk) nks=nks+1
  END DO

  ALLOCATE(wk_(nks))
  ALLOCATE(xk_(3,nks))

  nks=0
  fact=0.0d0
  DO nk=startk,lastk
     IF (equiv(nk)==nk) THEN
        nks=nks+1
        IF (nks>npk) CALL errore('kpoint_grid','too many k-points',1)
        wk_(nks) = wkk(nk)
        fact    = fact+wk_(nks)
        !  bring back into to the first BZ
        DO i=1,3
           xk_(i,nks) = xkg(i,nk)-nint(xkg(i,nk))
        ENDDO
     ENDIF
  ENDDO
  !  go to cartesian axis (in units 2pi/a0)
  CALL cryst_to_cart(nks,xk_,bg,1)

  CALL mp_sum(fact, comm)
  !  normalize weights to one
  DO nk=1,nks
     wk_(nk) = wk_(nk)/fact
  ENDDO
!
!  Now collect the nks and find the total number of points
!
  ALLOCATE(nks_save(nproc))
  nks_save=0
  DO iproc=1, nproc
     IF ((me_proc+1) ==iproc) nks_save(iproc)=nks
  ENDDO
  CALL mp_sum(nks, comm)
  CALL mp_sum(nks_save, comm)
!
!  put the points in the correct position in the list
!
  pos=0
  xk=0.0_DP
  wk=0.0_DP
  DO iproc=1, nproc
     IF ((me_proc+1)== iproc) THEN
        xk(:,pos+1:pos+nks_save(iproc))=xk_(:,1:nks_save(iproc))
        wk(pos+1:pos+nks_save(iproc))=wk_(1:nks_save(iproc))
     ENDIF
     pos=pos+nks_save(iproc)
  ENDDO
  CALL mp_sum(xk, comm)
  CALL mp_sum(wk, comm)

  DEALLOCATE(nks_save)
  DEALLOCATE(npos)
  DEALLOCATE(equiv)
  DEALLOCATE(xkg,wkk)
  DEALLOCATE(xk_,wk_)

  RETURN
END SUBROUTINE kpoint_grid_tpw
