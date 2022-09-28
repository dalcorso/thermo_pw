!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! Copyright (C) 2022 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE kpoint_grid_serial_tpw( nrot, time_reversal, skip_equivalence, &
                  s, t_rev, bg, npk, k1, k2, k3, nk1, nk2, nk3, nks, xk, wk )
  !-----------------------------------------------------------------------
  !  Automatic generation of a uniform grid of k-points.
  !  This routine uses less memory than the original routine and can
  !  be used for very large meshes.
  !
  USE kinds, ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nrot
  !! number of bravais lattice symmetries
  INTEGER, INTENT(IN) :: npk
  !! max number of k-points
  INTEGER, INTENT(IN) :: k1
  !! the offset from the origin, direction 1
  INTEGER, INTENT(IN) :: k2
  !! the offset from the origin, direction 2
  INTEGER, INTENT(IN) :: k3
  !! the offset from the origin, direction 3
  INTEGER, INTENT(IN) :: nk1
  !! the special-point grid, direction 1
  INTEGER, INTENT(IN) :: nk2
  !! the special-point grid, direction 2
  INTEGER, INTENT(IN) :: nk3
  !! the special-point grid, direction 3
  INTEGER, INTENT(IN) :: t_rev(48)
  !! time reversal flag, for noncolinear magnetism
  INTEGER, INTENT(IN) :: s(3,3,48)
  !! symmetry matrices, in crystal axis
  LOGICAL, INTENT(IN) :: time_reversal
  !! if .TRUE. the system has time reversal symmetry
  LOGICAL, INTENT(IN) :: skip_equivalence
  !! if .TRUE. skip check of k-points equivalence
  REAL(DP), INTENT(IN) :: bg(3,3)
  !! bg(:,i) are the reciprocal lattice vectors, b_i,
  !! in tpiba=2pi/alat units: b_i(:) = bg(:,i)/tpiba
  INTEGER,  INTENT(OUT) :: nks
  !! number of k points
  REAL(DP), INTENT(OUT) :: xk(3,npk)
  !! coordinates of k points
  REAL(DP), INTENT(OUT) :: wk(npk)
  !! weight of k points
  !
  ! ... local variables
  !
  REAL(DP) :: xkr(3), fact, xx, yy, zz, xkg(3)
  INTEGER :: nkr, i, j, k, il, jl, kl, ns, n, nk
  LOGICAL, ALLOCATABLE :: equiv(:)
  LOGICAL :: in_the_list
  REAL(DP), PARAMETER :: eps=1.0d-5
  !
  nkr=nk1*nk2*nk3
  IF ( skip_equivalence ) THEN
     CALL infomsg('kpoint_grid', 'ATTENTION: skip check of k-points &
                                                       &equivalence')
     IF (nkr>npk) CALL errore('kpoint_grid_serial','npk is too small',1)
     wk = 1.d0/DBLE(nkr)
     DO i=1,nk1
        DO j=1,nk2
           DO k=1,nk3
              !  this is nothing but consecutive ordering
              n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
              !  xkg are the components of the complete grid in crystal axis
              xk(1,n) = DBLE(i-1)/nk1 + DBLE(k1)/2/nk1
              xk(2,n) = DBLE(j-1)/nk2 + DBLE(k2)/2/nk2
              xk(3,n) = DBLE(k-1)/nk3 + DBLE(k3)/2/nk3
              xk(:,nk) = xk(:,nk)-NINT(xk(:,nk))
           ENDDO
        ENDDO
     ENDDO
     nks=nkr
!
!    Bring the k points in cartesian coordinates, in units 2 pi / a0
!
     CALL cryst_to_cart(nks,xk,bg,1)
     RETURN
  ENDIF 
!
!   Here we generate the k points mesh and count the inequivalent
!   points and their weights
!
  ALLOCATE (equiv(nkr))
  !
  DO nk=1,nkr
     equiv(nk)=.FALSE.
  ENDDO
  nks=0
  DO il=1,nk1
     DO jl=1,nk2
        DO kl=1,nk3
           !  this is nothing but consecutive ordering
           nk = (kl-1) + (jl-1)*nk3 + (il-1)*nk2*nk3 + 1
  !  equiv(nk) =.FALSE. : k-point nk is not equivalent to any previous k-point
  !  equiv(nk) =.TRUE. : k-point nk is equivalent to one previous k-point
  !  check if this k-point has already been found equivalent to another
           IF (.NOT.equiv(nk)) THEN
           !  xkg are the components of this new k point in crystal axis
              xkg(1) = DBLE(il-1)/nk1 + DBLE(k1)/2/nk1
              xkg(2) = DBLE(jl-1)/nk2 + DBLE(k2)/2/nk2
              xkg(3) = DBLE(kl-1)/nk3 + DBLE(k3)/2/nk3
              nks=nks+1
              IF (nks>npk) CALL errore('kpoint_grid_serial_tpw', &
                                                 'too many k-points',1)
              wk(nks)   = 1.0d0
              xk(:,nks) = xkg(:)
        !  check if there are equivalent k-point to this in the list
        !  (excepted those previously found to be equivalent to another)
        !  check both k and -k
              DO ns=1,nrot
                 DO i=1,3
                    xkr(i) = s(i,1,ns) * xkg(1) &
                           + s(i,2,ns) * xkg(2) &
                           + s(i,3,ns) * xkg(3)
                    xkr(i) = xkr(i) - NINT( xkr(i) )
                 ENDDO
                 IF (t_rev(ns)==1) xkr = -xkr
                 xx = xkr(1)*nk1 - 0.5d0*k1
                 yy = xkr(2)*nk2 - 0.5d0*k2
                 zz = xkr(3)*nk3 - 0.5d0*k3
                 in_the_list = ABS(xx-NINT(xx))<=eps .AND. &
                               ABS(yy-NINT(yy))<=eps .AND. &
                               ABS(zz-NINT(zz))<=eps
                 IF (in_the_list) THEN
                    i = MOD ( NINT ( xkr(1)*nk1 - 0.5d0*k1 + 2*nk1), nk1 ) + 1
                    j = MOD ( NINT ( xkr(2)*nk2 - 0.5d0*k2 + 2*nk2), nk2 ) + 1
                    k = MOD ( NINT ( xkr(3)*nk3 - 0.5d0*k3 + 2*nk3), nk3 ) + 1
                    n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                    IF (n>nk .AND. .NOT.equiv(n)) THEN
                       equiv(n) = .TRUE.
                       wk(nks)=wk(nks)+1.0d0
                    ELSE
                       IF ( n<nk ) CALL errore('kpoint_grid_serial_tpw', &
                              'something wrong in the checking algorithm',1)
                    ENDIF
                 ENDIF
                 IF ( time_reversal ) THEN
                    xx =-xkr(1)*nk1 - 0.5d0*k1
                    yy =-xkr(2)*nk2 - 0.5d0*k2
                    zz =-xkr(3)*nk3 - 0.5d0*k3
                    in_the_list=ABS(xx-NINT(xx))<=eps .AND. &
                                ABS(yy-NINT(yy))<=eps .AND. &
                                ABS(zz-NINT(zz))<=eps
                    IF (in_the_list) THEN
                       i = MOD ( NINT (-xkr(1)*nk1 - 0.5d0 * k1    &
                                                   + 2*nk1), nk1 ) + 1
                       j = MOD ( NINT (-xkr(2)*nk2 - 0.5d0 * k2    &
                                                   + 2*nk2), nk2 ) + 1
                       k = MOD ( NINT (-xkr(3)*nk3 - 0.5d0 * k3    &
                                                   + 2*nk3), nk3 ) + 1
                       n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                       IF (n>nk .AND. .NOT.equiv(n)) THEN
                          equiv(n) = .TRUE.
                          wk(nks)=wk(nks)+1.0d0
                       ELSE
                          IF (n<nk) CALL errore('kpoint_grid_serial_tpw', &
                              'something wrong in the checking algorithm',2)
                       ENDIF
                    ENDIF
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  !
  !  check that we have considered all k points, and bring k points
  !  inside the first Brillouin zone.
  !
  fact=0.0d0
  DO nk=1,nks
     fact = fact+wk(nk)
     DO i=1,3
        xk(i,nk) = xk(i,nk)-NINT(xk(i,nk))
     ENDDO
  ENDDO
  IF (NINT(fact) /= nkr) CALL errore('kpoint_grid_serial_tpw','Wrong weight',1)
  !  go to cartesian axis (in units 2pi/a0)
  CALL cryst_to_cart(nks,xk,bg,1)
  !  normalize weights to one
  DO nk=1,nks
     wk(nk) = wk(nk)/fact
  ENDDO

  DEALLOCATE(equiv)

  RETURN
END SUBROUTINE kpoint_grid_serial_tpw
!
