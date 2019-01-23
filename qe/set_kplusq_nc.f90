!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine set_kplusq_nc (xk, wk, xq, nks, npk)
  !-----------------------------------------------------------------------
  !     This routine sets the k and k+q points (with zero weight) used in
  !     the preparatory run for a linear response calculation.
  !
  !     on input: xk and wk contain k-points and corresponding weights
  !
  !     on output: the number of points is doubled and xk and wk in the
  !                odd  positions are the original ones; those in the
  !                even positions are the corresponding k+q values.
  !     the gamma point is treated in a special way. No change is done
  !     to the k-points
  !
  USE kinds, only : DP
  implicit none
  !
  !    First the dummy variables
  !

  integer :: npk, nks
  ! input-output: maximum allowed number of k
  ! input-output: starting and ending number of
  real(DP) :: xk (3, npk), wk (npk), eps, xq (3)
  ! input-output: coordinates of k points
  ! input-output: weights of k points
  ! the smallest xq
  ! input: coordinates of a q-point
  !
  !    And then the local variables
  !

  logical :: lgamma
  ! true if xq is the gamma point
  integer :: ik, j
  ! counter on k
  ! counter
  !
  eps = 1.d-12
  !
  ! shift the k points in the odd positions and fill the even ones with k+
  !

  lgamma = abs (xq (1) ) .lt.eps.and.abs (xq (2) ) .lt.eps.and.abs ( &
       xq (3) ) .lt.eps

  if (.not.lgamma) then

     if (4 * nks.gt.npk) call errore ('set_kplusq', 'too many k points', &
          & nks)
     do ik = nks, 1, - 1
        xk (:, 4 * ik - 3) = xk (:, ik)
        xk (:, 4 * ik - 2) = xk (:, ik) + xq (:)
        xk (:, 4 * ik - 1) = -xk (:, ik)
        xk (:, 4 * ik    ) = -xk (:, ik) - xq (:)
        wk (4 * ik - 3) = wk (ik) 
        wk (4 * ik - 2) = 0.0_DP
        wk (4 * ik - 1) = 0.0_DP
        wk (4 * ik    ) = 0.d0
     enddo
     nks = 4 * nks
  else
     if (2 * nks.gt.npk) call errore ('set_kplusq', 'too many k points', &
          & nks)
     do ik = nks, 1, - 1
        xk (:, 2 * ik - 1) = xk (:, ik)
        xk (:, 2 * ik    ) = -xk (:, ik) 
        wk (2 * ik - 1) = wk (ik) 
        wk (2 * ik    ) = 0.0_DP
     enddo
     nks = 2 * nks
  ENDIF
  return
end subroutine set_kplusq_nc

!-----------------------------------------------------------------------
subroutine set_kplusq_nc_tpw (xk, wk, xq, nks, npk, diago_bands, &
                                                    isym_bands, ik_origin)
  !-----------------------------------------------------------------------
  !     This routine sets the k and k+q points (with zero weight) used in
  !     the preparatory run for a linear response calculation.
  !
  !     on input: xk and wk contain k-points and corresponding weights
  !
  !     on output: the number of points is doubled and xk and wk in the
  !                odd  positions are the original ones; those in the
  !                even positions are the corresponding k+q values.
  !     the gamma point is treated in a special way. No change is done
  !     to the k-points
  !
  USE kinds, only : DP
  USE cell_base, ONLY : at
  USE symm_base, ONLY : nsym, s, sr, t_rev
  USE start_k,   ONLY : nk1, nk2, nk3
  USE band_computation, ONLY : nks0

  implicit none
  !
  !    First the dummy variables
  !
  integer :: npk, nks
  INTEGER :: isym_bands(npk), ik_origin(npk)
  LOGICAL :: diago_bands(npk)
  ! input-output: maximum allowed number of k
  ! input-output: starting and ending number of
  real(DP) :: xk (3, npk), wk (npk), eps, xq (3)
  ! input-output: coordinates of k points
  ! input-output: weights of k points
  ! the smallest xq
  ! input: coordinates of a q-point
  !
  !    And then the local variables
  !
  REAL(DP) :: xkq(3), xks(3)
  INTEGER :: invs(3,3,48), table(48,48)
  REAL(DP), ALLOCATABLE :: xk_save(:,:)
  LOGICAL :: diago_bands_save(nks)

  logical :: lgamma, samev
  ! true if xq is the gamma point
  integer :: ik, j, ika, jk, isym, jsym, kpol, iflag
  ! counter on k
  ! counter
  !
  eps = 1.d-12
  !
  ! shift the k points in the odd positions and fill the even ones with k+
  !
  lgamma = abs (xq (1) ) .lt.eps.and.abs (xq (2) ) .lt.eps.and.abs ( &
       xq (3) ) .lt.eps

  ALLOCATE(xk_save(3,nks))
  xk_save(:,1:nks)=xk(:,1:nks)
  diago_bands_save(1:nks)=diago_bands(1:nks)
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

  if (.not.lgamma) then

     if (4 * nks.gt.npk) call errore ('set_kplusq', 'too many k points', &
          & nks)
     do ik = nks, 1, - 1
        xk (:, 4 * ik - 3) = xk (:, ik)
        xk (:, 4 * ik - 2) = xk (:, ik) + xq (:)
        xk (:, 4 * ik - 1) = -xk (:, ik)
        xk (:, 4 * ik    ) = -xk (:, ik) - xq (:)
        wk (4 * ik - 3) = wk (ik) 
        wk (4 * ik - 2) = 0.0_DP
        wk (4 * ik - 1) = 0.0_DP
        wk (4 * ik    ) = 0.d0
        diago_bands( 4 * ik - 3 ) = diago_bands( ik )
        isym_bands( 4 * ik - 3 ) = isym_bands( ik )
        ik_origin( 4 * ik - 3) = 4 * ik_origin( ik )- 3

        DO ika=2, 0, -1
           xkq(:) = xk(:, 4*ik-ika)
           CALL cryst_to_cart(1, xkq, at, -1)
           IF (ika==2) THEN
!
!   k+q does not need to be time reversed, so we search operations that
!   have not time reversal, so this operator is not applied.
!
              iflag=0
           ELSEIF (ika==1.OR.ika==0) THEN
!
!   -k and -k-q need to be time reversed, so we search first operations that
!    have time reversa, so the operator T^2 is not applied.
!
              iflag=1
           ENDIF 
           CALL check_k_in_list(xkq, xk_save, diago_bands_save, nks, iflag, &
                                samev, isym, jk, invs)
           IF (.NOT.samev) THEN
              iflag=MOD(iflag+1,2)
              CALL check_k_in_list(xkq, xk_save, diago_bands_save, &
                                          nks, iflag, samev, isym, jk, invs)
           ENDIF
           IF (samev) THEN
              diago_bands( 4 * ik - ika) = .FALSE.
              isym_bands( 4 * ik - ika ) = isym
              ik_origin( 4 * ik - ika ) = 4 * jk - 3
           ELSE
              diago_bands( 4 * ik -ika  ) = .TRUE.
           ENDIF
        ENDDO
     enddo
     nks = 4 * nks

  else
     if (2 * nks.gt.npk) call errore ('set_kplusq', 'too many k points', &
          & nks)
     do ik = nks, 1, - 1
        xk (:, 2 * ik - 1) = xk (:, ik)
        xk (:, 2 * ik    ) = -xk (:, ik) 
        wk (2 * ik - 1) = wk (ik) 
        wk (2 * ik    ) = 0.0_DP
        diago_bands( 2 * ik - 1 ) = diago_bands( ik )
        isym_bands( 2 * ik - 1 ) = isym_bands( ik )
        ik_origin( 2 * ik - 1) = 2 * ik_origin( ik )- 1

        xkq(:) = xk(:, 2*ik)
        CALL cryst_to_cart(1, xkq, at, -1)
!
!   We search first operation with time reversal, since the wavefunction
!   at -k must be time reversed and the application of T^2 avoid to
!   apply it.
!
        iflag=1
        CALL check_k_in_list(xkq, xk_save, diago_bands_save, nks, iflag, & 
                                                   samev, isym, jk, invs)
        IF (.NOT.samev) THEN
           iflag=0
           CALL check_k_in_list(xkq, xk_save, diago_bands_save, nks, &
                                                iflag, samev, isym, jk, invs)
        ENDIF
        IF (samev) THEN
           diago_bands( 2 * ik ) = .FALSE.
           isym_bands( 2 * ik  ) = isym
           ik_origin( 2 * ik  ) = 2 * jk - 1
        ELSE
           diago_bands( 2 * ik  ) = .TRUE.
        ENDIF
     enddo
     nks = 2 * nks
  ENDIF

  DEALLOCATE(xk_save)

  WRITE(6,*) 'Exiting set_kplus_q', nks0, nks
  DO ik=1,nks
     WRITE(6,*) diago_bands(ik), ik_origin(ik), isym_bands(ik)
  ENDDO

  return
end subroutine set_kplusq_nc_tpw

SUBROUTINE check_k_in_list(xkq, xk_save, diago_bands, nkstot, iflag, &
                               samev, isym_out, indk, invs)
!
!  This routine receives a k point xkq, and a list of nks0 kpoints
!  all of them in crystal coordinates. It checks if xkq can be obtained
!  from one of the xk_save k points by rotation with a symmetry of
!  the point group of the crystal. On output samev is .TRUE. if this
!  point exists. isym is the number of the rotation to use and indk
!  the number of the k point to rotate.
!  iflag can be 0 only operation with t_rev=0
!               1 only operation with t_rev=1
!               2 all operations
!  are searched.
!
USE kinds, ONLY : DP
USE symm_base, ONLY : nsym, s, sr, t_rev

IMPLICIT NONE
INTEGER :: nkstot, iflag, invs(3,3,48)
LOGICAL :: diago_bands(nkstot)
REAL(DP) :: xkq(3), xk_save(3, nkstot) 
LOGICAL :: samev
INTEGER :: isym_out, indk

REAL(DP) :: xks(3)
INTEGER :: jk, isym, kpol

samev=.FALSE.
DO jk=1,nkstot
   IF (.NOT.diago_bands(jk)) CYCLE
   DO isym= 1, nsym
      IF (t_rev(isym)/=iflag.AND.(iflag/=2)) CYCLE
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
         isym_out=isym
         indk=jk
         RETURN
      ENDIF
   ENDDO
ENDDO

RETURN
END SUBROUTINE check_k_in_list
!
