!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine set_kup_and_kdw_tpw (xk, wk, isk, nkstot, npk, lgamma, &
                                        diago_bands, isym_bands, ik_origin)
  !-----------------------------------------------------------------------
  !     This routine sets the k vectors for the up and down spin wfc
  !
  !     on input: xk and wk contain k-points and corresponding weights
  !
  !     on output: the number of points is doubled and xk and wk in the
  !                first (nkstot/2) positions correspond to up spin
  !                those in the second (nkstot/2) ones correspond to down spin
  !
  USE kinds, ONLY : DP
  implicit none
  LOGICAL, INTENT(IN) :: lgamma
  !
  ! I/O variables first
  !
  integer :: npk, isk (npk), nkstot
  ! input: maximum allowed number of k-points
  ! output: spin associated to a given k-point
  ! input-output: starting and ending number of k-points 
  real(DP) :: xk (3, npk), wk (npk), xksave(3,npk), wksave(npk)
  ! input-output: coordinates of k points
  ! input-output: weights of k points
  INTEGER :: isym_bands(npk), ik_origin(npk)
  LOGICAL :: diago_bands(npk)
  INTEGER :: isym_bands_save(npk), ik_origin_save(npk)
  LOGICAL :: diago_bands_save(npk)
  !
  integer :: ik, iq, ikq
  !
  !
  if (2*nkstot > npk) call errore ('set_kup_and_kdw_tpw','too many k points',nkstot)
  xksave(:,1:nkstot) = xk(:,1:nkstot)
  wksave(1:nkstot) = wk(1:nkstot)
  diago_bands_save(1:nkstot)=diago_bands(1:nkstot)
  isym_bands_save(1:nkstot)=isym_bands(1:nkstot)
  ik_origin_save(1:nkstot)=ik_origin(1:nkstot)
!
!   The magnon calculation requires k and k+q of different spins on the same 
!   pool, so we put them close to each other and use kunit=4, or kunit=2 in
!   the gamma case
!
  IF (lgamma) THEN 
     do ik = 1, nkstot
        xk(:,2*ik-1) = xksave(:,ik)
        xk(:,2*ik) = xksave(:,ik)
        wk(2*ik-1) = wksave(ik)
        wk(2*ik) = wksave(ik)
        isk(2*ik-1) = 1
        isk(2*ik) = 2
        diago_bands(2*ik-1)=diago_bands_save(ik)
        diago_bands(2*ik)=diago_bands_save(ik)
        isym_bands(2*ik-1)=isym_bands_save(ik)
        isym_bands(2*ik)=isym_bands_save(ik)
        ik_origin(2*ik-1)=2*ik_origin_save(ik)-1
        ik_origin(2*ik)=2*ik_origin_save(ik)
     enddo
  ELSE
     do ik = 1, nkstot/2
        xk(:,4*ik-3) = xksave(:,2*ik-1)
        xk(:,4*ik-2) = xksave(:,2*ik)
        xk(:,4*ik-1) = xksave(:,2*ik-1)
        xk(:,4*ik) = xksave(:,2*ik)
        wk(4*ik-3) = wksave(2*ik-1)
        wk(4*ik-2) = wksave(2*ik)
        wk(4*ik-1) = wksave(2*ik-1)
        wk(4*ik) = wksave(2*ik)
        isk(4*ik-3) = 1
        isk(4*ik-2) = 1
        isk(4*ik-1) = 2
        isk(4*ik) = 2
        diago_bands(4*ik-3)=diago_bands_save(2*ik-1)
        diago_bands(4*ik-2)=diago_bands_save(2*ik)
        diago_bands(4*ik-1)=diago_bands_save(2*ik-1)
        diago_bands(4*ik)=diago_bands_save(2*ik)
        isym_bands(4*ik-3)=isym_bands_save(2*ik-1)
        isym_bands(4*ik-2)=isym_bands_save(2*ik)
        isym_bands(4*ik-1)=isym_bands_save(2*ik-1)
        isym_bands(4*ik)=isym_bands_save(2*ik)

        ik_origin(4*ik-3)=2*(ik_origin_save(2*ik-1)+1)-3
        ik_origin(4*ik-2)=2*(ik_origin_save(2*ik)+1)-3
        ik_origin(4*ik-1)=2*(ik_origin_save(2*ik-1)+1)-1
        ik_origin(4*ik)=2*(ik_origin_save(2*ik)+1)-1
     enddo
  ENDIF
  nkstot = 2 * nkstot
  !
  return
end subroutine set_kup_and_kdw_tpw
