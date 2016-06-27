!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine set_kup_and_kdw_tpw (xk, wk, isk, nkstot, npk)
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
  !
  integer :: ik, iq, ikq
  !
  !
  if (2*nkstot > npk) call errore ('set_kup_and_kdw_tpw','too many k points',nkstot)
  xksave(:,1:nkstot) = xk(:,1:nkstot)
  wksave(1:nkstot) = wk(1:nkstot)
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
  enddo
  nkstot = 2 * nkstot

  return
end subroutine set_kup_and_kdw_tpw
