!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine add_zstar_ue_tpw (imode0, npe)
  !-----------------------------------------------------------------------
  ! add the contribution of the modes imode0+1 -> imode+npe
  ! to the effective charges Z(Us,E) (Us=scf,E=bare)
  !
  ! trans =.true. is needed for this calculation to be meaningful
  !
  USE kinds, only : DP
  USE klist, ONLY : xk, wk, ngk, igk_k
  USE uspp,  ONLY : vkb
  USE wvfct, ONLY : npwx, npw
  USE wavefunctions,  ONLY: evc
  USE noncollin_module,      ONLY: noncolin
  USE buffers,  ONLY: get_buffer
  USE qpoint,   ONLY: npwq, nksq, ikks
  USE eqv,      ONLY: dpsi, dvpsi
  USE efield_mod, ONLY: zstarue0_rec, zstarue0
  USE control_lr, ONLY : nbnd_occ
  USE units_ph,   ONLY : iudwf, lrdwf
  USE units_lr,   ONLY : iuwfc, lrwfc

  USE mp_global,  ONLY : inter_pool_comm, intra_bgrp_comm
  USE mp,         ONLY : mp_sum

  implicit none

  integer, intent(in) :: imode0, npe

  integer :: ibnd, jpol, ipert, nrec, mode, ik, ikk
  ! counter on bands
  ! counter on polarization
  ! counter on pertubations
  ! counter on records
  ! counter on modes
  ! counter on k points

  real(DP) :: weight

  complex(DP), EXTERNAL :: zdotc


  call start_clock('add_zstar_ue')
  zstarue0_rec=(0.0_DP,0.0_DP)
  do ik = 1, nksq
     ikk=ikks(ik)
     npw=ngk(ikk)
     npwq = npw
     weight = wk (ikk)
     if (nksq.gt.1) call get_buffer (evc, lrwfc, iuwfc, ikk)
     call init_us_2 (npw, igk_k(1,ikk), xk (1, ikk), vkb)
     do jpol = 1, 3
        !
        ! read/compute DeltaV*psi(bare) for electric field
        !
        call dvpsi_e_tpw (ik, jpol)
        !
        do ipert = 1, npe
           mode = imode0 + ipert
           nrec = (ipert - 1) * nksq + ik
           !
           ! read dpsi(scf)/du for phonon mode # mode
           !
           call get_buffer (dpsi, lrdwf, iudwf, nrec)
           do ibnd = 1, nbnd_occ(ikk)
              zstarue0_rec(mode,jpol)=zstarue0_rec(mode,jpol)-weight * &
                   2.0_DP*zdotc (npw, dvpsi (1, ibnd), 1, dpsi (1, ibnd), 1)
              IF (noncolin) &
                 zstarue0_rec(mode,jpol)=zstarue0_rec(mode,jpol)-weight*&
                 2.0_DP*zdotc (npw, dvpsi (1+npwx, ibnd), 1, dpsi (1+npwx, ibnd), 1)
           enddo
        enddo
     enddo

  enddo
  call mp_sum ( zstarue0_rec, intra_bgrp_comm )
  call mp_sum ( zstarue0_rec, inter_pool_comm )

  zstarue0=zstarue0+zstarue0_rec

  call stop_clock('add_zstar_ue')
  return
end subroutine add_zstar_ue_tpw
