!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine compute_drhous_nc_tpw (drhous, dbecsum, wgg, becq, alpq)
  !-----------------------------------------------------------------------
  !
  !    This routine computes the part of the change of the charge density
  !    which is due to the orthogonalization constraint on wavefunctions
  !
  !
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat
  USE lsda_mod,   ONLY : lsda, nspin, current_spin, isk
  USE klist,      ONLY : xk, wk, ngk, igk_k
  USE buffers,    ONLY : get_buffer
  USE fft_base,   ONLY : dffts, dfftp
  USE fft_interfaces, ONLY : invfft
  USE gvecs,      ONLY : nls
  USE wvfct,      ONLY : npwx, nbnd
  USE noncollin_module, ONLY : npol, nspin_mag
  USE wavefunctions_module, ONLY : evc
  USE uspp,       ONLY : okvan, nkb, vkb
  USE uspp_param, ONLY : nhm

  USE qpoint,     ONLY : nksq, ikks, ikqs
  USE eqv,        ONLY : evq, dvpsi, dpsi
  USE control_lr, ONLY : lgamma, nbnd_occ
  USE control_ph, ONLY : zeu, zue
  USE zstar_add,  ONLY : done_start_zstar

  USE efield_mod, ONLY : zstarue0
  USE units_ph,   ONLY : lrwfc, iuwfc
  USE becmod,     ONLY : bec_type
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum

  implicit none
  !
  !     the dummy variables
  !

  complex(DP) :: dbecsum (nhm, nhm, nat, nspin, 3 * nat), &
         drhous (dfftp%nnr, nspin_mag, 3 * nat)
  !output:the derivative of becsum
  ! output: add the orthogonality term
  type (bec_type) :: becq(nksq), & ! (nkb, nbnd)
                     alpq (3, nksq)
  ! input: the becp with psi_{k+q}
  ! input: the alphap with psi_{k+q}

  real(DP) :: wgg (nbnd, nbnd, nksq)
  ! input: the weights

  integer :: npw, npwq, ik, ikq, ikk, ig, nu_i, ibnd, jpol, ios
  ! counter on k points
  ! the point k+q
  ! record for wfcs at k point
  ! counter on spin
  ! counter on g vectors
  ! counter on modes
  ! counter on the bands
  ! integer variable for I/O control

  real(DP) :: weight
  ! the weight of the k point

  complex(DP), allocatable :: evcr (:,:,:), dpsi_save(:,:)
  ! the wavefunctions in real space
  COMPLEX(DP) :: zdotc

  if (.not.okvan) return

  call start_clock ('com_drhous')

  allocate (evcr( dffts%nnr, npol, nbnd))
  IF (zeu.or.zue) allocate ( dpsi_save ( npwx*npol , nbnd))
  !
  zstarue0  = (0.d0, 0.d0)
  drhous(:,:,:) = (0.d0, 0.d0)
  dbecsum  = (0.d0, 0.d0)

  do ik = 1, nksq
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npw = ngk(ikk)
     npwq= ngk(ikq)
     weight = wk (ikk)
     if (lsda) current_spin = isk (ikk)
     !
     !   For each k point we construct the beta functions
     !
     call init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb)
     !
     !   Read the wavefunctions at k and transform to real space
     !

     call get_buffer(evc, lrwfc, iuwfc, ikk)
     evcr = (0.d0, 0.d0)
     do ibnd = 1, nbnd
        do ig = 1, npw
           evcr (nls (igk_k(ig,ikk) ), 1, ibnd) = evc (ig, ibnd)
           evcr (nls (igk_k(ig,ikk) ), 2, ibnd) = evc (ig+npwx, ibnd)
        enddo
        CALL invfft ('Wave', evcr (:, 1, ibnd), dffts)
        CALL invfft ('Wave', evcr (:, 2, ibnd), dffts)
     enddo
     !
     !   Read the wavefunctions at k+q
     !
     if (.not.lgamma.and.nksq.gt.1) call get_buffer(evq, lrwfc, iuwfc, ikq)
     !
     !   And compute the contribution of this k point to the change of
     !   the charge density
     !
     do nu_i = 1, 3 * nat
        call incdrhous_nc (drhous (1, 1, nu_i), weight, ik, &
                 dbecsum (1, 1, 1, 1, nu_i), evcr, wgg, becq, alpq, nu_i)
!
!   After this call dpsi contains the change of the wavefunctions due
!   to the change of the orthogonality constraint. We use this term
!   to calculate the part of the effective charges that depends
!   on this orthogonality term.
!
       IF (zeu.or.zue) THEN
          dpsi_save = dpsi
          DO jpol=1,3
             dvpsi=(0.0,0.0)
             call dvpsi_e(ik, jpol)
!
! NB: The minus sign is due to the fact that dpsi_save contains
!     -|psi_j><psi_j| dS/du |psi_i>
!
             DO ibnd=1,nbnd_occ(ikk)
                zstarue0(nu_i,jpol)=zstarue0(nu_i,jpol) - wk(ikk)* &
                          zdotc(npw,dpsi_save(1,ibnd),1,dvpsi(1,ibnd),1)
                zstarue0(nu_i,jpol)=zstarue0(nu_i,jpol) - wk(ikk)* &
                          zdotc(npw,dpsi_save(npwx+1,ibnd),1,dvpsi(npwx+1,ibnd),1)
             ENDDO
          ENDDO
       ENDIF
     enddo
  enddo

  IF (zeu.OR.zue.AND..NOT.done_start_zstar) &
         CALL mp_sum ( zstarue0, intra_bgrp_comm )

  deallocate(evcr)
  IF (zeu.or.zue) deallocate ( dpsi_save )

  call stop_clock ('com_drhous')
  return

end subroutine compute_drhous_nc_tpw
