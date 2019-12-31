!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE compute_drhous_nc_tpw (drhous, dbecsum, wgg, becq, alpq)
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
  USE wvfct,      ONLY : npwx, nbnd
  USE noncollin_module, ONLY : npol, nspin_mag
  USE wavefunctions, ONLY : evc
  USE uspp,       ONLY : okvan, nkb, vkb
  USE uspp_param, ONLY : nhm

  USE qpoint,     ONLY : nksq, ikks, ikqs
  USE eqv,        ONLY : evq, dvpsi, dpsi
  USE control_lr, ONLY : lgamma, nbnd_occ
  USE control_ph, ONLY : zeu, zue
  USE zstar_add,  ONLY : done_start_zstar

  USE efield_mod, ONLY : zstarue0
  USE units_lr,   ONLY : lrwfc, iuwfc
  USE becmod,     ONLY : bec_type
  USE partial,    ONLY : done_irr, comp_irr
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum

  IMPLICIT NONE
  !
  !     the dummy variables
  !

  COMPLEX(DP) :: dbecsum (nhm, nhm, nat, nspin, 3 * nat), &
         drhous (dfftp%nnr, nspin_mag, 3 * nat)
  !output:the derivative of becsum
  ! output: add the orthogonality term
  TYPE (bec_type) :: becq(nksq), & ! (nkb, nbnd)
                     alpq (3, nksq)
  ! input: the becp with psi_{k+q}
  ! input: the alphap with psi_{k+q}

  REAL(DP) :: wgg (nbnd, nbnd, nksq)
  ! input: the weights

  INTEGER :: npw, npwq, ik, ikq, ikk, ig, nu_i, ibnd, jpol, ios
  ! counter on k points
  ! the point k+q
  ! record for wfcs at k point
  ! counter on spin
  ! counter on g vectors
  ! counter on modes
  ! counter on the bands
  ! integer variable for I/O control

  REAL(DP) :: weight
  ! the weight of the k point

  COMPLEX(DP), ALLOCATABLE :: evcr (:,:,:), dpsi_save(:,:)
  ! the wavefunctions in real space
  COMPLEX(DP) :: zdotc
  LOGICAL :: add_zstar

  IF (.NOT.okvan) RETURN
  add_zstar= (zeu.OR.zue.AND.(.NOT.done_start_zstar)).AND.comp_irr(0).AND. &
             (.NOT. done_irr(0))

  CALL start_clock ('com_drhous')

  ALLOCATE (evcr( dffts%nnr, npol, nbnd))
  IF (add_zstar) allocate ( dpsi_save ( npwx*npol , nbnd))
  !
  zstarue0  = (0.d0, 0.d0)
  drhous(:,:,:) = (0.d0, 0.d0)
  dbecsum  = (0.d0, 0.d0)

  DO ik = 1, nksq
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npw = ngk(ikk)
     npwq= ngk(ikq)
     weight = wk (ikk)
     IF (lsda) current_spin = isk (ikk)
     !
     !   For each k point we construct the beta functions
     !
     CALL init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb)
     !
     !   Read the wavefunctions at k and transform to real space
     !

     CALL get_buffer(evc, lrwfc, iuwfc, ikk)
     evcr = (0.d0, 0.d0)
     DO ibnd = 1, nbnd
        DO ig = 1, npw
           evcr (dffts%nl (igk_k(ig,ikk) ), 1, ibnd) = evc (ig, ibnd)
           evcr (dffts%nl (igk_k(ig,ikk) ), 2, ibnd) = evc (ig+npwx, ibnd)
        ENDDO
        CALL invfft ('Wave', evcr (:, 1, ibnd), dffts)
        CALL invfft ('Wave', evcr (:, 2, ibnd), dffts)
     ENDDO
     !
     !   Read the wavefunctions at k+q
     !
     IF (.NOT.lgamma.AND.nksq>1) call get_buffer(evq, lrwfc, iuwfc, ikq)
     !
     !   And compute the contribution of this k point to the change of
     !   the charge density
     !
     DO nu_i = 1, 3 * nat
        CALL incdrhous_nc (drhous (1, 1, nu_i), weight, ik, &
                 dbecsum (1, 1, 1, 1, nu_i), evcr, wgg, becq, alpq, nu_i)
!
!   After this call dpsi contains the change of the wavefunctions due
!   to the change of the orthogonality constraint. We use this term
!   to calculate the part of the effective charges that depends
!   on this orthogonality term.
!
       
       IF (add_zstar) THEN
          dpsi_save = dpsi
          DO jpol=1,3
             dvpsi=(0.0,0.0)
             CALL dvpsi_e_tpw(ik, jpol)
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
     ENDDO
  ENDDO

  IF (add_zstar) CALL mp_sum ( zstarue0, intra_bgrp_comm )

  DEALLOCATE(evcr)
  IF (add_zstar) DEALLOCATE ( dpsi_save )

  CALL stop_clock ('com_drhous')
  RETURN

END SUBROUTINE compute_drhous_nc_tpw
