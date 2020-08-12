!
! Copyright (C) 2001-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE zstar_eu_tpw(drhoscf)
  !-----------------------------------------------------------------------
  ! calculate the effective charges Z(E,Us) (E=scf,Us=bare)
  ! This expression is obtained as the derivative of the forces with
  ! respect to the electric field.
  !
  ! epsil =.true. is needed for this calculation to be meaningful
  !
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : bg
  USE ions_base, ONLY : nat, zv, ityp, atm
  USE klist,     ONLY : wk, xk, ngk, igk_k
  USE symme,     ONLY : symtensor
  USE buffers,   ONLY : get_buffer
  USE wvfct,     ONLY : npw, npwx
  USE uspp,      ONLY : vkb
  USE fft_base,  ONLY : dffts
  use noncollin_module, ONLY : nspin_mag, noncolin
  USE wavefunctions,  ONLY: evc

  USE modes,     ONLY : u, nirr, npert
  USE qpoint,    ONLY : npwq, nksq, ikks
  USE eqv,       ONLY : dvpsi, dpsi
  USE efield_mod,   ONLY : zstareu0, zstareu
  USE zstar_add, ONLY : zstareu0_rec
  USE units_ph,  ONLY : iudwf, lrdwf
  USE units_lr,  ONLY : lrwfc, iuwfc
  USE control_ph,ONLY : done_zeu
  USE control_lr,ONLY : nbnd_occ
  USE phus,      ONLY : alphap
  USE lrus,      ONLY : becp1
  USE ph_restart, ONLY : ph_writefile
  USE io_global, ONLY : stdout

  USE mp_pools,  ONLY : inter_pool_comm
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE ldaU,      ONLY : lda_plus_u

  IMPLICIT NONE

  COMPLEX(DP) :: drhoscf (dffts%nnr, nspin_mag, 3)
  ! output: the change of the scf charge (smooth part only)

  INTEGER :: ipol, jpol, icart, na, nu, mu, imode0, irr, &
       ipert, nrec, mode, ik, ibnd, ikk, ierr
  ! counters
  REAL(DP) :: weight
  !  auxiliary space
  COMPLEX(DP), ALLOCATABLE :: zstareu0_wrk(:,:)
  !
  COMPLEX(DP), EXTERNAL :: zdotc
  !  scalar product
  !
  CALL start_clock ('zstar_eu')

  ALLOCATE( zstareu0_wrk( 3, 3 * nat ) )
  zstareu0_wrk(:,:)=(0.0_DP, 0.0_DP)

  DO ik = 1, nksq
     ikk=ikks(ik)
     npw=ngk(ikk)
     npwq = npw
     weight = wk (ikk)
     if (nksq > 1) CALL get_buffer(evc, lrwfc, iuwfc, ikk)
     CALL init_us_2 (npw, igk_k(1,ikk), xk (1, ikk), vkb)
     imode0 = 0
     DO irr = 1, nirr
        DO ipert = 1, npert (irr)
           mode = ipert+imode0
           dvpsi(:,:) = (0.d0, 0.d0)
           !
           ! recalculate  DeltaV*psi(ion) for mode nu
           !
           CALL dvqpsi_us_only (ik, u (1, mode), becp1, alphap)
           !
           ! DFPT+U: add the bare variation of the Hubbard potential 
           !
           IF (lda_plus_u) CALL dvqhub_barepsi_us (ik, u(:,mode))
           !
           DO jpol = 1, 3
              nrec = (jpol - 1) * nksq + ik
              !
              ! read dpsi(scf)/dE for electric field in jpol direction
              !
              CALL get_buffer(dpsi, lrdwf, iudwf, nrec)
              DO ibnd=1,nbnd_occ(ikk)
                 zstareu0_wrk(jpol,mode)=zstareu0_wrk(jpol,mode)-weight*&
                     ( zdotc(npw,dpsi(1,ibnd),1,dvpsi(1,ibnd),1) + &
                      zdotc(npw,dvpsi(1,ibnd),1,dpsi(1,ibnd),1) )
                 IF (noncolin) &
                    zstareu0_wrk(jpol,mode)=zstareu0_wrk(jpol, mode)- &
                      weight*(zdotc(npw,dpsi(npwx+1,ibnd),1, &
                              dvpsi(npwx+1,ibnd),1) +    &
                    zdotc(npw,dvpsi(npwx+1,ibnd),1,dpsi(npwx+1,ibnd),1) )
              ENDDO
           ENDDO
        ENDDO
        imode0 = imode0 + npert (irr)
     ENDDO
  ENDDO
  CALL zstar_eu_loc (drhoscf, zstareu0_wrk)

  CALL mp_sum ( zstareu0_wrk, intra_bgrp_comm )
  CALL mp_sum ( zstareu0_wrk, inter_pool_comm )

!  WRITE(6,*) ' term Z^{(1} wrk'
!  CALL tra_write_zstar(zstareu0_wrk, zstareu, .TRUE.)
!  WRITE(6,*) ' term Z^{(1} 0'
!  CALL tra_write_zstar(zstareu0, zstareu, .TRUE.)
!  WRITE(6,*) ' term Z^{(1} rec'
!  CALL tra_write_zstar(zstareu0_rec, zstareu, .TRUE.)

  zstareu0_wrk = zstareu0 + zstareu0_wrk + zstareu0_rec 
  !
  ! bring the mode index to cartesian coordinates
  !
  CALL tra_write_zstar(zstareu0_wrk, zstareu, .FALSE.)
  !
  !  symmetrize
  !
  CALL symtensor ( nat, zstareu )
  !
  ! add the diagonal part
  !
  DO ipol = 1, 3
     DO na = 1, nat
        zstareu (ipol, ipol, na) = zstareu (ipol, ipol, na) + zv (ityp ( na) )
     ENDDO
  ENDDO

  done_zeu=.TRUE.
  CALL summarize_zeu()
  CALL ph_writefile('tensors',0,0,ierr)
  DEALLOCATE ( zstareu0_wrk )

  CALL stop_clock ('zstar_eu')
  RETURN
END SUBROUTINE zstar_eu_tpw
