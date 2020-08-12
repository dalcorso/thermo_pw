!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!  
!
!------------------------------------------------
SUBROUTINE addnlcc_zstar_eu_us_tpw( drhoscf ) 
!------------------------------------------------

  USE kinds,     ONLY : DP
  USE ions_base, ONLY : nat
  USE funct,     ONLY : dft_is_gradient, dft_is_nonlocc
  USE scf,       ONLY : rho, rho_core
  USE gvect,     ONLY : g
  USE cell_base, ONLY : omega, alat
  USE fft_base,  ONLY : dfftp
  USE noncollin_module, ONLY : nspin_lsda, nspin_gga, nspin_mag
  USE uspp,      ONLY : nlcc_any
  USE fft_base,  ONLY : dfftp

  USE zstar_add, ONLY : zstareu0_rec
  USE qpoint,    ONLY : xq
  USE modes,     ONLY : npert, nirr
  USE eqv,       ONLY : dmuxc
  USE gc_lr,     ONLY : grho, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s

  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum

  IMPLICIT NONE
  
  COMPLEX(DP) :: drhoscf (dfftp%nnr,nspin_mag,3)

  INTEGER :: nrtot, ipert, is, is1, irr, ir, mode
  INTEGER :: imode0, npe, ipol

  REAL(DP) :: fac
  COMPLEX(DP), EXTERNAL :: zdotc
  
  COMPLEX(DP), ALLOCATABLE :: drhoc(:)
  COMPLEX(DP), ALLOCATABLE :: dvaux(:,:), zstareu0_wrk(:,:)

  zstareu0_rec=(0.0_DP,0.0_DP)
  IF (.NOT.nlcc_any) RETURN

  ALLOCATE( drhoc(dfftp%nnr) )
  ALLOCATE( dvaux(dfftp%nnr,nspin_mag) )
  ALLOCATE( zstareu0_wrk(3,3*nat)  )

  nrtot = dfftp%nr1 * dfftp%nr2 * dfftp%nr3
  zstareu0_wrk=(0.0_DP, 0.0_DP)
  fac = 1.d0 / DBLE (nspin_lsda)
  rho%of_r(:,1) = rho%of_r(:,1) + rho_core(:)
  DO ipol = 1, 3
     imode0 = 0
     DO irr = 1, nirr
        npe = npert(irr)
        !
        !  compute the exchange and correlation potential for this mode
        !
        DO ipert = 1, npe
           mode = imode0 + ipert
           CALL addcore (mode, drhoc)
           dvaux = (0.0_dp,0.0_dp)
           DO is = 1, nspin_lsda
              DO is1 = 1, nspin_mag
                 DO ir = 1, dfftp%nnr
                    dvaux (ir, is) = dvaux (ir, is) +     &
                         dmuxc (ir, is, is1) *            &
                         drhoscf (ir, is1, ipol)
                 ENDDO
              ENDDO
           END DO
           !
           ! add gradient correction to xc, NB: if nlcc is true we need to add here
           ! its contribution. grho contains already the core charge
           !

           IF ( dft_is_gradient() ) &
              CALL dgradcorr (dfftp, rho%of_r, grho, dvxc_rr, dvxc_sr, &
              dvxc_ss, dvxc_s, xq, drhoscf (1, 1, ipol),  &
              nspin_mag, nspin_gga, g, dvaux)
           IF (dft_is_nonlocc()) CALL dnonloccorr( rho%of_r, & 
                                      drhoscf(1, 1, ipol), xq, dvaux )

           DO is = 1, nspin_lsda
              zstareu0_wrk(ipol,mode) = zstareu0_wrk(ipol,mode) -             &
                   omega * fac / DBLE(nrtot) *         &
                   zdotc(dfftp%nnr,dvaux(1,is),1,drhoc,1) 
           END DO
        END DO
        imode0 = imode0 + npe
     END DO
  END DO
  rho%of_r(:,1) = rho%of_r(:,1) - rho_core
  call mp_sum ( zstareu0_wrk, intra_bgrp_comm )
!
! Note: drhoscf is already collected among pools.
!       All pools make the same calculation, so we do not need
!       to collect the results over all pools.
!
  zstareu0_rec=zstareu0_rec+zstareu0_wrk

  DEALLOCATE( zstareu0_wrk )
  DEALLOCATE( drhoc )
  DEALLOCATE( dvaux )

  RETURN
END SUBROUTINE addnlcc_zstar_eu_us_tpw
