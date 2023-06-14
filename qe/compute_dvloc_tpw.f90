!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE compute_dvloc_tpw (uact, dvloc)
  !----------------------------------------------------------------------
  !
  !! This routine calculates \(dV_\text{bare}/d\tau ) for one 
  !! perturbation with a given q. The displacements are described by a 
  !! vector uact.  
  !
  !
  USE kinds, only : DP
  USE ions_base, ONLY : nat, ityp
  USE cell_base, ONLY : tpiba
  USE fft_base,  ONLY : dffts
  USE fft_interfaces, ONLY: invfft
  USE gvect,     ONLY : eigts1, eigts2, eigts3, mill, g, &
                        ngm
  USE gvecs,     ONLY : ngms
  USE eqv,       ONLY : vlocq
  USE qpoint,    ONLY : xq, eigqts

  USE Coul_cut_2D, ONLY: do_cutoff_2D  
  USE Coul_cut_2D_ph, ONLY : cutoff_localq
  implicit none
  !
  !   The dummy variables
  !
  COMPLEX(DP) :: uact(3*nat)
  ! ... local variables
  COMPLEX(DP) :: dvloc(dffts%nnr)
  !
  ! local variables
  !
  INTEGER ::  na  
  !! counter on atoms
  INTEGER :: mu
  !! counter on modes
  INTEGER :: ig
  !! counter on G vectors
  INTEGER :: nt
  !! the type of atom

  COMPLEX(DP) :: gtau, gu, fact, u1, u2, u3, gu0
  REAL(DP) :: fac

  call start_clock ('dvloc')
  !
  !    We start by computing the contribution of the local potential.
  !    The computation of the derivative of the local potential is done in
  !    reciprocal space while the product with the wavefunction is done in
  !    real space
  !
  dvloc(:) = (0.d0, 0.d0)
  do na = 1, nat
     fact = tpiba * (0.d0, -1.d0) * eigqts (na)
     mu = 3 * (na - 1)
     if (abs (uact (mu + 1) ) + abs (uact (mu + 2) ) + abs (uact (mu + &
          3) ) .gt.1.0d-12) then
        nt = ityp (na)
        u1 = uact (mu + 1)
        u2 = uact (mu + 2)
        u3 = uact (mu + 3)
        gu0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3
        do ig = 1, ngms
           gtau = eigts1 (mill(1,ig), na) * eigts2 (mill(2,ig), na) * &
                  eigts3 (mill(3,ig), na)
           gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) * u3
           dvloc(dffts%nl(ig))=dvloc(dffts%nl(ig))+vlocq(ig,nt)*gu*&
                fact * gtau
        enddo
        IF (do_cutoff_2D) then  
           call cutoff_localq( dvloc, fact, u1, u2, u3, gu0, nt, na) 
        ENDIF
        !
     endif
  enddo
  !
  ! Now we compute dV_loc/dtau in real space
  !
  CALL invfft ('Rho', dvloc, dffts)

  call stop_clock ('dvloc')
  return
END SUBROUTINE compute_dvloc_tpw
