!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dveqpsi_us (ik)
  !----------------------------------------------------------------------
  !
  ! This routine calculates applies e^iqr to a wavefunction.
  ! In the US and PAW case apply e^iqr K(r,r',r'')
  !
  USE kinds, only : DP
  USE ions_base, ONLY : nat, ityp
  USE cell_base, ONLY : tpiba
  USE fft_base,   ONLY: dfftp, dffts
  USE fft_interfaces, ONLY: fwfft, invfft
  USE gvecs,     ONLY : nls
  USE noncollin_module, ONLY : noncolin
  USE wvfct,     ONLY : nbnd, npwx

  USE wavefunctions_module,  ONLY: evc
  USE eqv,        ONLY : dvpsi
  USE qpoint,     ONLY : ikks, ikqs
  USE klist,      ONLY : ngk, igk_k 
  implicit none
  !
  !   The dummy variables
  !

  integer :: ik, npw, npwq, ikk, ikq
  ! input: the k point

  integer :: ibnd, ig
  ! counter on bands
  ! counter on G vectors
  ! counter on polarizations
  ! counter on reciprocal mesh

  complex(DP) , allocatable ::  aux2 (:)
  ! work space
  logical :: htg

  call start_clock ('dveqpsi_us')
  allocate (aux2(dffts%nnr))
  !
  !  The unperturbed wavefunctions must be multiplied by e^{iqr}.
  !  This means that we have to order the coefficients with the mesh
  !  of k+q. We do this by distributing the vectors on the FFT mesh
  !  with the indices of k+G (igk) and then saving them with the mesh 
  !  of k+q+G (igkq)
  !
  dvpsi(:,:) = (0.d0, 0.d0)

  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)

  do ibnd = 1, nbnd
     aux2(:) = (0.d0, 0.d0)
     do ig = 1, npw
        aux2 (nls (igk_k (ig,ikk) ) ) = evc (ig, ibnd)
     enddo
     do ig = 1, npwq
        dvpsi (ig, ibnd) = aux2 (nls (igk_k (ig,ikq) ) )
     enddo
     IF (noncolin) THEN
        aux2(:) = (0.d0, 0.d0)
        do ig = 1, npw
           aux2 (nls (igk_k (ig,ikk) ) ) = evc (ig+npwx, ibnd)
        enddo
        do ig = 1, npwq
           dvpsi (ig+npwx, ibnd) = aux2 (nls (igk_k (ig,ikq) ) )
        enddo
     END IF
  enddo
  !
  deallocate (aux2)
  !
  !   We add the contribution of the nonlocal potential in the US form
  !   First a term similar to the KB case.
  !   Then a term due to the change of the D coefficients.
  !
  call dveqpsi_us_only (ik)

  call stop_clock ('dveqpsi_us')
  return
end subroutine dveqpsi_us
