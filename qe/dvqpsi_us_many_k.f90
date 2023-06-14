!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvqpsi_us_many_k (ik, id, uact, becp1, alphap, aux1)
  !----------------------------------------------------------------------
  !
  !! This routine calculates \(dV_\text{bare}/d\tau \cdot \psi\) for one 
  !! perturbation with a given q. The displacements are described by a 
  !! vector u.  
  !! The result is stored in \(\text{dvpsi}\). The routine is called for
  !! each k-point and for each pattern u. It computes simultaneously all 
  !! the bands. It implements Eq. (B29) of PRB 64, 235118 (2001). The 
  !! contribution of the local pseudopotential is calculated here, that 
  !! of the nonlocal pseudopotential in \(\texttt{dvqpsi_us_only}\).
  !
  !
  USE kinds, only : DP
  USE funct,     ONLY : dft_is_nonlocc
  USE xc_lib,    ONLY : xclib_dft_is
  USE ions_base, ONLY : nat, ityp
  USE cell_base, ONLY : tpiba
  USE fft_base,  ONLY : dfftp, dffts
  USE fft_interfaces, ONLY: fwfft, invfft
  USE gvect,     ONLY : eigts1, eigts2, eigts3, mill, g, &
                        ngm
  USE gvecs,     ONLY : ngms, doublegrid
  USE lsda_mod,  ONLY : lsda, isk
  USE scf,       ONLY : rho, rho_core
  USE noncollin_module, ONLY : nspin_gga, nspin_mag, npol
  use uspp_param,ONLY : upf
  USE wvfct,     ONLY : nbnd, npwx
  USE wavefunctions,  ONLY: evc
  USE nlcc_ph,    ONLY : drc
  USE uspp,       ONLY : nlcc_any
  USE eqv,        ONLY : dvpsi, dmuxc, vlocq
  USE qpoint,     ONLY : xq, eigqts, ikqs, ikks
  USE klist,      ONLY : ngk, igk_k
  USE gc_lr,      ONLY: grho, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s
  USE qpoint,     ONLY : nksq
  USE becmod,     ONLY : bec_type

  USE Coul_cut_2D, ONLY: do_cutoff_2D  
  USE Coul_cut_2D_ph, ONLY : cutoff_localq
  implicit none
  !
  !   The dummy variables
  !
  INTEGER, INTENT(in) :: ik, id
  !! input: the k point
  COMPLEX(DP) :: uact(3*nat)
  !! input: the pattern of displacements
  COMPLEX(DP)::aux1(dffts%nnr)
  LOGICAL :: addnlcc
  !!
  !
  TYPE(bec_type) :: becp1(nksq), alphap(3,nksq)
  !
  ! ... local variables
  !
  INTEGER ::  na  
  !! counter on atoms
  INTEGER :: mu
  !! counter on modes
  INTEGER :: npw
  !! Number of pw
  INTEGER :: ikk
  !! the point k
  INTEGER :: npwq
  !! Number of q
  INTEGER :: ikq
  !! k-q index
  INTEGER :: iks
  !!
  INTEGER :: ig
  !! counter on G vectors
  INTEGER :: nt
  !! the type of atom
  INTEGER :: ibnd
  !! counter on bands
  INTEGER :: ir 
  !! counter on real mesh
  INTEGER :: is
  !! 
  INTEGER :: ip

  complex(DP) :: gtau, gu, fact, u1, u2, u3, gu0
  complex(DP) , allocatable, target :: aux (:)
!  complex(DP) , allocatable :: aux1 (:), aux2 (:)
  complex(DP) , allocatable ::  aux2 (:)
  complex(DP) , pointer :: auxs (:)
  REAL(DP) :: fac
  COMPLEX(DP), ALLOCATABLE :: drhoc(:)

  call start_clock ('dvqpsi_us')
  if (nlcc_any.and.addnlcc) then
     allocate (drhoc( dfftp%nnr))
     allocate (aux( dfftp%nnr))
     allocate (auxs(dffts%nnr))
  endif
!  allocate (aux1(dffts%nnr))
  allocate (aux2(dffts%nnr))
  !
  !    We start by computing the contribution of the local potential.
  !
  dvpsi(:,:) = (0.d0, 0.d0)
  !
  ! Now we compute dV_loc/dtau in real space
  !
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  do ibnd = 1, nbnd
     do ip=1,npol
        aux2(:) = (0.d0, 0.d0)
        if (ip==1) then
           do ig = 1, npw
              aux2 (dffts%nl (igk_k (ig,ikk) ) ) = evc (ig, ibnd)
           enddo
        else
           do ig = 1, npw
              aux2 (dffts%nl (igk_k (ig,ikk) ) ) = evc (ig+npwx, ibnd)
           enddo
        end if
        !
        !  This wavefunction is computed in real space
        !
        CALL invfft ('Wave', aux2, dffts)
        do ir = 1, dffts%nnr
           aux2 (ir) = aux2 (ir) * aux1 (ir)
        enddo
        !
        ! and finally dV_loc/dtau * psi is transformed in reciprocal space
        !
        CALL fwfft ('Wave', aux2, dffts)
        if (ip==1) then
           do ig = 1, npwq
              dvpsi (ig, ibnd) = aux2 (dffts%nl (igk_k (ig,ikq) ) ) 
           enddo
        else
           do ig = 1, npwq
              dvpsi (ig+npwx, ibnd) = aux2 (dffts%nl (igk_k (ig,ikq) ) ) 
           enddo
        end if
     enddo
  enddo
  !
  deallocate (aux2)
!  deallocate (aux1)
  if (nlcc_any.and.addnlcc) then
     deallocate (drhoc)
     deallocate (aux)
     deallocate (auxs)
  endif
  !
  !   We add the contribution of the nonlocal potential in the US form
  !   First a term similar to the KB case.
  !   Then a term due to the change of the D coefficients.
  !
  call dvqpsi_us_only_tpw (ik, id, uact, becp1, alphap)

  call stop_clock ('dvqpsi_us')
  return
end subroutine dvqpsi_us_many_k
