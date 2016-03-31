!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
subroutine symeq (dvsym)
  !---------------------------------------------------------------------
  !
  !     This routine symmetrize the change of the potential due to an
  !     periodic perturbation of wavevector q.  
  !
  USE kinds, only : DP
  USE cell_base, ONLY : at
  USE fft_base,  only : dfftp
  USE symm_base, only : s, ftau
  USE noncollin_module, only : nspin_lsda, nspin_mag
  USE modes, ONLY : gi, nsymq
  USE constants, ONLY : tpi
  implicit none

  complex(DP) :: dvsym (dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, nspin_mag)
  complex(DP), allocatable ::  aux (:,:,:)
  ! the potential to symmetrize
  ! auxiliary quantity
  real(DP) :: gf(3), n(3)
  !  temp variables
  complex(DP) :: term (3, 48), phase (48)

  integer :: is, ri, rj, rk, i, j, k, irot, ipol, jpol
  ! counter on spin polarization
  ! the rotated points
  ! the point
  ! counter on symmetries
  ! counter on polarizations

  if (nsymq == 1) return
  allocate (aux(dfftp%nr1x , dfftp%nr2x , dfftp%nr3x))
  !
  ! Here we symmetrize with respect to the small group of q
  !
  n(1) = tpi / DBLE (dfftp%nr1)
  n(2) = tpi / DBLE (dfftp%nr2)
  n(3) = tpi / DBLE (dfftp%nr3)
  do irot = 1, nsymq
     gf(:) =  gi (1,irot) * at (1, :) * n(:) + &
              gi (2,irot) * at (2, :) * n(:) + &
              gi (3,irot) * at (3, :) * n(:)
     term (:, irot) = CMPLX(cos (gf (:) ), sin (gf (:) ) ,kind=DP)
  enddo

  do is = 1, nspin_lsda
     aux(:,:,:) = dvsym(:,:,:,is)
     dvsym(:,:,:,is) = (0.d0, 0.d0)
     !
     !  symmmetrize
     !
     phase(:) = (1.d0, 0.d0)
     do k = 1, dfftp%nr3
        do j = 1, dfftp%nr2
           do i = 1, dfftp%nr1
              do irot = 1, nsymq
                 call ruotaijk (s(1,1,irot), ftau(1,irot), i, j, k, &
                                dfftp%nr1, dfftp%nr2, dfftp%nr3, ri, rj, rk)
                 !
                 !    ruotaijk find the rotated of i,j,k with the inverse of S
                 !
                 dvsym(i,j,k,is) = dvsym(i,j,k,is) + aux(ri,rj,rk)*phase(irot)
              enddo
              do irot = 1, nsymq
                 phase (irot) = phase (irot) * term (1, irot)
              enddo
           enddo
           do irot = 1, nsymq
              phase (irot) = phase (irot) * term (2, irot)
           enddo
        enddo
        do irot = 1, nsymq
           phase (irot) = phase (irot) * term (3, irot)
        enddo
     enddo
     dvsym(:,:,:,is) = dvsym(:,:,:,is) / DBLE(nsymq)
  enddo
  deallocate (aux)
  return
end subroutine symeq

