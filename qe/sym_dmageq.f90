!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine sym_dmageq (dvsym)
  !---------------------------------------------------------------------
  !
  !     This routine symmetrize the change of the potential due to an
  !     electric field perturbation. It is assumed that the perturbations
  !     are on the basis of the crystal
  !
  !
  USE kinds, only : DP
  USE constants, ONLY : tpi
  USE cell_base,only : at, bg
  USE fft_base, only : dfftp
  USE symm_base,only : sname, s, ftau, t_rev, invs
  USE modes,    ONLY : gi, nsymq
  USE lsda_mod, only : nspin
  implicit none

  complex(DP) :: dvsym (dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, nspin)
  complex(DP), allocatable ::  aux (:,:,:,:)
  complex(DP) ::  dmags(3), mag(3), magrot(3)
  ! the potential to symmetrize
  ! auxiliary quantity
  real(DP) :: gf(3), n(3)
  !  temp variables
  complex(DP) :: term (3, 48), phase (48)


  integer :: is, ri, rj, rk, i, j, k, irot, ipol, jpol, kpol
  ! counter on spin polarization
  ! the rotated points
  ! the point
  ! counter on symmetries
  ! counter on polarizations

  if (nsymq == 1) return
  allocate (aux(dfftp%nr1x , dfftp%nr2x , dfftp%nr3x , 3))

  do is = 2, 4
     aux(:,:,:,is-1) = dvsym(:,:,:,is)
     dvsym(:,:,:,is) = (0.d0, 0.d0)
  enddo
  !
  !  symmmetrize
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
              dmags=(0.d0,0.d0)
              do is=1,3
                 dmags(is)=aux(ri,rj,rk,is)*phase(irot)
              enddo
              do kpol = 1, 3
                 mag(kpol)=bg(1,kpol)*dmags(1) + &
                           bg(2,kpol)*dmags(2) + &
                           bg(3,kpol)*dmags(3)
              enddo
! rotate the magnetic moment
              do kpol = 1, 3
                 magrot(kpol) = s(1,kpol,invs(irot))*mag(1) + &
                                s(2,kpol,invs(irot))*mag(2) + &
                                s(3,kpol,invs(irot))*mag(3)
              enddo
              if (sname(irot)(1:3)=='inv') magrot=-magrot
              if(t_rev(irot).eq.1) magrot=-magrot
! go back to cartesian coordinates
              do kpol = 1, 3
                 mag(kpol)=at(kpol,1)*magrot(1) + &
                           at(kpol,2)*magrot(2) + &
                           at(kpol,3)*magrot(3)
              enddo
              dvsym(i,j,k,2) = dvsym(i,j,k,2) + mag(1)
              dvsym(i,j,k,3) = dvsym(i,j,k,3) + mag(2)
              dvsym(i,j,k,4) = dvsym(i,j,k,4) + mag(3)
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

  do is=2,4
     dvsym(:,:,:,is) = dvsym(:,:,:,is) / DBLE(nsymq)
  enddo
  deallocate (aux)
  return
end subroutine sym_dmageq
