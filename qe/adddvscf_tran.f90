!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine adddvscf_tran (spin_psi0, ik)
  !----------------------------------------------------------------------
  !
  !     This routine computes the contribution of the selfconsistent
  !     change of the potential to the known part of the linear
  !     system and adds it to dvpsi.
  !     It implements the second term in Eq. B30 of PRB 64, 235118 (2001).
  !

  USE kinds,      ONLY : DP
  USE uspp_param, ONLY : upf, nh
  USE uspp,       ONLY : vkb, okvan
! modules from pwcom
  USE lsda_mod,   ONLY : lsda, current_spin, isk
  USE ions_base,  ONLY : ntyp => nsp, nat, ityp
  USE wvfct,      ONLY : nbnd, npwx
  USE noncollin_module, ONLY : noncolin, npol
! modules from phcom
  USE qpoint,     ONLY : npwq
  USE lrus,       ONLY : int3, int3_nc, becp1
  USE eqv,        ONLY : dvpsi
  implicit none
  !
  !   The dummy variables
  !
  integer :: ik, spin_psi0
  ! input: the k point
  ! input: the spin of the unpertubed wavefunctions psi0
  !
  !   And the local variables
  !
  integer :: na, nt, ibnd, ih, jh, ijkb0, ikb, jkb, is, js, ijs
  ! counter on atoms
  ! counter on atomic types
  ! counter on bands
  ! counter on beta functions
  ! counter on beta functions
  ! auxiliary variable for indexing
  ! counter on the k points
  ! counter on vkb
  ! counter on vkb
  complex(DP) :: sump
  ! auxiliary variable

  if (.not.okvan) return
  call start_clock ('adddvscf')
  ijkb0 = 0
  do nt = 1, ntyp
     if (upf(nt)%tvanp  ) then
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              !
              !   we multiply the integral for the becp term and the beta_n
              !
              do ibnd = 1, nbnd
                 do ih = 1, nh (nt)
                    ikb = ijkb0 + ih
                    sump = (0.d0, 0.d0)
                    do jh = 1, nh (nt)
                       jkb = ijkb0 + jh
                       sump = sump + int3 (ih, jh, 1, na, 1)*&
                                   becp1(ik)%k(jkb, ibnd)
                    enddo
                    call zaxpy(npwq,sump,vkb(1,ikb),1,dvpsi(1,ibnd),1)
                 enddo
              enddo
              ijkb0 = ijkb0 + nh (nt)
           endif
        enddo
     else
        do na = 1, nat
           if (ityp (na) .eq.nt) ijkb0 = ijkb0 + nh (nt)
        enddo
     endif
  enddo

  call stop_clock ('adddvscf')
  return
end subroutine adddvscf_tran
