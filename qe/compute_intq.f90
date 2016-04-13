!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine compute_intq
  !----------------------------------------------------------------------
  !
  !     This routine computes the contribution of the selfconsistent
  !     change of the potential to the known part of the linear
  !     system and adds it to dvpsi.
  !
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE noncollin_module,     ONLY : noncolin
  USE cell_base,            ONLY : omega
  USE uspp,                 ONLY : okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm

  USE optical,              ONLY : intq
  USE qpoint,               ONLY : xq, eigqts
  USE control_lr,           ONLY : lgamma

  USE mp_bands,  ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum

  implicit none

  integer :: na, ig, nt, ir, ih, jh
  ! countera

  real(DP), allocatable ::  ylmk0 (:,:)
  ! the modulus of q+G
  ! the values of q+G
  ! the spherical harmonics


  ! work space
  complex(DP) :: qg(3), qgm, aux1
  real(DP) :: qmod, zero(3)

  if (.not.okvan) return
  call start_clock ('compute_intq')

  intq (:,:,:) = (0.d0, 0.0d0)
  allocate (ylmk0(1 , lmaxq * lmaxq))
  !
  !    first compute the spherical harmonics
  !
  zero=0.0_DP
  call setqmod (1, xq, zero, qmod, qg)
  call ylmr2 (lmaxq * lmaxq, 1, qg, qmod, ylmk0)
  qmod = sqrt (qmod  )

  do nt = 1, ntyp
     if (upf(nt)%tvanp ) then
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
              call qvan2 (1, ih, jh, nt, qmod, qgm, ylmk0)
              do na = 1, nat
                 if (ityp (na) == nt) then
                    aux1 = qgm * eigqts(na)
                    intq(ih,jh,na) = omega * CONJG(aux1)
                 endif
              enddo
           enddo
        enddo
        do na = 1, nat
           if (ityp(na) == nt) then
              !
              !    We use the symmetry properties of the ps factor
              !
              do ih = 1, nh (nt)
                 do jh = ih, nh (nt)
                    intq(jh,ih,na) = intq(ih,jh,na)
                 enddo
              enddo
           endif
        enddo
     endif
  enddo

  IF (noncolin) CALL set_intq_nc()

  deallocate (ylmk0)

  call stop_clock ('compute_intq')
  return
end subroutine compute_intq
