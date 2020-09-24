!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------
SUBROUTINE ccg_psi_tpw (lda, n, m, psi, h_diag, flag)
  !-----------------------------------------------------------------
  !
  !    This routine gives a preconditioning to the linear system solver.
  !    The preconditioning is diagonal in reciprocal space
  !
  !
  USE kinds, only : DP
  USE noncollin_module, ONLY : noncolin, npol
  implicit none

  integer :: lda, n, m, flag
  ! input: the leading dimension of the psi vector
  ! input: the real dimension of the vector
  ! input: the number of vectors
  ! input: flag=1 use h_diag, flag=-1 use conjg(h_diag)
  complex(kind=DP) :: psi (lda, m)
  ! inp/out: the vector to be preconditioned

  complex(kind=DP) :: h_diag (lda, m)
  ! input: the preconditioning vector

  integer :: k, i, shift
  ! counter on bands
  ! counter on the elements of the vector
  !
  do k = 1, m
     do i = 1, n
       if (flag .eq. 1) then
         psi (i, k) = psi (i, k) * h_diag (i, k)
       else if (flag .eq. -1) then
         psi (i, k) = psi (i, k) * CONJG(h_diag (i, k))
       else
         print*, 'flag is neither 1 nor -1. Stop'
       endif  
     enddo
     IF (noncolin) THEN
        shift=lda/npol
        do i = 1, n
           if (flag .eq. 1) then
              psi (i+shift, k) = psi (i+shift, k) * h_diag (i+shift, k)
           else if (flag .eq. -1) then
              psi (i+shift, k) = psi (i+shift, k) * CONJG(h_diag (i+shift, k))
           else
              print*, 'flag is neither 1 nor -1. Stop'
           endif  
        end do
     END IF
  enddo
  return
END SUBROUTINE ccg_psi_tpw
