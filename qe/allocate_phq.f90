!
! Copyright (C) 2017 Dal Corso Andrea
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine allocate_phq_tpw
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: quantities needed for the linear
  ! response problem
  !

  USE kinds, only : DP
  USE ions_base, ONLY : nat
  USE zstar_add, ONLY : zstareu0_rec 

  implicit none

  ALLOCATE (zstareu0_rec (3, 3 * nat))

  RETURN
END SUBROUTINE allocate_phq_tpw

!-----------------------------------------------------------------------
SUBROUTINE deallocate_phq_tpw()
  !-----------------------------------------------------------------------

USE zstar_add,    ONLY : zstareu0_rec
IMPLICIT NONE

IF (ALLOCATED(zstareu0_rec)) DEALLOCATE(zstareu0_rec)

RETURN
END SUBROUTINE deallocate_phq_tpw
