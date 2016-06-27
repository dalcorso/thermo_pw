!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE dv_of_drho_tran (dvscf)
  !-----------------------------------------------------------------------
  !
  !   This routine computes the change of the self consistent potential
  !   due to a transverse magnetic field perturbation. In input it
  !   receives m+ or m- (in dvscf) and the term 1/|m| d B_xc / d |m| 
  !   in the array dmuxc_tran
  !
  USE kinds,      ONLY : DP
  USE fft_base,   ONLY : dfftp
  USE optical,    ONLY : dmuxc_tran

  IMPLICIT NONE

  COMPLEX(DP), INTENT(INOUT) :: dvscf (dfftp%nnr)
  ! input: the change of the charge,
  ! output: change of the potential

  INTEGER :: ir
  ! counter on r vectors

  CALL start_clock ('dv_of_drho_tran')
!
!  In input the m+ or m- components of the exchange and
!  correlation magnetization. In output the exchange and correlation fields.
!
!
! the exchange-correlation contribution is computed in real space
! Note that gradient correction, if any is already contained in dmuxc_tran
! Probably not yet working, there should be a term due to the gradient
! of dvscf
!
  DO ir = 1, dfftp%nnr
     dvscf(ir) = dmuxc_tran(ir) * dvscf(ir)
  ENDDO

  CALL stop_clock ('dv_of_drho_tran')
  RETURN
END SUBROUTINE dv_of_drho_tran
