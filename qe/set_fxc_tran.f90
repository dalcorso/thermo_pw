!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_fxc_tran( )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes 1/|m| dE_xc / d |m| in the lsda case and
  !     set it to dmuxc_tran
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE lsda_mod,         ONLY : nspin
  USE scf,              ONLY : scf_type, rho, rho_core, rhog_core
  USE optical,          ONLY : dmuxc_tran
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: ir
    ! counter on mesh points
    ! counter on nspin
  !
  REAL(DP) :: etxc, vtxc
  REAL(DP), ALLOCATABLE :: v(:,:)
  !
  CALL start_clock( 'set_fxc_tran' )
  !
  etxc=0.0_DP
  vtxc=0.0_DP
  dmuxc_tran(:) = 0.0_DP
  ALLOCATE(v(dfftp%nnr,nspin))
  v=0.0_DP

  CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, v )
  DO ir = 1, dfftp%nnr
     IF (ABS(rho%of_r(ir,1) - rho%of_r(ir,2)) > 1.D-11 ) THEN
        dmuxc_tran(ir) = 0.5_DP*( v(ir, 1) - v(ir,2) ) / &
                         (rho%of_r(ir,1) - rho%of_r(ir,2))
     ENDIF
  ENDDO

  DEALLOCATE(v)

  CALL stop_clock( 'set_fxc_tran' )
  !
  RETURN
  !
END SUBROUTINE set_fxc_tran
