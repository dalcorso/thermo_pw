!
! Copyright (C) 2017 Dal Corso Andrea
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE allocate_phq_tpw
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: quantities needed for the linear
  ! response problem
  !  
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat
  USE zstar_add,        ONLY : zstareu0_rec 
  USE magnetic_charges, ONLY : mag_charge_mode, mag_charge, alpha_me
  USE qpoint_aux,       ONLY : becpt, alphapt
  USE becmod,           ONLY : allocate_bec_type
  USE noncollin_module, ONLY : noncolin, npol, domag
  USE lsda_mod,         ONLY : nspin
  USE qpoint,           ONLY : nksq
  USE wvfct,            ONLY : nbnd
  USE uspp,             ONLY : nkb, okvan
  USE uspp_param,       ONLY : nhm
  USE control_qe,       ONLY : this_pcxpsi_is_on_file_tpw
  USE control_lr,       ONLY : lgamma
  USE control_ph,       ONLY : epsil

  IMPLICIT NONE

  INTEGER :: ik, ipol

  ALLOCATE (zstareu0_rec (3, 3 * nat))

  IF (noncolin.AND.domag) THEN
     IF (lgamma) THEN 
        ALLOCATE (mag_charge_mode(3*nat,3))
        ALLOCATE (mag_charge(3, nat, 3))
        IF (epsil) ALLOCATE (alpha_me(3,3))
     END IF
     ALLOCATE (this_pcxpsi_is_on_file_tpw(nksq,3,2))
  ELSE
     ALLOCATE (this_pcxpsi_is_on_file_tpw(nksq,3,1))
  ENDIF
  this_pcxpsi_is_on_file_tpw=.FALSE.

  RETURN
END SUBROUTINE allocate_phq_tpw

!-----------------------------------------------------------------------
SUBROUTINE deallocate_phq_tpw()
  !-----------------------------------------------------------------------

USE zstar_add,        ONLY : zstareu0_rec
USE qpoint,           ONLY : nksq
USE becmod,           ONLY : deallocate_bec_type
USE control_qe,       ONLY : this_pcxpsi_is_on_file_tpw
USE magnetic_charges, ONLY : mag_charge_mode, mag_charge, alpha_me

IMPLICIT NONE

INTEGER :: ik, ipol

IF (ALLOCATED(zstareu0_rec)) DEALLOCATE(zstareu0_rec)
IF (ALLOCATED(mag_charge_mode)) DEALLOCATE(mag_charge_mode)
IF (ALLOCATED(mag_charge)) DEALLOCATE(mag_charge)
IF (ALLOCATED(alpha_me)) DEALLOCATE(alpha_me)

IF (ALLOCATED(this_pcxpsi_is_on_file_tpw)) DEALLOCATE(this_pcxpsi_is_on_file_tpw)

RETURN
END SUBROUTINE deallocate_phq_tpw
