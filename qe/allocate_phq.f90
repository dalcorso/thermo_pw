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
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat
  USE zstar_add,        ONLY : zstareu0_rec 
  USE magnetic_charges, ONLY : mag_charge_mode, mag_charge, alpha_me
  USE qpoint_aux,       ONLY : becpt, alphapt
  USE becmod,           ONLY : allocate_bec_type
  USE noncollin_module, ONLY : noncolin, npol
  USE spin_orb,         ONLY : domag
  USE lsda_mod,         ONLY : nspin
  USE qpoint,           ONLY : nksq
  USE wvfct,            ONLY : nbnd
  USE uspp,             ONLY : nkb, okvan
  USE uspp_param,       ONLY : nhm
  USE nc_mag_aux,       ONLY : int1_nc_save, deeq_nc_save
  USE control_qe,       ONLY : this_pcxpsi_is_on_file_tpw
  USE control_lr,       ONLY : lgamma
  USE control_ph,       ONLY : epsil

  IMPLICIT NONE

  INTEGER :: ik, ipol

  ALLOCATE (zstareu0_rec (3, 3 * nat))

  IF (noncolin.AND.domag) THEN
     ALLOCATE (becpt(nksq))
     ALLOCATE (alphapt(3,nksq))
     IF (lgamma) THEN 
        ALLOCATE (mag_charge_mode(3*nat,3))
        ALLOCATE (mag_charge(3, nat, 3))
        IF (epsil) ALLOCATE (alpha_me(3,3))
     END IF
     DO ik=1,nksq
        CALL allocate_bec_type ( nkb, nbnd, becpt(ik) )
        DO ipol=1,3
           CALL allocate_bec_type ( nkb, nbnd, alphapt(ipol,ik) )
        ENDDO
     ENDDO
     IF (okvan) THEN
        ALLOCATE (int1_nc_save( nhm, nhm, 3, nat, nspin, 2))
        ALLOCATE (deeq_nc_save( nhm, nhm, nat, nspin, 2))
     ENDIF
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
USE qpoint_aux,       ONLY : ikmks, ikmkmqs, becpt, alphapt
USE qpoint,           ONLY : nksq
USE becmod,           ONLY : deallocate_bec_type
USE nc_mag_aux,       ONLY : int1_nc_save, deeq_nc_save
USE control_qe,       ONLY : this_pcxpsi_is_on_file_tpw
USE magnetic_charges, ONLY : mag_charge_mode, mag_charge, alpha_me

IMPLICIT NONE

INTEGER :: ik, ipol

IF (ALLOCATED(zstareu0_rec)) DEALLOCATE(zstareu0_rec)
IF (ALLOCATED(mag_charge_mode)) DEALLOCATE(mag_charge_mode)
IF (ALLOCATED(mag_charge)) DEALLOCATE(mag_charge)
IF (ALLOCATED(alpha_me)) DEALLOCATE(alpha_me)
IF (ALLOCATED(ikmks)) DEALLOCATE(ikmks)
IF (ALLOCATED(ikmkmqs)) DEALLOCATE(ikmkmqs)

IF (ALLOCATED(alphapt)) THEN
   DO ik=1,nksq
      DO ipol=1,3
         CALL deallocate_bec_type ( alphapt(ipol,ik) )
      ENDDO
   ENDDO
   DEALLOCATE (alphapt)
ENDIF
IF (ALLOCATED(int1_nc_save)) DEALLOCATE (int1_nc_save)
IF (ALLOCATED(deeq_nc_save)) DEALLOCATE (deeq_nc_save)

IF (ALLOCATED(becpt))  THEN
   DO ik=1, nksq
      CALL deallocate_bec_type ( becpt(ik) )
   ENDDO
   DEALLOCATE(becpt)
ENDIF
IF (ALLOCATED(this_pcxpsi_is_on_file_tpw)) DEALLOCATE(this_pcxpsi_is_on_file_tpw)

RETURN
END SUBROUTINE deallocate_phq_tpw
