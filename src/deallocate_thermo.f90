!
! Copyright (C) 2013 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE deallocate_thermo()
  !-----------------------------------------------------------------------
  !
  !  This routine deallocates the variables that control the thermo calculation
  !
  USE kinds,          ONLY : DP
  USE thermo_mod,     ONLY : alat_geo, energy_geo, vmin_t, b0_t, free_e_min_t
  USE control_paths,  ONLY : xqaux, wqaux, letter, label_list, letter_path, &
                             label_disp_q
  USE control_conv,   ONLY : ke, keden, nk_test, sigma_test
  USE anharmonic,     ONLY : alpha_t, beta_t, gamma_t, cv_t, cp_t, b0_s
  USE thermodynamics, ONLY : omegav
  USE ifc,            ONLY : disp_q, disp_wq

  IMPLICIT NONE
  !
  IF ( ALLOCATED (alat_geo) )        DEALLOCATE(alat_geo)
  IF ( ALLOCATED (energy_geo) )      DEALLOCATE(energy_geo)
  IF ( ALLOCATED (omegav) )          DEALLOCATE(omegav) 
  IF ( ALLOCATED (vmin_t) )          DEALLOCATE(vmin_t) 
  IF ( ALLOCATED (ke) )              DEALLOCATE(ke) 
  IF ( ALLOCATED (keden) )           DEALLOCATE(keden) 
  IF ( ALLOCATED (nk_test) )         DEALLOCATE(nk_test) 
  IF ( ALLOCATED (sigma_test) )      DEALLOCATE(sigma_test) 
  IF ( ALLOCATED (b0_t) )            DEALLOCATE(b0_t) 
  IF ( ALLOCATED (b0_s) )            DEALLOCATE(b0_s) 
  IF ( ALLOCATED (cv_t) )            DEALLOCATE(cv_t) 
  IF ( ALLOCATED (cp_t) )            DEALLOCATE(cp_t) 
  IF ( ALLOCATED (alpha_t) )         DEALLOCATE(alpha_t) 
  IF ( ALLOCATED (beta_t) )          DEALLOCATE(beta_t) 
  IF ( ALLOCATED (gamma_t) )         DEALLOCATE(gamma_t) 
  IF ( ALLOCATED (free_e_min_t) )    DEALLOCATE(free_e_min_t) 
  IF ( ALLOCATED (xqaux) )           DEALLOCATE(xqaux)
  IF ( ALLOCATED (wqaux) )           DEALLOCATE(wqaux)
  IF ( ALLOCATED (letter) )          DEALLOCATE(letter)
  IF ( ALLOCATED (letter_path) )     DEALLOCATE(letter_path)
  IF ( ALLOCATED (label_list) )      DEALLOCATE(label_list)
  IF ( ALLOCATED (label_disp_q) )    DEALLOCATE(label_disp_q)
  IF ( ALLOCATED (disp_q) )          DEALLOCATE(disp_q)
  IF ( ALLOCATED (disp_wq) )         DEALLOCATE(disp_wq)
  ! 
  RETURN
  !
END SUBROUTINE deallocate_thermo
