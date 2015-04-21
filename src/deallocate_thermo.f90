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
  USE thermo_mod,     ONLY : alat_geo, energy_geo, omega_geo
  USE thermodynamics, ONLY : ph_free_ener, ph_ener, ph_entropy, ph_cv
  USE ph_freq_thermodynamics, ONLY : phf_free_ener, phf_ener, phf_entropy, phf_cv
  USE anharmonic,     ONLY : vmin_t, b0_t, free_e_min_t, &
                             alpha_t, beta_t, gamma_t, cv_t, cp_t, b0_s
  USE ph_freq_anharmonic, ONLY : vminf_t, b0f_t, free_e_minf_t, &
                             alphaf_t, betaf_t, gammaf_t, cvf_t, cpf_t, b0f_s
  USE grun_anharmonic,  ONLY : betab, grun_gamma_t
  USE control_paths,    ONLY : xqaux, wqaux, letter, label_list, letter_path, &
                               label_disp_q, disp_q, disp_wq, nrap_plot_in,   &
                               rap_plot_in, nrap_plot, rap_plot
  USE control_2d_bands, ONLY : averag, vacuum, aux_ind_sur

  USE control_conv,     ONLY : ke, keden, nk_test, sigma_test
  USE elastic_constants, ONLY : epsilon_geo, sigma_geo, epsilon_voigt

  IMPLICIT NONE
  !
  IF ( ALLOCATED (alat_geo) )        DEALLOCATE(alat_geo)
  IF ( ALLOCATED (energy_geo) )      DEALLOCATE(energy_geo)
  IF ( ALLOCATED (omega_geo) )       DEALLOCATE(omega_geo) 
  IF ( ALLOCATED (vmin_t) )          DEALLOCATE(vmin_t) 
  IF ( ALLOCATED (ke) )              DEALLOCATE(ke) 
  IF ( ALLOCATED (keden) )           DEALLOCATE(keden) 
  IF ( ALLOCATED (nk_test) )         DEALLOCATE(nk_test) 
  IF ( ALLOCATED (sigma_test) )      DEALLOCATE(sigma_test) 

  IF ( ALLOCATED(ph_free_ener) )     DEALLOCATE(ph_free_ener)
  IF ( ALLOCATED(ph_ener) )          DEALLOCATE(ph_ener)
  IF ( ALLOCATED(ph_entropy) )       DEALLOCATE(ph_entropy)
  IF ( ALLOCATED(ph_cv) )            DEALLOCATE(ph_cv)
  IF ( ALLOCATED(phf_free_ener) )    DEALLOCATE(phf_free_ener)
  IF ( ALLOCATED(phf_ener) )         DEALLOCATE(phf_ener)
  IF ( ALLOCATED(phf_entropy) )      DEALLOCATE(phf_entropy)
  IF ( ALLOCATED(phf_cv) )           DEALLOCATE(phf_cv)

  IF ( ALLOCATED (b0_t) )            DEALLOCATE(b0_t) 
  IF ( ALLOCATED (b0_s) )            DEALLOCATE(b0_s) 
  IF ( ALLOCATED (cv_t) )            DEALLOCATE(cv_t) 
  IF ( ALLOCATED (cp_t) )            DEALLOCATE(cp_t) 
  IF ( ALLOCATED (alpha_t) )         DEALLOCATE(alpha_t) 
  IF ( ALLOCATED (beta_t) )          DEALLOCATE(beta_t) 
  IF ( ALLOCATED (gamma_t) )         DEALLOCATE(gamma_t) 
  IF ( ALLOCATED (free_e_min_t) )    DEALLOCATE(free_e_min_t) 
  IF ( ALLOCATED (b0f_t) )           DEALLOCATE(b0f_t) 
  IF ( ALLOCATED (b0f_s) )           DEALLOCATE(b0f_s) 
  IF ( ALLOCATED (cvf_t) )           DEALLOCATE(cvf_t) 
  IF ( ALLOCATED (cpf_t) )           DEALLOCATE(cpf_t) 
  IF ( ALLOCATED (alphaf_t) )        DEALLOCATE(alphaf_t) 
  IF ( ALLOCATED (betaf_t) )         DEALLOCATE(betaf_t) 
  IF ( ALLOCATED (gammaf_t) )        DEALLOCATE(gammaf_t) 
  IF ( ALLOCATED (betab) )           DEALLOCATE(betab) 
  IF ( ALLOCATED (grun_gamma_t) )    DEALLOCATE(grun_gamma_t) 
  IF ( ALLOCATED (free_e_minf_t) )   DEALLOCATE(free_e_minf_t) 
  IF ( ALLOCATED (xqaux) )           DEALLOCATE(xqaux)
  IF ( ALLOCATED (wqaux) )           DEALLOCATE(wqaux)
  IF ( ALLOCATED (letter) )          DEALLOCATE(letter)
  IF ( ALLOCATED (letter_path) )     DEALLOCATE(letter_path)
  IF ( ALLOCATED (label_list) )      DEALLOCATE(label_list)
  IF ( ALLOCATED (label_disp_q) )    DEALLOCATE(label_disp_q)
  IF ( ALLOCATED (disp_q) )          DEALLOCATE(disp_q)
  IF ( ALLOCATED (disp_wq) )         DEALLOCATE(disp_wq)

  IF ( ALLOCATED (nrap_plot_in) )    DEALLOCATE(nrap_plot_in)
  IF ( ALLOCATED (rap_plot_in) )     DEALLOCATE(rap_plot_in)
  IF ( ALLOCATED (nrap_plot) )       DEALLOCATE(nrap_plot)
  IF ( ALLOCATED (rap_plot) )        DEALLOCATE(rap_plot)

  IF ( ALLOCATED (averag) )          DEALLOCATE(averag)
  IF ( ALLOCATED (vacuum) )          DEALLOCATE(vacuum)
  IF ( ALLOCATED (aux_ind_sur) )     DEALLOCATE(aux_ind_sur)

  IF ( ALLOCATED (epsilon_voigt) )   DEALLOCATE(epsilon_voigt)
  IF ( ALLOCATED (epsilon_geo) )     DEALLOCATE(epsilon_geo)
  IF ( ALLOCATED (sigma_geo) )       DEALLOCATE(sigma_geo)
  ! 
  RETURN
  !
END SUBROUTINE deallocate_thermo
