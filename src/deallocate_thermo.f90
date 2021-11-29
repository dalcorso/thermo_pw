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
  USE thermo_mod,     ONLY : celldm_geo, energy_geo, omega_geo, ibrav_geo, &
                             no_ph, tot_ngeo, in_degree, dynmat_on_file,   &
                             ef_geo
  USE thermodynamics, ONLY : ph_free_ener, ph_ener, ph_entropy, ph_ce, &
                             ph_e0, ph_b_fact
  USE ph_freq_thermodynamics, ONLY : phf_free_ener, phf_ener, phf_entropy, &
                             phf_b_fact, phf_e0, phf_ce
  USE el_thermodynamics, ONLY : el_ener, el_free_ener, el_entr, el_mu, &
                           el_ce
  USE anharmonic,     ONLY : vmin_t, b0_t, b01_t, b02_t, free_e_min_t, a_t,      &
                             alpha_t, beta_t, gamma_t, cv_t, ce_t, cp_t,  &
                             b0_s, ener_t, free_ener_t, entropy_t,        &
                             celldm_t, alpha_anis_t, cpmce_anis, el_cons_t, &
                             el_comp_t, macro_el_t, bths_t, ggamma_t,     &
                             el_cons_s, el_comp_s, macro_el_s, v_t, v_s,  &
                             bfact_t, el_con_geo_t, betap,                &
                             vmin_pt, b0_pt, b01_pt, b02_pt, emin_pt,     &
                             beta_pt, ce_pt, &
                             cp_pt, gamma_pt, b0_s_pt, vmin_ptt, b0_ptt,  &
                             b01_ptt, b02_ptt, gamma_ptt, ce_ptt
  USE el_anharmonic, ONLY :  vmine_t, b0e_t, b01e_t, b02e_t, free_e_mine_t
  USE ph_freq_anharmonic, ONLY : vminf_t, b0f_t, b01f_t, b02f_t, free_e_minf_t, &
                             alphaf_t, betaf_t, gammaf_t, cvf_t, cef_t, &
                             cpf_t, b0f_s, enerf_t, free_enerf_t, entropyf_t, &
                             celldmf_t, alphaf_anis_t, cpmcef_anis, &
                             el_consf_t, el_compf_t, macro_elf_t, bthsf_t, &
                             ggammaf_t, el_consf_s, el_compf_s, macro_elf_s, &
                             vf_t, vf_s, bfactf_t, el_conf_geo_t
  USE grun_anharmonic,  ONLY : betab, alpha_an_g, cp_grun_t, cv_grun_t, &
                             ce_grun_t, b0_grun_s, &
                             grun_gamma_t, poly_grun, poly_grun_red, &
                             grun_cpmce_anis, el_cons_grun_t, el_comp_grun_t
  USE control_paths,    ONLY : xqaux, wqaux, letter, label_list, letter_path, &
                               label_disp_q, disp_q, disp_wq, nrap_plot_in,   &
                               rap_plot_in, nrap_plot, rap_plot, high_sym_path
  USE control_2d_bands, ONLY : averag, vacuum, aux_ind_sur
  USE initial_conf,     ONLY : ityp_save, tau_save, tau_save_crys, &
                               collect_info_save
  USE equilibrium_conf, ONLY : tau0, tau0_crys
  USE control_thermo,   ONLY : all_geometries_together
  USE control_grun,     ONLY : vgrun_t, b0_grun_t, celldm_grun_t
  USE temperature,      ONLY : temp, temp_plot, itemp_plot

  USE control_conv,     ONLY : ke, keden, nk_test, sigma_test
  USE control_mur_p,    ONLY : vmin_p, b0_p, b01_p, b02_p, emin_p
  USE control_mur,      ONLY : p0
  USE control_ev,       ONLY : e0, v0
  USE elastic_constants, ONLY : epsilon_geo, sigma_geo, epsilon_voigt
  USE control_elastic_constants, ONLY : rot_mat, el_con_geo, &
                             el_con_ibrav_geo, el_con_celldm_geo, &
                             el_con_tau_crys_geo, el_con_omega_geo, &
                             epsil_geo
  USE control_pressure, ONLY : press_plot, ipress_plot
  USE control_vol,      ONLY : ivol_plot
  USE control_debye,    ONLY : deb_energy, deb_free_energy, deb_entropy, &
                               deb_cv, deb_b_fact, deb_bfact
  USE control_quadratic_energy, ONLY : p2, hessian_v, hessian_e, x_pos_min
  USE control_quartic_energy, ONLY :  p4, x_min_4
  USE el_anharmonic, ONLY : el_energy_t, el_free_energy_t, el_entropy_t, &
                            el_ce_t, el_energyf_t, el_free_energyf_t,    &
                            el_entropyf_t, el_cef_t
  USE collect_info,  ONLY : destroy_collect_info_type
  USE control_eldos, ONLY : dos_k, dos_wk
  USE polynomial,    ONLY : clean_poly

  IMPLICIT NONE
  INTEGER :: igeom
  !
  IF ( ALLOCATED (energy_geo) )      DEALLOCATE(energy_geo)
  IF ( ALLOCATED (ef_geo) )          DEALLOCATE(ef_geo)
  IF ( ALLOCATED (celldm_geo) )      DEALLOCATE(celldm_geo) 
  IF ( ALLOCATED (omega_geo) )       DEALLOCATE(omega_geo) 
  IF ( ALLOCATED (ibrav_geo) )       DEALLOCATE(ibrav_geo) 
  IF ( ALLOCATED (no_ph) )           DEALLOCATE(no_ph) 
  IF ( ALLOCATED (in_degree) )       DEALLOCATE(in_degree) 
  IF ( ALLOCATED (dynmat_on_file) )  DEALLOCATE(dynmat_on_file) 
  IF ( ALLOCATED (vmin_t) )          DEALLOCATE(vmin_t) 
  IF ( ALLOCATED (ke) )              DEALLOCATE(ke) 
  IF ( ALLOCATED (keden) )           DEALLOCATE(keden) 
  IF ( ALLOCATED (nk_test) )         DEALLOCATE(nk_test) 
  IF ( ALLOCATED (sigma_test) )      DEALLOCATE(sigma_test) 

  IF ( ALLOCATED(ph_free_ener) )     DEALLOCATE(ph_free_ener)
  IF ( ALLOCATED(ph_ener) )          DEALLOCATE(ph_ener)
  IF ( ALLOCATED(ph_entropy) )       DEALLOCATE(ph_entropy)
  IF ( ALLOCATED(ph_e0) )            DEALLOCATE(ph_e0)
  IF ( ALLOCATED(ph_ce) )            DEALLOCATE(ph_ce)
  IF ( ALLOCATED(ph_b_fact) )        DEALLOCATE(ph_b_fact)
  IF ( ALLOCATED(phf_free_ener) )    DEALLOCATE(phf_free_ener)
  IF ( ALLOCATED(phf_ener) )         DEALLOCATE(phf_ener)
  IF ( ALLOCATED(phf_entropy) )      DEALLOCATE(phf_entropy)
  IF ( ALLOCATED(phf_e0) )           DEALLOCATE(phf_e0)
  IF ( ALLOCATED(phf_ce) )           DEALLOCATE(phf_ce)
  IF ( ALLOCATED(phf_b_fact) )       DEALLOCATE(phf_b_fact)

  IF (ALLOCATED(el_free_ener))       DEALLOCATE(el_free_ener)
  IF (ALLOCATED(el_ener))            DEALLOCATE(el_ener)
  IF (ALLOCATED(el_entr))            DEALLOCATE(el_entr)
  IF (ALLOCATED(el_mu))              DEALLOCATE(el_mu)
  IF (ALLOCATED(el_ce))              DEALLOCATE(el_ce)


  IF ( ALLOCATED (b0_t) )            DEALLOCATE(b0_t) 
  IF ( ALLOCATED (b01_t) )           DEALLOCATE(b01_t) 
  IF ( ALLOCATED (b02_t) )           DEALLOCATE(b02_t) 
  IF ( ALLOCATED (b0_s) )            DEALLOCATE(b0_s) 
  IF ( ALLOCATED (cv_t) )            DEALLOCATE(cv_t) 
  IF ( ALLOCATED (ce_t) )            DEALLOCATE(ce_t) 
  IF ( ALLOCATED (ener_t) )          DEALLOCATE(ener_t) 
  IF ( ALLOCATED (free_ener_t) )     DEALLOCATE(free_ener_t)
  IF ( ALLOCATED (a_t) )             DEALLOCATE(a_t)
  IF ( ALLOCATED (entropy_t) )       DEALLOCATE(entropy_t) 
  IF ( ALLOCATED (cp_t) )            DEALLOCATE(cp_t) 
  IF ( ALLOCATED (alpha_t) )         DEALLOCATE(alpha_t) 
  IF ( ALLOCATED (beta_t) )          DEALLOCATE(beta_t) 
  IF ( ALLOCATED (gamma_t) )         DEALLOCATE(gamma_t) 
  IF ( ALLOCATED (bths_t) )          DEALLOCATE(bths_t) 
  IF ( ALLOCATED (ggamma_t) )        DEALLOCATE(ggamma_t) 
  IF ( ALLOCATED (bfact_t) )         DEALLOCATE(bfact_t) 
  IF ( ALLOCATED (celldm_t) )        DEALLOCATE(celldm_t) 
  IF ( ALLOCATED (alpha_anis_t) )    DEALLOCATE(alpha_anis_t) 
  IF ( ALLOCATED (cpmce_anis) )      DEALLOCATE(cpmce_anis) 
  IF ( ALLOCATED (free_e_min_t) )    DEALLOCATE(free_e_min_t) 
  IF ( ALLOCATED (el_cons_t) )       DEALLOCATE(el_cons_t)
  IF ( ALLOCATED (el_comp_t) )       DEALLOCATE(el_comp_t)
  IF ( ALLOCATED (el_cons_s) )       DEALLOCATE(el_cons_s)
  IF ( ALLOCATED (el_comp_s) )       DEALLOCATE(el_comp_s)
  IF ( ALLOCATED (macro_el_t) )      DEALLOCATE(macro_el_t)
  IF ( ALLOCATED (macro_el_s) )      DEALLOCATE(macro_el_s)
  IF ( ALLOCATED (v_t) )             DEALLOCATE(v_t)
  IF ( ALLOCATED (v_s) )             DEALLOCATE(v_s)
  IF ( ALLOCATED (el_con_geo_t) )    DEALLOCATE(el_con_geo_t) 

  IF ( ALLOCATED (vmin_p) )          DEALLOCATE(vmin_p) 
  IF ( ALLOCATED (p0) )              DEALLOCATE(p0) 
  IF ( ALLOCATED (b0_p) )            DEALLOCATE(b0_p) 
  IF ( ALLOCATED (b01_p) )           DEALLOCATE(b01_p) 
  IF ( ALLOCATED (b02_p) )           DEALLOCATE(b02_p) 
  IF ( ALLOCATED (emin_p) )          DEALLOCATE(emin_p) 
  IF ( ALLOCATED (betap) )           DEALLOCATE(betap) 
  IF ( ALLOCATED (vmin_ptt) )        DEALLOCATE(vmin_ptt) 
  IF ( ALLOCATED (b0_ptt) )          DEALLOCATE(b0_ptt) 
  IF ( ALLOCATED (b01_ptt) )         DEALLOCATE(b01_ptt) 
  IF ( ALLOCATED (b02_ptt) )         DEALLOCATE(b02_ptt) 
  IF ( ALLOCATED (ce_ptt) )          DEALLOCATE(ce_ptt) 
  IF ( ALLOCATED (gamma_ptt) )       DEALLOCATE(gamma_ptt) 
  IF ( ALLOCATED (vmin_pt) )         DEALLOCATE(vmin_pt) 
  IF ( ALLOCATED (b0_pt) )           DEALLOCATE(b0_pt) 
  IF ( ALLOCATED (b01_pt) )          DEALLOCATE(b01_pt) 
  IF ( ALLOCATED (b02_pt) )          DEALLOCATE(b02_pt) 
  IF ( ALLOCATED (emin_pt) )         DEALLOCATE(emin_pt) 
  IF ( ALLOCATED (beta_pt) )         DEALLOCATE(beta_pt) 
  IF ( ALLOCATED (ce_pt) )           DEALLOCATE(ce_pt) 
  IF ( ALLOCATED (cp_pt) )           DEALLOCATE(cp_pt) 
  IF ( ALLOCATED (gamma_pt) )        DEALLOCATE(gamma_pt) 
  IF ( ALLOCATED (b0_s_pt) )         DEALLOCATE(b0_s_pt) 
  IF ( ALLOCATED (e0) )              DEALLOCATE(e0)
  IF ( ALLOCATED (v0) )              DEALLOCATE(v0)

  IF ( ALLOCATED (b0f_t) )           DEALLOCATE(b0f_t) 
  IF ( ALLOCATED (b01f_t) )          DEALLOCATE(b01f_t) 
  IF ( ALLOCATED (b02f_t) )          DEALLOCATE(b02f_t) 
  IF ( ALLOCATED (b0f_s) )           DEALLOCATE(b0f_s) 
  IF ( ALLOCATED (cvf_t) )           DEALLOCATE(cvf_t) 
  IF ( ALLOCATED (enerf_t) )         DEALLOCATE(enerf_t)
  IF ( ALLOCATED (free_enerf_t) )    DEALLOCATE(free_enerf_t) 
  IF ( ALLOCATED (entropyf_t) )      DEALLOCATE(entropyf_t) 
  IF ( ALLOCATED (cef_t) )           DEALLOCATE(cef_t) 
  IF ( ALLOCATED (cpf_t) )           DEALLOCATE(cpf_t) 
  IF ( ALLOCATED (alphaf_t) )        DEALLOCATE(alphaf_t) 
  IF ( ALLOCATED (betaf_t) )         DEALLOCATE(betaf_t) 
  IF ( ALLOCATED (gammaf_t) )        DEALLOCATE(gammaf_t) 
  IF ( ALLOCATED (bthsf_t) )         DEALLOCATE(bthsf_t) 
  IF ( ALLOCATED (ggammaf_t) )       DEALLOCATE(ggammaf_t) 
  IF ( ALLOCATED (bfactf_t) )        DEALLOCATE(bfactf_t) 
  IF ( ALLOCATED (celldmf_t) )       DEALLOCATE(celldmf_t) 
  IF ( ALLOCATED (alphaf_anis_t) )   DEALLOCATE(alphaf_anis_t) 
  IF ( ALLOCATED (cpmcef_anis) )     DEALLOCATE(cpmcef_anis) 
  IF ( ALLOCATED (free_e_minf_t) )   DEALLOCATE(free_e_minf_t) 
  IF ( ALLOCATED (el_consf_t) )      DEALLOCATE(el_consf_t)
  IF ( ALLOCATED (el_compf_t) )      DEALLOCATE(el_compf_t)
  IF ( ALLOCATED (el_consf_s) )      DEALLOCATE(el_consf_s)
  IF ( ALLOCATED (el_compf_s) )      DEALLOCATE(el_compf_s)
  IF ( ALLOCATED (macro_elf_s) )     DEALLOCATE(macro_elf_s)
  IF ( ALLOCATED (macro_elf_t) )     DEALLOCATE(macro_elf_t)
  IF ( ALLOCATED (vf_t) )            DEALLOCATE(vf_t)
  IF ( ALLOCATED (vf_s) )            DEALLOCATE(vf_s)
  IF ( ALLOCATED (el_conf_geo_t) )   DEALLOCATE(el_conf_geo_t) 

  IF ( ALLOCATED (vmine_t) )  DEALLOCATE(vmine_t) 
  IF ( ALLOCATED (b0e_t) )    DEALLOCATE(b0e_t) 
  IF ( ALLOCATED (b01e_t) )   DEALLOCATE(b01e_t) 
  IF ( ALLOCATED (b02e_t) )   DEALLOCATE(b02e_t) 
  IF ( ALLOCATED (free_e_mine_t) ) DEALLOCATE(free_e_mine_t) 

  IF ( ALLOCATED (betab) )           DEALLOCATE(betab) 
  IF ( ALLOCATED (alpha_an_g) )      DEALLOCATE(alpha_an_g) 
  IF ( ALLOCATED (poly_grun) )       DEALLOCATE(poly_grun) 
  IF ( ALLOCATED (poly_grun_red) )   DEALLOCATE(poly_grun_red) 
  IF ( ALLOCATED (vgrun_t) )         DEALLOCATE(vgrun_t) 
  IF ( ALLOCATED (b0_grun_t) )       DEALLOCATE(b0_grun_t) 
  IF ( ALLOCATED (el_cons_grun_t) )  DEALLOCATE(el_cons_grun_t) 
  IF ( ALLOCATED (el_comp_grun_t) )  DEALLOCATE(el_comp_grun_t) 
  IF ( ALLOCATED (celldm_grun_t) )   DEALLOCATE(celldm_grun_t) 
  IF ( ALLOCATED (cp_grun_t) )       DEALLOCATE(cp_grun_t) 
  IF ( ALLOCATED (cv_grun_t) )       DEALLOCATE(cv_grun_t) 
  IF ( ALLOCATED (ce_grun_t) )       DEALLOCATE(ce_grun_t) 
  IF ( ALLOCATED (b0_grun_s) )       DEALLOCATE(b0_grun_s) 
  IF ( ALLOCATED (grun_cpmce_anis) ) DEALLOCATE(grun_cpmce_anis) 
  IF ( ALLOCATED (grun_gamma_t) )    DEALLOCATE(grun_gamma_t) 

  IF (ALLOCATED (el_energy_t) )      DEALLOCATE(el_energy_t)
  IF (ALLOCATED (el_free_energy_t) ) DEALLOCATE(el_free_energy_t)
  IF (ALLOCATED (el_entropy_t) )     DEALLOCATE(el_entropy_t)
  IF (ALLOCATED (el_ce_t) )          DEALLOCATE(el_ce_t)

  IF (ALLOCATED (el_energyf_t) )      DEALLOCATE(el_energyf_t)
  IF (ALLOCATED (el_free_energyf_t) ) DEALLOCATE(el_free_energyf_t)
  IF (ALLOCATED (el_entropyf_t) )     DEALLOCATE(el_entropyf_t)
  IF (ALLOCATED (el_cef_t) )          DEALLOCATE(el_cef_t)

  IF ( ALLOCATED (xqaux) )           DEALLOCATE(xqaux)
  IF ( ALLOCATED (wqaux) )           DEALLOCATE(wqaux)
  IF ( ALLOCATED (letter) )          DEALLOCATE(letter)
  IF ( ALLOCATED (letter_path) )     DEALLOCATE(letter_path)
  IF ( ALLOCATED (label_list) )      DEALLOCATE(label_list)
  IF ( ALLOCATED (label_disp_q) )    DEALLOCATE(label_disp_q)
  IF ( ALLOCATED (disp_q) )          DEALLOCATE(disp_q)
  IF ( ALLOCATED (disp_wq) )         DEALLOCATE(disp_wq)

  CALL clean_poly(p2)
  IF ( ALLOCATED (hessian_v) )       DEALLOCATE( hessian_v )
  IF ( ALLOCATED (hessian_e) )       DEALLOCATE( hessian_e )
  IF ( ALLOCATED (x_pos_min) )       DEALLOCATE( x_pos_min )
  IF ( ALLOCATED (x_min_4) )         DEALLOCATE( x_min_4 )
  CALL clean_poly(p4)

  IF ( ALLOCATED (temp) )            DEALLOCATE(temp)
  IF ( ALLOCATED (temp_plot) )       DEALLOCATE(temp_plot)
  IF ( ALLOCATED (itemp_plot) )      DEALLOCATE(itemp_plot)
  IF ( ALLOCATED (press_plot) )      DEALLOCATE(press_plot)
  IF ( ALLOCATED (ipress_plot) )     DEALLOCATE(ipress_plot)
  IF ( ALLOCATED (ivol_plot) )       DEALLOCATE(ivol_plot)

  IF ( ALLOCATED (nrap_plot_in) )    DEALLOCATE(nrap_plot_in)
  IF ( ALLOCATED (rap_plot_in) )     DEALLOCATE(rap_plot_in)
  IF ( ALLOCATED (nrap_plot) )       DEALLOCATE(nrap_plot)
  IF ( ALLOCATED (rap_plot) )        DEALLOCATE(rap_plot)
  IF ( ALLOCATED (high_sym_path) )   DEALLOCATE(high_sym_path)

  IF ( ALLOCATED (averag) )          DEALLOCATE(averag)
  IF ( ALLOCATED (vacuum) )          DEALLOCATE(vacuum)
  IF ( ALLOCATED (aux_ind_sur) )     DEALLOCATE(aux_ind_sur)

  IF ( ALLOCATED (epsilon_voigt) )   DEALLOCATE(epsilon_voigt)
  IF ( ALLOCATED (epsilon_geo) )     DEALLOCATE(epsilon_geo)
  IF ( ALLOCATED (epsil_geo) )       DEALLOCATE(epsil_geo)
  IF ( ALLOCATED (sigma_geo) )       DEALLOCATE(sigma_geo)
  IF ( ALLOCATED (rot_mat) )         DEALLOCATE(rot_mat)
  IF ( ALLOCATED (el_con_geo) )      DEALLOCATE(el_con_geo)

  IF ( ALLOCATED (tau_save) )        DEALLOCATE(tau_save)
  IF ( ALLOCATED (tau_save_crys) )   DEALLOCATE(tau_save_crys)
  IF ( ALLOCATED (ityp_save) )       DEALLOCATE(ityp_save)

  IF ( ALLOCATED (tau0) )            DEALLOCATE(tau0)
  IF ( ALLOCATED (tau0_crys) )       DEALLOCATE(tau0_crys)

  IF ( ALLOCATED (deb_energy) )      DEALLOCATE( deb_energy )
  IF ( ALLOCATED (deb_free_energy) ) DEALLOCATE( deb_free_energy )
  IF ( ALLOCATED (deb_entropy) )     DEALLOCATE( deb_entropy )
  IF ( ALLOCATED (deb_cv) )          DEALLOCATE( deb_cv )
  IF ( ALLOCATED (deb_b_fact) )      DEALLOCATE( deb_b_fact )
  IF ( ALLOCATED (deb_bfact) )       DEALLOCATE( deb_bfact )

  IF ( ALLOCATED (dos_k) )           DEALLOCATE( dos_k )
  IF ( ALLOCATED (dos_wk) )          DEALLOCATE( dos_wk )

  IF ( ALLOCATED (el_con_ibrav_geo) )  DEALLOCATE( el_con_ibrav_geo )  
  IF ( ALLOCATED (el_con_celldm_geo) ) DEALLOCATE( el_con_celldm_geo )  
  IF ( ALLOCATED (el_con_tau_crys_geo) ) DEALLOCATE( el_con_tau_crys_geo )  
  IF ( ALLOCATED (el_con_omega_geo) ) DEALLOCATE( el_con_omega_geo )  

  CALL deallocate_q2r()
  IF (ALLOCATED(collect_info_save)) THEN
     IF (all_geometries_together) THEN
        DO igeom=1, tot_ngeo
           CALL destroy_collect_info_type(collect_info_save(igeom))
        ENDDO
     ELSE
        CALL destroy_collect_info_type(collect_info_save(1))
     ENDIF
     DEALLOCATE(collect_info_save)
  ENDIF

  ! 
  RETURN
  !
END SUBROUTINE deallocate_thermo
