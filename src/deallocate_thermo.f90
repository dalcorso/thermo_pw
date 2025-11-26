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
  USE thermo_mod,     ONLY : celldm_geo, energy_geo, omega_geo, ibrav_geo,  &
                             no_ph, tot_ngeo, dynmat_on_file,    &
                             ef_geo, tau_geo, no_ph_eos
  USE thermodynamics, ONLY : ph_free_ener, ph_ener, ph_entropy, ph_ce,      &
                             ph_t_debye, ph_e0, ph_b_fact, ph_e0_eos,       &
                             ph_free_ener_eos, ph_ener_eos, ph_entropy_eos, &
                             ph_ce_eos, uint_geo_eos_t, tau_geo_eos_t

  USE ph_freq_thermodynamics, ONLY : phf_free_ener, phf_ener, phf_entropy,  &
                             phf_b_fact, phf_e0, phf_ce, phf_t_debye,       &
                             phf_free_ener_eos, phf_ener_eos, phf_e0_eos,   &
                             phf_entropy_eos, phf_ce_eos, uintf_geo_eos_t,  &
                             tauf_geo_eos_t
  USE el_thermodynamics, ONLY : el_ener, el_free_ener, el_entr, el_mu,      &
                                el_ce, el_ener_eos, el_free_ener_eos,       &
                                el_entr_eos, el_mu_eos, el_ce_eos,          &
                                elf_ener_eos, elf_free_ener_eos,            &
                                elf_entr_eos, elf_mu_eos, elf_ce_eos
  USE anharmonic,     ONLY : vmin_t, b0_t, b01_t, b02_t, free_e_min_t, a_t, &
                             ener_t, free_ener_t, entropy_t, ce_t, cv_t,    &
                             cp_t, b0_s, alpha_t, beta_t, gamma_t,          &
                             celldm_t, alpha_anis_t, cpmce_anis, bfact_t,   &
                             bths_t, ggamma_t, el_cons_t, el_comp_t,        &
                             el_cons_s, el_comp_s, macro_el_t, macro_el_s,  &
                             v_t, v_s, el_con_geo_t, e_piezo_tensor_t,      &
                             d_piezo_tensor_t,                              &
                             vmin_noe_t, b0_noe_t, b01_noe_t, b02_noe_t,    &
                             free_e_min_noe_t, a_noe_t, ce_noe_t, cv_noe_t, &
                             cp_noe_t, b0_noe_s, beta_noe_t, gamma_noe_t,   &
                             p1t_t, p2t_t, p3t_t, p4t_t, density_t, csmct_t, &
                             free_ener_noe_t, ener_noe_t, entropy_noe_t,    &
                             celldm_noe_t, p1t_noe_t, p2t_noe_t,            &
                             p3t_noe_t, p4t_noe_t, density_noe_t,           &
                             celldm_t_p1, celldm_t_m1,                      &
                             celldm_noe_t_p1, celldm_noe_t_m1, b0_ec_t,     &
                             b0_ec_s, debye_macro_el_t, debye_macro_el_s,   &
                             dyde_t, uint_t, tau_t, uint_zsisa_t,           &
                             tau_zsisa_t, alpha_int_t, alpha_int_zsisa_t
  USE anharmonic_pt,  ONLY : vmin_pt, b0_pt, b01_pt, b02_pt, emin_pt,       &
                             ce_pt, cv_pt, cp_pt, b0_s_pt, beta_pt,         &
                             gamma_pt, free_ener_pt, ener_pt, entr_pt,      &
                             celldm_pt, alpha_anis_pt, density_pt,          &
                             cpmce_anis_pt, bths_pt, ggamma_pt,             &
                             csmct_pt, el_cons_pt, el_cons_s_pt,            &
                             el_comp_pt, el_comp_s_pt, macro_el_pt,         &
                             macro_el_s_pt, v_pt, v_s_pt,                   &
                             e_piezo_tensor_pt, d_piezo_tensor_pt,          &
                             celldm_pt_p1, celldm_pt_m1,                    &
                             debye_macro_el_pt, debye_macro_el_s_pt, uint_pt,&
                             tau_pt, uint_zsisa_pt, tau_zsisa_pt,           &
                             alpha_int_pt, alpha_int_zsisa_pt

  USE ph_freq_anharmonic_pt, ONLY : vminf_pt, b0f_pt, b01f_pt, b02f_pt,     &
                             eminf_pt, densityf_pt, celldmf_pt, enerf_pt,   &  
                             free_enerf_pt, entrf_pt, cef_pt, cvf_pt,       &
                             cpf_pt, b0f_s_pt, betaf_pt, gammaf_pt,         &
                             alphaf_anis_pt, el_consf_pt, el_compf_pt,      &
                             cpmcef_anis_pt, el_consf_s_pt, el_compf_s_pt,  &
                             e_piezo_tensorf_pt, d_piezo_tensorf_pt,        &
                             bthsf_pt, ggammaf_pt, csmctf_pt, macro_elf_pt, &
                             macro_elf_s_pt, vf_pt, vf_s_pt, b0f_ec_pt,     &
                             debye_macro_elf_pt, debye_macro_elf_s_pt,      &
                             uintf_pt, tauf_pt, uintf_zsisa_pt,             &
                             tauf_zsisa_pt, alphaf_int_pt, alphaf_int_zsisa_pt
  USE anharmonic_ptt, ONLY : vmin_ptt, b0_ptt, b01_ptt, b02_ptt, emin_ptt,  &
                             ce_ptt, cv_ptt, cp_ptt, beta_ptt, b0_s_ptt,    &
                             celldm_ptt, gamma_ptt, el_cons_ptt,            &
                             el_cons_s_ptt, el_comp_ptt, el_comp_s_ptt,     &
                             macro_el_ptt, macro_el_s_ptt, v_ptt, v_s_ptt,  &
                             e_piezo_tensor_ptt, d_piezo_tensor_ptt,        &
                             celldm_ptt_p1, celldm_ptt_m1, emin_ptt_p1,     &
                             emin_ptt_m1, vmin_ptt_p1, vmin_ptt_m1,         &
                             alpha_anis_ptt, density_ptt, cpmce_anis_ptt,   &
                             bths_ptt, ggamma_ptt, csmct_ptt, ener_ptt,     &
                             entr_ptt, b0_ec_ptt, debye_macro_el_ptt,       &
                             debye_macro_el_s_ptt, uint_ptt, tau_ptt,       &
                             uint_zsisa_ptt, tau_zsisa_ptt, alpha_int_ptt,  &
                             alpha_int_zsisa_ptt, uint_ptt_p1, uint_ptt_m1, &
                             tau_ptt_p1, tau_ptt_m1
  USE ph_freq_anharmonic_ptt,  ONLY : vminf_ptt, b0f_ptt, b01f_ptt, b02f_ptt,&
                             eminf_ptt, celldmf_ptt, eminf_ptt_p1,           &
                             eminf_ptt_m1, celldmf_ptt_p1, celldmf_ptt_m1,  &
                             el_consf_ptt, e_piezo_tensorf_ptt,             &
                             d_piezo_tensorf_ptt,                           &
                             vminf_ptt_p1, vminf_ptt_m1, densityf_ptt,      &
                             enerf_ptt, entrf_ptt, cef_ptt,                 &
                             cvf_ptt, cpf_ptt, b0f_s_ptt, betaf_ptt,        &
                             gammaf_ptt, b0f_ec_ptt, vf_ptt, vf_s_ptt,      &
                             debye_macro_elf_ptt, debye_macro_elf_s_ptt,    &
                             uintf_ptt, tauf_ptt, uintf_zsisa_ptt,          &
                             tauf_zsisa_ptt, alphaf_int_ptt,                &
                             alphaf_int_zsisa_ptt, uintf_ptt_p1, uintf_ptt_m1,&
                             tauf_ptt_p1, tauf_ptt_m1
  USE anharmonic_vt,  ONLY : press_vt, press_vtt
  USE ph_freq_anharmonic_vt,  ONLY : pressf_vt, pressf_vtt
  USE el_anharmonic, ONLY :  vmine_t, b0e_t, b01e_t, b02e_t, free_e_mine_t, &
                             el_free_ener_pt, el_ener_pt, el_ce_pt,         &
                             el_entr_pt, el_ce_ptt, el_b0_t, el_beta_t,     &
                             el_betaf_t, el_gamma_t, el_gammaf_t,           &
                             el_cp_t, el_cpf_t, el_enerf_pt, el_cef_pt,     &
                             el_entrf_pt, el_cef_ptt
  USE ph_freq_anharmonic, ONLY : vminf_t, b0f_t, b01f_t, b02f_t,            &
                             free_e_minf_t, af_t, enerf_t, free_enerf_t,    &
                             entropyf_t, cef_t, cvf_t, cpf_t, b0f_s,        &
                             alphaf_t, betaf_t, gammaf_t, celldmf_t,        &
                             alphaf_anis_t, cpmcef_anis, bfactf_t,          &
                             bthsf_t, ggammaf_t, el_consf_t, el_compf_t,    &
                             el_consf_s, el_compf_s, macro_elf_s,           &
                             macro_elf_t, vf_t, vf_s, el_conf_geo_t,        &
                             densityf_t, csmctf_t, p1tf_t, p2tf_t, p3tf_t,  &
                             p4tf_t, p1tf_noe_t, p2tf_noe_t, p3tf_noe_t,    &
                             p4tf_noe_t, celldmf_noe_t, free_e_minf_noe_t,  &
                             vminf_noe_t, densityf_noe_t, b0f_noe_t,        &
                             b01f_noe_t, b02f_noe_t, e_piezo_tensorf_t,     &
                             d_piezo_tensorf_t,                             &
                             free_enerf_noe_t, enerf_noe_t, entropyf_noe_t, &
                             cef_noe_t, cvf_noe_t, b0f_ec_t, b0f_ec_s,      &
                             debye_macro_elf_t, debye_macro_elf_s, dydef_t, &
                             uintf_t, tauf_t, uintf_zsisa_t, tauf_zsisa_t,  &
                             alphaf_int_t, alphaf_int_zsisa_t
  USE grun_anharmonic,  ONLY : betab, alpha_an_g, cp_grun_t, cv_grun_t,     &
                             ce_grun_t, b0_grun_s,                          &
                             grun_gamma_t, poly_grun, poly_grun_red,        &
                             grun_cpmce_anis, el_cons_grun_t, el_comp_grun_t
  USE control_paths,    ONLY : xqaux, wqaux, letter, label_list, letter_path, &
                               label_disp_q, disp_q, disp_wq, nrap_plot_in,   &
                               rap_plot_in, nrap_plot, rap_plot, high_sym_path
  USE control_2d_bands, ONLY : averag, vacuum, aux_ind_sur
  USE control_gen_gruneisen, ONLY : xngeo, ind_rec3
  USE initial_conf,     ONLY : ityp_save, tau_save, tau_save_crys, &
                               collect_info_save, atm_save
  USE equilibrium_conf, ONLY : tau0, tau0_crys
  USE control_thermo,   ONLY : all_geometries_together
  USE control_grun,     ONLY : vgrun_t, b0_grun_t, celldm_grun_t
  USE temperature,      ONLY : temp, temp_plot, itemp_plot, sigma_ry
  USE control_pressure, ONLY : npress

  USE control_conv,     ONLY : ke, keden, nk_test, sigma_test
  USE control_mur_p,    ONLY : vmin_p, b0_p, b01_p, b02_p, emin_p
  USE control_mur,      ONLY : p0
  USE control_ev,       ONLY : e0, v0
  USE elastic_constants, ONLY : epsilon_geo, sigma_geo, epsilon_voigt
  USE piezoelectric_tensor, ONLY : deallocate_piezo
  USE control_elastic_constants, ONLY : rot_mat, el_con_geo, &
                             el_con_ibrav_geo, el_con_celldm_geo, &
                             el_con_tau_crys_geo, el_con_omega_geo, &
                             el_con_tau_geo,                       &
                             epsil_geo, all_geometry_done_geo,     &
                             found_dos_ec, found_ph_ec, tau_acc, min_y, &
                             epsil_y, min_y_t, dyde, tau_save_ec, &
                             el_con_at_geo
  USE control_pressure, ONLY : press_plot, ipress_plot
  USE uniform_pressure, ONLY : omega_p, density_p, celldm_p, p2_p, p4_p
  USE control_atomic_pos, ONLY : uint0, uint_eq, tau_eq, p2_eq, p4_eq, &
                                 uint_p, tau_p, nint_var
  USE control_vol,      ONLY : ivol_plot
  USE control_debye,    ONLY : deb_energy, deb_free_energy, deb_entropy, &
                               deb_cv, deb_b_fact, deb_bfact
  USE control_quadratic_energy, ONLY : p2, hessian_v, hessian_e, x_pos_min
  USE control_quartic_energy, ONLY :  p4, x_min_4, hessian4_v, hessian4_e, &
                            lquartic
  USE el_anharmonic, ONLY : el_energy_t, el_free_energy_t, el_entropy_t, &
                            el_ce_t, el_energyf_t, el_free_energyf_t,    &
                            el_entropyf_t, el_cef_t
  USE control_emp_free_ener, ONLY : emp_ener, emp_free_ener, emp_entr, &
                           emp_ce
  USE emp_anharmonic, ONLY : emp_energy_t, emp_free_energy_t, emp_entropy_t, &
                            emp_ce_t, emp_energyf_t, emp_free_energyf_t, &
                            emp_entropyf_t, emp_cef_t, emp_free_ener_pt,  &
                            emp_ener_pt, emp_entr_pt, emp_ce_pt, emp_ce_ptt, &
                            emp_free_enerf_pt, emp_enerf_pt, emp_entrf_pt, &
                            emp_cef_pt, emp_cef_ptt
  USE collect_info,  ONLY : destroy_collect_info_type
  USE control_eldos, ONLY : dos_k, dos_wk
  USE geometry_file, ONLY : deallocate_geometry_file
  USE polynomial,    ONLY : clean_poly

  IMPLICIT NONE
  INTEGER :: igeom, ipress, ivar
  !
  IF ( ALLOCATED (energy_geo) )      DEALLOCATE(energy_geo)
  IF ( ALLOCATED (ef_geo) )          DEALLOCATE(ef_geo)
  IF ( ALLOCATED (celldm_geo) )      DEALLOCATE(celldm_geo) 
  IF ( ALLOCATED (omega_geo) )       DEALLOCATE(omega_geo) 
  IF ( ALLOCATED (ibrav_geo) )       DEALLOCATE(ibrav_geo) 
  IF ( ALLOCATED (tau_geo) )         DEALLOCATE(tau_geo)
  IF ( ALLOCATED (no_ph) )           DEALLOCATE(no_ph) 
  IF ( ALLOCATED (no_ph_eos) )       DEALLOCATE(no_ph_eos) 
  IF ( ALLOCATED (found_dos_ec) )    DEALLOCATE(found_dos_ec) 
  IF ( ALLOCATED (found_ph_ec) )     DEALLOCATE(found_ph_ec) 
  IF ( ALLOCATED (dynmat_on_file) )  DEALLOCATE(dynmat_on_file) 
  IF ( ALLOCATED (vmin_t) )          DEALLOCATE(vmin_t) 
  IF ( ALLOCATED (ke) )              DEALLOCATE(ke) 
  IF ( ALLOCATED (keden) )           DEALLOCATE(keden) 
  IF ( ALLOCATED (nk_test) )         DEALLOCATE(nk_test) 
  IF ( ALLOCATED (sigma_test) )      DEALLOCATE(sigma_test) 
  IF ( ALLOCATED (all_geometry_done_geo)) DEALLOCATE(all_geometry_done_geo)
  IF ( ALLOCATED (uint_p))           DEALLOCATE(uint_p)
  IF ( ALLOCATED (tau_p))            DEALLOCATE(tau_p)

  IF ( ALLOCATED (ph_free_ener) )    DEALLOCATE(ph_free_ener)
  IF ( ALLOCATED (ph_ener) )         DEALLOCATE(ph_ener)
  IF ( ALLOCATED (ph_entropy) )      DEALLOCATE(ph_entropy)
  IF ( ALLOCATED (ph_t_debye) )      DEALLOCATE(ph_t_debye)
  IF ( ALLOCATED (ph_e0) )           DEALLOCATE(ph_e0)
  IF ( ALLOCATED (ph_ce) )           DEALLOCATE(ph_ce)
  IF ( ALLOCATED (ph_b_fact) )       DEALLOCATE(ph_b_fact)
  IF ( ALLOCATED (ph_free_ener_eos)) DEALLOCATE(ph_free_ener_eos)
  IF ( ALLOCATED (ph_e0_eos))        DEALLOCATE(ph_e0_eos)
  IF ( ALLOCATED (ph_ener_eos))      DEALLOCATE(ph_ener_eos)
  IF ( ALLOCATED (ph_entropy_eos))   DEALLOCATE(ph_entropy_eos)
  IF ( ALLOCATED (ph_ce_eos))        DEALLOCATE(ph_ce_eos)
  IF ( ALLOCATED (uint_geo_eos_t))   DEALLOCATE(uint_geo_eos_t)
  IF ( ALLOCATED (tau_geo_eos_t))    DEALLOCATE(tau_geo_eos_t)
  IF ( ALLOCATED (phf_free_ener) )   DEALLOCATE(phf_free_ener)
  IF ( ALLOCATED (phf_e0_eos))       DEALLOCATE(phf_e0_eos)
  IF ( ALLOCATED (phf_ener) )        DEALLOCATE(phf_ener)
  IF ( ALLOCATED (phf_entropy) )     DEALLOCATE(phf_entropy)
  IF ( ALLOCATED (phf_t_debye) )     DEALLOCATE(phf_t_debye)
  IF ( ALLOCATED (phf_e0) )          DEALLOCATE(phf_e0)
  IF ( ALLOCATED (phf_ce) )          DEALLOCATE(phf_ce)
  IF ( ALLOCATED (phf_b_fact) )      DEALLOCATE(phf_b_fact)
  IF ( ALLOCATED (phf_free_ener_eos)) DEALLOCATE(phf_free_ener_eos)
  IF ( ALLOCATED (phf_ener_eos))     DEALLOCATE(phf_ener_eos)
  IF ( ALLOCATED (phf_entropy_eos))  DEALLOCATE(phf_entropy_eos)
  IF ( ALLOCATED (phf_ce_eos))       DEALLOCATE(phf_ce_eos)
  IF ( ALLOCATED (uintf_geo_eos_t))  DEALLOCATE(uintf_geo_eos_t)
  IF ( ALLOCATED (tauf_geo_eos_t))   DEALLOCATE(tauf_geo_eos_t)

  IF ( ALLOCATED (el_free_ener) )    DEALLOCATE(el_free_ener)
  IF ( ALLOCATED (el_ener) )         DEALLOCATE(el_ener)
  IF ( ALLOCATED (el_entr) )         DEALLOCATE(el_entr)
  IF ( ALLOCATED (el_mu) )           DEALLOCATE(el_mu)
  IF ( ALLOCATED (el_ce) )           DEALLOCATE(el_ce)

  IF ( ALLOCATED (el_free_ener_eos) )DEALLOCATE(el_free_ener_eos)
  IF ( ALLOCATED (el_ener_eos) )     DEALLOCATE(el_ener_eos)
  IF ( ALLOCATED (el_entr_eos) )     DEALLOCATE(el_entr_eos)
  IF ( ALLOCATED (el_mu_eos) )       DEALLOCATE(el_mu_eos)
  IF ( ALLOCATED (el_ce_eos) )       DEALLOCATE(el_ce_eos)

  IF ( ALLOCATED (elf_free_ener_eos) )DEALLOCATE(elf_free_ener_eos)
  IF ( ALLOCATED (elf_ener_eos) )     DEALLOCATE(elf_ener_eos)
  IF ( ALLOCATED (elf_entr_eos) )     DEALLOCATE(elf_entr_eos)
  IF ( ALLOCATED (elf_mu_eos) )       DEALLOCATE(elf_mu_eos)
  IF ( ALLOCATED (elf_ce_eos) )       DEALLOCATE(elf_ce_eos)

  IF ( ALLOCATED (emp_free_ener) )   DEALLOCATE(emp_free_ener)
  IF ( ALLOCATED (emp_ener) )        DEALLOCATE(emp_ener)
  IF ( ALLOCATED (emp_entr) )        DEALLOCATE(emp_entr)
  IF ( ALLOCATED (emp_ce) )          DEALLOCATE(emp_ce)

  IF ( ALLOCATED (vmin_t) )          DEALLOCATE(vmin_t) 
  IF ( ALLOCATED (b0_t) )            DEALLOCATE(b0_t) 
  IF ( ALLOCATED (b0_ec_t) )         DEALLOCATE(b0_ec_t) 
  IF ( ALLOCATED (b01_t) )           DEALLOCATE(b01_t) 
  IF ( ALLOCATED (b02_t) )           DEALLOCATE(b02_t) 
  IF ( ALLOCATED (free_e_min_t) )    DEALLOCATE(free_e_min_t) 
  IF ( ALLOCATED (p1t_t) )           DEALLOCATE(p1t_t)
  IF ( ALLOCATED (p2t_t) )           DEALLOCATE(p2t_t)
  IF ( ALLOCATED (p3t_t) )           DEALLOCATE(p3t_t)
  IF ( ALLOCATED (p4t_t) )           DEALLOCATE(p4t_t)

  IF ( ALLOCATED (ener_t) )          DEALLOCATE(ener_t) 
  IF ( ALLOCATED (free_ener_t) )     DEALLOCATE(free_ener_t)
  IF ( ALLOCATED (entropy_t) )       DEALLOCATE(entropy_t) 
  IF ( ALLOCATED (ce_t) )            DEALLOCATE(ce_t) 
  IF ( ALLOCATED (cv_t) )            DEALLOCATE(cv_t) 

  IF ( ALLOCATED (cp_t) )            DEALLOCATE(cp_t) 
  IF ( ALLOCATED (b0_s) )            DEALLOCATE(b0_s) 
  IF ( ALLOCATED (b0_ec_s) )         DEALLOCATE(b0_ec_s) 
  IF ( ALLOCATED (alpha_t) )         DEALLOCATE(alpha_t) 
  IF ( ALLOCATED (beta_t) )          DEALLOCATE(beta_t) 
  IF ( ALLOCATED (gamma_t) )         DEALLOCATE(gamma_t) 
  IF ( ALLOCATED (uint_t) )          DEALLOCATE(uint_t) 
  IF ( ALLOCATED (tau_t) )           DEALLOCATE(tau_t) 
  IF ( ALLOCATED (alpha_int_t) )     DEALLOCATE(alpha_int_t) 
  IF ( ALLOCATED (uint_zsisa_t) )    DEALLOCATE(uint_zsisa_t) 
  IF ( ALLOCATED (tau_zsisa_t) )     DEALLOCATE(tau_zsisa_t) 
  IF ( ALLOCATED (alpha_int_zsisa_t) ) DEALLOCATE(alpha_int_zsisa_t) 

  IF ( ALLOCATED (celldm_t) )        DEALLOCATE(celldm_t) 
  IF ( ALLOCATED (celldm_t_p1) )     DEALLOCATE(celldm_t_p1) 
  IF ( ALLOCATED (celldm_t_m1) )     DEALLOCATE(celldm_t_m1) 
  IF ( ALLOCATED (density_t) )       DEALLOCATE(density_t) 
  IF ( ALLOCATED (alpha_anis_t) )    DEALLOCATE(alpha_anis_t) 
  IF ( ALLOCATED (cpmce_anis) )      DEALLOCATE(cpmce_anis) 
  IF ( ALLOCATED (bfact_t) )         DEALLOCATE(bfact_t) 
  IF ( ALLOCATED (csmct_t) )         DEALLOCATE(csmct_t) 
  IF ( ALLOCATED (bths_t) )          DEALLOCATE(bths_t) 
  IF ( ALLOCATED (ggamma_t) )        DEALLOCATE(ggamma_t) 

  IF ( ALLOCATED (el_cons_t) )       DEALLOCATE(el_cons_t)
  IF ( ALLOCATED (el_comp_t) )       DEALLOCATE(el_comp_t)
  IF ( ALLOCATED (el_cons_s) )       DEALLOCATE(el_cons_s)
  IF ( ALLOCATED (el_comp_s) )       DEALLOCATE(el_comp_s)
  IF ( ALLOCATED (macro_el_t) )      DEALLOCATE(macro_el_t)
  IF ( ALLOCATED (macro_el_s) )      DEALLOCATE(macro_el_s)
  IF ( ALLOCATED (debye_macro_el_t) ) DEALLOCATE(debye_macro_el_t)
  IF ( ALLOCATED (debye_macro_el_s) ) DEALLOCATE(debye_macro_el_s)
  IF ( ALLOCATED (v_t) )             DEALLOCATE(v_t)
  IF ( ALLOCATED (v_s) )             DEALLOCATE(v_s)
  IF ( ALLOCATED (el_con_geo_t) )    DEALLOCATE(el_con_geo_t) 
  IF ( ALLOCATED (dyde_t) )          DEALLOCATE(dyde_t)
  IF ( ALLOCATED (tau_save_ec) )     DEALLOCATE(tau_save_ec)
  IF ( ALLOCATED (e_piezo_tensor_t) ) DEALLOCATE(e_piezo_tensor_t)
  IF ( ALLOCATED (d_piezo_tensor_t) ) DEALLOCATE(d_piezo_tensor_t)

  IF ( ALLOCATED (vmin_noe_t) )       DEALLOCATE(vmin_noe_t) 
  IF ( ALLOCATED (density_noe_t) )    DEALLOCATE(density_noe_t) 
  IF ( ALLOCATED (celldm_noe_t) )     DEALLOCATE(celldm_noe_t) 
  IF ( ALLOCATED (celldm_noe_t_p1) )  DEALLOCATE(celldm_noe_t_p1) 
  IF ( ALLOCATED (celldm_noe_t_m1) )  DEALLOCATE(celldm_noe_t_m1) 
  IF ( ALLOCATED (b0_noe_t) )         DEALLOCATE(b0_noe_t) 
  IF ( ALLOCATED (b01_noe_t) )        DEALLOCATE(b01_noe_t) 
  IF ( ALLOCATED (b02_noe_t) )        DEALLOCATE(b02_noe_t) 
  IF ( ALLOCATED (free_e_min_noe_t) ) DEALLOCATE(free_e_min_noe_t)
  IF ( ALLOCATED (a_noe_t) )          DEALLOCATE(a_noe_t)
  IF ( ALLOCATED (p1t_noe_t) )        DEALLOCATE(p1t_noe_t)
  IF ( ALLOCATED (p2t_noe_t) )        DEALLOCATE(p2t_noe_t)
  IF ( ALLOCATED (p3t_noe_t) )        DEALLOCATE(p3t_noe_t)
  IF ( ALLOCATED (p4t_noe_t) )        DEALLOCATE(p4t_noe_t)

  IF ( ALLOCATED (free_ener_noe_t) )  DEALLOCATE(free_ener_noe_t) 
  IF ( ALLOCATED (ener_noe_t) )       DEALLOCATE(ener_noe_t) 
  IF ( ALLOCATED (entropy_noe_t) )    DEALLOCATE(entropy_noe_t) 
  IF ( ALLOCATED (ce_noe_t) )         DEALLOCATE(ce_noe_t) 
  IF ( ALLOCATED (cv_noe_t) )         DEALLOCATE(cv_noe_t) 

  IF ( ALLOCATED (cp_noe_t) )         DEALLOCATE(cp_noe_t) 
  IF ( ALLOCATED (b0_noe_s) )         DEALLOCATE(b0_noe_s) 
  IF ( ALLOCATED (beta_noe_t) )       DEALLOCATE(beta_noe_t)
  IF ( ALLOCATED (gamma_noe_t) )      DEALLOCATE(gamma_noe_t) 

  IF (ALLOCATED (emp_free_energy_t) ) DEALLOCATE(emp_free_energy_t)
  IF (ALLOCATED (emp_energy_t) )      DEALLOCATE(emp_energy_t)
  IF (ALLOCATED (emp_entropy_t) )     DEALLOCATE(emp_entropy_t)
  IF (ALLOCATED (emp_ce_t) )          DEALLOCATE(emp_ce_t)

  IF ( ALLOCATED (emp_free_energyf_t) ) DEALLOCATE(emp_free_energyf_t)
  IF ( ALLOCATED (emp_energyf_t) )      DEALLOCATE(emp_energyf_t)
  IF ( ALLOCATED (emp_entropyf_t) )     DEALLOCATE(emp_entropyf_t)
  IF ( ALLOCATED (emp_cef_t) )          DEALLOCATE(emp_cef_t)


  IF ( ALLOCATED (vmin_ptt) )        DEALLOCATE(vmin_ptt) 
  IF ( ALLOCATED (vmin_ptt_p1) )     DEALLOCATE(vmin_ptt_p1) 
  IF ( ALLOCATED (vmin_ptt_m1) )     DEALLOCATE(vmin_ptt_m1) 
  IF ( ALLOCATED (b0_ptt) )          DEALLOCATE(b0_ptt) 
  IF ( ALLOCATED (b0_ec_ptt) )       DEALLOCATE(b0_ec_ptt) 
  IF ( ALLOCATED (b01_ptt) )         DEALLOCATE(b01_ptt) 
  IF ( ALLOCATED (b02_ptt) )         DEALLOCATE(b02_ptt) 
  IF ( ALLOCATED (emin_ptt) )        DEALLOCATE(emin_ptt) 
  IF ( ALLOCATED (emin_ptt_p1) )     DEALLOCATE(emin_ptt_p1) 
  IF ( ALLOCATED (emin_ptt_m1) )     DEALLOCATE(emin_ptt_m1) 
  IF ( ALLOCATED (press_vtt) )       DEALLOCATE(press_vtt) 
  IF ( ALLOCATED (pressf_vtt) )      DEALLOCATE(pressf_vtt) 

  IF ( ALLOCATED (ener_ptt) )        DEALLOCATE(ener_ptt) 
  IF ( ALLOCATED (entr_ptt) )        DEALLOCATE(entr_ptt) 
  IF ( ALLOCATED (ce_ptt) )          DEALLOCATE(ce_ptt) 
  IF ( ALLOCATED (cv_ptt) )          DEALLOCATE(cv_ptt) 

  IF ( ALLOCATED (cp_ptt) )          DEALLOCATE(cp_ptt) 
  IF ( ALLOCATED (beta_ptt) )        DEALLOCATE(beta_ptt) 
  IF ( ALLOCATED (b0_s_ptt) )        DEALLOCATE(b0_s_ptt) 
  IF ( ALLOCATED (gamma_ptt) )       DEALLOCATE(gamma_ptt) 
  IF ( ALLOCATED (uint_ptt) )        DEALLOCATE(uint_ptt) 
  IF ( ALLOCATED (tau_ptt) )         DEALLOCATE(tau_ptt) 
  IF ( ALLOCATED (uint_ptt_p1) )     DEALLOCATE(uint_ptt_p1) 
  IF ( ALLOCATED (tau_ptt_p1) )      DEALLOCATE(tau_ptt_p1) 
  IF ( ALLOCATED (uint_ptt_m1) )     DEALLOCATE(uint_ptt_m1) 
  IF ( ALLOCATED (tau_ptt_m1) )      DEALLOCATE(tau_ptt_m1) 
  IF ( ALLOCATED (alpha_int_ptt) )   DEALLOCATE(alpha_int_ptt) 
  IF ( ALLOCATED (uint_zsisa_ptt) )  DEALLOCATE(uint_zsisa_ptt) 
  IF ( ALLOCATED (tau_zsisa_ptt) )   DEALLOCATE(tau_zsisa_ptt) 
  IF ( ALLOCATED (alpha_int_zsisa_ptt) ) DEALLOCATE(alpha_int_zsisa_ptt) 

  IF ( ALLOCATED (celldm_ptt) )      DEALLOCATE(celldm_ptt) 
  IF ( ALLOCATED (celldm_ptt_p1) )   DEALLOCATE(celldm_ptt_p1) 
  IF ( ALLOCATED (celldm_ptt_m1) )   DEALLOCATE(celldm_ptt_m1) 
  IF ( ALLOCATED (alpha_anis_ptt) )  DEALLOCATE(alpha_anis_ptt) 
  IF ( ALLOCATED (density_ptt) )     DEALLOCATE(density_ptt) 
  IF ( ALLOCATED (cpmce_anis_ptt) )  DEALLOCATE(cpmce_anis_ptt) 
  IF ( ALLOCATED (csmct_ptt) )       DEALLOCATE(csmct_ptt) 
  IF ( ALLOCATED (bths_ptt) )        DEALLOCATE(bths_ptt) 
  IF ( ALLOCATED (ggamma_ptt) )      DEALLOCATE(ggamma_ptt) 
  IF ( ALLOCATED (el_cons_ptt) )     DEALLOCATE(el_cons_ptt)
  IF ( ALLOCATED (el_comp_ptt) )     DEALLOCATE(el_comp_ptt)
  IF ( ALLOCATED (el_cons_s_ptt) )   DEALLOCATE(el_cons_s_ptt)
  IF ( ALLOCATED (el_comp_s_ptt) )   DEALLOCATE(el_comp_s_ptt)
  IF ( ALLOCATED (macro_el_ptt) )    DEALLOCATE(macro_el_ptt)
  IF ( ALLOCATED (macro_el_s_ptt) )  DEALLOCATE(macro_el_s_ptt)
  IF ( ALLOCATED (v_ptt) )           DEALLOCATE(v_ptt)
  IF ( ALLOCATED (v_s_ptt) )         DEALLOCATE(v_s_ptt)
  IF ( ALLOCATED (debye_macro_el_ptt) ) DEALLOCATE(debye_macro_el_ptt)
  IF ( ALLOCATED (debye_macro_el_s_ptt) ) DEALLOCATE(debye_macro_el_s_ptt)
  IF ( ALLOCATED (e_piezo_tensor_ptt) ) DEALLOCATE(e_piezo_tensor_ptt)
  IF ( ALLOCATED (d_piezo_tensor_ptt) ) DEALLOCATE(d_piezo_tensor_ptt)

  IF ( ALLOCATED (vminf_ptt) )        DEALLOCATE(vminf_ptt) 
  IF ( ALLOCATED (vminf_ptt_p1) )     DEALLOCATE(vminf_ptt_p1) 
  IF ( ALLOCATED (vminf_ptt_m1) )     DEALLOCATE(vminf_ptt_m1) 
  IF ( ALLOCATED (b0f_ptt) )          DEALLOCATE(b0f_ptt) 
  IF ( ALLOCATED (b0f_ec_ptt) )       DEALLOCATE(b0f_ec_ptt) 
  IF ( ALLOCATED (b01f_ptt) )         DEALLOCATE(b01f_ptt) 
  IF ( ALLOCATED (b02f_ptt) )         DEALLOCATE(b02f_ptt) 
  IF ( ALLOCATED (eminf_ptt) )        DEALLOCATE(eminf_ptt) 
  IF ( ALLOCATED (eminf_ptt_p1) )     DEALLOCATE(eminf_ptt_p1) 
  IF ( ALLOCATED (eminf_ptt_m1) )     DEALLOCATE(eminf_ptt_m1) 
  IF ( ALLOCATED (celldmf_ptt) )      DEALLOCATE(celldmf_ptt) 
  IF ( ALLOCATED (celldmf_ptt_p1) )   DEALLOCATE(celldmf_ptt_p1) 
  IF ( ALLOCATED (celldmf_ptt_m1) )   DEALLOCATE(celldmf_ptt_m1) 
  IF ( ALLOCATED (el_consf_ptt) )     DEALLOCATE(el_consf_ptt) 
  IF ( ALLOCATED (e_piezo_tensorf_ptt) ) DEALLOCATE(e_piezo_tensorf_ptt) 
  IF ( ALLOCATED (d_piezo_tensorf_ptt) ) DEALLOCATE(d_piezo_tensorf_ptt) 

  IF ( ALLOCATED (enerf_ptt) )        DEALLOCATE(enerf_ptt) 
  IF ( ALLOCATED (entrf_ptt) )        DEALLOCATE(entrf_ptt) 
  IF ( ALLOCATED (cef_ptt) )          DEALLOCATE(cef_ptt) 
  IF ( ALLOCATED (cvf_ptt) )          DEALLOCATE(cvf_ptt) 

  IF ( ALLOCATED (cpf_ptt) )          DEALLOCATE(cpf_ptt) 
  IF ( ALLOCATED (betaf_ptt) )        DEALLOCATE(betaf_ptt) 
  IF ( ALLOCATED (b0f_s_ptt) )        DEALLOCATE(b0f_s_ptt) 
  IF ( ALLOCATED (gammaf_ptt) )       DEALLOCATE(gammaf_ptt) 
  IF ( ALLOCATED (uintf_ptt) )        DEALLOCATE(uintf_ptt) 
  IF ( ALLOCATED (tauf_ptt) )         DEALLOCATE(tauf_ptt) 
  IF ( ALLOCATED (uintf_ptt_p1) )     DEALLOCATE(uintf_ptt_p1) 
  IF ( ALLOCATED (tauf_ptt_p1) )      DEALLOCATE(tauf_ptt_p1) 
  IF ( ALLOCATED (uintf_ptt_m1) )     DEALLOCATE(uintf_ptt_m1) 
  IF ( ALLOCATED (tauf_ptt_m1) )      DEALLOCATE(tauf_ptt_m1) 
  IF ( ALLOCATED (alphaf_int_ptt) )   DEALLOCATE(alphaf_int_ptt) 
  IF ( ALLOCATED (uintf_zsisa_ptt) )  DEALLOCATE(uintf_zsisa_ptt) 
  IF ( ALLOCATED (tauf_zsisa_ptt) )   DEALLOCATE(tauf_zsisa_ptt) 
  IF ( ALLOCATED (alphaf_int_zsisa_ptt) ) DEALLOCATE(alphaf_int_zsisa_ptt) 

  IF ( ALLOCATED (vf_ptt) )           DEALLOCATE(vf_ptt) 
  IF ( ALLOCATED (vf_s_ptt) )         DEALLOCATE(vf_s_ptt) 
  IF ( ALLOCATED (debye_macro_elf_ptt) ) DEALLOCATE(debye_macro_elf_ptt) 
  IF ( ALLOCATED (debye_macro_elf_s_ptt) ) DEALLOCATE(debye_macro_elf_s_ptt) 

  IF ( ALLOCATED (vmin_pt) )         DEALLOCATE(vmin_pt) 
  IF ( ALLOCATED (b0_pt) )           DEALLOCATE(b0_pt) 
  IF ( ALLOCATED (b01_pt) )          DEALLOCATE(b01_pt) 
  IF ( ALLOCATED (b02_pt) )          DEALLOCATE(b02_pt) 
  IF ( ALLOCATED (emin_pt) )         DEALLOCATE(emin_pt) 
  IF ( ALLOCATED (press_vt) )        DEALLOCATE(press_vt) 
  IF ( ALLOCATED (pressf_vt) )       DEALLOCATE(pressf_vt) 

  IF ( ALLOCATED (ener_pt) )         DEALLOCATE(ener_pt) 
  IF ( ALLOCATED (entr_pt) )         DEALLOCATE(entr_pt) 
  IF ( ALLOCATED (ce_pt) )           DEALLOCATE(ce_pt) 
  IF ( ALLOCATED (cv_pt) )           DEALLOCATE(cv_pt) 

  IF ( ALLOCATED (cp_pt) )           DEALLOCATE(cp_pt) 
  IF ( ALLOCATED (beta_pt) )         DEALLOCATE(beta_pt) 
  IF ( ALLOCATED (b0_s_pt) )         DEALLOCATE(b0_s_pt) 
  IF ( ALLOCATED (gamma_pt) )        DEALLOCATE(gamma_pt) 
  IF ( ALLOCATED (uint_pt) )         DEALLOCATE(uint_pt) 
  IF ( ALLOCATED (tau_pt) )          DEALLOCATE(tau_pt) 
  IF ( ALLOCATED (alpha_int_pt) )    DEALLOCATE(alpha_int_pt) 
  IF ( ALLOCATED (uint_zsisa_pt) )   DEALLOCATE(uint_zsisa_pt) 
  IF ( ALLOCATED (tau_zsisa_pt) )    DEALLOCATE(tau_zsisa_pt) 
  IF ( ALLOCATED (alpha_int_zsisa_pt) ) DEALLOCATE(alpha_int_zsisa_pt) 

  IF ( ALLOCATED (celldm_pt) )        DEALLOCATE(celldm_pt) 
  IF ( ALLOCATED (celldm_pt_p1) )     DEALLOCATE(celldm_pt_p1) 
  IF ( ALLOCATED (celldm_pt_m1) )     DEALLOCATE(celldm_pt_m1) 
  IF ( ALLOCATED (alpha_anis_pt) )    DEALLOCATE(alpha_anis_pt)
  IF ( ALLOCATED (density_pt) )       DEALLOCATE(density_pt)
  IF ( ALLOCATED (cpmce_anis_pt) )    DEALLOCATE(cpmce_anis_pt)
  IF ( ALLOCATED (csmct_pt) )         DEALLOCATE(csmct_pt)
  IF ( ALLOCATED (bths_pt) )          DEALLOCATE(bths_pt)
  IF ( ALLOCATED (ggamma_pt) )        DEALLOCATE(ggamma_pt)
  IF ( ALLOCATED (el_cons_pt) )       DEALLOCATE(el_cons_pt)
  IF ( ALLOCATED (el_comp_pt) )       DEALLOCATE(el_comp_pt)
  IF ( ALLOCATED (el_cons_s_pt) )     DEALLOCATE(el_cons_s_pt)
  IF ( ALLOCATED (el_comp_s_pt) )     DEALLOCATE(el_comp_s_pt)
  IF ( ALLOCATED (macro_el_pt) )      DEALLOCATE(macro_el_pt)
  IF ( ALLOCATED (macro_el_s_pt) )    DEALLOCATE(macro_el_s_pt)
  IF ( ALLOCATED (v_pt) )             DEALLOCATE(v_pt)
  IF ( ALLOCATED (v_s_pt) )           DEALLOCATE(v_s_pt)
  IF ( ALLOCATED (debye_macro_el_pt) ) DEALLOCATE(debye_macro_el_pt)
  IF ( ALLOCATED (debye_macro_el_s_pt) ) DEALLOCATE(debye_macro_el_s_pt)
  IF ( ALLOCATED (e_piezo_tensor_pt) ) DEALLOCATE(e_piezo_tensor_pt)
  IF ( ALLOCATED (d_piezo_tensor_pt) ) DEALLOCATE(d_piezo_tensor_pt)


  IF ( ALLOCATED (emp_free_ener_pt) ) DEALLOCATE(emp_free_ener_pt)
  IF ( ALLOCATED (emp_ener_pt) )      DEALLOCATE(emp_ener_pt)
  IF ( ALLOCATED (emp_entr_pt) )      DEALLOCATE(emp_entr_pt)
  IF ( ALLOCATED (emp_ce_pt) )        DEALLOCATE(emp_ce_pt)
  IF ( ALLOCATED (emp_ce_ptt) )       DEALLOCATE(emp_ce_ptt)

  IF ( ALLOCATED (vminf_pt) )         DEALLOCATE(vminf_pt) 
  IF ( ALLOCATED (b0f_pt) )           DEALLOCATE(b0f_pt) 
  IF ( ALLOCATED (b0f_ec_pt) )        DEALLOCATE(b0f_ec_pt) 
  IF ( ALLOCATED (b01f_pt) )          DEALLOCATE(b01f_pt) 
  IF ( ALLOCATED (b02f_pt) )          DEALLOCATE(b02f_pt) 
  IF ( ALLOCATED (eminf_pt) )         DEALLOCATE(eminf_pt) 

  IF ( ALLOCATED (celldmf_pt) )       DEALLOCATE(celldmf_pt) 
  IF ( ALLOCATED (densityf_pt) )      DEALLOCATE(densityf_pt) 
  IF ( ALLOCATED (enerf_pt) )         DEALLOCATE(enerf_pt) 
  IF ( ALLOCATED (free_enerf_pt) )    DEALLOCATE(free_enerf_pt) 
  IF ( ALLOCATED (entrf_pt) )         DEALLOCATE(entrf_pt) 
  IF ( ALLOCATED (cef_pt) )           DEALLOCATE(cef_pt) 
  IF ( ALLOCATED (cvf_pt) )           DEALLOCATE(cvf_pt) 
  IF ( ALLOCATED (cpf_pt) )           DEALLOCATE(cpf_pt) 
  IF ( ALLOCATED (b0f_s_pt) )         DEALLOCATE(b0f_s_pt) 
  IF ( ALLOCATED (betaf_pt) )         DEALLOCATE(betaf_pt) 
  IF ( ALLOCATED (gammaf_pt) )        DEALLOCATE(gammaf_pt) 
  IF ( ALLOCATED (densityf_pt) )      DEALLOCATE(densityf_pt) 
  IF ( ALLOCATED (uintf_pt) )         DEALLOCATE(uintf_pt) 
  IF ( ALLOCATED (tauf_pt) )          DEALLOCATE(tauf_pt) 
  IF ( ALLOCATED (alphaf_int_pt) )    DEALLOCATE(alphaf_int_pt) 
  IF ( ALLOCATED (uintf_zsisa_pt) )   DEALLOCATE(uintf_zsisa_pt) 
  IF ( ALLOCATED (tauf_zsisa_pt) )    DEALLOCATE(tauf_zsisa_pt) 
  IF ( ALLOCATED (alphaf_int_zsisa_pt) ) DEALLOCATE(alphaf_int_zsisa_pt) 

  IF ( ALLOCATED (alphaf_anis_pt) )   DEALLOCATE(alphaf_anis_pt)
  IF ( ALLOCATED (cpmcef_anis_pt) )   DEALLOCATE(cpmcef_anis_pt)
  IF ( ALLOCATED (csmctf_pt) )        DEALLOCATE(csmctf_pt)
  IF ( ALLOCATED (bthsf_pt) )         DEALLOCATE(bthsf_pt)
  IF ( ALLOCATED (ggammaf_pt) )       DEALLOCATE(ggammaf_pt)
  IF ( ALLOCATED (el_consf_pt) )      DEALLOCATE(el_consf_pt)
  IF ( ALLOCATED (el_compf_pt) )      DEALLOCATE(el_compf_pt)
  IF ( ALLOCATED (el_consf_s_pt) )    DEALLOCATE(el_consf_s_pt)
  IF ( ALLOCATED (el_compf_s_pt) )    DEALLOCATE(el_compf_s_pt)
  IF ( ALLOCATED (macro_elf_pt) )     DEALLOCATE(macro_elf_pt)
  IF ( ALLOCATED (macro_elf_s_pt) )   DEALLOCATE(macro_elf_s_pt)
  IF ( ALLOCATED (vf_pt) )            DEALLOCATE(vf_pt)
  IF ( ALLOCATED (vf_s_pt) )          DEALLOCATE(vf_s_pt)
  IF ( ALLOCATED (debye_macro_elf_pt) ) DEALLOCATE(debye_macro_elf_pt)
  IF ( ALLOCATED (debye_macro_elf_s_pt) ) DEALLOCATE(debye_macro_elf_s_pt)
  IF ( ALLOCATED (e_piezo_tensorf_pt) ) DEALLOCATE(e_piezo_tensorf_pt)
  IF ( ALLOCATED (d_piezo_tensorf_pt) ) DEALLOCATE(d_piezo_tensorf_pt)

  IF ( ALLOCATED (emp_free_enerf_pt) ) DEALLOCATE(emp_free_enerf_pt)
  IF ( ALLOCATED (emp_enerf_pt) )      DEALLOCATE(emp_enerf_pt)
  IF ( ALLOCATED (emp_entrf_pt) )      DEALLOCATE(emp_entrf_pt)
  IF ( ALLOCATED (emp_cef_pt) )        DEALLOCATE(emp_cef_pt)
  IF ( ALLOCATED (emp_cef_ptt) )       DEALLOCATE(emp_cef_ptt)

  IF (ALLOCATED(omega_p)) DEALLOCATE(omega_p)
  IF (ALLOCATED(density_p)) DEALLOCATE(density_p)
  IF (ALLOCATED(celldm_p)) DEALLOCATE(celldm_p)

  IF (ALLOCATED(p2_p)) THEN
     DO ipress=1,npress
        CALL clean_poly(p2_p(ipress))
     ENDDO
     DEALLOCATE(p2_p)
  ENDIF
  IF (ALLOCATED(p4_p)) THEN
     DO ipress=1,npress
        IF (lquartic) CALL clean_poly(p4_p(ipress))
     ENDDO
     DEALLOCATE(p4_p)
  ENDIF
  IF (ALLOCATED(p2_eq)) THEN
     DO ivar=1,nint_var
        CALL clean_poly(p2_eq(ivar))
     ENDDO
     DEALLOCATE(p2_eq)
  ENDIF
  IF (ALLOCATED(p4_eq)) THEN
     DO ivar=1,nint_var
        CALL clean_poly(p4_eq(ivar))
     ENDDO
     DEALLOCATE(p4_eq)
  ENDIF

  IF ( ALLOCATED (vmin_p) )           DEALLOCATE(vmin_p) 
  IF ( ALLOCATED (p0) )               DEALLOCATE(p0) 
  IF ( ALLOCATED (b0_p) )             DEALLOCATE(b0_p) 
  IF ( ALLOCATED (b01_p) )            DEALLOCATE(b01_p) 
  IF ( ALLOCATED (b02_p) )            DEALLOCATE(b02_p) 
  IF ( ALLOCATED (emin_p) )           DEALLOCATE(emin_p) 

  IF ( ALLOCATED (e0) )               DEALLOCATE(e0)
  IF ( ALLOCATED (v0) )               DEALLOCATE(v0)

  IF ( ALLOCATED (vminf_t) )          DEALLOCATE(vminf_t) 
  IF ( ALLOCATED (b0f_t) )            DEALLOCATE(b0f_t) 
  IF ( ALLOCATED (b0f_ec_t) )         DEALLOCATE(b0f_ec_t) 
  IF ( ALLOCATED (b01f_t) )           DEALLOCATE(b01f_t) 
  IF ( ALLOCATED (b02f_t) )           DEALLOCATE(b02f_t) 
  IF ( ALLOCATED (free_e_minf_t) )    DEALLOCATE(free_e_minf_t) 
  IF ( ALLOCATED (af_t) )             DEALLOCATE(af_t) 
  IF ( ALLOCATED (p1tf_t) )           DEALLOCATE(p1tf_t)
  IF ( ALLOCATED (p2tf_t) )           DEALLOCATE(p2tf_t)
  IF ( ALLOCATED (p3tf_t) )           DEALLOCATE(p3tf_t)
  IF ( ALLOCATED (p4tf_t) )           DEALLOCATE(p4tf_t)

  IF ( ALLOCATED (enerf_t) )          DEALLOCATE(enerf_t)
  IF ( ALLOCATED (free_enerf_t) )     DEALLOCATE(free_enerf_t) 
  IF ( ALLOCATED (entropyf_t) )       DEALLOCATE(entropyf_t) 
  IF ( ALLOCATED (cef_t) )            DEALLOCATE(cef_t) 
  IF ( ALLOCATED (cvf_t) )            DEALLOCATE(cvf_t) 

  IF ( ALLOCATED (cpf_t) )            DEALLOCATE(cpf_t) 
  IF ( ALLOCATED (b0f_s) )            DEALLOCATE(b0f_s) 
  IF ( ALLOCATED (b0f_ec_s) )         DEALLOCATE(b0f_ec_s) 
  IF ( ALLOCATED (alphaf_t) )         DEALLOCATE(alphaf_t) 
  IF ( ALLOCATED (betaf_t) )          DEALLOCATE(betaf_t) 
  IF ( ALLOCATED (gammaf_t) )         DEALLOCATE(gammaf_t) 
  IF ( ALLOCATED (uintf_t) )          DEALLOCATE(uintf_t)
  IF ( ALLOCATED (tauf_t) )           DEALLOCATE(tauf_t)
  IF ( ALLOCATED (alphaf_int_t) )     DEALLOCATE(alphaf_int_t)
  IF ( ALLOCATED (uintf_zsisa_t) )    DEALLOCATE(uintf_zsisa_t)
  IF ( ALLOCATED (tauf_zsisa_t) )     DEALLOCATE(tauf_zsisa_t)
  IF ( ALLOCATED (alphaf_int_zsisa_t) ) DEALLOCATE(alphaf_int_zsisa_t)

  IF ( ALLOCATED (celldmf_t) )        DEALLOCATE(celldmf_t) 
  IF ( ALLOCATED (densityf_t) )       DEALLOCATE(densityf_t) 
  IF ( ALLOCATED (alphaf_anis_t) )    DEALLOCATE(alphaf_anis_t) 
  IF ( ALLOCATED (cpmcef_anis) )      DEALLOCATE(cpmcef_anis) 
  IF ( ALLOCATED (bfactf_t) )         DEALLOCATE(bfactf_t) 
  IF ( ALLOCATED (csmctf_t) )         DEALLOCATE(csmctf_t) 
  IF ( ALLOCATED (bthsf_t) )          DEALLOCATE(bthsf_t) 
  IF ( ALLOCATED (ggammaf_t) )        DEALLOCATE(ggammaf_t) 

  IF ( ALLOCATED (el_consf_t) )       DEALLOCATE(el_consf_t)
  IF ( ALLOCATED (el_compf_t) )       DEALLOCATE(el_compf_t)
  IF ( ALLOCATED (el_consf_s) )       DEALLOCATE(el_consf_s)
  IF ( ALLOCATED (el_compf_s) )       DEALLOCATE(el_compf_s)
  IF ( ALLOCATED (macro_elf_s) )      DEALLOCATE(macro_elf_s)
  IF ( ALLOCATED (macro_elf_t) )      DEALLOCATE(macro_elf_t)
  IF ( ALLOCATED (vf_t) )             DEALLOCATE(vf_t)
  IF ( ALLOCATED (vf_s) )             DEALLOCATE(vf_s)
  IF ( ALLOCATED (debye_macro_elf_t) ) DEALLOCATE(debye_macro_elf_t)
  IF ( ALLOCATED (debye_macro_elf_s) ) DEALLOCATE(debye_macro_elf_s)
  IF ( ALLOCATED (el_conf_geo_t) )    DEALLOCATE(el_conf_geo_t) 
  IF ( ALLOCATED (dydef_t) )          DEALLOCATE(dydef_t)

  IF ( ALLOCATED (vminf_noe_t) )      DEALLOCATE(vminf_noe_t) 
  IF ( ALLOCATED (b0f_noe_t) )        DEALLOCATE(b0f_noe_t) 
  IF ( ALLOCATED (b01f_noe_t) )       DEALLOCATE(b01f_noe_t) 
  IF ( ALLOCATED (b02f_noe_t) )       DEALLOCATE(b02f_noe_t) 
  IF ( ALLOCATED (densityf_noe_t) )   DEALLOCATE(densityf_noe_t) 
  IF ( ALLOCATED (celldmf_noe_t) )    DEALLOCATE(celldmf_noe_t) 
  IF ( ALLOCATED (free_e_minf_noe_t) )DEALLOCATE(free_e_minf_noe_t) 
  IF ( ALLOCATED (free_enerf_noe_t) ) DEALLOCATE(free_enerf_noe_t) 
  IF ( ALLOCATED (enerf_noe_t) )      DEALLOCATE(enerf_noe_t) 
  IF ( ALLOCATED (entropyf_noe_t) )   DEALLOCATE(entropyf_noe_t) 
  IF ( ALLOCATED (cef_noe_t) )        DEALLOCATE(cef_noe_t) 
  IF ( ALLOCATED (cvf_noe_t) )        DEALLOCATE(cvf_noe_t) 

  IF ( ALLOCATED (vmine_t) )          DEALLOCATE(vmine_t) 
  IF ( ALLOCATED (b0e_t) )            DEALLOCATE(b0e_t) 
  IF ( ALLOCATED (b01e_t) )           DEALLOCATE(b01e_t) 
  IF ( ALLOCATED (b02e_t) )           DEALLOCATE(b02e_t) 
  IF ( ALLOCATED (free_e_mine_t) )    DEALLOCATE(free_e_mine_t) 

  IF ( ALLOCATED (betab) )            DEALLOCATE(betab) 
  IF ( ALLOCATED (alpha_an_g) )       DEALLOCATE(alpha_an_g) 
  IF ( ALLOCATED (poly_grun) )        DEALLOCATE(poly_grun) 
  IF ( ALLOCATED (poly_grun_red) )    DEALLOCATE(poly_grun_red) 
  IF ( ALLOCATED (vgrun_t) )          DEALLOCATE(vgrun_t) 
  IF ( ALLOCATED (b0_grun_t) )        DEALLOCATE(b0_grun_t) 
  IF ( ALLOCATED (el_cons_grun_t) )   DEALLOCATE(el_cons_grun_t) 
  IF ( ALLOCATED (el_comp_grun_t) )   DEALLOCATE(el_comp_grun_t) 
  IF ( ALLOCATED (celldm_grun_t) )    DEALLOCATE(celldm_grun_t) 
  IF ( ALLOCATED (cp_grun_t) )        DEALLOCATE(cp_grun_t) 
  IF ( ALLOCATED (cv_grun_t) )        DEALLOCATE(cv_grun_t) 
  IF ( ALLOCATED (ce_grun_t) )        DEALLOCATE(ce_grun_t) 
  IF ( ALLOCATED (b0_grun_s) )        DEALLOCATE(b0_grun_s) 
  IF ( ALLOCATED (grun_cpmce_anis) )  DEALLOCATE(grun_cpmce_anis) 
  IF ( ALLOCATED (grun_gamma_t) )     DEALLOCATE(grun_gamma_t) 
 
  IF ( ALLOCATED (el_ener_pt) )       DEALLOCATE(el_ener_pt)
  IF ( ALLOCATED (el_entr_pt) )       DEALLOCATE(el_entr_pt)
  IF ( ALLOCATED (el_ce_pt) )         DEALLOCATE(el_ce_pt)

  IF ( ALLOCATED (el_energy_t) )      DEALLOCATE(el_energy_t)
  IF ( ALLOCATED (el_free_energy_t) ) DEALLOCATE(el_free_energy_t)
  IF ( ALLOCATED (el_entropy_t) )     DEALLOCATE(el_entropy_t)
  IF ( ALLOCATED (el_ce_t) )          DEALLOCATE(el_ce_t)

  IF ( ALLOCATED (el_ce_ptt) )        DEALLOCATE(el_ce_ptt)
  IF ( ALLOCATED (el_beta_t ) )       DEALLOCATE(el_beta_t)
  IF ( ALLOCATED (el_cp_t ) )         DEALLOCATE(el_cp_t)
  IF ( ALLOCATED (el_gamma_t ) )      DEALLOCATE(el_gamma_t)

  IF ( ALLOCATED (el_energyf_t) )     DEALLOCATE(el_energyf_t)
  IF ( ALLOCATED (el_free_energyf_t) ) DEALLOCATE(el_free_energyf_t)
  IF ( ALLOCATED (el_entropyf_t) )    DEALLOCATE(el_entropyf_t)
  IF ( ALLOCATED (el_cef_t) )         DEALLOCATE(el_cef_t)
 
  IF ( ALLOCATED (el_betaf_t ) )      DEALLOCATE(el_betaf_t)
  IF ( ALLOCATED (el_cpf_t ) )        DEALLOCATE(el_cpf_t)
  IF ( ALLOCATED (el_gammaf_t ) )     DEALLOCATE(el_gammaf_t)

  IF ( ALLOCATED (el_enerf_pt) )      DEALLOCATE(el_enerf_pt)
  IF ( ALLOCATED (el_entrf_pt) )      DEALLOCATE(el_entrf_pt)
  IF ( ALLOCATED (el_cef_pt) )        DEALLOCATE(el_cef_pt)
  IF ( ALLOCATED (el_cef_ptt) )       DEALLOCATE(el_cef_ptt)

  IF ( ALLOCATED (xqaux) )            DEALLOCATE(xqaux)
  IF ( ALLOCATED (wqaux) )            DEALLOCATE(wqaux)
  IF ( ALLOCATED (letter) )           DEALLOCATE(letter)
  IF ( ALLOCATED (letter_path) )      DEALLOCATE(letter_path)
  IF ( ALLOCATED (label_list) )       DEALLOCATE(label_list)
  IF ( ALLOCATED (label_disp_q) )     DEALLOCATE(label_disp_q)
  IF ( ALLOCATED (disp_q) )           DEALLOCATE(disp_q)
  IF ( ALLOCATED (disp_wq) )          DEALLOCATE(disp_wq)

  CALL clean_poly(p2)
  IF ( ALLOCATED (hessian_v) )       DEALLOCATE( hessian_v )
  IF ( ALLOCATED (hessian_e) )       DEALLOCATE( hessian_e )
  IF ( ALLOCATED (x_pos_min) )       DEALLOCATE( x_pos_min )
  IF ( ALLOCATED (hessian4_v) )      DEALLOCATE( hessian4_v )
  IF ( ALLOCATED (hessian4_e) )      DEALLOCATE( hessian4_e )
  IF ( ALLOCATED (x_min_4) )         DEALLOCATE( x_min_4 )
  CALL clean_poly(p4)

  IF ( ALLOCATED (temp) )            DEALLOCATE(temp)
  IF ( ALLOCATED (temp_plot) )       DEALLOCATE(temp_plot)
  IF ( ALLOCATED (itemp_plot) )      DEALLOCATE(itemp_plot)
  IF ( ALLOCATED (uint0) )           DEALLOCATE(uint0)
  IF ( ALLOCATED (uint_eq) )         DEALLOCATE(uint_eq)
  IF ( ALLOCATED (tau_eq) )          DEALLOCATE(tau_eq)
  IF ( ALLOCATED (press_plot) )      DEALLOCATE(press_plot)
  IF ( ALLOCATED (ipress_plot) )     DEALLOCATE(ipress_plot)
  IF ( ALLOCATED (ivol_plot) )       DEALLOCATE(ivol_plot)
  IF ( ALLOCATED (sigma_ry) )        DEALLOCATE(sigma_ry)

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
  IF ( ALLOCATED (tau_acc) )         DEALLOCATE(tau_acc)
  IF ( ALLOCATED (min_y) )           DEALLOCATE(min_y)
  IF ( ALLOCATED (min_y_t) )         DEALLOCATE(min_y_t)
  IF ( ALLOCATED (dyde) )            DEALLOCATE(dyde)
  IF ( ALLOCATED (epsil_y) )         DEALLOCATE(epsil_y)
  IF ( ALLOCATED (el_con_geo) )      DEALLOCATE(el_con_geo)

  IF ( ALLOCATED (tau_save) )        DEALLOCATE(tau_save)
  IF ( ALLOCATED (tau_save_crys) )   DEALLOCATE(tau_save_crys)
  IF ( ALLOCATED (ityp_save) )       DEALLOCATE(ityp_save)
  IF ( ALLOCATED (atm_save) )        DEALLOCATE(atm_save)

  IF ( ALLOCATED (tau0) )            DEALLOCATE(tau0)
  IF ( ALLOCATED (tau0_crys) )       DEALLOCATE(tau0_crys)

  IF (ALLOCATED (xngeo))             DEALLOCATE(xngeo)
  IF (ALLOCATED (ind_rec3))          DEALLOCATE(ind_rec3)

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
  IF ( ALLOCATED (el_con_at_geo) ) DEALLOCATE( el_con_at_geo ) 
  IF ( ALLOCATED (el_con_tau_crys_geo) ) DEALLOCATE( el_con_tau_crys_geo )  
  IF ( ALLOCATED (el_con_tau_geo) )   DEALLOCATE( el_con_tau_geo )  
  IF ( ALLOCATED (el_con_omega_geo) ) DEALLOCATE( el_con_omega_geo )  

  IF ( ALLOCATED (p1tf_noe_t) )        DEALLOCATE(p1tf_noe_t)
  IF ( ALLOCATED (p2tf_noe_t) )        DEALLOCATE(p2tf_noe_t)
  IF ( ALLOCATED (p3tf_noe_t) )        DEALLOCATE(p3tf_noe_t)
  IF ( ALLOCATED (p4tf_noe_t) )        DEALLOCATE(p4tf_noe_t)

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
  CALL deallocate_geometry_file()
  CALL deallocate_piezo()

  RETURN
  !
END SUBROUTINE deallocate_thermo
