!
! Copyright (C) 2013 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE allocate_thermodynamics()
  !-----------------------------------------------------------------------
  !
  !  This routine deallocates the variables that control the thermo calculation
  !
  USE kinds,          ONLY : DP
  USE ions_base,      ONLY : nat
  USE thermo_mod,     ONLY : tot_ngeo, tot_ngeo_eos
  USE temperature,    ONLY : ntemp
  USE thermodynamics, ONLY : ph_free_ener, ph_ener, ph_entropy, ph_e0, &
                             ph_t_debye, ph_ce, ph_b_fact, ph_e0_eos,       &
                             ph_free_ener_eos, ph_ener_eos, ph_entropy_eos, &
                             ph_ce_eos, uint_geo_eos_t, tau_geo_eos_t,      &
                             e_piezo_tensor_eos_t
  USE ph_freq_thermodynamics, ONLY : phf_free_ener, phf_ener, phf_entropy, &
                                     phf_e0, phf_ce, phf_b_fact,           &
                                     ph_freq_save, phf_t_debye,            &
                                     phf_free_ener_eos, phf_ener_eos,      &
                                     phf_entropy_eos, phf_ce_eos,          &  
                                     phf_e0_eos,uintf_geo_eos_t,           &
                                     tauf_geo_eos_t, e_piezo_tensorf_eos_t
  USE el_thermodynamics, ONLY : el_ener, el_free_ener, el_entr, el_mu,     &
                           el_ce, el_ener_eos, el_free_ener_eos, el_entr_eos, &
                           el_mu_eos, el_ce_eos, elf_ener_eos, &
                           elf_free_ener_eos, elf_entr_eos, elf_mu_eos, &
                           elf_ce_eos
  USE el_anharmonic,  ONLY : vmine_t, b0e_t, b01e_t, b02e_t, free_e_mine_t

  USE control_emp_free_ener, ONLY : emp_ener, emp_free_ener, emp_entr, &
                           emp_ce
  USE control_atomic_pos, ONLY : nint_var

  IMPLICIT NONE

  IF (.NOT.ALLOCATED(ph_free_ener))  ALLOCATE(ph_free_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(ph_ener))       ALLOCATE(ph_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(ph_entropy))    ALLOCATE(ph_entropy(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(ph_t_debye))    ALLOCATE(ph_t_debye(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(ph_e0))         ALLOCATE(ph_e0(tot_ngeo))
  IF (.NOT.ALLOCATED(ph_ce))         ALLOCATE(ph_ce(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(ph_b_fact))     ALLOCATE(ph_b_fact(3,3,nat,ntemp,tot_ngeo))

  IF (.NOT.ALLOCATED(ph_free_ener_eos)) &
                                ALLOCATE(ph_free_ener_eos(ntemp,tot_ngeo_eos))
  IF (.NOT.ALLOCATED(ph_e0_eos)) ALLOCATE(ph_e0_eos(tot_ngeo_eos))
  IF (.NOT.ALLOCATED(ph_ener_eos)) ALLOCATE(ph_ener_eos(ntemp,tot_ngeo_eos))
  IF (.NOT.ALLOCATED(ph_entropy_eos)) &
                                 ALLOCATE(ph_entropy_eos(ntemp,tot_ngeo_eos))
  IF (.NOT.ALLOCATED(ph_ce_eos)) ALLOCATE(ph_ce_eos(ntemp,tot_ngeo_eos))
  IF (.NOT.ALLOCATED(uint_geo_eos_t)) &
                       ALLOCATE(uint_geo_eos_t(nint_var,ntemp,tot_ngeo_eos))
  IF (.NOT.ALLOCATED(tau_geo_eos_t)) &
                       ALLOCATE(tau_geo_eos_t(3,nat,ntemp,tot_ngeo_eos))
  IF (.NOT.ALLOCATED(e_piezo_tensor_eos_t)) &
                       ALLOCATE(e_piezo_tensor_eos_t(3,6,ntemp,&
                                                         tot_ngeo_eos))

  IF (.NOT.ALLOCATED(phf_free_ener)) ALLOCATE(phf_free_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(phf_ener))      ALLOCATE(phf_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(phf_entropy))   ALLOCATE(phf_entropy(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(phf_t_debye))   ALLOCATE(phf_t_debye(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(phf_e0))        ALLOCATE(phf_e0(tot_ngeo))
  IF (.NOT.ALLOCATED(phf_ce))        ALLOCATE(phf_ce(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(phf_b_fact))    ALLOCATE(phf_b_fact(3,3,nat,ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(phf_free_ener_eos)) &
                             ALLOCATE(phf_free_ener_eos(ntemp,tot_ngeo_eos))
  IF (.NOT.ALLOCATED(phf_ener_eos)) ALLOCATE(phf_ener_eos(ntemp,tot_ngeo_eos))
  IF (.NOT.ALLOCATED(phf_entropy_eos)) &
                                  ALLOCATE(phf_entropy_eos(ntemp,tot_ngeo_eos))
  IF (.NOT.ALLOCATED(phf_e0_eos)) ALLOCATE(phf_e0_eos(tot_ngeo_eos))
  IF (.NOT.ALLOCATED(phf_ce_eos)) ALLOCATE(phf_ce_eos(ntemp,tot_ngeo_eos))
  
  IF (.NOT.ALLOCATED(uintf_geo_eos_t)) &
                        ALLOCATE(uintf_geo_eos_t(nint_var,ntemp,tot_ngeo_eos))
  IF (.NOT.ALLOCATED(tauf_geo_eos_t)) &
                        ALLOCATE(tauf_geo_eos_t(3,nat,ntemp,tot_ngeo_eos))
  IF (.NOT.ALLOCATED(e_piezo_tensorf_eos_t)) &
                       ALLOCATE(e_piezo_tensorf_eos_t(3,6,ntemp,&
                                                         tot_ngeo_eos))
!
!   allocate the structures needed to save the frequencies
!
  IF (.NOT.ALLOCATED(ph_freq_save))   ALLOCATE(ph_freq_save(tot_ngeo))
!
!   allocate the variables for the electronic thermodynamic
!
  IF (.NOT.ALLOCATED(el_free_ener))      ALLOCATE(el_free_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(el_ener))           ALLOCATE(el_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(el_entr))           ALLOCATE(el_entr(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(el_mu))             ALLOCATE(el_mu(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(el_ce))             ALLOCATE(el_ce(ntemp,tot_ngeo))

  IF (.NOT.ALLOCATED(el_free_ener_eos))  ALLOCATE(el_free_ener_eos(ntemp,&
                                                          tot_ngeo_eos))
  IF (.NOT.ALLOCATED(el_ener_eos))       ALLOCATE(el_ener_eos(ntemp, &
                                                          tot_ngeo_eos))
  IF (.NOT.ALLOCATED(el_entr_eos))       ALLOCATE(el_entr_eos(ntemp,&
                                                          tot_ngeo_eos))
  IF (.NOT.ALLOCATED(el_mu_eos))         ALLOCATE(el_mu_eos(ntemp,&
                                                          tot_ngeo_eos))
  IF (.NOT.ALLOCATED(el_ce_eos))         ALLOCATE(el_ce_eos(ntemp,&
                                                            tot_ngeo_eos))
  IF (.NOT.ALLOCATED(elf_free_ener_eos))  ALLOCATE(elf_free_ener_eos(ntemp,&
                                                          tot_ngeo_eos))
  IF (.NOT.ALLOCATED(elf_ener_eos))       ALLOCATE(elf_ener_eos(ntemp, &
                                                          tot_ngeo_eos))
  IF (.NOT.ALLOCATED(elf_entr_eos))       ALLOCATE(elf_entr_eos(ntemp,&
                                                          tot_ngeo_eos))
  IF (.NOT.ALLOCATED(elf_mu_eos))         ALLOCATE(elf_mu_eos(ntemp,&
                                                          tot_ngeo_eos))
  IF (.NOT.ALLOCATED(elf_ce_eos))         ALLOCATE(elf_ce_eos(ntemp,&
                                                            tot_ngeo_eos))
!
!  allocate the variables for the anharmonic properties calculated adding
!  the electronic free energy
!
  IF (.NOT. ALLOCATED (vmine_t) )        ALLOCATE(vmine_t(ntemp))
  IF (.NOT. ALLOCATED (b0e_t) )          ALLOCATE(b0e_t(ntemp))
  IF (.NOT. ALLOCATED (b01e_t) )         ALLOCATE(b01e_t(ntemp))
  IF (.NOT. ALLOCATED (b02e_t) )         ALLOCATE(b02e_t(ntemp))
  IF (.NOT. ALLOCATED (free_e_mine_t) )  ALLOCATE(free_e_mine_t(ntemp))
!
!  Allocate the variables for the empirical free energy
!
  IF (.NOT.ALLOCATED(emp_free_ener))   ALLOCATE(emp_free_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(emp_ener))        ALLOCATE(emp_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(emp_entr))        ALLOCATE(emp_entr(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(emp_ce))          ALLOCATE(emp_ce(ntemp,tot_ngeo))


  RETURN
  !
END SUBROUTINE allocate_thermodynamics
!
!-------------------------------------------------------------------------
SUBROUTINE allocate_el_thermodynamics(tot_ngeo)
  !-----------------------------------------------------------------------
  !
  !  This routine deallocates the variables that control the thermo calculation
  !
  USE kinds,          ONLY : DP
  USE temperature,    ONLY : ntemp
  USE el_thermodynamics, ONLY : el_ener, el_free_ener, el_entr, el_mu, &
                           el_ce
  USE el_anharmonic,  ONLY : vmine_t, b0e_t, b01e_t, b02e_t, free_e_mine_t
  IMPLICIT NONE

INTEGER, INTENT(IN) :: tot_ngeo
!
!   allocate the variables for the electronic thermodynamic
!
  IF (.NOT.ALLOCATED(el_free_ener))      ALLOCATE(el_free_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(el_ener))           ALLOCATE(el_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(el_entr))           ALLOCATE(el_entr(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(el_mu))             ALLOCATE(el_mu(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(el_ce))             ALLOCATE(el_ce(ntemp,tot_ngeo))
!
!  allocate the variables for the anharmonic properties calculated adding
!  the electronic free energy
!
  IF (.NOT. ALLOCATED (vmine_t) )        ALLOCATE(vmine_t(ntemp))
  IF (.NOT. ALLOCATED (b0e_t) )          ALLOCATE(b0e_t(ntemp))
  IF (.NOT. ALLOCATED (b01e_t) )         ALLOCATE(b01e_t(ntemp))
  IF (.NOT. ALLOCATED (b02e_t) )         ALLOCATE(b02e_t(ntemp))
  IF (.NOT. ALLOCATED (free_e_mine_t) )  ALLOCATE(free_e_mine_t(ntemp))

  RETURN
  !
END SUBROUTINE allocate_el_thermodynamics
!-------------------------------------------------------------------------
SUBROUTINE allocate_anharmonic()
!-------------------------------------------------------------------------

  USE thermo_mod,          ONLY : tot_ngeo
  USE temperature,         ONLY : ntemp, ntemp_plot
  USE ions_base,           ONLY : nat
  USE anharmonic,          ONLY : vmin_t, b0_t, b01_t, b02_t, free_e_min_t,  &
                                  a_t, ener_t, free_ener_t, entropy_t, ce_t, &
                                  cv_t, cp_t, b0_s, alpha_t, beta_t, gamma_t,&
                                  celldm_t, alpha_anis_t, cpmce_anis,        &
                                  bfact_t, bths_t, ggamma_t, el_cons_t,      &
                                  e_piezo_tensor_t, d_piezo_tensor_t,        &
                                  epsilon_infty_t, zeu_t,                    &
                                  el_comp_t, el_cons_s, el_comp_s,           &
                                  macro_el_t, macro_el_s, v_s, v_t,          &
                                  el_con_geo_t, vmin_noe_t, b0_noe_t,        &
                                  b01_noe_t, b02_noe_t, free_e_min_noe_t,    &
                                  a_noe_t, beta_noe_t, cv_noe_t, ce_noe_t,   &
                                  cp_noe_t, b0_noe_s, gamma_noe_t, p1t_t,    &
                                  p2t_t, p3t_t, p4t_t, density_t, csmct_t,   &
                                  celldm_t_p1, celldm_t_m1, free_ener_noe_t, &
                                  ener_noe_t, entropy_noe_t, celldm_noe_t,   &
                                  celldm_noe_t_p1, celldm_noe_t_m1,          &
                                  p1t_noe_t, p2t_noe_t, p3t_noe_t, p4t_noe_t, &
                                  density_noe_t, b0_ec_t, b0_ec_s, &
                                  debye_macro_el_t, debye_macro_el_s, &
                                  dyde_t, uint_t, tau_t, uint_zsisa_t, &
                                  tau_zsisa_t, alpha_int_t, alpha_int_zsisa_t
  USE anharmonic_pt,       ONLY : vmin_pt, b0_pt, b01_pt, b02_pt, emin_pt,   &
                                  ce_pt, cv_pt, cp_pt, beta_pt, b0_s_pt,     &
                                  gamma_pt, free_ener_pt, ener_pt, entr_pt,  &
                                  celldm_pt, alpha_anis_pt, density_pt,      &
                                  cpmce_anis_pt, bths_pt, ggamma_pt,         &
                                  csmct_pt, el_cons_pt, el_cons_s_pt,        &
                                  el_comp_pt, el_comp_s_pt,                  &
                                  macro_el_pt, macro_el_s_pt, v_pt,          &
                                  e_piezo_tensor_pt, d_piezo_tensor_pt,      &
                                  epsilon_infty_pt, zeu_pt,                  &
                                  v_s_pt, celldm_pt_p1, celldm_pt_m1,        &
                                  b0_ec_pt, debye_macro_el_pt,               &
                                  debye_macro_el_s_pt, uint_pt, tau_pt,      &
                                  uint_zsisa_pt, tau_zsisa_pt, alpha_int_pt, &
                                  alpha_int_zsisa_pt
  USE anharmonic_ptt,      ONLY : vmin_ptt, b0_ptt, b01_ptt, b02_ptt,        &
                                  emin_ptt, ce_ptt, cv_ptt, cp_ptt,          &
                                  beta_ptt, b0_s_ptt, gamma_ptt, celldm_ptt, &
                                  el_cons_ptt, el_cons_s_ptt,                &
                                  el_comp_ptt, el_comp_s_ptt,                &
                                  macro_el_ptt, macro_el_s_ptt, v_ptt,       &
                                  v_s_ptt, e_piezo_tensor_ptt,               &
                                  d_piezo_tensor_ptt,                        &
                                  epsilon_infty_ptt, zeu_ptt,                &
                                  celldm_ptt_p1, celldm_ptt_m1,              &
                                  emin_ptt_p1, emin_ptt_m1, vmin_ptt_p1,     &
                                  vmin_ptt_m1, alpha_anis_ptt, density_ptt,  &
                                  cpmce_anis_ptt, bths_ptt, ggamma_ptt,      &
                                  csmct_ptt, ener_ptt, entr_ptt, b0_ec_ptt,  &
                                  debye_macro_el_ptt, debye_macro_el_s_ptt,  &
                                  uint_ptt, tau_ptt, uint_zsisa_ptt,         &
                                  tau_zsisa_ptt, alpha_int_ptt,              &
                                  alpha_int_zsisa_ptt, uint_ptt_p1,          &
                                  uint_ptt_m1, tau_ptt_p1, tau_ptt_m1
  USE ph_freq_anharmonic,  ONLY : vminf_t, b0f_t, b01f_t, b02f_t,            &
                                  free_e_minf_t, af_t,                       &
                                  enerf_t, free_enerf_t, entropyf_t,         &
                                  cef_t, cvf_t, cpf_t, b0f_s, alphaf_t,      &
                                  betaf_t, gammaf_t, celldmf_t,              &
                                  alphaf_anis_t, cpmcef_anis,                &
                                  bfactf_t, bthsf_t, ggammaf_t,              &
                                  el_consf_t, el_compf_t, el_consf_s,        &
                                  el_compf_s, macro_elf_t, macro_elf_s,      &
                                  e_piezo_tensorf_t, d_piezo_tensorf_t,      &
                                  epsilon_inftyf_t,  zeuf_t,                  &
                                  vf_t, vf_s, el_conf_geo_t, densityf_t,     &
                                  csmctf_t, p1tf_t, p2tf_t, p3tf_t, p4tf_t,  &
                                  p1tf_noe_t, p2tf_noe_t, p3tf_noe_t,        &
                                  p4tf_noe_t, celldmf_noe_t,                 &
                                  free_e_minf_noe_t, densityf_noe_t,         &
                                  vminf_noe_t, celldmf_t_p1, celldmf_t_m1,   &
                                  celldmf_noe_t_p1, celldmf_noe_t_m1,        &
                                  b0f_noe_t, b01f_noe_t, b02f_noe_t,         &
                                  free_enerf_noe_t, enerf_noe_t,             &
                                  entropyf_noe_t, cef_noe_t, cvf_noe_t,      &
                                  betaf_noe_t, cpf_noe_t, b0f_noe_s,         &
                                  gammaf_noe_t, af_noe_t, b0f_ec_t,          &
                                  b0f_ec_s, debye_macro_elf_t,               &
                                  debye_macro_elf_s, dydef_t, uintf_t,       &
                                  tauf_t, uintf_zsisa_t, tauf_zsisa_t,       & 
                                  alphaf_int_t, alphaf_int_zsisa_t
  USE ph_freq_anharmonic_pt, ONLY : vminf_pt, b0f_pt, b01f_pt, b02f_pt,      &
                                  eminf_pt, densityf_pt, celldmf_pt,         &
                                  enerf_pt, free_enerf_pt, entrf_pt,         &
                                  cef_pt, cvf_pt, cpf_pt, b0f_s_pt,          &
                                  betaf_pt, gammaf_pt, celldmf_pt_p1,        &
                                  celldmf_pt_m1, alphaf_anis_pt,             &
                                  el_consf_pt, el_compf_pt, cpmcef_anis_pt,  &
                                  el_consf_s_pt, el_compf_s_pt,              &
                                  e_piezo_tensorf_pt, d_piezo_tensorf_pt,    &
                                  epsilon_inftyf_pt, zeuf_pt,                &
                                  bthsf_pt, ggammaf_pt, csmctf_pt,           &
                                  macro_elf_pt, macro_elf_s_pt, vf_pt,       &
                                  vf_s_pt, b0f_ec_pt,                        &
                                  debye_macro_elf_pt, debye_macro_elf_s_pt,  &
                                  uintf_pt, tauf_pt, uintf_zsisa_pt,         &
                                  tauf_zsisa_pt, alphaf_int_pt,              &
                                  alphaf_int_zsisa_pt
  USE ph_freq_anharmonic_ptt,  ONLY : vminf_ptt, b0f_ptt, b01f_ptt,          &
                                  b02f_ptt,eminf_ptt, celldmf_ptt,           & 
                                  celldmf_ptt_p1, celldmf_ptt_m1,            &
                                  eminf_ptt_p1, eminf_ptt_m1, vminf_ptt_p1,  &
                                  vminf_ptt_m1, densityf_ptt,                &
                                  enerf_ptt, entrf_ptt, cef_ptt,             &
                                  cvf_ptt, cpf_ptt, b0f_s_ptt, betaf_ptt,    &
                                  gammaf_ptt, alphaf_anis_ptt,                &
                                  el_consf_ptt, el_compf_ptt, cpmcef_anis_ptt,&
                                  el_consf_s_ptt, el_compf_s_ptt,            &
                                  bthsf_ptt, ggammaf_ptt, csmctf_ptt,        &
                                  macro_elf_ptt, macro_elf_s_ptt, vf_ptt,    &
                                  vf_s_ptt, b0f_ec_ptt,                      &
                                  e_piezo_tensorf_ptt, d_piezo_tensorf_ptt,  &
                                  epsilon_inftyf_ptt, zeuf_ptt,              &
                                  debye_macro_elf_ptt, debye_macro_elf_s_ptt,&
                                  uintf_ptt, tauf_ptt, uintf_zsisa_ptt,      &
                                  tauf_zsisa_ptt, alphaf_int_ptt,            &
                                  alphaf_int_zsisa_ptt, uintf_ptt_p1,        &
                                  uintf_ptt_m1, tauf_ptt_p1, tauf_ptt_m1
  USE anharmonic_vt,       ONLY : press_vt, press_vtt
  USE ph_freq_anharmonic_vt,       ONLY : pressf_vt, pressf_vtt
  USE grun_anharmonic,     ONLY : betab, alpha_an_g, cp_grun_t, cv_grun_t,   &
                                  ce_grun_t, b0_grun_s, grun_gamma_t,        &
                                  grun_cpmce_anis, el_cons_grun_t,           &
                                  el_comp_grun_t
  USE el_anharmonic, ONLY : el_energy_t, el_free_energy_t, el_entropy_t,     &
                            el_ce_t, el_energyf_t, el_free_energyf_t,        &
                            el_entropyf_t, el_cef_t, el_ce_pt, el_ce_ptt,    &
                            el_free_ener_pt, el_ener_pt, el_entr_pt,         &
                            el_b0_t, el_beta_t, el_betaf_t, el_cp_t,         &
                            el_cpf_t, el_gamma_t, el_gammaf_t,               &
                            el_enerf_pt, el_free_enerf_pt, el_entrf_pt,      &
                            el_cef_pt, el_cef_ptt
  USE emp_anharmonic, ONLY : emp_energy_t, emp_free_energy_t, emp_entropy_t, &
                            emp_ce_t, emp_energyf_t, emp_free_energyf_t,     &
                            emp_entropyf_t, emp_cef_t, emp_free_ener_pt,     &
                            emp_ener_pt, emp_entr_pt, emp_ce_pt, emp_ce_ptt, &
                            emp_free_ener_pt, emp_enerf_pt, emp_entrf_pt,    &
                            emp_cef_pt, emp_cef_ptt

  USE control_atomic_pos,  ONLY : nint_var
  USE control_quartic_energy, ONLY : poly_degree_ph
  USE control_vol,         ONLY : nvol_plot, nvol
  USE control_grun,        ONLY : vgrun_t, celldm_grun_t, b0_grun_t
  USE control_pressure,    ONLY : npress, npress_plot

  IMPLICIT NONE
  INTEGER :: m1

  m1=poly_degree_ph+1

  IF (.NOT. ALLOCATED (vmin_t) )        ALLOCATE(vmin_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0_t) )          ALLOCATE(b0_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0_ec_t) )       ALLOCATE(b0_ec_t(ntemp)) 
  IF (.NOT. ALLOCATED (b01_t) )         ALLOCATE(b01_t(ntemp)) 
  IF (.NOT. ALLOCATED (b02_t) )         ALLOCATE(b02_t(ntemp)) 
  IF (.NOT. ALLOCATED (free_e_min_t) )  ALLOCATE(free_e_min_t(ntemp)) 
  IF (.NOT. ALLOCATED (a_t) )           ALLOCATE(a_t(m1,ntemp)) 
  IF (.NOT. ALLOCATED (p2t_t) )         ALLOCATE(p2t_t(ntemp)) 
  IF (.NOT. ALLOCATED (p1t_t) )         ALLOCATE(p1t_t(ntemp)) 
  IF (.NOT. ALLOCATED (p3t_t) )         ALLOCATE(p3t_t(ntemp)) 
  IF (.NOT. ALLOCATED (p4t_t) )         ALLOCATE(p4t_t(ntemp)) 

  IF (.NOT. ALLOCATED (ener_t) )        ALLOCATE(ener_t(ntemp))
  IF (.NOT. ALLOCATED (free_ener_t) )   ALLOCATE(free_ener_t(ntemp)) 
  IF (.NOT. ALLOCATED (entropy_t) )     ALLOCATE(entropy_t(ntemp)) 
  IF (.NOT. ALLOCATED (ce_t) )          ALLOCATE(ce_t(ntemp)) 
  IF (.NOT. ALLOCATED (cv_t) )          ALLOCATE(cv_t(ntemp)) 

  IF (.NOT. ALLOCATED (cp_t) )          ALLOCATE(cp_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0_s) )          ALLOCATE(b0_s(ntemp)) 
  IF (.NOT. ALLOCATED (b0_ec_s) )       ALLOCATE(b0_ec_s(ntemp)) 
  IF (.NOT. ALLOCATED (alpha_t) )       ALLOCATE(alpha_t(ntemp)) 
  IF (.NOT. ALLOCATED (beta_t) )        ALLOCATE(beta_t(ntemp)) 
  IF (.NOT. ALLOCATED (gamma_t) )       ALLOCATE(gamma_t(ntemp)) 
  IF (.NOT. ALLOCATED (uint_t))         ALLOCATE(uint_t(nint_var,ntemp))
  IF (.NOT. ALLOCATED (tau_t))          ALLOCATE(tau_t(3,nat,ntemp))
  IF (.NOT. ALLOCATED (alpha_int_t))    ALLOCATE(alpha_int_t(nint_var,ntemp))
  IF (.NOT. ALLOCATED (uint_zsisa_t))   ALLOCATE(uint_zsisa_t(nint_var,ntemp))
  IF (.NOT. ALLOCATED (tau_zsisa_t))    ALLOCATE(tau_zsisa_t(3,nat,ntemp))
  IF (.NOT. ALLOCATED (alpha_int_zsisa_t)) ALLOCATE(alpha_int_zsisa_t&
                                                         (nint_var,ntemp))
  IF (.NOT. ALLOCATED (celldm_t) )      ALLOCATE(celldm_t(6,ntemp)) 
  IF (.NOT. ALLOCATED (density_t) )     ALLOCATE(density_t(ntemp)) 
  IF (.NOT. ALLOCATED (alpha_anis_t) )  ALLOCATE(alpha_anis_t(6,ntemp)) 
  IF (.NOT. ALLOCATED (cpmce_anis) )    ALLOCATE(cpmce_anis(ntemp)) 
  IF (.NOT. ALLOCATED (bfact_t))        ALLOCATE(bfact_t(6,nat,ntemp))
  IF (.NOT. ALLOCATED (csmct_t) )       ALLOCATE(csmct_t(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (bths_t) )        ALLOCATE(bths_t(3,3,ntemp)) 
  IF (.NOT. ALLOCATED (ggamma_t) )      ALLOCATE(ggamma_t(3,3,ntemp)) 

  IF (.NOT. ALLOCATED (el_cons_t) )     ALLOCATE(el_cons_t(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (el_comp_t) )     ALLOCATE(el_comp_t(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (el_cons_s) )     ALLOCATE(el_cons_s(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (el_comp_s) )     ALLOCATE(el_comp_s(6,6,ntemp))
  IF (.NOT. ALLOCATED (macro_el_s) )    ALLOCATE(macro_el_s(8,ntemp))
  IF (.NOT. ALLOCATED (macro_el_t) )    ALLOCATE(macro_el_t(8,ntemp)) 
  IF (.NOT. ALLOCATED (debye_macro_el_t) ) ALLOCATE(debye_macro_el_t(ntemp)) 
  IF (.NOT. ALLOCATED (debye_macro_el_s) ) ALLOCATE(debye_macro_el_s(ntemp)) 
  IF (.NOT. ALLOCATED (v_s) )           ALLOCATE(v_s(3,ntemp))
  IF (.NOT. ALLOCATED (v_t) )           ALLOCATE(v_t(3,ntemp))
  IF (.NOT. ALLOCATED (el_con_geo_t) )  ALLOCATE(el_con_geo_t(6,6,ntemp,&
                                                                   tot_ngeo)) 
  IF (.NOT. ALLOCATED (e_piezo_tensor_t) ) ALLOCATE(e_piezo_tensor_t(3,6,ntemp)) 
  IF (.NOT. ALLOCATED (d_piezo_tensor_t) ) ALLOCATE(d_piezo_tensor_t(3,6,ntemp)) 
  IF (.NOT. ALLOCATED (epsilon_infty_t) ) ALLOCATE(epsilon_infty_t(3,3,ntemp)) 
  IF (.NOT. ALLOCATED (zeu_t) )           ALLOCATE(zeu_t(3,3,nat,ntemp)) 

  IF (.NOT. ALLOCATED (dyde_t) )        ALLOCATE(dyde_t(21,ntemp)) 
  IF (.NOT. ALLOCATED (celldm_t_p1) )   ALLOCATE(celldm_t_p1(6,ntemp)) 
  IF (.NOT. ALLOCATED (celldm_t_m1) )   ALLOCATE(celldm_t_m1(6,ntemp)) 
  IF (.NOT. ALLOCATED (celldm_noe_t_p1) ) ALLOCATE(celldm_noe_t_p1(6,ntemp)) 
  IF (.NOT. ALLOCATED (celldm_noe_t_m1) ) ALLOCATE(celldm_noe_t_m1(6,ntemp)) 

  IF (.NOT. ALLOCATED (vmin_noe_t) )     ALLOCATE(vmin_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (density_noe_t) )  ALLOCATE(density_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (celldm_noe_t) )   ALLOCATE(celldm_noe_t(6,ntemp)) 
  IF (.NOT. ALLOCATED (b0_noe_t) )       ALLOCATE(b0_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (b01_noe_t) )      ALLOCATE(b01_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (b02_noe_t) )      ALLOCATE(b02_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (free_e_min_noe_t) ) ALLOCATE(free_e_min_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (a_noe_t) )        ALLOCATE(a_noe_t(m1,ntemp)) 
  IF (.NOT. ALLOCATED (p2t_noe_t) )      ALLOCATE(p2t_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (p1t_noe_t) )      ALLOCATE(p1t_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (p3t_noe_t) )      ALLOCATE(p3t_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (p4t_noe_t) )      ALLOCATE(p4t_noe_t(ntemp)) 

  IF (.NOT. ALLOCATED (free_ener_noe_t) ) ALLOCATE(free_ener_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (ener_noe_t) )     ALLOCATE(ener_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (entropy_noe_t) )  ALLOCATE(entropy_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (ce_noe_t) )       ALLOCATE(ce_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (cv_noe_t) )       ALLOCATE(cv_noe_t(ntemp)) 

  IF (.NOT. ALLOCATED (cp_noe_t) )       ALLOCATE(cp_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0_noe_s) )       ALLOCATE(b0_noe_s(ntemp)) 
  IF (.NOT. ALLOCATED (beta_noe_t) )     ALLOCATE(beta_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (gamma_noe_t) )    ALLOCATE(gamma_noe_t(ntemp)) 

  IF (ntemp_plot>0) THEN
     IF (.NOT. ALLOCATED (vmin_ptt) )   ALLOCATE(vmin_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (vmin_ptt_p1) ) ALLOCATE(vmin_ptt_p1&
                                                         (npress,ntemp_plot))
     IF (.NOT. ALLOCATED (vmin_ptt_m1) ) ALLOCATE(vmin_ptt_m1&
                                                         (npress,ntemp_plot))
     IF (.NOT. ALLOCATED (b0_ptt) )     ALLOCATE(b0_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (b0_ec_ptt) )  ALLOCATE(b0_ec_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (b01_ptt) )    ALLOCATE(b01_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (b02_ptt) )    ALLOCATE(b02_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (emin_ptt) )   ALLOCATE(emin_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (emin_ptt_p1) ) ALLOCATE(emin_ptt_p1(npress,&
                                                                 ntemp_plot)) 
     IF (.NOT. ALLOCATED (emin_ptt_m1) )   ALLOCATE(emin_ptt_m1(npress,&
                                                                 ntemp_plot)) 
     IF (.NOT. ALLOCATED (celldm_ptt) ) ALLOCATE(celldm_ptt(6,npress,&
                                                                 ntemp_plot)) 
     IF (.NOT. ALLOCATED (celldm_ptt_p1) ) ALLOCATE(celldm_ptt_p1(6,npress,&
                                                                 ntemp_plot)) 
     IF (.NOT. ALLOCATED (celldm_ptt_m1) ) ALLOCATE(celldm_ptt_m1(6,npress,&
                                                                 ntemp_plot)) 
     IF (.NOT. ALLOCATED(uint_ptt)) &
                       ALLOCATE(uint_ptt(nint_var,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(tau_ptt)) &
                       ALLOCATE(tau_ptt(3,nat,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(uint_ptt_p1)) &
                       ALLOCATE(uint_ptt_p1(nint_var,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(tau_ptt_p1)) &
                       ALLOCATE(tau_ptt_p1(3,nat,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(uint_ptt_m1)) &
                       ALLOCATE(uint_ptt_m1(nint_var,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(tau_ptt_m1)) &
                       ALLOCATE(tau_ptt_m1(3,nat,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(alpha_int_ptt)) &
                       ALLOCATE(alpha_int_ptt(nint_var,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(uint_zsisa_ptt)) &
                       ALLOCATE(uint_zsisa_ptt(nint_var,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(tau_zsisa_ptt)) &
                       ALLOCATE(tau_zsisa_ptt(3,nat,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(alpha_int_zsisa_ptt)) &
                     ALLOCATE(alpha_int_zsisa_ptt(nint_var,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (ener_ptt) )   ALLOCATE(ener_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (entr_ptt) )   ALLOCATE(entr_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (ce_ptt) )     ALLOCATE(ce_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (cv_ptt) )     ALLOCATE(cv_ptt(npress,ntemp_plot)) 

     IF (.NOT. ALLOCATED (cp_ptt) )     ALLOCATE(cp_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (b0_s_ptt) )   ALLOCATE(b0_s_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (beta_ptt) )   ALLOCATE(beta_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (gamma_ptt) )  ALLOCATE(gamma_ptt(npress,ntemp_plot)) 

     IF (.NOT. ALLOCATED (press_vt) )   ALLOCATE(press_vt(nvol,ntemp_plot)) 

     IF (.NOT. ALLOCATED (vminf_ptt) )   ALLOCATE(vminf_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (vminf_ptt_p1) )   &
                                      ALLOCATE(vminf_ptt_p1(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (vminf_ptt_m1) )  &
                                      ALLOCATE(vminf_ptt_m1(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (b0f_ptt) )  ALLOCATE(b0f_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (b0f_ec_ptt) ) ALLOCATE(b0f_ec_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (b01f_ptt) )  ALLOCATE(b01f_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (b02f_ptt) )  ALLOCATE(b02f_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (celldmf_ptt) ) ALLOCATE(celldmf_ptt(6,npress, &
                                                                 ntemp_plot)) 
     IF (.NOT. ALLOCATED (eminf_ptt) )   ALLOCATE(eminf_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (eminf_ptt_p1) ) ALLOCATE(eminf_ptt_p1(npress,&
                                                                 ntemp_plot)) 
     IF (.NOT. ALLOCATED (eminf_ptt_m1) )   ALLOCATE(eminf_ptt_m1(npress,&
                                                                 ntemp_plot)) 
     IF (.NOT. ALLOCATED (celldmf_ptt_p1) ) ALLOCATE(celldmf_ptt_p1(6,npress,&
                                                                 ntemp_plot)) 
     IF (.NOT. ALLOCATED (celldmf_ptt_m1) ) ALLOCATE(celldmf_ptt_m1(6,npress,&
                                                                 ntemp_plot)) 
     IF (.NOT. ALLOCATED (enerf_ptt) ) ALLOCATE(enerf_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (entrf_ptt) ) ALLOCATE(entrf_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (cef_ptt) )   ALLOCATE(cef_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (cvf_ptt) )   ALLOCATE(cvf_ptt(npress,ntemp_plot)) 

     IF (.NOT. ALLOCATED (cpf_ptt) )   ALLOCATE(cpf_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (b0f_s_ptt) ) ALLOCATE(b0f_s_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (betaf_ptt) ) ALLOCATE(betaf_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (gammaf_ptt) ) ALLOCATE(gammaf_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (pressf_vt) )   ALLOCATE(pressf_vt(nvol,ntemp_plot)) 
  ENDIF

  IF (npress_plot>0) THEN
     IF (.NOT. ALLOCATED (vmin_pt) )    ALLOCATE(vmin_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (b0_pt) )      ALLOCATE(b0_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (b0_ec_pt) )   ALLOCATE(b0_ec_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (b01_pt) )     ALLOCATE(b01_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (b02_pt) )     ALLOCATE(b02_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (emin_pt) )    ALLOCATE(emin_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (celldm_pt) ) ALLOCATE(celldm_pt(6,ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (celldm_pt_p1) ) ALLOCATE(celldm_pt_p1(6,ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (celldm_pt_m1) ) ALLOCATE(celldm_pt_m1(6,ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED(uint_pt)) &
                       ALLOCATE(uint_pt(nint_var,ntemp,npress_plot))
     IF (.NOT. ALLOCATED(tau_pt)) &
                       ALLOCATE(tau_pt(3,nat,ntemp,npress_plot))
     IF (.NOT. ALLOCATED(alpha_int_pt)) &
                       ALLOCATE(alpha_int_pt(nint_var,ntemp,npress_plot))
     IF (.NOT. ALLOCATED(uint_zsisa_pt)) &
                       ALLOCATE(uint_zsisa_pt(nint_var,ntemp,npress_plot))
     IF (.NOT. ALLOCATED(tau_zsisa_pt)) &
                       ALLOCATE(tau_zsisa_pt(3,nat,ntemp,npress_plot))
     IF (.NOT. ALLOCATED(alpha_int_zsisa_pt)) &
                     ALLOCATE(alpha_int_zsisa_pt(nint_var,ntemp,npress_plot))

     IF (.NOT. ALLOCATED (free_ener_pt) ) ALLOCATE(free_ener_pt(ntemp,&
                                                               npress_plot)) 
     IF (.NOT. ALLOCATED (ener_pt) )    ALLOCATE(ener_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (entr_pt) )    ALLOCATE(entr_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (ce_pt) )      ALLOCATE(ce_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (cv_pt) )      ALLOCATE(cv_pt(ntemp,npress_plot)) 

     IF (.NOT. ALLOCATED (cp_pt) )      ALLOCATE(cp_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (beta_pt) )    ALLOCATE(beta_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (b0_s_pt) )    ALLOCATE(b0_s_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (gamma_pt) )   ALLOCATE(gamma_pt(ntemp,npress_plot)) 

     IF (.NOT. ALLOCATED (vminf_pt) )   ALLOCATE(vminf_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (b0f_pt) )     ALLOCATE(b0f_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (b0f_ec_pt) )  ALLOCATE(b0f_ec_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (b01f_pt) )    ALLOCATE(b01f_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (b02f_pt) )    ALLOCATE(b02f_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (celldmf_pt) ) &
                                    ALLOCATE(celldmf_pt(6,ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (eminf_pt) ) ALLOCATE(eminf_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (densityf_pt) ) ALLOCATE(densityf_pt(ntemp,npress_plot)) 

     IF (.NOT. ALLOCATED (enerf_pt) ) ALLOCATE(enerf_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (free_enerf_pt) ) ALLOCATE(free_enerf_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (entrf_pt) ) ALLOCATE(entrf_pt(ntemp, npress_plot)) 
     IF (.NOT. ALLOCATED (cef_pt) )   ALLOCATE(cef_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (cvf_pt) )   ALLOCATE(cvf_pt(ntemp,npress_plot)) 

     IF (.NOT. ALLOCATED (cpf_pt) )    ALLOCATE(cpf_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (b0f_s_pt) )  ALLOCATE(b0f_s_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (betaf_pt) )  ALLOCATE(betaf_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (gammaf_pt) ) ALLOCATE(gammaf_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (celldmf_pt_p1) ) &
                                 ALLOCATE(celldmf_pt_p1(6,ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (celldmf_pt_m1) ) &
                                 ALLOCATE(celldmf_pt_m1(6,ntemp,npress_plot)) 
  ENDIF

  IF (nvol_plot>0) THEN
     IF (.NOT. ALLOCATED (press_vtt) ) ALLOCATE(press_vtt(ntemp,nvol_plot)) 
     IF (.NOT. ALLOCATED (pressf_vtt) ) ALLOCATE(pressf_vtt(ntemp,nvol_plot)) 
  ENDIF

  IF (.NOT. ALLOCATED (vminf_t) )       ALLOCATE(vminf_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0f_t) )         ALLOCATE(b0f_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0f_ec_t) )      ALLOCATE(b0f_ec_t(ntemp)) 
  IF (.NOT. ALLOCATED (b01f_t) )        ALLOCATE(b01f_t(ntemp)) 
  IF (.NOT. ALLOCATED (b02f_t) )        ALLOCATE(b02f_t(ntemp)) 
  IF (.NOT. ALLOCATED (free_e_minf_t) ) ALLOCATE(free_e_minf_t(ntemp)) 
  IF (.NOT. ALLOCATED (af_t) )          ALLOCATE(af_t(m1,ntemp)) 
  IF (.NOT. ALLOCATED (p2tf_t) )        ALLOCATE(p2tf_t(ntemp)) 
  IF (.NOT. ALLOCATED (p1tf_t) )        ALLOCATE(p1tf_t(ntemp)) 
  IF (.NOT. ALLOCATED (p3tf_t) )        ALLOCATE(p3tf_t(ntemp)) 
  IF (.NOT. ALLOCATED (p4tf_t) )        ALLOCATE(p4tf_t(ntemp)) 

  IF (.NOT. ALLOCATED (enerf_t) )       ALLOCATE(enerf_t(ntemp)) 
  IF (.NOT. ALLOCATED (free_enerf_t) )  ALLOCATE(free_enerf_t(ntemp))
  IF (.NOT. ALLOCATED (entropyf_t) )    ALLOCATE(entropyf_t(ntemp)) 
  IF (.NOT. ALLOCATED (cef_t) )         ALLOCATE(cef_t(ntemp)) 
  IF (.NOT. ALLOCATED (cvf_t) )         ALLOCATE(cvf_t(ntemp)) 

  IF (.NOT. ALLOCATED (cpf_t) )         ALLOCATE(cpf_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0f_s) )         ALLOCATE(b0f_s(ntemp)) 
  IF (.NOT. ALLOCATED (b0f_ec_s) )      ALLOCATE(b0f_ec_s(ntemp)) 
  IF (.NOT. ALLOCATED (alphaf_t) )      ALLOCATE(alphaf_t(ntemp)) 
  IF (.NOT. ALLOCATED (betaf_t) )       ALLOCATE(betaf_t(ntemp)) 
  IF (.NOT. ALLOCATED (gammaf_t) )      ALLOCATE(gammaf_t(ntemp)) 
  IF (.NOT. ALLOCATED (uintf_t))        ALLOCATE(uintf_t(nint_var,ntemp))
  IF (.NOT. ALLOCATED (tauf_t))         ALLOCATE(tauf_t(3,nat,ntemp))
  IF (.NOT. ALLOCATED (alphaf_int_t))   ALLOCATE(alphaf_int_t(nint_var,ntemp))
  IF (.NOT. ALLOCATED (uintf_zsisa_t))  ALLOCATE(uintf_zsisa_t(nint_var,ntemp))
  IF (.NOT. ALLOCATED (tauf_zsisa_t))   ALLOCATE(tauf_zsisa_t(3,nat,ntemp))
  IF (.NOT. ALLOCATED (alphaf_int_zsisa_t)) ALLOCATE(alphaf_int_zsisa_t&
                                                          (nint_var,ntemp))

  IF (.NOT. ALLOCATED (celldmf_t) )     ALLOCATE(celldmf_t(6,ntemp)) 
  IF (.NOT. ALLOCATED (densityf_t) )    ALLOCATE(densityf_t(ntemp)) 
  IF (.NOT. ALLOCATED (alphaf_anis_t) ) ALLOCATE(alphaf_anis_t(6,ntemp)) 
  IF (.NOT. ALLOCATED (cpmcef_anis) )   ALLOCATE(cpmcef_anis(ntemp)) 
  IF (.NOT. ALLOCATED (bfactf_t))       ALLOCATE(bfactf_t(6,nat,ntemp))
  IF (.NOT. ALLOCATED (csmctf_t) )      ALLOCATE(csmctf_t(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (bthsf_t) )       ALLOCATE(bthsf_t(3,3,ntemp)) 
  IF (.NOT. ALLOCATED (ggammaf_t) )     ALLOCATE(ggammaf_t(3,3,ntemp)) 

  IF (.NOT. ALLOCATED (el_consf_t) )    ALLOCATE(el_consf_t(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (el_compf_t) )    ALLOCATE(el_compf_t(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (el_consf_s) )    ALLOCATE(el_consf_s(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (el_compf_s) )    ALLOCATE(el_compf_s(6,6,ntemp))
  IF (.NOT. ALLOCATED (macro_elf_s) )   ALLOCATE(macro_elf_s(8,ntemp)) 
  IF (.NOT. ALLOCATED (macro_elf_t) )   ALLOCATE(macro_elf_t(8,ntemp))
  IF (.NOT. ALLOCATED (vf_s) )          ALLOCATE(vf_s(3,ntemp))
  IF (.NOT. ALLOCATED (vf_t) )          ALLOCATE(vf_t(3,ntemp))
  IF (.NOT. ALLOCATED (debye_macro_elf_t) ) ALLOCATE(debye_macro_elf_t(ntemp)) 
  IF (.NOT. ALLOCATED (debye_macro_elf_s) ) ALLOCATE(debye_macro_elf_s(ntemp)) 
  IF (.NOT. ALLOCATED (el_conf_geo_t) ) ALLOCATE(el_conf_geo_t(6,6,&
                                                           ntemp,tot_ngeo)) 
  IF (.NOT. ALLOCATED (dydef_t) )        ALLOCATE(dydef_t(21,ntemp)) 

  IF (.NOT. ALLOCATED (e_piezo_tensorf_t) ) ALLOCATE(e_piezo_tensorf_t(3,6,ntemp)) 
  IF (.NOT. ALLOCATED (d_piezo_tensorf_t) ) ALLOCATE(d_piezo_tensorf_t(3,6,ntemp)) 
  IF (.NOT. ALLOCATED (epsilon_inftyf_t) ) ALLOCATE(epsilon_inftyf_t(3,3,ntemp)) 
  IF (.NOT. ALLOCATED (zeuf_t) )           ALLOCATE(zeuf_t(3,3,nat,ntemp)) 

  IF (.NOT. ALLOCATED (vminf_noe_t) )       ALLOCATE(vminf_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (densityf_noe_t) )    ALLOCATE(densityf_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (celldmf_noe_t) )     ALLOCATE(celldmf_noe_t(6,ntemp)) 
  IF (.NOT. ALLOCATED (free_e_minf_noe_t) ) ALLOCATE(free_e_minf_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0f_noe_t) )         ALLOCATE(b0f_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (b01f_noe_t) )        ALLOCATE(b01f_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (b02f_noe_t) )        ALLOCATE(b02f_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (af_noe_t) )          ALLOCATE(af_noe_t(m1,ntemp)) 

  IF (.NOT. ALLOCATED (free_enerf_noe_t) ) ALLOCATE(free_enerf_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (enerf_noe_t) )     ALLOCATE(enerf_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (entropyf_noe_t) )  ALLOCATE(entropyf_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (cef_noe_t) )       ALLOCATE(cef_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (cvf_noe_t) )       ALLOCATE(cvf_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (betaf_noe_t) )     ALLOCATE(betaf_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (cpf_noe_t) )       ALLOCATE(cpf_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0f_noe_s) )       ALLOCATE(b0f_noe_s(ntemp)) 
  IF (.NOT. ALLOCATED (gammaf_noe_t) )    ALLOCATE(gammaf_noe_t(ntemp)) 

  IF (.NOT. ALLOCATED (vgrun_t) )       ALLOCATE(vgrun_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0_grun_t) )     ALLOCATE(b0_grun_t(ntemp)) 
  IF (.NOT. ALLOCATED (celldm_grun_t) ) ALLOCATE(celldm_grun_t(6,ntemp)) 
  IF (.NOT. ALLOCATED (el_cons_grun_t) )ALLOCATE(el_cons_grun_t(6,6,ntemp))
  IF (.NOT. ALLOCATED (el_comp_grun_t) )ALLOCATE(el_comp_grun_t(6,6,ntemp))
  IF (.NOT. ALLOCATED (b0_grun_s) )     ALLOCATE(b0_grun_s(ntemp)) 
  IF (.NOT. ALLOCATED (cp_grun_t) )     ALLOCATE(cp_grun_t(ntemp)) 
  IF (.NOT. ALLOCATED (ce_grun_t) )     ALLOCATE(ce_grun_t(ntemp)) 
  IF (.NOT. ALLOCATED (cv_grun_t) )     ALLOCATE(cv_grun_t(ntemp)) 
  IF (.NOT. ALLOCATED (grun_cpmce_anis) ) ALLOCATE(grun_cpmce_anis(ntemp)) 
  IF (.NOT. ALLOCATED (alpha_an_g) )    ALLOCATE(alpha_an_g(6,ntemp)) 
  IF (.NOT. ALLOCATED (betab) )         ALLOCATE(betab(ntemp))
  IF (.NOT. ALLOCATED (grun_gamma_t) )  ALLOCATE(grun_gamma_t(ntemp)) 

  IF (.NOT. ALLOCATED (el_energy_t) )   ALLOCATE(el_energy_t(ntemp))
  IF (.NOT. ALLOCATED (el_free_energy_t) ) ALLOCATE(el_free_energy_t(ntemp))
  IF (.NOT. ALLOCATED (el_entropy_t) )  ALLOCATE(el_entropy_t(ntemp))
  IF (.NOT. ALLOCATED (el_ce_t) )       ALLOCATE(el_ce_t(ntemp))
  IF (.NOT. ALLOCATED (el_b0_t) )       ALLOCATE(el_b0_t(ntemp))

  IF (.NOT. ALLOCATED (emp_energy_t) )   ALLOCATE(emp_energy_t(ntemp))
  IF (.NOT. ALLOCATED (emp_free_energy_t) ) ALLOCATE(emp_free_energy_t(ntemp))
  IF (.NOT. ALLOCATED (emp_entropy_t) )  ALLOCATE(emp_entropy_t(ntemp))
  IF (.NOT. ALLOCATED (emp_ce_t) )       ALLOCATE(emp_ce_t(ntemp))

  IF (.NOT. ALLOCATED (emp_energyf_t) )   ALLOCATE(emp_energyf_t(ntemp))
  IF (.NOT. ALLOCATED (emp_free_energyf_t) ) ALLOCATE(emp_free_energyf_t(ntemp))
  IF (.NOT. ALLOCATED (emp_entropyf_t) )  ALLOCATE(emp_entropyf_t(ntemp))
  IF (.NOT. ALLOCATED (emp_cef_t) )       ALLOCATE(emp_cef_t(ntemp))

  IF (npress_plot>0) THEN
     IF (.NOT. ALLOCATED (el_free_ener_pt) ) ALLOCATE(el_free_ener_pt(ntemp,&
                                                                 npress_plot))
     IF (.NOT. ALLOCATED (el_ener_pt) ) ALLOCATE(el_ener_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (el_entr_pt) ) ALLOCATE(el_entr_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (el_ce_pt) )   ALLOCATE(el_ce_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (emp_ener_pt) )  &
                                      ALLOCATE(emp_ener_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (emp_entr_pt) )  &
                                      ALLOCATE(emp_entr_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (emp_ce_pt) )  ALLOCATE(emp_ce_pt(ntemp,npress_plot))

     IF (.NOT. ALLOCATED (alpha_anis_pt) ) &
                                ALLOCATE(alpha_anis_pt(6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (density_pt) ) ALLOCATE(density_pt(ntemp,npress_plot))

     IF (.NOT. ALLOCATED(uintf_pt)) &
                       ALLOCATE(uintf_pt(nint_var,ntemp,npress_plot))
     IF (.NOT.ALLOCATED(tauf_pt)) &
                       ALLOCATE(tauf_pt(3,nat,ntemp,npress_plot))
     IF (.NOT. ALLOCATED(alphaf_int_pt)) &
                       ALLOCATE(alphaf_int_pt(nint_var,ntemp,npress_plot))
     IF (.NOT. ALLOCATED(uintf_zsisa_pt)) &
                       ALLOCATE(uintf_zsisa_pt(nint_var,ntemp,npress_plot))
     IF (.NOT.ALLOCATED(tauf_zsisa_pt)) &
                       ALLOCATE(tauf_zsisa_pt(3,nat,ntemp,npress_plot))
     IF (.NOT. ALLOCATED(alphaf_int_zsisa_pt)) &
                      ALLOCATE(alphaf_int_zsisa_pt(nint_var,ntemp,npress_plot))

     IF (.NOT. ALLOCATED (cpmce_anis_pt) ) &
                                ALLOCATE(cpmce_anis_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (bths_pt) ) &
                                ALLOCATE(bths_pt(3,3,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (ggamma_pt) ) &
                                ALLOCATE(ggamma_pt(3,3,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (csmct_pt) ) ALLOCATE(csmct_pt(6,6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (el_cons_pt) ) &
                          ALLOCATE(el_cons_pt(6,6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (el_comp_pt) ) &
                          ALLOCATE(el_comp_pt(6,6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (el_cons_s_pt) ) &
                          ALLOCATE(el_cons_s_pt(6,6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (el_comp_s_pt) ) &
                          ALLOCATE(el_comp_s_pt(6,6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (macro_el_pt) ) &
                          ALLOCATE(macro_el_pt(8,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (macro_el_s_pt) ) &
                          ALLOCATE(macro_el_s_pt(8,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (v_pt) ) &
                          ALLOCATE(v_pt(3,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (v_s_pt) ) &
                          ALLOCATE(v_s_pt(3,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (debye_macro_el_pt) ) &
                          ALLOCATE(debye_macro_el_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (debye_macro_el_s_pt) ) &
                          ALLOCATE(debye_macro_el_s_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (e_piezo_tensor_pt) ) &
                          ALLOCATE(e_piezo_tensor_pt(3,6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (d_piezo_tensor_pt) ) &
                          ALLOCATE(d_piezo_tensor_pt(3,6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (epsilon_infty_pt) ) &
                          ALLOCATE(epsilon_infty_pt(3,3,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (zeu_pt) ) &
                          ALLOCATE(zeu_pt(3,3,nat,ntemp,npress_plot))

     IF (.NOT. ALLOCATED (el_free_enerf_pt) ) ALLOCATE(el_free_enerf_pt(ntemp,&
                                                                 npress_plot))
     IF (.NOT. ALLOCATED (el_enerf_pt) ) ALLOCATE(el_enerf_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (el_entrf_pt) ) ALLOCATE(el_entrf_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (el_cef_pt) )   ALLOCATE(el_cef_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (emp_enerf_pt) )  &
                                    ALLOCATE(emp_enerf_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (emp_entrf_pt) )  &
                                    ALLOCATE(emp_entrf_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (emp_cef_pt) ) ALLOCATE(emp_cef_pt(ntemp,npress_plot))

     IF (.NOT. ALLOCATED (alphaf_anis_pt) ) &
                                ALLOCATE(alphaf_anis_pt(6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (cpmcef_anis_pt) ) &
                                ALLOCATE(cpmcef_anis_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (bthsf_pt) ) &
                                ALLOCATE(bthsf_pt(3,3,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (ggammaf_pt) ) &
                                ALLOCATE(ggammaf_pt(3,3,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (csmctf_pt) ) ALLOCATE(csmctf_pt(6,6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (el_consf_pt) ) &
                          ALLOCATE(el_consf_pt(6,6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (el_compf_pt) ) &
                          ALLOCATE(el_compf_pt(6,6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (el_consf_s_pt) ) &
                          ALLOCATE(el_consf_s_pt(6,6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (el_compf_s_pt) ) &
                          ALLOCATE(el_compf_s_pt(6,6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (macro_elf_pt) ) &
                          ALLOCATE(macro_elf_pt(8,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (macro_elf_s_pt) ) &
                          ALLOCATE(macro_elf_s_pt(8,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (vf_pt) ) &
                          ALLOCATE(vf_pt(3,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (vf_s_pt) ) &
                          ALLOCATE(vf_s_pt(3,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (debye_macro_elf_pt) ) &
                          ALLOCATE(debye_macro_elf_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (debye_macro_elf_s_pt) ) &
                          ALLOCATE(debye_macro_elf_s_pt(ntemp,npress_plot))
     IF (.NOT. ALLOCATED (e_piezo_tensorf_pt) ) &
                          ALLOCATE(e_piezo_tensorf_pt(3,6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (d_piezo_tensorf_pt) ) &
                          ALLOCATE(d_piezo_tensorf_pt(3,6,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (epsilon_inftyf_pt) ) &
                          ALLOCATE(epsilon_inftyf_pt(3,3,ntemp,npress_plot))
     IF (.NOT. ALLOCATED (zeuf_pt) ) &
                          ALLOCATE(zeuf_pt(3,3,nat,ntemp,npress_plot))
  ENDIF
  IF (ntemp_plot>0) THEN
     IF (.NOT. ALLOCATED (alpha_anis_ptt) ) ALLOCATE(alpha_anis_ptt&
                                                    (6,npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (density_ptt) )  ALLOCATE(density_ptt&
                                                          (npress,ntemp_plot))
     IF (.NOT. ALLOCATED (cpmce_anis_ptt) )  ALLOCATE(cpmce_anis_ptt&
                                                          (npress,ntemp_plot))
     IF (.NOT. ALLOCATED (bths_ptt) )  ALLOCATE(bths_ptt&
                                                    (3,3,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (ggamma_ptt) )  ALLOCATE(ggamma_ptt&
                                                    (3,3,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(uintf_ptt)) &
                       ALLOCATE(uintf_ptt(nint_var,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(tauf_ptt)) &
                       ALLOCATE(tauf_ptt(3,nat,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(uintf_ptt_p1)) &
                       ALLOCATE(uintf_ptt_p1(nint_var,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(tauf_ptt_p1)) &
                       ALLOCATE(tauf_ptt_p1(3,nat,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(uintf_ptt_m1)) &
                       ALLOCATE(uintf_ptt_m1(nint_var,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(tauf_ptt_m1)) &
                       ALLOCATE(tauf_ptt_m1(3,nat,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(alphaf_int_ptt)) &
                       ALLOCATE(alphaf_int_ptt(nint_var,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(uintf_zsisa_ptt)) &
                       ALLOCATE(uintf_zsisa_ptt(nint_var,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(tauf_zsisa_ptt)) &
                       ALLOCATE(tauf_zsisa_ptt(3,nat,npress,ntemp_plot))
     IF (.NOT. ALLOCATED(alphaf_int_zsisa_ptt)) &
                       ALLOCATE(alphaf_int_zsisa_ptt(nint_var,npress,&
                                                                ntemp_plot))
     IF (.NOT. ALLOCATED (csmct_ptt) )  ALLOCATE(csmct_ptt&
                                                    (6,6,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (el_ce_ptt) )  ALLOCATE(el_ce_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (emp_cef_ptt) ) ALLOCATE(emp_cef_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (el_cef_ptt) ) ALLOCATE(el_cef_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (emp_ce_ptt) ) ALLOCATE(emp_ce_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (el_cons_ptt) ) &
                          ALLOCATE(el_cons_ptt(6,6,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (el_comp_ptt) ) &
                          ALLOCATE(el_comp_ptt(6,6,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (el_cons_s_ptt) ) &
                          ALLOCATE(el_cons_s_ptt(6,6,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (el_comp_s_ptt) ) &
                          ALLOCATE(el_comp_s_ptt(6,6,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (macro_el_ptt) ) &
                          ALLOCATE(macro_el_ptt(8,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (macro_el_s_ptt) ) &
                          ALLOCATE(macro_el_s_ptt(8,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (v_ptt) ) &
                          ALLOCATE(v_ptt(3,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (v_s_ptt) ) &
                          ALLOCATE(v_s_ptt(3,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (debye_macro_el_ptt) ) &
                          ALLOCATE(debye_macro_el_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (debye_macro_el_s_ptt) ) &
                          ALLOCATE(debye_macro_el_s_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (e_piezo_tensor_ptt) ) &
                          ALLOCATE(e_piezo_tensor_ptt(3,6,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (d_piezo_tensor_ptt) ) &
                          ALLOCATE(d_piezo_tensor_ptt(3,6,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (epsilon_infty_ptt) ) &
                          ALLOCATE(epsilon_infty_ptt(3,3,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (zeu_ptt) ) &
                          ALLOCATE(zeu_ptt(3,3,nat,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (alphaf_anis_ptt) ) ALLOCATE(alphaf_anis_ptt&
                                                    (6,npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (densityf_ptt) )  ALLOCATE(densityf_ptt&
                                                          (npress,ntemp_plot))
     IF (.NOT. ALLOCATED (cpmcef_anis_ptt) )  ALLOCATE(cpmcef_anis_ptt&
                                                          (npress,ntemp_plot))
     IF (.NOT. ALLOCATED (bthsf_ptt) )  ALLOCATE(bthsf_ptt&
                                                    (3,3,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (ggammaf_ptt) )  ALLOCATE(ggammaf_ptt&
                                                    (3,3,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (csmctf_ptt) )  ALLOCATE(csmctf_ptt&
                                                    (6,6,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (el_consf_ptt) ) &
                          ALLOCATE(el_consf_ptt(6,6,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (el_compf_ptt) ) &
                          ALLOCATE(el_compf_ptt(6,6,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (el_consf_s_ptt) ) &
                          ALLOCATE(el_consf_s_ptt(6,6,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (el_compf_s_ptt) ) &
                          ALLOCATE(el_compf_s_ptt(6,6,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (macro_elf_ptt) ) &
                          ALLOCATE(macro_elf_ptt(8,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (macro_elf_s_ptt) ) &
                          ALLOCATE(macro_elf_s_ptt(8,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (vf_ptt) ) &
                          ALLOCATE(vf_ptt(3,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (vf_s_ptt) ) &
                          ALLOCATE(vf_s_ptt(3,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (debye_macro_elf_ptt) ) &
                          ALLOCATE(debye_macro_elf_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (debye_macro_elf_s_ptt) ) &
                          ALLOCATE(debye_macro_elf_s_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (e_piezo_tensorf_ptt) ) &
                          ALLOCATE(e_piezo_tensorf_ptt(3,6,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (d_piezo_tensorf_ptt) ) &
                          ALLOCATE(d_piezo_tensorf_ptt(3,6,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (epsilon_inftyf_ptt) ) &
                          ALLOCATE(epsilon_inftyf_ptt(3,3,npress,ntemp_plot))
     IF (.NOT. ALLOCATED (zeuf_ptt) ) &
                          ALLOCATE(zeuf_ptt(3,3,nat,npress,ntemp_plot))
  ENDIF

  IF (.NOT. ALLOCATED (el_energyf_t) )  ALLOCATE(el_energyf_t(ntemp))
  IF (.NOT. ALLOCATED (el_free_energyf_t) ) ALLOCATE(el_free_energyf_t(ntemp))
  IF (.NOT. ALLOCATED (el_entropyf_t) ) ALLOCATE(el_entropyf_t(ntemp))
  IF (.NOT. ALLOCATED (el_cef_t) )      ALLOCATE(el_cef_t(ntemp))
  IF (.NOT. ALLOCATED (el_beta_t) )     ALLOCATE(el_beta_t(ntemp))
  IF (.NOT. ALLOCATED (el_betaf_t) )    ALLOCATE(el_betaf_t(ntemp))
  IF (.NOT. ALLOCATED (el_cp_t) )       ALLOCATE(el_cp_t(ntemp))
  IF (.NOT. ALLOCATED (el_cpf_t) )      ALLOCATE(el_cpf_t(ntemp))
  IF (.NOT. ALLOCATED (el_gamma_t) )    ALLOCATE(el_gamma_t(ntemp))
  IF (.NOT. ALLOCATED (el_gammaf_t) )   ALLOCATE(el_gammaf_t(ntemp))

  IF (.NOT. ALLOCATED (celldmf_t_p1) ) ALLOCATE(celldmf_t_p1(6,ntemp)) 
  IF (.NOT. ALLOCATED (celldmf_t_m1) ) ALLOCATE(celldmf_t_m1(6,ntemp)) 
  IF (.NOT. ALLOCATED (celldmf_noe_t_p1) ) ALLOCATE(celldmf_noe_t_p1(6,ntemp)) 
  IF (.NOT. ALLOCATED (celldmf_noe_t_m1) ) ALLOCATE(celldmf_noe_t_m1(6,ntemp)) 
  IF (.NOT. ALLOCATED (p2tf_noe_t) )    ALLOCATE(p2tf_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (p1tf_noe_t) )    ALLOCATE(p1tf_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (p3tf_noe_t) )    ALLOCATE(p3tf_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (p4tf_noe_t) )    ALLOCATE(p4tf_noe_t(ntemp)) 

  RETURN
  !
END SUBROUTINE allocate_anharmonic

!---------------------------------------------------------------------
SUBROUTINE allocate_debye()
!---------------------------------------------------------------------

USE ions_base,     ONLY : nat
USE temperature,   ONLY : ntemp
USE control_debye, ONLY : deb_energy, deb_free_energy, deb_entropy, deb_cv, &
                          deb_b_fact, deb_bfact

IMPLICIT NONE

IF (.NOT. ALLOCATED (deb_energy) )       ALLOCATE( deb_energy(ntemp) )
IF (.NOT. ALLOCATED (deb_free_energy) )  ALLOCATE( deb_free_energy(ntemp) )
IF (.NOT. ALLOCATED (deb_entropy) )      ALLOCATE( deb_entropy(ntemp) )
IF (.NOT. ALLOCATED (deb_cv) )           ALLOCATE( deb_cv(ntemp) )
IF (.NOT. ALLOCATED (deb_bfact) )        ALLOCATE( deb_bfact(ntemp) )
IF (.NOT. ALLOCATED (deb_b_fact) )       ALLOCATE( deb_b_fact(3,3,nat,ntemp) )

RETURN
END SUBROUTINE allocate_debye
