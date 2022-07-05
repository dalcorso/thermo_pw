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
  USE thermo_mod,     ONLY : tot_ngeo
  USE temperature,    ONLY : ntemp
  USE thermodynamics, ONLY : ph_free_ener, ph_ener, ph_entropy, ph_e0, &
                             ph_ce, ph_b_fact
  USE ph_freq_thermodynamics, ONLY : phf_free_ener, phf_ener, phf_entropy, &
                                     phf_e0, phf_ce, phf_b_fact, ph_freq_save
  USE el_thermodynamics, ONLY : el_ener, el_free_ener, el_entr, el_mu, &
                           el_ce
  USE el_anharmonic,  ONLY : vmine_t, b0e_t, b01e_t, b02e_t, free_e_mine_t

  USE control_emp_free_ener, ONLY : emp_ener, emp_free_ener, emp_entr, &
                           emp_ce

  IMPLICIT NONE

  IF (.NOT.ALLOCATED(ph_free_ener))  ALLOCATE(ph_free_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(ph_ener))       ALLOCATE(ph_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(ph_entropy))    ALLOCATE(ph_entropy(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(ph_e0))         ALLOCATE(ph_e0(tot_ngeo))
  IF (.NOT.ALLOCATED(ph_ce))         ALLOCATE(ph_ce(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(ph_b_fact))     ALLOCATE(ph_b_fact(3,3,nat,ntemp,tot_ngeo))

  IF (.NOT.ALLOCATED(phf_free_ener)) ALLOCATE(phf_free_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(phf_ener))      ALLOCATE(phf_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(phf_entropy))   ALLOCATE(phf_entropy(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(phf_e0))        ALLOCATE(phf_e0(tot_ngeo))
  IF (.NOT.ALLOCATED(phf_ce))        ALLOCATE(phf_ce(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(phf_b_fact))    ALLOCATE(phf_b_fact(3,3,nat,ntemp,tot_ngeo))
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
                                  el_comp_t, el_cons_s, el_comp_s,           &
                                  macro_el_t, macro_el_s, v_s, v_t,          &
                                  el_con_geo_t, vmin_noe_t, b0_noe_t,        &
                                  b01_noe_t, b02_noe_t, free_e_min_noe_t,    &
                                  a_noe_t, beta_noe_t, cv_noe_t, ce_noe_t,   &
                                  cp_noe_t, b0_noe_s, gamma_noe_t           
  USE anharmonic_pt,       ONLY : vmin_pt, b0_pt, b01_pt, b02_pt, emin_pt,   &
                                  ce_pt, cv_pt, cp_pt, beta_pt, b0_s_pt,     &
                                  gamma_pt, free_ener_pt, ener_pt, entr_pt
  USE anharmonic_ptt,      ONLY : vmin_ptt, b0_ptt, b01_ptt, b02_ptt,        &
                                  emin_ptt, ce_ptt, cv_ptt, cp_ptt,          &
                                  beta_ptt, b0_s_ptt, gamma_ptt
  USE anharmonic_vt,       ONLY : press_vt, press_vtt
  USE ph_freq_anharmonic,  ONLY : vminf_t, b0f_t, b01f_t, b02f_t,            &
                                  free_e_minf_t, af_t,                       &
                                  enerf_t, free_enerf_t, entropyf_t,         &
                                  cef_t, cvf_t, cpf_t, b0f_s, alphaf_t,      &
                                  betaf_t, gammaf_t, celldmf_t,              &
                                  alphaf_anis_t, cpmcef_anis,                &
                                  bfactf_t, bthsf_t, ggammaf_t,              &
                                  el_consf_t, el_compf_t, el_consf_s,        &
                                  el_compf_s, macro_elf_t, macro_elf_s,      &
                                  vf_t, vf_s, el_conf_geo_t
  USE grun_anharmonic,     ONLY : betab, alpha_an_g, cp_grun_t, cv_grun_t,   &
                                  ce_grun_t, b0_grun_s, grun_gamma_t,        &
                                  grun_cpmce_anis, el_cons_grun_t,           &
                                  el_comp_grun_t
  USE el_anharmonic, ONLY : el_energy_t, el_free_energy_t, el_entropy_t,     &
                            el_ce_t, el_energyf_t, el_free_energyf_t,        &
                            el_entropyf_t, el_cef_t, el_ce_pt, el_ce_ptt,    &
                            el_free_ener_pt, el_ener_pt, el_entr_pt,         &
                            el_b0_t, el_beta_t, el_betaf_t, el_cp_t,         &
                            el_cpf_t, el_gamma_t, el_gammaf_t
  USE emp_anharmonic, ONLY : emp_energy_t, emp_free_energy_t, emp_entropy_t, &
                            emp_ce_t, emp_energyf_t, emp_free_energyf_t,     &
                            emp_entropyf_t, emp_cef_t, emp_free_ener_pt,     &
                            emp_ener_pt, emp_entr_pt, emp_ce_pt, emp_ce_ptt

  USE control_quartic_energy, ONLY : poly_degree_ph
  USE control_vol,         ONLY : nvol_plot, nvol
  USE control_grun,        ONLY : vgrun_t, celldm_grun_t, b0_grun_t
  USE control_pressure,    ONLY : npress, npress_plot

  IMPLICIT NONE
  INTEGER :: m1

  m1=poly_degree_ph+1

  IF (.NOT. ALLOCATED (vmin_t) )        ALLOCATE(vmin_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0_t) )          ALLOCATE(b0_t(ntemp)) 
  IF (.NOT. ALLOCATED (b01_t) )         ALLOCATE(b01_t(ntemp)) 
  IF (.NOT. ALLOCATED (b02_t) )         ALLOCATE(b02_t(ntemp)) 
  IF (.NOT. ALLOCATED (free_e_min_t) )  ALLOCATE(free_e_min_t(ntemp)) 
  IF (.NOT. ALLOCATED (a_t) )           ALLOCATE(a_t(m1,ntemp)) 

  IF (.NOT. ALLOCATED (ener_t) )        ALLOCATE(ener_t(ntemp))
  IF (.NOT. ALLOCATED (free_ener_t) )   ALLOCATE(free_ener_t(ntemp)) 
  IF (.NOT. ALLOCATED (entropy_t) )     ALLOCATE(entropy_t(ntemp)) 
  IF (.NOT. ALLOCATED (ce_t) )          ALLOCATE(ce_t(ntemp)) 
  IF (.NOT. ALLOCATED (cv_t) )          ALLOCATE(cv_t(ntemp)) 

  IF (.NOT. ALLOCATED (cp_t) )          ALLOCATE(cp_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0_s) )          ALLOCATE(b0_s(ntemp)) 
  IF (.NOT. ALLOCATED (alpha_t) )       ALLOCATE(alpha_t(ntemp)) 
  IF (.NOT. ALLOCATED (beta_t) )        ALLOCATE(beta_t(ntemp)) 
  IF (.NOT. ALLOCATED (gamma_t) )       ALLOCATE(gamma_t(ntemp)) 

  IF (.NOT. ALLOCATED (celldm_t) )      ALLOCATE(celldm_t(6,ntemp)) 
  IF (.NOT. ALLOCATED (alpha_anis_t) )  ALLOCATE(alpha_anis_t(6,ntemp)) 
  IF (.NOT. ALLOCATED (cpmce_anis) )    ALLOCATE(cpmce_anis(ntemp)) 
  IF (.NOT. ALLOCATED (bfact_t))        ALLOCATE(bfact_t(6,nat,ntemp))
  IF (.NOT. ALLOCATED (bths_t) )        ALLOCATE(bths_t(3,3,ntemp)) 
  IF (.NOT. ALLOCATED (ggamma_t) )      ALLOCATE(ggamma_t(3,3,ntemp)) 

  IF (.NOT. ALLOCATED (el_cons_t) )     ALLOCATE(el_cons_t(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (el_comp_t) )     ALLOCATE(el_comp_t(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (el_cons_s) )     ALLOCATE(el_cons_s(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (el_comp_s) )     ALLOCATE(el_comp_s(6,6,ntemp))
  IF (.NOT. ALLOCATED (macro_el_s) )    ALLOCATE(macro_el_s(8,ntemp))
  IF (.NOT. ALLOCATED (macro_el_t) )    ALLOCATE(macro_el_t(8,ntemp)) 
  IF (.NOT. ALLOCATED (v_s) )           ALLOCATE(v_s(3,ntemp))
  IF (.NOT. ALLOCATED (v_t) )           ALLOCATE(v_t(3,ntemp))
  IF (.NOT. ALLOCATED (el_con_geo_t) )  ALLOCATE(el_con_geo_t(6,6,ntemp,&
                                                                   tot_ngeo)) 

  IF (.NOT. ALLOCATED (vmin_noe_t) )     ALLOCATE(vmin_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0_noe_t) )       ALLOCATE(b0_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (b01_noe_t) )      ALLOCATE(b01_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (b02_noe_t) )      ALLOCATE(b02_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (free_e_min_noe_t) ) ALLOCATE(free_e_min_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (a_noe_t) )        ALLOCATE(a_noe_t(m1,ntemp)) 

  IF (.NOT. ALLOCATED (ce_noe_t) )       ALLOCATE(ce_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (cv_noe_t) )       ALLOCATE(cv_noe_t(ntemp)) 

  IF (.NOT. ALLOCATED (cp_noe_t) )       ALLOCATE(cp_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0_noe_s) )       ALLOCATE(b0_noe_s(ntemp)) 
  IF (.NOT. ALLOCATED (beta_noe_t) )     ALLOCATE(beta_noe_t(ntemp)) 
  IF (.NOT. ALLOCATED (gamma_noe_t) )    ALLOCATE(gamma_noe_t(ntemp)) 

  IF (ntemp_plot>0) THEN
     IF (.NOT. ALLOCATED (vmin_ptt) )   ALLOCATE(vmin_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (b0_ptt) )     ALLOCATE(b0_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (b01_ptt) )    ALLOCATE(b01_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (b02_ptt) )    ALLOCATE(b02_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (emin_ptt) )   ALLOCATE(emin_ptt(npress,ntemp_plot)) 

     IF (.NOT. ALLOCATED (ce_ptt) )     ALLOCATE(ce_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (cv_ptt) )     ALLOCATE(cv_ptt(npress,ntemp_plot)) 

     IF (.NOT. ALLOCATED (cp_ptt) )     ALLOCATE(cp_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (b0_s_ptt) )   ALLOCATE(b0_s_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (beta_ptt) )   ALLOCATE(beta_ptt(npress,ntemp_plot)) 
     IF (.NOT. ALLOCATED (gamma_ptt) )  ALLOCATE(gamma_ptt(npress,ntemp_plot)) 

     IF (.NOT. ALLOCATED (press_vt) )   ALLOCATE(press_vt(nvol,ntemp_plot)) 
  ENDIF

  IF (npress_plot>0) THEN
     IF (.NOT. ALLOCATED (vmin_pt) )    ALLOCATE(vmin_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (b0_pt) )      ALLOCATE(b0_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (b01_pt) )     ALLOCATE(b01_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (b02_pt) )     ALLOCATE(b02_pt(ntemp,npress_plot)) 
     IF (.NOT. ALLOCATED (emin_pt) )    ALLOCATE(emin_pt(ntemp,npress_plot)) 

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
  ENDIF

  IF (nvol_plot>0) THEN
     IF (.NOT. ALLOCATED (press_vtt) )   ALLOCATE(press_vtt(ntemp,nvol_plot)) 
  ENDIF

  IF (.NOT. ALLOCATED (vminf_t) )       ALLOCATE(vminf_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0f_t) )         ALLOCATE(b0f_t(ntemp)) 
  IF (.NOT. ALLOCATED (b01f_t) )        ALLOCATE(b01f_t(ntemp)) 
  IF (.NOT. ALLOCATED (b02f_t) )        ALLOCATE(b02f_t(ntemp)) 
  IF (.NOT. ALLOCATED (free_e_minf_t) ) ALLOCATE(free_e_minf_t(ntemp)) 
  IF (.NOT. ALLOCATED (af_t) )          ALLOCATE(af_t(m1,ntemp)) 

  IF (.NOT. ALLOCATED (enerf_t) )       ALLOCATE(enerf_t(ntemp)) 
  IF (.NOT. ALLOCATED (free_enerf_t) )  ALLOCATE(free_enerf_t(ntemp))
  IF (.NOT. ALLOCATED (entropyf_t) )    ALLOCATE(entropyf_t(ntemp)) 
  IF (.NOT. ALLOCATED (cef_t) )         ALLOCATE(cef_t(ntemp)) 
  IF (.NOT. ALLOCATED (cvf_t) )         ALLOCATE(cvf_t(ntemp)) 

  IF (.NOT. ALLOCATED (cpf_t) )         ALLOCATE(cpf_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0f_s) )         ALLOCATE(b0f_s(ntemp)) 
  IF (.NOT. ALLOCATED (alphaf_t) )      ALLOCATE(alphaf_t(ntemp)) 
  IF (.NOT. ALLOCATED (betaf_t) )       ALLOCATE(betaf_t(ntemp)) 
  IF (.NOT. ALLOCATED (gammaf_t) )      ALLOCATE(gammaf_t(ntemp)) 

  IF (.NOT. ALLOCATED (celldmf_t) )     ALLOCATE(celldmf_t(6,ntemp)) 
  IF (.NOT. ALLOCATED (alphaf_anis_t) ) ALLOCATE(alphaf_anis_t(6,ntemp)) 
  IF (.NOT. ALLOCATED (cpmcef_anis) )   ALLOCATE(cpmcef_anis(ntemp)) 
  IF (.NOT. ALLOCATED (bfactf_t))       ALLOCATE(bfactf_t(6,nat,ntemp))
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
  IF (.NOT. ALLOCATED (el_conf_geo_t) ) ALLOCATE(el_conf_geo_t(6,6,&
                                                           ntemp,tot_ngeo)) 

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
  ENDIF
  IF (ntemp_plot>0) THEN
     IF (.NOT. ALLOCATED (el_ce_ptt) )  ALLOCATE(el_ce_ptt(npress,ntemp_plot))
     IF (.NOT. ALLOCATED (emp_ce_ptt) ) ALLOCATE(emp_ce_ptt(npress,ntemp_plot))
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
