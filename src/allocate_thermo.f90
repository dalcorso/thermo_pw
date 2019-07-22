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
  USE thermodynamics, ONLY : ph_free_ener, ph_ener, ph_entropy, ph_ce, ph_b_fact
  USE ph_freq_thermodynamics, ONLY : phf_free_ener, phf_ener, phf_entropy, &
                                     phf_ce, phf_b_fact

  IMPLICIT NONE

  IF (.NOT.ALLOCATED(ph_free_ener))  ALLOCATE(ph_free_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(ph_ener))       ALLOCATE(ph_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(ph_entropy))    ALLOCATE(ph_entropy(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(ph_ce))         ALLOCATE(ph_ce(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(ph_b_fact))     ALLOCATE(ph_b_fact(3,3,nat,ntemp,tot_ngeo))

  IF (.NOT.ALLOCATED(phf_free_ener)) ALLOCATE(phf_free_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(phf_ener))      ALLOCATE(phf_ener(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(phf_entropy))   ALLOCATE(phf_entropy(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(phf_ce))        ALLOCATE(phf_ce(ntemp,tot_ngeo))
  IF (.NOT.ALLOCATED(phf_b_fact))    ALLOCATE(phf_b_fact(3,3,nat,ntemp,tot_ngeo))

  RETURN
  !
END SUBROUTINE allocate_thermodynamics

SUBROUTINE allocate_anharmonic()

  USE thermo_mod,          ONLY : tot_ngeo
  USE temperature,         ONLY : ntemp
  USE ions_base,           ONLY : nat
  USE anharmonic,          ONLY : vmin_t, b0_t, b01_t, free_e_min_t,         &
                                  alpha_t, beta_t, gamma_t, cv_t, ce_t,      &
                                  cp_t, b0_s, &
                                  celldm_t, alpha_anis_t, cpmce_anis,        &
                                  el_cons_t, el_comp_t, macro_el_t,          &
                                  el_cons_s, el_comp_s,                      &
                                  bths_t, ggamma_t, bfact_t, el_con_geo_t
  USE ph_freq_anharmonic,  ONLY : vminf_t, b0f_t, b01f_t, free_e_minf_t,     &
                                  alphaf_t, betaf_t, gammaf_t, cvf_t, cef_t, &
                                  cpf_t, b0f_s, celldmf_t, alphaf_anis_t,    & 
                                  cpmcef_anis, el_consf_t, el_compf_t,       &
                                  el_consf_s, el_compf_s,                    &
                                  macro_elf_t, bthsf_t, ggammaf_t, bfactf_t, &
                                  el_conf_geo_t
  USE grun_anharmonic,     ONLY : betab, alpha_an_g, cp_grun_t, cv_grun_t,   &
                                  ce_grun_t, b0_grun_s, grun_gamma_t,        &
                                  grun_cpmce_anis, el_cons_grun_t,           &
                                  el_comp_grun_t
  USE control_grun,        ONLY : vgrun_t, celldm_grun_t, b0_grun_t

  IMPLICIT NONE

  IF (.NOT. ALLOCATED (vmin_t) )        ALLOCATE(vmin_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0_t) )          ALLOCATE(b0_t(ntemp)) 
  IF (.NOT. ALLOCATED (b01_t) )         ALLOCATE(b01_t(ntemp)) 
  IF (.NOT. ALLOCATED (free_e_min_t) )  ALLOCATE(free_e_min_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0_s) )          ALLOCATE(b0_s(ntemp)) 
  IF (.NOT. ALLOCATED (cv_t) )          ALLOCATE(cv_t(ntemp)) 
  IF (.NOT. ALLOCATED (ce_t) )          ALLOCATE(ce_t(ntemp)) 
  IF (.NOT. ALLOCATED (cp_t) )          ALLOCATE(cp_t(ntemp)) 
  IF (.NOT. ALLOCATED (alpha_t) )       ALLOCATE(alpha_t(ntemp)) 
  IF (.NOT. ALLOCATED (beta_t) )        ALLOCATE(beta_t(ntemp)) 
  IF (.NOT. ALLOCATED (gamma_t) )       ALLOCATE(gamma_t(ntemp)) 
  IF (.NOT. ALLOCATED (bths_t) )        ALLOCATE(bths_t(3,3,ntemp)) 
  IF (.NOT. ALLOCATED (ggamma_t) )      ALLOCATE(ggamma_t(3,3,ntemp)) 
  IF (.NOT. ALLOCATED (bfact_t))        ALLOCATE(bfact_t(6,nat,ntemp))
  IF (.NOT. ALLOCATED (celldm_t) )      ALLOCATE(celldm_t(6,ntemp)) 
  IF (.NOT. ALLOCATED (alpha_anis_t) )  ALLOCATE(alpha_anis_t(6,ntemp)) 
  IF (.NOT. ALLOCATED (cpmce_anis) )    ALLOCATE(cpmce_anis(ntemp)) 
  IF (.NOT. ALLOCATED (el_cons_t) )     ALLOCATE(el_cons_t(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (el_comp_t) )     ALLOCATE(el_comp_t(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (el_cons_s) )     ALLOCATE(el_cons_s(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (el_comp_s) )     ALLOCATE(el_comp_s(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (macro_el_t) )    ALLOCATE(macro_el_t(8,ntemp)) 
  IF (.NOT. ALLOCATED (el_con_geo_t) )  ALLOCATE(el_con_geo_t(6,6,ntemp,&
                                                                   tot_ngeo)) 

  IF (.NOT. ALLOCATED (vminf_t) )       ALLOCATE(vminf_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0f_t) )         ALLOCATE(b0f_t(ntemp)) 
  IF (.NOT. ALLOCATED (b01f_t) )        ALLOCATE(b01f_t(ntemp)) 
  IF (.NOT. ALLOCATED (free_e_minf_t) ) ALLOCATE(free_e_minf_t(ntemp)) 
  IF (.NOT. ALLOCATED (b0f_s) )         ALLOCATE(b0f_s(ntemp)) 
  IF (.NOT. ALLOCATED (cvf_t) )         ALLOCATE(cvf_t(ntemp)) 
  IF (.NOT. ALLOCATED (cef_t) )         ALLOCATE(cef_t(ntemp)) 
  IF (.NOT. ALLOCATED (cpf_t) )         ALLOCATE(cpf_t(ntemp)) 
  IF (.NOT. ALLOCATED (alphaf_t) )      ALLOCATE(alphaf_t(ntemp)) 
  IF (.NOT. ALLOCATED (betaf_t) )       ALLOCATE(betaf_t(ntemp)) 
  IF (.NOT. ALLOCATED (gammaf_t) )      ALLOCATE(gammaf_t(ntemp)) 
  IF (.NOT. ALLOCATED (bthsf_t) )       ALLOCATE(bthsf_t(3,3,ntemp)) 
  IF (.NOT. ALLOCATED (ggammaf_t) )     ALLOCATE(ggammaf_t(3,3,ntemp)) 
  IF (.NOT. ALLOCATED (celldmf_t) )     ALLOCATE(celldmf_t(6,ntemp)) 
  IF (.NOT. ALLOCATED (bfactf_t))       ALLOCATE(bfactf_t(6,nat,ntemp))
  IF (.NOT. ALLOCATED (alphaf_anis_t) ) ALLOCATE(alphaf_anis_t(6,ntemp)) 
  IF (.NOT. ALLOCATED (cpmcef_anis) )   ALLOCATE(cpmcef_anis(ntemp)) 
  IF (.NOT. ALLOCATED (el_consf_t) )    ALLOCATE(el_consf_t(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (el_compf_t) )    ALLOCATE(el_compf_t(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (el_consf_s) )    ALLOCATE(el_consf_s(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (el_compf_s) )    ALLOCATE(el_compf_s(6,6,ntemp)) 
  IF (.NOT. ALLOCATED (macro_elf_t) )   ALLOCATE(macro_elf_t(8,ntemp)) 
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

  RETURN
  !
END SUBROUTINE allocate_anharmonic

SUBROUTINE allocate_debye()

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
