!
! Copyright (C) 2025 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE remove_internal_coordinates()
!---------------------------------------------------------------------------
!
USE thermo_mod,            ONLY : tot_ngeo, celldm_geo_eos, uint_geo, &
                                  tot_ngeo_eos, energy_geo
USE ions_base,             ONLY : nat
USE temperature,           ONLY : ntemp, temp
USE control_atomic_pos,    ONLY : linternal_thermo, nint_var, iconstr_internal
USE thermodynamics,        ONLY : ph_free_ener, ph_e0, ph_ener, ph_entropy, &
                                  ph_ce, uint_geo_eos_t, tau_geo_eos_t,     &
                                  ph_free_ener_eos, ph_e0_eos, ph_ener_eos, &
                                  ph_entropy_eos, ph_ce_eos
USE control_eldos,         ONLY : lel_free_energy
USE el_thermodynamics,     ONLY : el_free_ener, el_ener, el_entr, el_ce,      &
                                  el_free_ener_eos, el_ener_eos, el_entr_eos, &
                                  el_ce_eos
IMPLICIT NONE

INTEGER :: igeom, itemp
!
!  find the minimum internal parameters for each set of external 
!  parameters
!
   IF (linternal_thermo) THEN
      CALL redefine_u_min_t(energy_geo, ph_free_ener, uint_geo, &
                        tot_ngeo, uint_geo_eos_t, tot_ngeo_eos)
!
!  interpolate the vibrational free energy, energy, entropy and heat capacity
!  at the minimum of the internal parameters.
!
      CALL redefine_thermo_t(ph_free_ener, ph_e0, ph_ener, ph_entropy, &
              ph_ce,  uint_geo, tot_ngeo, ph_free_ener_eos, ph_e0_eos, &
              ph_ener_eos, ph_entropy_eos, ph_ce_eos, uint_geo_eos_t,  &
              tot_ngeo_eos)
!
!  interpolate the electronic free energy, energy, entropy and heat capacity
!  at the minimum of the internal parameters.
!
      IF (lel_free_energy) &
         CALL redefine_thermo_t(el_free_ener, el_ener, el_entr, el_ce, &
              uint_geo, tot_ngeo, el_free_ener_eos, el_ener_eos,       &
              el_entr_eos, el_ce_eos, uint_geo_eos_t, tot_ngeo_eos)
!
!  for each temperature find the atomic positions from the internal 
!  parameters at each set of external parameters.
!
      DO igeom = 1, tot_ngeo_eos
         DO itemp=1, ntemp
            CALL internal_to_tau(celldm_geo_eos(1,igeom),              &
                          tau_geo_eos_t(1,1,itemp,igeom),              &
                          uint_geo_eos_t(1,itemp, igeom),              &
                          nat, nint_var, iconstr_internal, 1)
         ENDDO
      ENDDO
   ELSE
!
!   If there is no set of internal parameters just copy the thermodynamic
!   variables in those suited to compute eos.
!
      ph_free_ener_eos(1:ntemp,1:tot_ngeo_eos) = &
                                           ph_free_ener(1:ntemp,1:tot_ngeo)
      ph_ener_eos(1:ntemp,1:tot_ngeo_eos) = ph_ener(1:ntemp,1:tot_ngeo)
      ph_entropy_eos(1:ntemp,1:tot_ngeo_eos) = ph_entropy(1:ntemp,1:tot_ngeo)
      ph_ce_eos(1:ntemp,1:tot_ngeo_eos) = ph_ce(1:ntemp,1:tot_ngeo)
      ph_e0_eos(1:tot_ngeo_eos) = ph_e0(1:tot_ngeo)
      IF (lel_free_energy) THEN
         el_free_ener_eos(1:ntemp,1:tot_ngeo_eos) = &
                                           el_free_ener(1:ntemp,1:tot_ngeo)
         el_ener_eos(1:ntemp,1:tot_ngeo_eos) = el_ener(1:ntemp,1:tot_ngeo)
         el_entr_eos(1:ntemp,1:tot_ngeo_eos) = el_entr(1:ntemp,1:tot_ngeo)
         el_ce_eos(1:ntemp,1:tot_ngeo_eos) = el_ce(1:ntemp,1:tot_ngeo)
      ENDIF
   ENDIF

   RETURN
END SUBROUTINE remove_internal_coordinates

!---------------------------------------------------------------------------
SUBROUTINE remove_internal_coordinates_ph()
!---------------------------------------------------------------------------
!
USE thermo_mod,            ONLY : tot_ngeo, celldm_geo_eos, uint_geo, &
                                  tot_ngeo_eos, energy_geo
USE ions_base,             ONLY : nat
USE temperature,           ONLY : ntemp
USE control_atomic_pos,    ONLY : linternal_thermo, nint_var, iconstr_internal
USE ph_freq_thermodynamics,ONLY : phf_free_ener, phf_e0, phf_ener,      &
                                  phf_entropy, phf_ce, uintf_geo_eos_t, &
                                  tauf_geo_eos_t, phf_free_ener_eos,    &
                                  phf_e0_eos, phf_ener_eos,             &
                                  phf_entropy_eos, phf_ce_eos
USE control_eldos,         ONLY : lel_free_energy
USE el_thermodynamics,     ONLY : el_free_ener, el_ener, el_entr, el_ce, &
                                  elf_free_ener_eos, elf_ener_eos,       &
                                  elf_entr_eos, elf_ce_eos
IMPLICIT NONE

INTEGER :: igeom, itemp
!
!  find the minimum internal parameters for each set of external 
!  parameters
!
   IF (linternal_thermo) THEN
      CALL redefine_u_min_t(energy_geo, phf_free_ener, uint_geo, &
                        tot_ngeo, uintf_geo_eos_t, tot_ngeo_eos)
!
!  interpolate the vibrational free energy, energy, entropy and heat capacity
!  at the minimum of the internal parameters.
!
      CALL redefine_thermo_t(phf_free_ener, phf_e0, phf_ener, phf_entropy, &
              phf_ce,  uint_geo, tot_ngeo, phf_free_ener_eos, phf_e0_eos,  &
              phf_ener_eos, phf_entropy_eos, phf_ce_eos, uintf_geo_eos_t,  &
              tot_ngeo_eos)
!
!  interpolate the electronic free energy, energy, entropy and heat capacity
!  at the minimum of the internal parameters.
!
      IF (lel_free_energy) &
         CALL redefine_thermo_t(el_free_ener, el_ener, el_entr, el_ce,  &
              uint_geo, tot_ngeo, elf_free_ener_eos, elf_ener_eos,      &
              elf_entr_eos, elf_ce_eos, uintf_geo_eos_t, tot_ngeo_eos)
!
!  for each temperature find the atomic positions from the internal 
!  parameters at each set of external parameters.
!
      DO igeom = 1, tot_ngeo_eos
         DO itemp=1, ntemp
            CALL internal_to_tau(celldm_geo_eos(1,igeom),               &
                          tauf_geo_eos_t(1,1,itemp,igeom),              &
                          uintf_geo_eos_t(1,itemp, igeom),              &
                          nat, nint_var, iconstr_internal, 1)
         ENDDO
      ENDDO
   ELSE
!
!   If there is no set of internal parameters just copy the thermodynamic
!   variables in those suited to compute eos.
!
      phf_free_ener_eos(1:ntemp,1:tot_ngeo_eos) = &
                                           phf_free_ener(1:ntemp,1:tot_ngeo)
      phf_ener_eos(1:ntemp,1:tot_ngeo_eos) = phf_ener(1:ntemp,1:tot_ngeo)
      phf_entropy_eos(1:ntemp,1:tot_ngeo_eos) = phf_entropy(1:ntemp,1:tot_ngeo)
      phf_ce_eos(1:ntemp,1:tot_ngeo_eos) = phf_ce(1:ntemp,1:tot_ngeo)
      IF (lel_free_energy) THEN
         elf_free_ener_eos(1:ntemp,1:tot_ngeo_eos) = &
                                           el_free_ener(1:ntemp,1:tot_ngeo)
         elf_ener_eos(1:ntemp,1:tot_ngeo_eos) = el_ener(1:ntemp,1:tot_ngeo)
         elf_entr_eos(1:ntemp,1:tot_ngeo_eos) = el_entr(1:ntemp,1:tot_ngeo)
         elf_ce_eos(1:ntemp,1:tot_ngeo_eos) = el_ce(1:ntemp,1:tot_ngeo)
      ENDIF
   ENDIF

   RETURN
END SUBROUTINE remove_internal_coordinates_ph
