!
! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE write_anhar_anis()
!-----------------------------------------------------------------------
!
!   This routine computes the thermal expansion tensor and other
!   related quantities for anisotropic solids. (phdos solids)
!
USE kinds,          ONLY : DP
USE temperature,    ONLY : ntemp, temp
USE thermo_mod,     ONLY : ibrav_geo
USE thermo_sym,     ONLY : laue
USE thermodynamics, ONLY : ph_e0, ph_ce, ph_b_fact, ph_ener, ph_free_ener, &
                           ph_entropy
USE el_thermodynamics, ONLY : el_ener, el_free_ener, el_entr, el_ce
USE anharmonic,     ONLY : alpha_anis_t, vmin_t, b0_t, celldm_t, beta_t, &
                           gamma_t, cv_t, ce_t, cp_t, ener_t, free_ener_t, &
                           entropy_t, b0_s, cpmce_anis, el_cons_t, el_comp_t, &
                           free_e_min_t, bths_t, ggamma_t, &
                           bfact_t, el_cons_s, el_comp_s, macro_el_t, &
                           macro_el_s, v_t, v_s, density_t, csmct_t
USE initial_conf,   ONLY : ibrav_save
USE control_thermo, ONLY : with_eigen
USE control_eldos,  ONLY : lel_free_energy
USE control_elastic_constants, ONLY : lelastic
USE elastic_constants, ONLY : compute_elastic_compliances, &
                              write_el_cons_on_file, print_macro_elasticity, &
                              write_macro_el_on_file, print_sound_velocities,& 
                              write_sound_on_file
USE isoentropic,    ONLY : isostress_heat_capacity, thermal_stress,      &
                           gen_average_gruneisen, isoentropic_elastic_constants
USE el_anharmonic,  ONLY : el_energy_t, el_free_energy_t, el_entropy_t, &
                           el_ce_t
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : meta_ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm, ibrav
INTEGER :: find_free_unit
REAL(DP) :: compute_omega_geo

REAL(DP) :: e0

CALL compute_alpha_anis(celldm_t, alpha_anis_t, temp, ntemp, ibrav_save)
!
!  compute the volume as a function of temperature
!
DO itemp=1,ntemp
   vmin_t(itemp)=compute_omega_geo(ibrav_save, celldm_t(1,itemp))
   CALL compute_density(vmin_t(itemp),density_t(itemp),.TRUE.)
ENDDO
!
!  compute the volume thermal expansion as the derivative of the volume
!  with respect to temperature
!
CALL compute_beta(vmin_t, beta_t, temp, ntemp)

!
!   here the heat capacity, vibrational energy, entropy, zero point energy 
!

CALL interpolate_thermo(vmin_t, celldm_t, ph_ce, ce_t)

CALL interpolate_thermo(vmin_t, celldm_t, ph_ener, ener_t)

CALL interpolate_thermo(vmin_t, celldm_t, ph_free_ener, free_ener_t)

CALL interpolate_thermo(vmin_t, celldm_t, ph_entropy, entropy_t)
                                                       
CALL interpolate_e0(vmin_t, celldm_t, ph_e0, e0)

IF (lel_free_energy) THEN
   CALL interpolate_thermo(vmin_t, celldm_t, el_ce, el_ce_t)
   ce_t = ce_t + el_ce_t
   CALL interpolate_thermo(vmin_t, celldm_t, el_ener, el_energy_t)
   ener_t = ener_t + el_energy_t
   CALL interpolate_thermo(vmin_t, celldm_t, el_free_ener, el_free_energy_t)
   free_ener_t = free_ener_t + el_free_energy_t
   CALL interpolate_thermo(vmin_t, celldm_t, el_entr, el_entropy_t)
   entropy_t = entropy_t + el_entropy_t
ENDIF


IF (lelastic) THEN
   CALL isostress_heat_capacity(vmin_t,el_cons_t,alpha_anis_t,temp, &
                                                         cpmce_anis,ntemp)
   cp_t=ce_t+cpmce_anis
   CALL compute_cv_bs_g(beta_t, vmin_t, b0_t, cv_t, cp_t, b0_s, gamma_t)

   CALL thermal_stress(el_cons_t,alpha_anis_t,bths_t,ntemp)
   CALL gen_average_gruneisen(vmin_t,bths_t,cv_t,ggamma_t,ntemp)

   CALL isoentropic_elastic_constants(vmin_t,bths_t,cv_t,temp,csmct_t,ntemp)
   el_cons_s=el_cons_t + csmct_t
   DO itemp=2, ntemp-1
      CALL compute_elastic_compliances(el_cons_s(:,:,itemp), &
                                                   el_comp_s(:,:,itemp))
     !
     ! Compute macro-elasticity variables and sound velocities using adiabatic 
     ! elastic constants end elastic compliances just computed vs. temperature 
     !
      CALL print_macro_elasticity(ibrav_save, el_cons_s(:,:,itemp), &
                       el_comp_s(:,:,itemp), macro_el_s(:,itemp), .FALSE.)

      CALL print_sound_velocities(ibrav_save, el_cons_s(:,:,itemp), &
           el_comp_s(:,:,itemp), density_t(itemp), v_s(1,itemp), v_s(2,itemp),&
                                                      v_s(3,itemp),.TRUE.)
      !
      ! Compute macro-elasticity variables and sound velocities using 
      ! isothermal 
      ! elastic constants and elastic compliances vs. temperature computed in 
      ! previous routines (set_elastic_constant_t and write_elastic_t_qha).
      !
      CALL print_macro_elasticity(ibrav_save, el_cons_t(:,:,itemp), &
                       el_comp_t(:,:,itemp), macro_el_t(:,itemp), .FALSE.)

      CALL print_sound_velocities(ibrav_save, el_cons_t(:,:,itemp), &
           el_comp_t(:,:,itemp), density_t(itemp), v_t(1,itemp), v_t(2,itemp),&
                                                   v_t(3,itemp),.TRUE.)
   ENDDO
ENDIF

ibrav=ibrav_geo(1)
IF (meta_ionode) THEN
!
!   here we write the anharmonic quantities calculated from the phonon dos
!   on disk
!
   filename='anhar_files/'//TRIM(flanhar)
   CALL add_pressure(filename)

   CALL write_ener_beta(temp, vmin_t, free_e_min_t, beta_t, ntemp, filename)
!
!  Here we write on output the celldm parameters and their derivative
!  with respect to temperature. 
!
   filename='anhar_files/'//TRIM(flanhar)//'.celldm'
   CALL add_pressure(filename)

   CALL write_alpha_anis(ibrav_save, celldm_t, alpha_anis_t, temp, ntemp, &
                                                              filename, 0 )

   filename="anhar_files/"//TRIM(flanhar)//'.therm'
   CALL add_pressure(filename)
   CALL write_thermo_anharm(temp, ntemp, e0, ener_t, free_ener_t, & 
                                                  entropy_t, ce_t, filename )
!
!   here auxiliary quantities calculated from the phonon dos
!
   IF (lelastic) THEN
      !
      !   here the bulk modulus and the gruneisen parameter
      !
      filename="anhar_files/"//TRIM(flanhar)//'.bulk'
      CALL add_pressure(filename)

      CALL write_bulk_anharm(temp, b0_t, b0_s, ntemp, filename)
!
!   the heat capacities
!
      filename="anhar_files/"//TRIM(flanhar)//'.heat'
      CALL add_pressure(filename)

      CALL write_heat_anhar(temp, ce_t, cv_t, cp_t, ntemp, filename)

!  Here we write on output the anharmonic properties computed for
!  anisotropic solids, using the thermal expansion tensor, as opposed
!  to the volume thermal expansion used in the file heat
!
      filename='anhar_files/'//TRIM(flanhar)//'.heat_anis'
      CALL add_pressure(filename)
      CALL write_heat_anhar_anis(temp, ce_t, cv_t, cp_t, ntemp, filename)
!
!   Here the thermal stresses
!
      filename='anhar_files/'//TRIM(flanhar)//'.tstress'
      CALL add_pressure(filename)
      CALL write_thermal_stress(temp, bths_t, ntemp, filename)
!
!   Here the average Gruneisen paramater and the quantities that contribute
!   to it
!
      filename='anhar_files/'//TRIM(flanhar)//'.gamma'
      CALL add_pressure(filename)
      CALL write_gamma_anharm(temp, gamma_t, cv_t, beta_t, b0_t, ntemp, &
                                                                 filename)
!
!   Here the generalized Gruneisen parameters
!
      filename='anhar_files/'//TRIM(flanhar)//'.ggamma'
      CALL add_pressure(filename)
      CALL write_generalized_gamma(temp, ggamma_t, ntemp, filename)
   ELSE
!
!   only the interpolated heat capacity is available
!
      filename="anhar_files/"//TRIM(flanhar)//'.heat'
      CALL add_pressure(filename)
      CALL write_heat_anhar_small(temp, ce_t, ntemp, filename)
   ENDIF
ENDIF
IF (lelastic) THEN
!
!   Here the elastic constants at constant entropy
!
   filename='anhar_files/'//TRIM(flanhar)//'.el_cons_s'
   CALL add_pressure(filename)
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_cons_s, &
                                                  b0_s, filename, 0)
!
!  and here the elastic compliances
!
   filename='anhar_files/'//TRIM(flanhar)//'.el_comp_s'
   CALL add_pressure(filename)
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_comp_s, &
                                                     b0_s, filename, 1)
!
!   Isothermal macro-elasticity variables
!
   filename='anhar_files/'//TRIM(flanhar)//'.macro_el'
   CALL add_pressure(filename)
   CALL write_macro_el_on_file(temp, ntemp, macro_el_t, filename, 0)
!
!   Adiabatic macro-elasticity variables
!
   filename='anhar_files/'//TRIM(flanhar)//'.macro_el_s'
   CALL add_pressure(filename)
   CALL write_macro_el_on_file(temp, ntemp, macro_el_s, filename, 0)
!
!   Isothermal sound velocities
!
   filename='anhar_files/'//TRIM(flanhar)//'.sound_vel'
   CALL add_pressure(filename)
   CALL write_sound_on_file(temp, ntemp, v_t, filename, 0)
!
!   Adiabatic sound velocities
!
   filename='anhar_files/'//TRIM(flanhar)//'.sound_vel_s'
   CALL add_pressure(filename)
   CALL write_sound_on_file(temp, ntemp, v_s, filename, 0)
ENDIF

IF (with_eigen) THEN
   CALL interpolate_b_fact_anis(celldm_t, ph_b_fact, bfact_t)
   filename="anhar_files/"//TRIM(flanhar)//'.anis'
   CALL add_pressure(filename)
   CALL write_anharm_bfact(temp, bfact_t, ntemp, filename)
END IF

RETURN
END SUBROUTINE write_anhar_anis
!
!-----------------------------------------------------------------------
SUBROUTINE write_ph_freq_anhar_anis()
!-----------------------------------------------------------------------
!
!   This routine computes the thermal expansion tensor and other related
!   quantities for anisotropic solids using the free energy calculated 
!   from direct integration of the phonon frequencies.
!
USE kinds,          ONLY : DP
USE temperature,    ONLY : ntemp, temp
USE thermo_mod,     ONLY : ibrav_geo
USE thermo_sym,     ONLY : laue
USE ph_freq_thermodynamics, ONLY : phf_e0, phf_ce, phf_b_fact, phf_ener, &
                                   phf_free_ener, phf_entropy
USE el_thermodynamics,  ONLY : el_ener, el_free_ener, el_entr, el_ce
USE ph_freq_anharmonic, ONLY : alphaf_anis_t, vminf_t, b0f_t, celldmf_t, &
                               betaf_t, gammaf_t, cvf_t, cef_t, cpf_t,   &
                               enerf_t, free_enerf_t, entropyf_t, b0f_s, &
                               cpmcef_anis, el_consf_t, el_compf_t, &
                               free_e_minf_t, bthsf_t, ggammaf_t, bfactf_t, &
                               el_consf_s, el_compf_s, macro_elf_t, &
                               macro_elf_s, vf_s, vf_t, densityf_t, csmctf_t
USE control_elastic_constants, ONLY : lelasticf
USE elastic_constants, ONLY : compute_elastic_compliances, &
                              write_el_cons_on_file, print_macro_elasticity, &
                              write_macro_el_on_file, print_sound_velocities, &  
                              write_sound_on_file
USE el_anharmonic,  ONLY : el_energyf_t, el_free_energyf_t, el_entropyf_t, &
                           el_cef_t
USE initial_conf,   ONLY : ibrav_save
USE control_thermo, ONLY : with_eigen
USE control_eldos,  ONLY : lel_free_energy
USE isoentropic,    ONLY : isostress_heat_capacity, thermal_stress,      &
                           gen_average_gruneisen, isoentropic_elastic_constants
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : meta_ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit, ibrav
REAL(DP) :: compute_omega_geo, el_con_t(6,6,ntemp)

REAL(DP) :: e0

CALL compute_alpha_anis(celldmf_t, alphaf_anis_t, temp, ntemp, ibrav_save)
!
!  compute the volume as a function of temperature
!
DO itemp=1,ntemp
   vminf_t(itemp)=compute_omega_geo(ibrav_save, celldmf_t(1,itemp))
   CALL compute_density(vminf_t(itemp),densityf_t(itemp),.TRUE.)
ENDDO
!
!  compute the volume thermal expansion as the derivative of the volume
!  with respect to temperature
!
CALL compute_beta(vminf_t, betaf_t, temp, ntemp)

!
!   here the heat capacity, vibrational energy, entropy, zero point energy 
!

CALL interpolate_thermo(vminf_t, celldmf_t, phf_ce, cef_t)

CALL interpolate_thermo(vminf_t, celldmf_t, phf_ener, enerf_t)

CALL interpolate_thermo(vminf_t, celldmf_t, phf_free_ener, free_enerf_t)

CALL interpolate_thermo(vminf_t, celldmf_t, phf_entropy, entropyf_t)

CALL interpolate_e0(vminf_t, celldmf_t, phf_e0, e0)

IF (lel_free_energy) THEN
   CALL interpolate_thermo(vminf_t, celldmf_t, el_ce, el_cef_t)
   cef_t = cef_t + el_cef_t
   CALL interpolate_thermo(vminf_t, celldmf_t, el_ener, el_energyf_t)
   enerf_t = enerf_t + el_energyf_t
   CALL interpolate_thermo(vminf_t, celldmf_t, el_free_ener, el_free_energyf_t)
   free_enerf_t = free_enerf_t + el_free_energyf_t
   CALL interpolate_thermo(vminf_t, celldmf_t, el_entr, el_entropyf_t)
   entropyf_t = entropyf_t + el_entropyf_t
ENDIF


IF (lelasticf) THEN
   CALL isostress_heat_capacity(vminf_t,el_consf_t,alphaf_anis_t,temp,&
                                                         cpmcef_anis,ntemp)
   cpf_t=cef_t+cpmcef_anis
   CALL compute_cv_bs_g(betaf_t, vminf_t, b0f_t, cvf_t, cpf_t, b0f_s, gammaf_t)
   CALL thermal_stress(el_consf_t,alphaf_anis_t,bthsf_t,ntemp)
   CALL gen_average_gruneisen(vminf_t,bthsf_t,cvf_t,ggammaf_t,ntemp)

   CALL isoentropic_elastic_constants(vminf_t,bthsf_t,cvf_t,temp,csmctf_t,&
                                                                     ntemp)
   el_consf_s=el_consf_t + csmctf_t
   DO itemp=2, ntemp-1
      CALL compute_elastic_compliances(el_consf_s(:,:,itemp), &
                                                        el_compf_s(:,:,itemp))

      ! Compute macro-elasticity variables and sound velocities using adiabatic 
      ! elastic constants end elastic compliances just computed vs. temperature 
      CALL print_macro_elasticity(ibrav_save, el_consf_s(:,:,itemp), &
                         el_compf_s(:,:,itemp), macro_elf_s(:,itemp), .FALSE.)

      CALL print_sound_velocities(ibrav_save, el_consf_s(:,:,itemp), &
           el_compf_s(:,:,itemp), densityf_t(itemp), vf_s(1,itemp),    &
           vf_s(2,itemp), vf_s(3,itemp), .TRUE.)


      ! Compute macro-elasticity variables and sound velocities using isothermal 
      ! elastic constants and elastic compliances vs. temperature computed in previous 
      ! routines (set_elastic_constant_t and write_elastic_t_qha).
      CALL print_macro_elasticity(ibrav_save, el_consf_t(:,:,itemp), &
                          el_compf_t(:,:,itemp), macro_elf_t(:,itemp), .FALSE.)

      CALL print_sound_velocities(ibrav_save, el_consf_t(:,:,itemp), &
           el_compf_t(:,:,itemp), densityf_t(itemp), vf_t(1,itemp),    &
           vf_t(2,itemp), vf_t(3,itemp), .TRUE.)

   ENDDO
ENDIF

ibrav=ibrav_geo(1)
IF (meta_ionode) THEN
!
!   here we plot the anharmonic quantities calculated from the phonon dos
!
   iu_therm=find_free_unit()
   filename='anhar_files/'//TRIM(flanhar)//'_ph'
   CALL add_pressure(filename)

   CALL write_ener_beta(temp, vminf_t, free_e_minf_t, betaf_t, ntemp, filename)

!
!  Here we write on output the celldm parameters and their derivative
!  with respect to temperature. 
!
   filename='anhar_files/'//TRIM(flanhar)//'.celldm_ph'
   CALL add_pressure(filename)

   CALL write_alpha_anis(ibrav_save, celldmf_t, alphaf_anis_t, temp, ntemp, &
                                                               filename, 0)

   filename="anhar_files/"//TRIM(flanhar)//'.therm_ph'
   CALL add_pressure(filename)
   CALL write_thermo_anharm(temp, ntemp, e0, enerf_t, free_enerf_t, & 
                                               entropyf_t, cef_t, filename)

!
!   here auxiliary quantities calculated from the phonon dos
!
   IF (lelasticf) THEN
      !
      !   here the bulk modulus and the gruneisen parameter
      !
      filename="anhar_files/"//TRIM(flanhar)//'.bulk_ph'
      CALL add_pressure(filename)

      CALL write_bulk_anharm(temp, b0f_t, b0f_s, ntemp, filename)
!
!   the heat capacities
!
      filename="anhar_files/"//TRIM(flanhar)//'.heat_ph'
      CALL add_pressure(filename)

      CALL write_heat_anhar(temp, cef_t, cvf_t, cpf_t, ntemp, filename)
!
!   Here the thermal stresses
!
      filename='anhar_files/'//TRIM(flanhar)//'.tstress_ph'
      CALL add_pressure(filename)
      CALL write_thermal_stress(temp, bthsf_t, ntemp, filename)
!
!   Here the average Gruneisen paramater and the quantities that contribute
!   to it
!
      filename='anhar_files/'//TRIM(flanhar)//'.gamma_ph'
      CALL add_pressure(filename)
      CALL write_gamma_anharm(temp, gammaf_t, cvf_t, betaf_t, b0f_t, ntemp, &
                                                                 filename)
!
!  Here we write on output the anharmonic properties computed for
!  anisotropic solids, using the thermal expansion tensor, as opposed
!  to the volume thermal expansion used in the file aux
!
      filename='anhar_files/'//TRIM(flanhar)//'.heat_anis_ph'
      CALL add_pressure(filename)
      CALL write_heat_anhar_anis(temp, cef_t, cvf_t, cpf_t, ntemp, filename)
!
!   Here the generalized Gruneisen parameters
!
      filename='anhar_files/'//TRIM(flanhar)//'.ggamma_ph'
      CALL add_pressure(filename)
      CALL write_generalized_gamma(temp, ggammaf_t, ntemp, filename)
   ELSE
!
!   only the interpolated heat capacity is available
!
      filename="anhar_files/"//TRIM(flanhar)//'.heat_ph'
      CALL add_pressure(filename)
      CALL write_heat_anhar_small(temp, cef_t, ntemp, filename)
   ENDIF
ENDIF

IF (lelasticf) THEN
!
!   Here the elastic constants at constant entropy
!
   filename='anhar_files/'//TRIM(flanhar)//'.el_cons_s_ph'
   CALL add_pressure(filename)
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_consf_s, b0f_s, &
                                                              filename, 0)
!
!   and here the elastic compliances at constant entropy
!
   filename='anhar_files/'//TRIM(flanhar)//'.el_comp_s_ph'
   CALL add_pressure(filename)
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_compf_s, b0f_s, &
                                                              filename, 1)
!
!   Isothermal macro-elasticity variables
!
   filename='anhar_files/'//TRIM(flanhar)//'.macro_el_ph'
   CALL add_pressure(filename)
   CALL write_macro_el_on_file(temp, ntemp, macro_elf_t, filename,0)
!
!   Adiabatic macro-elasticity variables
!
   filename='anhar_files/'//TRIM(flanhar)//'.macro_el_s_ph'
   CALL add_pressure(filename)
   CALL write_macro_el_on_file(temp, ntemp, macro_elf_s, filename,0)
!
!   Isothermal sound velocities
!
   filename='anhar_files/'//TRIM(flanhar)//'.sound_vel_ph'
   CALL add_pressure(filename)
   CALL write_sound_on_file(temp, ntemp, vf_t, filename, 0)
!
!   Adiabatic sound velocities
!
   filename='anhar_files/'//TRIM(flanhar)//'.sound_vel_s_ph'
   CALL add_pressure(filename)
   CALL write_sound_on_file(temp, ntemp, vf_s, filename, 0)
ENDIF

IF (with_eigen) THEN
   CALL interpolate_b_fact_anis(celldmf_t, phf_b_fact, bfactf_t)
   filename="anhar_files/"//TRIM(flanhar)//'.anis_ph'
   CALL add_pressure(filename)
   CALL write_anharm_bfact(temp, bfactf_t, ntemp, filename)
END IF

RETURN
END SUBROUTINE write_ph_freq_anhar_anis
!
!-----------------------------------------------------------------------
SUBROUTINE write_grun_anhar_anis()
!-----------------------------------------------------------------------
!
!  This routine computes the thermal expansion and other related quantities
!  using the Gruneisen parameters.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE ions_base,      ONLY : nat
USE cell_base,      ONLY : ibrav
USE thermo_mod,     ONLY : ngeo, no_ph, reduced_grid
USE temperature,    ONLY : ntemp, temp
USE control_grun,   ONLY : vgrun_t, celldm_grun_t, b0_grun_t, lb0_t
USE ph_freq_thermodynamics, ONLY : ph_freq_save
USE grun_anharmonic, ONLY : alpha_an_g, grun_gamma_t, done_grun,          &
                           cp_grun_t, b0_grun_s, betab, grun_cpmce_anis,  &
                           cv_grun_t, ce_grun_t, lelastic_grun,           &
                           el_cons_grun_t, el_comp_grun_t, lelastic_grun, &
                           poly_degree_grun, p_grun_p2, poly_grun_red
USE polyfit_mod,    ONLY : compute_poly, compute_poly_deriv
USE ph_freq_module, ONLY : thermal_expansion_ph, ph_freq_type,  &
                           destroy_ph_freq, init_ph_freq
USE lattices,       ONLY : compress_celldm, crystal_parameters
USE quadratic_surfaces, ONLY : evaluate_fit_quadratic,      &
                               evaluate_quadratic_grad
USE isoentropic,    ONLY : isostress_heat_capacity
USE control_dosq,   ONLY : nq1_d, nq2_d, nq3_d
USE polynomial,     ONLY : clean_poly
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : meta_ionode, stdout
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm, i, nq, imode, iq, nvar, nwork
INTEGER :: itens, jtens, startq, lastq, iq_eff
TYPE(ph_freq_type) :: ph_freq    ! the frequencies at the volumes at
                                 ! which the gruneisen parameters are 
                                 ! calculated
TYPE(ph_freq_type), ALLOCATABLE :: ph_grun(:)  ! the gruneisen parameters 
                                 ! recomputed at each temperature at the 
                                 ! geometry corresponding to that temperature
REAL(DP) :: cm(6), aux(6), alpha_aux(6), alpha(6), f, g, vm
REAL(DP), ALLOCATABLE :: grad(:), x(:)
INTEGER :: compute_nwork, central_geo, find_free_unit

done_grun=.FALSE.
!
!  Not implemented cases
!
IF ( ibrav<1 .OR. ibrav>11 ) THEN
   WRITE(stdout,'(5x,"Thermal expansions from Gruneisen parameters &
                     & not available")' )
   RETURN
END IF
!
!  If the elastic constants are not available, this calculation cannot be done
!
IF (.NOT.lelastic_grun ) THEN
   WRITE(stdout,'(5x,"The elastic constants are needed to compute ")')
   WRITE(stdout,'(5x,"thermal expansions from Gruneisen parameters")')
   RETURN
ENDIF

WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(5x,"Computing the anharmonic properties from &
                                                   &Gruneisen parameters")')
WRITE(stdout,'(5x,"Writing on file anhar_files/",a)') TRIM(flanhar)// &
                                                                  '.aux_grun'
WRITE(stdout,'(2x,76("+"),/)')
!
! divide the q points among the processors. Each processor has only a part
! of the q points and computes the contribution of these points to the
! anharmonic properties
!
CALL find_central_geo(ngeo, no_ph, central_geo)
nq=ph_freq_save(central_geo)%nq
startq=ph_freq_save(central_geo)%startq
lastq=ph_freq_save(central_geo)%lastq
CALL init_ph_freq(ph_freq, nat, nq1_d, nq2_d, nq3_d, startq, lastq,  &
                                                               nq, .FALSE.)
ph_freq%wg=ph_freq_save(central_geo)%wg
!
! now allocate space for each set of gruneisen parameters
!
nvar=crystal_parameters(ibrav)
nwork=compute_nwork()
ALLOCATE(ph_grun(nvar))
ALLOCATE(grad(nvar))
ALLOCATE(x(nvar))
DO i=1, nvar
   CALL init_ph_freq(ph_grun(i), nat, nq1_d, nq2_d, nq3_d, startq, &
                                                    lastq, nq, .FALSE.)
END DO
!
!  computes the anharmonic quantities at each temperature
!
alpha_an_g=0.0_DP
DO itemp = 1, ntemp
   IF (MOD(itemp,30)==0) &
             WRITE(6,'(5x,"Computing temperature ", i5, " / ",&
       & i5, 4x," T=",f12.2," K")') itemp, ntemp, temp(itemp)
!
!  Computes the volume at which the Gruneisen parameters and the frequencies
!  are interpolated
!
   cm(:)=celldm_grun_t(:,itemp)
   vm = vgrun_t(itemp)

   CALL compress_celldm(cm,x,nvar,ibrav)
!
!  compute the frequencies and the gruneisen parameters from the interpolating 
!  polynomial
!
   ph_freq%nu= 0.0_DP
   DO i=1,nvar
      ph_grun(i)%nu= 0.0_DP
   END DO
   iq_eff=0
   DO iq=startq, lastq
      iq_eff=iq_eff+1
      IF (reduced_grid) THEN
         DO i=1,nvar
            DO imode=1,3*nat
               CALL compute_poly(x(i), poly_degree_grun, &
                                            poly_grun_red(1,imode,i,iq),f)
!
!  this function gives the derivative with respect to x(i) multiplied by x(i)
!
               CALL compute_poly_deriv(x(i), poly_degree_grun, &
                                             poly_grun_red(1,imode,i,iq),g)

               ph_freq%nu(imode,iq_eff) = f
               IF (f > 0.0_DP ) THEN
                  ph_grun(i)%nu(imode,iq_eff)= - g / f 
               ELSE
                  ph_grun(i)%nu(imode,iq_eff)=0.0_DP
               END IF
            ENDDO
         ENDDO
      ELSE
         DO imode=1,3*nat
            CALL evaluate_fit_quadratic(nvar,x,f,p_grun_p2(imode,iq))
            CALL evaluate_quadratic_grad(nvar,x,grad,p_grun_p2(imode,iq))
            ph_freq%nu(imode,iq_eff) = f 
            IF (f > 0.0_DP ) THEN
               DO i=1,nvar
                  ph_grun(i)%nu(imode,iq_eff)=-grad(i) / f
               END DO
            ELSE
               DO i=1,nvar
                  ph_grun(i)%nu(imode,iq_eff)=0.0_DP
               END DO
            END IF
         END DO
      ENDIF
   END DO
!
!  Compute thermal expansion from Gruneisen parameters. Loop over the number 
!  of independent crystal parameters for this Bravais lattice
!  alpha calculated by thermal_expansion_ph is not multiplied by the elastic 
!  compliances
!
   DO i=1,nvar
      CALL thermal_expansion_ph(ph_freq, ph_grun(i), temp(itemp), &
                                alpha_aux(i), ce_grun_t(itemp))
   END DO
!
!  Here convert from derivatives with respect to crystal parameters to
!  derivatives with respect to strain. 
!
   CALL convert_ac_alpha(alpha_aux, alpha, cm, ibrav)
!
!  To get the thermal expansion we need to multiply by the elastic compliances
!
   aux=0.0_DP
   DO itens=1,6
      DO jtens=1,6
         aux(itens)=aux(itens) + el_comp_grun_t(itens,jtens,itemp)*alpha(jtens)
      END DO
   END DO
   alpha_an_g(:,itemp) = aux(:) * ry_kbar / vm
END DO
CALL mp_sum(alpha_an_g, world_comm)
CALL mp_sum(ce_grun_t, world_comm)
!
!  compute the volume thermal expansion as the trace of the thermal expansion 
!  tensor
!
betab(:)=alpha_an_g(1,:)+alpha_an_g(2,:)+alpha_an_g(3,:)
!
!  computes the other anharmonic quantities

CALL isostress_heat_capacity(vgrun_t,el_cons_grun_t,alpha_an_g,temp,&
                                                    grun_cpmce_anis,ntemp)
cp_grun_t = ce_grun_t + grun_cpmce_anis
CALL compute_cv_bs_g(betab, vgrun_t, b0_grun_t, cv_grun_t, &
                                     cp_grun_t, b0_grun_s, grun_gamma_t)

IF (meta_ionode) THEN
!
!   here quantities calculated from the gruneisen parameters
!
   filename='anhar_files/'//TRIM(flanhar)//'.celldm_grun'
   CALL add_pressure(filename)

   CALL write_alpha_anis(ibrav, celldm_grun_t, alpha_an_g, temp, ntemp, &
                                                                filename, 0 )
!
!  Here the average Gruneisen parameter and the quantities that form it
!
   filename="anhar_files/"//TRIM(flanhar)//'.aux_grun'
   CALL add_pressure(filename)
   CALL write_aux_grun(temp, betab, cp_grun_t, cv_grun_t, b0_grun_s, &
                                               b0_grun_t, ntemp, filename)
!
!  Here the average Gruneisen paramater and the quantities that contribute
!  to it
!
   filename='anhar_files/'//TRIM(flanhar)//'.gamma_grun'
   CALL add_pressure(filename)
   CALL write_gamma_anharm(temp, grun_gamma_t, cv_grun_t, betab, &
                                               b0_grun_t, ntemp, filename)
!
!  Here we write on output the anharmonic properties computed for
!  anisotropic solids, using the thermal expansion tensor, as opposed
!  to the volume thermal expansion used in the file aux
!
   filename='anhar_files/'//TRIM(flanhar)//'.heat_anis_grun'
   CALL add_pressure(filename)
   CALL write_heat_anhar_anis(temp, ce_grun_t, cv_grun_t, cp_grun_t, &
                                                           ntemp, filename)
ENDIF

done_grun=.TRUE.

CALL destroy_ph_freq(ph_freq)
DO i=1,nvar
   CALL destroy_ph_freq(ph_grun(i))
END DO
DEALLOCATE(ph_grun)

IF (.NOT.reduced_grid) THEN
   DO imode=1,3*nat
      DO iq=startq, lastq
         CALL clean_poly(p_grun_p2(imode,iq))
      ENDDO
   ENDDO
   DEALLOCATE(p_grun_p2)
ENDIF

DEALLOCATE(grad)
DEALLOCATE(x)

RETURN
END SUBROUTINE write_grun_anhar_anis
!
!-----------------------------------------------------------------------
SUBROUTINE write_anhar_anis_pt()
!-----------------------------------------------------------------------
!
!   This routine computes the thermal expansion tensor and other related
!   quantities for anisotropic solids as a function of temperature 
!   for several pressures.
!
USE kinds,          ONLY : DP
USE temperature,    ONLY : ntemp, temp, itemp300
USE thermo_sym,     ONLY : laue
USE thermodynamics, ONLY : ph_e0
USE anharmonic_pt,  ONLY : alpha_anis_pt, vmin_pt, b0_pt, celldm_pt, beta_pt, &
                           gamma_pt, cv_pt, ce_pt, cp_pt, ener_pt, &
                           free_ener_pt, entr_pt, b0_s_pt, &
                           cpmce_anis_pt, el_cons_pt, el_comp_pt, &
                           emin_pt, bths_pt, ggamma_pt, &
                           el_cons_pt, el_cons_s_pt, el_comp_pt, &
                           el_comp_s_pt, macro_el_pt, macro_el_s_pt, &
                           v_pt, v_s_pt, density_pt, csmct_pt
USE initial_conf,   ONLY : ibrav_save
USE control_eldos,  ONLY : lel_free_energy
USE control_elastic_constants, ONLY : lelastic
USE control_pressure, ONLY : press, npress_plot, ipress_plot
USE elastic_constants, ONLY : compute_elastic_compliances, &
                              write_el_cons_on_file, print_macro_elasticity, &
                              write_macro_el_on_file, print_sound_velocities,& 
                              write_sound_on_file
USE isoentropic,    ONLY : isostress_heat_capacity, thermal_stress,      &
                           gen_average_gruneisen, isoentropic_elastic_constants
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : meta_ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, ipress, ipressp, iu_therm
INTEGER :: find_free_unit
REAL(DP) :: compute_omega_geo, aux(ntemp)

REAL(DP) :: e0

DO ipressp=1, npress_plot
   ipress=ipress_plot(ipressp)
   CALL compute_alpha_anis(celldm_pt(:,:,ipressp), &
                     alpha_anis_pt(:,:,:,ipressp), temp, ntemp, ibrav_save)
!
!  compute the volume as a function of temperature
!
   DO itemp=1,ntemp
      vmin_pt(itemp,ipressp)=compute_omega_geo(ibrav_save, &
                             celldm_pt(1,itemp,ipressp))
      CALL compute_density(vmin_pt(itemp,ipressp), density_pt(itemp,ipressp),&
                                                                     .FALSE.)
      CALL interpolate_e0(vmin_pt(:,ipressp), celldm_pt(:,:,ipressp), &
                                                                   ph_e0, e0)
   ENDDO
!
!  compute the volume thermal expansion as the derivative of the volume
!  with respect to temperature
!
   CALL compute_beta(vmin_pt(:,ipressp), beta_pt(:,ipressp), temp, ntemp)


   IF (lelastic) THEN
      CALL isostress_heat_capacity(vmin_pt(:,ipressp), &
            el_cons_pt(:,:,:,ipressp), alpha_anis_pt(:,:,:,ipressp), &
            temp, cpmce_anis_pt(:,ipressp),ntemp)
      cp_pt(:,ipressp)=ce_pt(:,ipressp)+cpmce_anis_pt(:,ipressp)
      CALL compute_cv_bs_g(beta_pt(:,ipressp), vmin_pt(:,ipressp), &
          b0_pt(:,ipressp), cv_pt(:,ipressp), cp_pt(:,ipressp), &
          b0_s_pt(:,ipressp), gamma_pt(:,ipressp))

      CALL thermal_stress(el_cons_pt(:,:,:,ipressp),&
                   alpha_anis_pt(:,:,:,ipressp), bths_pt(:,:,:,ipressp),ntemp)
      CALL gen_average_gruneisen(vmin_pt(:,ipressp),bths_pt(:,:,:,ipressp), &
                         cv_pt(:,ipressp),ggamma_pt(:,:,:,ipressp),ntemp)
      CALL isoentropic_elastic_constants(vmin_pt(:,ipressp), &
                    bths_pt(:,:,:,ipressp),cv_pt(:,ipressp),temp,&
                                           csmct_pt(:,:,:,ipressp),ntemp)
      el_cons_s_pt(:,:,:,ipressp)=el_cons_pt(:,:,:,ipressp) + &
                                    csmct_pt(:,:,:,ipressp)
      DO itemp=2, ntemp-1
         CALL compute_elastic_compliances(el_cons_s_pt(:,:,itemp,ipressp), &
                                         el_comp_s_pt(:,:,itemp,ipressp))
     !
     ! Compute macro-elasticity variables and sound velocities using adiabatic 
     ! elastic constants end elastic compliances just computed vs. temperature 
     !
         CALL print_macro_elasticity(ibrav_save, el_cons_s_pt(:,:,itemp,&
              ipressp),el_comp_s_pt(:,:,itemp,ipressp), &
              macro_el_s_pt(:,itemp,ipressp), .FALSE.)

         CALL print_sound_velocities(ibrav_save, &
              el_cons_s_pt(:,:,itemp,ipressp),   &
              el_comp_s_pt(:,:,itemp,ipressp), density_pt(itemp,ipressp), &
              v_s_pt(1,itemp,ipressp), v_s_pt(2,itemp,ipressp), &
                                           v_s_pt(3,itemp,ipressp),.TRUE.)
      !
      ! Compute macro-elasticity variables and sound velocities using 
      ! isothermal 
      ! elastic constants and elastic compliances vs. temperature computed in 
      ! previous routines (set_elastic_constant_t and write_elastic_t_qha).
      !
         CALL print_macro_elasticity(ibrav_save, &
           el_cons_pt(:,:,itemp,ipressp), el_comp_pt(:,:,itemp,ipressp), &
           macro_el_pt(:,itemp,ipressp),.FALSE.)

         CALL print_sound_velocities(ibrav_save, &
           el_cons_pt(:,:,itemp,ipressp), el_comp_pt(:,:,itemp,ipressp), &
           density_pt(itemp,ipressp),                    &
           v_pt(1,itemp,ipressp), v_pt(2,itemp,ipressp), &
           v_pt(3,itemp,ipressp), .TRUE.)
      ENDDO
   ENDIF

   IF (itemp300 > 0) THEN
      aux(:)=vmin_pt(:,ipressp)/vmin_pt(itemp300,ipressp)
   ELSE
      aux(:)=vmin_pt(:,ipressp)
   ENDIF

   IF (meta_ionode) THEN
!
!   here we write the anharmonic quantities calculated from the phonon dos
!   on disk
!
      filename='anhar_files/'//TRIM(flanhar)//'.press'
      CALL add_value(filename, press(ipress))

      CALL write_ener_beta_vol(temp, vmin_pt(:,ipressp), aux(:), &
                  emin_pt(:,ipressp), beta_pt(:,ipressp), ntemp, filename)

!
!  Here we write on output the celldm parameters and their derivative
!  with respect to temperature. 
!
      filename='anhar_files/'//TRIM(flanhar)//'.celldm_press'
      CALL add_value(filename,press(ipress))

      CALL write_alpha_anis(ibrav_save, celldm_pt(:,:,ipressp), &
                 alpha_anis_pt(:,:,:,ipressp), temp, ntemp, filename, 0 )

      filename="anhar_files/"//TRIM(flanhar)//'.therm_press'
      CALL add_value(filename,press(ipress))
      CALL write_thermo_anharm(temp, ntemp, e0, ener_pt(:,ipressp), &
           free_ener_pt(:,ipressp), entr_pt(:,ipressp), ce_pt(:,ipressp), &
                                                          filename )
!
!   here auxiliary quantities calculated from the phonon dos
!
      IF (lelastic) THEN
      !
      !   here the bulk modulus and the gruneisen parameter
      !
         filename="anhar_files/"//TRIM(flanhar)//'.bulk_press'
         CALL add_value(filename,press(ipress))

         CALL write_bulk_anharm(temp, b0_pt(:,ipressp), b0_s_pt(:,ipressp), &
                                                        ntemp, filename)
!
!   the heat capacities
!
         filename="anhar_files/"//TRIM(flanhar)//'.heat_press'
         CALL add_value(filename,press(ipress))

         CALL write_heat_anhar(temp, ce_pt(:,ipressp), cv_pt(:,ipressp), &
                     cp_pt(:,ipressp), ntemp, filename)

!  Here we write on output the anharmonic properties computed for
!  anisotropic solids, using the thermal expansion tensor, as opposed
!  to the volume thermal expansion used in the file heat
!
         filename='anhar_files/'//TRIM(flanhar)//'.heat_anis.press'
         CALL add_value(filename,press(ipress))
         CALL write_heat_anhar_anis(temp, ce_pt(:,ipressp), cv_pt(:,ipressp),&
                                          cp_pt(:,ipressp), ntemp, filename)
!
!   Here the thermal stresses
!
         filename='anhar_files/'//TRIM(flanhar)//'.tstress_press'
         CALL add_value(filename,press(ipress))
         CALL write_thermal_stress(temp, bths_pt(:,:,:,ipressp), ntemp, &
                                                                  filename)
!
!   Here the average Gruneisen paramater and the quantities that contribute
!   to it
!
         filename='anhar_files/'//TRIM(flanhar)//'.gamma_press'
         CALL add_value(filename,press(ipress))
         CALL write_gamma_anharm(temp, gamma_pt(:,ipressp), cv_pt(:,ipressp),&
                   beta_pt(:,ipressp), b0_pt(:,ipressp), ntemp, filename)
!
!   Here the generalized Gruneisen parameters
!
         filename='anhar_files/'//TRIM(flanhar)//'.ggamma_press'
         CALL add_value(filename,press(ipress))
         CALL write_generalized_gamma(temp, ggamma_pt(:,:,:,ipressp), ntemp, &
                                                                  filename)
      ELSE
!
!   only the interpolated heat capacity is available
!
         filename="anhar_files/"//TRIM(flanhar)//'.heat_press'
         CALL add_value(filename,press(ipress))
         CALL write_heat_anhar_small(temp, ce_pt(:,ipressp), ntemp, filename)
      ENDIF
   ENDIF
   IF (lelastic) THEN
!
!   Here the elastic constants at constant entropy
!
      filename='anhar_files/'//TRIM(flanhar)//'.el_cons_s_press'
      CALL add_value(filename,press(ipress))
      CALL write_el_cons_on_file(temp, ntemp, ibrav_save, laue,          &
                 el_cons_s_pt(:,:,:,ipressp), b0_s_pt(:,ipressp), &
                 filename, 0)
!
!  and here the elastic compliances
!
      filename='anhar_files/'//TRIM(flanhar)//'.el_comp_s_press'
      CALL add_value(filename,press(ipress))
      CALL write_el_cons_on_file(temp, ntemp, ibrav_save, laue, &
           el_comp_s_pt(:,:,:,ipressp), b0_s_pt(:,ipressp), filename, 1)
!
!   Isothermal macro-elasticity variables
!
      filename='anhar_files/'//TRIM(flanhar)//'.macro_el_press'
      CALL add_value(filename,press(ipress))
      CALL write_macro_el_on_file(temp, ntemp, macro_el_pt(:,:,ipressp), &
                                                               filename, 0)
!
!   Adiabatic macro-elasticity variables
!
      filename='anhar_files/'//TRIM(flanhar)//'.macro_el_s_press'
      CALL add_value(filename,press(ipress))
      CALL write_macro_el_on_file(temp, ntemp, macro_el_s_pt(:,:,ipressp), &
                                                                 filename, 0)
!
!   Isothermal sound velocities
!
      filename='anhar_files/'//TRIM(flanhar)//'.sound_vel_press'
      CALL add_value(filename,press(ipress))
      CALL write_sound_on_file(temp, ntemp, v_pt(:,:,ipressp), filename, 0)
!
!   Adiabatic sound velocities
!
      filename='anhar_files/'//TRIM(flanhar)//'.sound_vel_s_press'
      CALL add_value(filename,press(ipress))
      CALL write_sound_on_file(temp, ntemp, v_s_pt(:,:,ipressp), filename, 0)
   ENDIF
ENDDO

RETURN
END SUBROUTINE write_anhar_anis_pt

!-----------------------------------------------------------------------
SUBROUTINE write_anhar_anis_ptt()
!-----------------------------------------------------------------------
!
!   This routine computes the thermal expansion tensor and other 
!   related quantities for anisotropic solids as a function of pressure
!   for selected pressures.
!
USE kinds,          ONLY : DP
USE temperature,    ONLY : ntemp, temp, itemp300
USE thermo_sym,     ONLY : laue
USE thermodynamics, ONLY : ph_e0
USE anharmonic,     ONLY : noelcvg, vmin_t
USE anharmonic_ptt, ONLY : vmin_ptt, vmin_ptt_p1, vmin_ptt_m1, &
                           alpha_anis_ptt, celldm_ptt, celldm_ptt_p1, & 
                           celldm_ptt_m1, beta_ptt, el_cons_ptt,      &
                           cpmce_anis_ptt, cp_ptt, ce_ptt, cv_ptt, b0_ptt, &
                           b0_s_ptt, gamma_ptt, bths_ptt, ggamma_ptt, &
                           el_cons_s_ptt, el_comp_s_ptt, macro_el_s_ptt, &
                           el_cons_ptt, el_comp_ptt, macro_el_ptt, emin_ptt, &
                           v_ptt, v_s_ptt, density_ptt, ener_ptt, &
                           entr_ptt, csmct_ptt
USE el_anharmonic,  ONLY : el_ce_ptt
USE control_eldos,  ONLY : lel_free_energy
USE initial_conf,   ONLY : ibrav_save
USE control_eldos,  ONLY : lel_free_energy
USE control_elastic_constants, ONLY : lelastic
USE control_pressure, ONLY : press, npress
USE temperature,    ONLY : temp, ntemp_plot, itemp_plot
USE elastic_constants, ONLY : compute_elastic_compliances, &
                              write_el_cons_on_file, print_macro_elasticity, &
                              write_macro_el_on_file, print_sound_velocities,& 
                              write_sound_on_file
USE isoentropic,    ONLY : isostress_heat_capacity, thermal_stress_p,      &
                           gen_average_gruneisen_p,                        &
                           isoentropic_elastic_constants_p
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : meta_ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, ipress, ipressp, iu_therm
INTEGER :: find_free_unit
REAL(DP) :: compute_omega_geo, aux(npress)

INTEGER  :: itempp
REAL(DP) :: e0_p(npress)
LOGICAL  :: subtract_el

DO itempp=1, ntemp_plot
   itemp=itemp_plot(itempp)
   CALL compute_alpha_anis_p(celldm_ptt(:,:,itempp), &
                     celldm_ptt_p1(:,:,itempp), celldm_ptt_m1(:,:,itempp), &
                     alpha_anis_ptt(:,:,:,itempp), press, npress, ibrav_save)
!
!  compute the volume as a function of temperature
!
   DO ipress=1,npress
      vmin_ptt(ipress,itempp)=compute_omega_geo(ibrav_save, &
                             celldm_ptt(1,ipress,itempp))
      CALL compute_density(vmin_ptt(ipress,itempp), &
                           density_ptt(ipress,itempp),.FALSE.)
      CALL interpolate_e0_p(vmin_ptt(:,itempp), celldm_ptt(:,:,itempp), &
                                                            ph_e0, e0_p)
      vmin_ptt_p1(ipress,itempp)=compute_omega_geo(ibrav_save, &
                             celldm_ptt_p1(1,ipress,itempp))
      vmin_ptt_m1(ipress,itempp)=compute_omega_geo(ibrav_save, &
                             celldm_ptt_m1(1,ipress,itempp))
   ENDDO
!
!  compute the volume thermal expansion as the derivative of the volume
!  with respect to temperature
!
   CALL compute_beta_p(vmin_ptt(:,itempp), vmin_ptt_p1(:,itempp), &
                       vmin_ptt_m1(:,itempp), beta_ptt(:,itempp), &
                       press, npress)

   IF (lelastic) THEN
      CALL isostress_heat_capacity(vmin_ptt(:,itempp),               &
            el_cons_ptt(:,:,:,itempp), alpha_anis_ptt(:,:,:,itempp), &
            temp, cpmce_anis_ptt(:,itempp),npress)
      cp_ptt(:,itempp)=ce_ptt(:,itempp)+cpmce_anis_ptt(:,itempp)
      subtract_el=(lel_free_energy.AND.noelcvg)
      CALL compute_cp_bs_gp(beta_ptt(:,itempp), vmin_ptt(:,itempp),   &
          b0_ptt(:,itempp), cv_ptt(:,itempp), cp_ptt(:,itempp),       &
          b0_s_ptt(:,itempp), gamma_ptt(:,itempp), el_ce_ptt, itemp,  &
          subtract_el)

      CALL thermal_stress_p(el_cons_ptt(:,:,:,itempp),                &
                   alpha_anis_ptt(:,:,:,itempp), bths_ptt(:,:,:,itempp),npress)
      CALL gen_average_gruneisen_p(vmin_ptt(:,itempp),bths_ptt(:,:,:,itempp), &
                       cv_ptt(:,itempp),ggamma_ptt(:,:,:,itempp),npress)
      CALL isoentropic_elastic_constants_p(vmin_ptt(:,itempp),        &
            bths_ptt(:,:,:,itempp),cv_ptt(:,itempp),temp,             &
                            csmct_ptt(:,:,:,itempp),ntemp, npress,itemp)
      el_cons_s_ptt(:,:,:,itempp)=el_cons_ptt(:,:,:,itempp) + &
                            csmct_ptt(:,:,:,itempp)
      DO ipress=1, npress
         CALL compute_elastic_compliances(el_cons_s_ptt(:,:,ipress,itempp), &
                                          el_comp_s_ptt(:,:,ipress,itempp))
     !
     ! Compute macro-elasticity variables and sound velocities using adiabatic 
     ! elastic constants end elastic compliances just computed vs. temperature 
     !
         CALL print_macro_elasticity(ibrav_save, el_cons_s_ptt(:,:,ipress,&
              itempp),el_comp_s_ptt(:,:,ipress,itempp), &
              macro_el_s_ptt(:,ipress,itempp), .FALSE.)

         CALL print_sound_velocities(ibrav_save, &
              el_cons_s_ptt(:,:,ipress,itempp),   &
              el_comp_s_ptt(:,:,ipress,itempp), density_ptt(ipress,itempp), &
              v_s_ptt(1,ipress,itempp), v_s_ptt(2,ipress,itempp), &
                                           v_s_ptt(3,ipress,itempp),.TRUE.)
     !
     ! Compute macro-elasticity variables and sound velocities using 
     ! isothermal 
     ! elastic constants and elastic compliances vs. temperature computed in 
     ! previous routines (set_elastic_constant_t and write_elastic_t_qha).
     !
         CALL print_macro_elasticity(ibrav_save, &
           el_cons_ptt(:,:,ipress,itempp), el_comp_ptt(:,:,ipress,itempp), &
           macro_el_ptt(:,ipress,itempp),.FALSE.)

         CALL print_sound_velocities(ibrav_save, &
           el_cons_ptt(:,:,ipress,itempp), el_comp_ptt(:,:,ipress,itempp), &
           density_ptt(ipress,itempp),                    &
           v_ptt(1,ipress,itempp), v_ptt(2,ipress,itempp), &
           v_ptt(3,ipress,itempp), .TRUE.)
      ENDDO
   ENDIF

   IF (itemp300 > 0) THEN
      aux(:)=vmin_ptt(:,itempp)/vmin_t(itemp300)
   ELSE
      aux(:)=vmin_ptt(:,itempp)
   ENDIF


   IF (meta_ionode) THEN
!
!   here we write the anharmonic quantities calculated from the phonon dos
!   on disk
!
      filename='anhar_files/'//TRIM(flanhar)//'.temp'
      CALL add_value(filename, temp(itemp))

      CALL write_ener_beta_ptt(press, vmin_ptt(:,itempp),  aux(:),         &
                  emin_ptt(:,itempp), beta_ptt(:,itempp), npress, itemp, &
                  filename)

!
!  Here we write on output the celldm parameters and their derivative
!  with respect to temperature. 
!
      filename='anhar_files/'//TRIM(flanhar)//'.celldm_temp'
      CALL add_value(filename,temp(itemp))

      CALL write_alpha_anis(ibrav_save, celldm_ptt(:,:,itempp), &
                 alpha_anis_ptt(:,:,:,itempp), press, npress, filename, 1)

      filename="anhar_files/"//TRIM(flanhar)//'.therm_temp'
      CALL add_value(filename,temp(itemp))
      CALL write_thermo_anharm_p(press, npress, e0_p, ener_ptt(:,itempp), &
           emin_ptt(:,itempp), entr_ptt(:,itempp), ce_ptt(:,itempp), &
                                         temp, ntemp, itemp, filename )
!
!   here auxiliary quantities calculated from the phonon dos
!
      IF (lelastic) THEN
      !
      !   here the bulk modulus and the gruneisen parameter
      !
         filename="anhar_files/"//TRIM(flanhar)//'.bulk_temp'
         CALL add_value(filename,temp(itemp))

         CALL write_bulk_anharm_ptt(press, b0_ptt(:,itempp), &
                               b0_s_ptt(:,itempp), npress, itemp, filename)
!
!   the heat capacities
!
         filename="anhar_files/"//TRIM(flanhar)//'.heat_temp'
         CALL add_value(filename,temp(itemp))

         CALL write_heat_anhar_t(press, ce_ptt(:,itempp), ce_ptt(:,itempp), &
              cp_ptt(:,itempp), npress, itemp, filename)


!  Here we write on output the anharmonic properties computed for
!  anisotropic solids, using the thermal expansion tensor, as opposed
!  to the volume thermal expansion used in the file heat
!
         filename='anhar_files/'//TRIM(flanhar)//'.heat_anis.temp'
         CALL add_value(filename,temp(itemp))
         CALL write_heat_anhar_anis_p(press, ce_ptt(:,itempp), &
               cv_ptt(:,itempp), cp_ptt(:,itempp), npress, filename)
!
!   Here the thermal stresses
!
         filename='anhar_files/'//TRIM(flanhar)//'.tstress_temp'
         CALL add_value(filename,temp(itemp))
         CALL write_thermal_stress_p(press, bths_ptt(:,:,:,ipressp), npress, &
                                                                  filename)
!
!   Here the average Gruneisen paramater and the quantities that contribute
!   to it
!
         filename='anhar_files/'//TRIM(flanhar)//'.gamma_temp'
         CALL add_value(filename,temp(itemp))
         CALL write_gamma_anharm_t(press, gamma_ptt(:,itempp), &
            ce_ptt(:,itempp), beta_ptt(:,itempp), b0_ptt(:,itempp), &
                                           npress, itemp, filename)
!
!   Here the generalized Gruneisen parameters
!
         filename='anhar_files/'//TRIM(flanhar)//'.ggamma_temp'
         CALL add_value(filename,temp(itemp))
         CALL write_generalized_gamma_p(press, ggamma_ptt(:,:,:,itempp), &
                                                          ntemp, filename)
      ELSE
!
!   only the interpolated heat capacity is available
!
         filename="anhar_files/"//TRIM(flanhar)//'.heat_temp'
         CALL add_value(filename,temp(itemp))
         CALL write_heat_anhar_small_p(press, ce_ptt(:,itempp), &
                                                        npress, filename)
      ENDIF
   ENDIF
   IF (lelastic) THEN
!
!   Here the elastic constants at constant entropy
!
      filename='anhar_files/'//TRIM(flanhar)//'.el_cons_s_temp'
      CALL add_value(filename,temp(itemp))
      CALL write_el_cons_on_file(press, npress, ibrav_save, laue,          &
                 el_cons_s_ptt(:,:,:,itempp), b0_s_ptt(:,itempp), &
                 filename, 2)
!
!  and here the elastic compliances
!
      filename='anhar_files/'//TRIM(flanhar)//'.el_comp_s_temp'
      CALL add_value(filename,temp(itemp))
      CALL write_el_cons_on_file(press, npress, ibrav_save, laue, &
           el_comp_s_ptt(:,:,:,itempp), b0_s_ptt(:,itempp), filename, 3)
!
!   Isothermal macro-elasticity variables
!
      filename='anhar_files/'//TRIM(flanhar)//'.macro_el_temp'
      CALL add_value(filename,temp(itemp))
      CALL write_macro_el_on_file(press, npress, macro_el_ptt(:,:,itempp), &
                                                               filename, 1)
!
!   Adiabatic macro-elasticity variables
!
      filename='anhar_files/'//TRIM(flanhar)//'.macro_el_s_temp'
      CALL add_value(filename,temp(itemp))
      CALL write_macro_el_on_file(press, npress, macro_el_s_ptt(:,:,itempp), &
                                                                 filename, 1)
!
!   Isothermal sound velocities
!
      filename='anhar_files/'//TRIM(flanhar)//'.sound_vel_temp'
      CALL add_value(filename,temp(itemp))
      CALL write_sound_on_file(press, npress, v_ptt(:,:,itempp), filename, 1)
!
!   Adiabatic sound velocities
!
      filename='anhar_files/'//TRIM(flanhar)//'.sound_vel_s_temp'
      CALL add_value(filename,temp(itemp))
      CALL write_sound_on_file(press, npress, v_s_ptt(:,:,itempp), filename, 1)
   ENDIF
ENDDO

RETURN
END SUBROUTINE write_anhar_anis_ptt

!-----------------------------------------------------------------------
SUBROUTINE write_alpha_anis(ibrav, celldmf_t, alpha_t, temp, ntemp, &
                                                         filename, flag)
!-----------------------------------------------------------------------
!
!  This routine writes on file the thermal expansion tensor, as a function
!  temperature (flag=0) or of pressure (flag=1)
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav, ntemp, flag
REAL(DP), INTENT(IN) :: celldmf_t(6,ntemp), alpha_t(6,ntemp), temp(ntemp)
CHARACTER(LEN=*), INTENT(IN) :: filename
CHARACTER(LEN=256) :: label
INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

iu_therm=find_free_unit()

IF (flag==0) THEN
   label='T (K)   '
ELSE
   label='p (kbar)'
ENDIF

OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

IF (ibrav==1 .OR. ibrav==2 .OR. ibrav==3 ) THEN
   WRITE(iu_therm,'("# ",a,"        celldm(1)         alpha_xx(x10^6)")' ) &
                   TRIM(label)
   DO itemp = 1, ntemp-1
      WRITE(iu_therm, '(e12.5,4e20.9)') temp(itemp), celldmf_t(1,itemp), &
                                            alpha_t(1,itemp)*1.D6
   END DO
ELSEIF (ibrav==4 .OR. ibrav==6 .OR. ibrav==7 ) THEN
   WRITE(iu_therm,'("#  ",a,"      celldm(1)          celldm(3)          alpha_xx(x10^6)         alpha_zz (x10^6)")' ) TRIM(label)
   DO itemp = 1, ntemp-1
      WRITE(iu_therm, '(e12.5,4e20.9)') temp(itemp), celldmf_t(1,itemp), &
                                                     celldmf_t(3,itemp), &
                                                     alpha_t(1,itemp)*1.D6, &
                                                     alpha_t(3,itemp)*1.D6
   END DO
ELSEIF ( ibrav==5 ) THEN
   WRITE(iu_therm,'("#  ",a,"  celldm(1)   celldm(4)    alpha_xx(x10^6)   alpha_zz (x10^6)")' ) TRIM(label)
   DO itemp = 1, ntemp-1
      WRITE(iu_therm, '(e12.5,4e20.9)') temp(itemp), celldmf_t(1,itemp), &
                                                     celldmf_t(4,itemp), &
                                                     alpha_t(1,itemp)*1.D6, &
                                                     alpha_t(3,itemp)*1.D6
   END DO
ELSEIF (ibrav==8 .OR. ibrav==9 .OR. ibrav==10 .OR. ibrav==11) THEN
   WRITE(iu_therm,'("# ",a,"     celldm(1)     celldm(2)    celldm(3) &
             &alpha_xx(x10^6)  alpha_yy(x10^6)   alpha_zz(x10^6)")' ) &
                               TRIM(label)
   DO itemp = 1, ntemp-1
      WRITE(iu_therm, '(e12.5,6e20.9)') temp(itemp), celldmf_t(1,itemp), &
                                                     celldmf_t(2,itemp), &
                                                     celldmf_t(3,itemp), &
                                                     alpha_t(1,itemp)*1.D6, &
                                                     alpha_t(2,itemp)*1.D6, &
                                                     alpha_t(3,itemp)*1.D6  
   END DO
ELSEIF (ibrav==12 .OR. ibrav==13) THEN
   WRITE(iu_therm,'("#  ",a,"       celldm(1)         celldm(2)        celldm(3)        celldm(4)")' ) TRIM(label)
   DO itemp = 1, ntemp
      WRITE(iu_therm, '(e12.5,4e17.9)') temp(itemp), celldmf_t(1,itemp), &
                                                     celldmf_t(2,itemp), &
                                                     celldmf_t(3,itemp), &
                                                     celldmf_t(4,itemp)
   END DO
ELSEIF (ibrav==-12 .OR. ibrav==-13) THEN
   WRITE(iu_therm,'("#  ",a,"      celldm(1)         celldm(2)        celldm(3)        celldm(5)")' ) TRIM(label)
   DO itemp = 1, ntemp
      WRITE(iu_therm, '(e12.5,4e17.9)') temp(itemp), celldmf_t(1,itemp), &
                                                     celldmf_t(2,itemp), &
                                                     celldmf_t(3,itemp), &
                                                     celldmf_t(5,itemp)
   END DO
ELSEIF (ibrav==14) THEN
   WRITE(iu_therm,'("#   ",a,"       celldm(1)         celldm(2)        &
               &celldm(3)        celldm(4)        celldm(5)        celldm(6)")' ) TRIM(label)
   DO itemp = 1, ntemp
      WRITE(iu_therm, '(e12.5,6e15.7)') temp(itemp), celldmf_t(1,itemp), &
                                                     celldmf_t(2,itemp), &
                                                     celldmf_t(3,itemp), &
                                                     celldmf_t(4,itemp), &
                                                     celldmf_t(5,itemp), &
                                                     celldmf_t(6,itemp)
   END DO
ELSE IF (ibrav==0) THEN
!
!  In this case we write nothing but do not stop
!
ELSE
   CALL errore('write_alpha_anis','ibrav not programmed',1)
END IF

CLOSE(iu_therm)

RETURN
END SUBROUTINE write_alpha_anis
!
!-----------------------------------------------------------------------
SUBROUTINE compute_alpha_anis(celldm_t, alpha_anis_t, temp, ntemp, ibrav)
!-----------------------------------------------------------------------
!
!  This routine computes the thermal expansion tensor from temperature
!  derivatives of the crystal parameters. The input and output 
!  quantities are function of temperature.
!
USE kinds, ONLY : DP

IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav, ntemp
REAL(DP), INTENT(IN) :: celldm_t(6,ntemp), temp(ntemp)
REAL(DP), INTENT(INOUT) :: alpha_anis_t(6,ntemp)

INTEGER :: itemp
REAL(DP) :: fact1, fact2, deriv1, deriv2

alpha_anis_t=0.0_DP
SELECT CASE (ibrav) 
   CASE(1,2,3) 
      DO itemp = 2, ntemp-1
         alpha_anis_t(1,itemp) = (celldm_t(1,itemp+1)-celldm_t(1,itemp-1)) / &
                         (temp(itemp+1)-temp(itemp-1)) / celldm_t(1,itemp)
         alpha_anis_t(2,itemp) = alpha_anis_t(1,itemp)
         alpha_anis_t(3,itemp) = alpha_anis_t(1,itemp)
      END DO
   CASE(4,6,7)
      DO itemp = 2, ntemp-1
         alpha_anis_t(1,itemp) = (celldm_t(1,itemp+1)-celldm_t(1,itemp-1)) / &
                         (temp(itemp+1)-temp(itemp-1)) / celldm_t(1,itemp)
         alpha_anis_t(2,itemp) = alpha_anis_t(1,itemp)
         alpha_anis_t(3,itemp) = ( celldm_t(3,itemp+1)*celldm_t(1,itemp+1)-   &
                               celldm_t(3,itemp-1)*celldm_t(1,itemp-1) )/ &
                         (temp(itemp+1)-temp(itemp-1))/(celldm_t(3,itemp)*  &
                          celldm_t(1,itemp)) 
      END DO
   CASE(5)
      DO itemp = 2, ntemp-1
         fact1 = 2.0_DP * ( 1.0_DP - celldm_t(4,itemp) )
         fact2 = 1.0_DP + 2.0_DP * celldm_t(4,itemp)
         deriv1 = ( celldm_t(1,itemp+1) - celldm_t(1,itemp-1) )/      &
                                     ( temp(itemp+1) - temp(itemp-1) )
         deriv2 = ( celldm_t(4,itemp+1) - celldm_t(4,itemp-1)) /      &
                                     ( temp(itemp+1) - temp(itemp-1) )
         alpha_anis_t(1,itemp) = deriv1 / celldm_t(1, itemp) - deriv2 / fact1 
         alpha_anis_t(2,itemp) = alpha_anis_t(1,itemp)
         alpha_anis_t(3,itemp) = deriv1 / celldm_t(1, itemp) + deriv2 / fact2
      END DO
   CASE (8,9,10,11)
      DO itemp = 2, ntemp-1
         alpha_anis_t(1,itemp) = (celldm_t(1,itemp+1)-celldm_t(1,itemp-1)) / &
                        (temp(itemp+1)-temp(itemp-1)) / celldm_t(1,itemp)
         alpha_anis_t(2,itemp) = ( celldm_t(2,itemp+1)*celldm_t(1,itemp+1)-   &
                                   celldm_t(2,itemp-1)*celldm_t(1,itemp-1) )/ &
                            (temp(itemp+1)-temp(itemp-1))/(celldm_t(2,itemp)*  &
                             celldm_t(1,itemp)) 
         alpha_anis_t(3,itemp) = ( celldm_t(3,itemp+1)*celldm_t(1,itemp+1)-   &
                                   celldm_t(3,itemp-1)*celldm_t(1,itemp-1) )/ &
                         (temp(itemp+1)-temp(itemp-1))/(celldm_t(3,itemp)*  &
                          celldm_t(1,itemp)) 
      END DO
   CASE DEFAULT
! 
!   In this case do nothing. The thermal expansion tensor is not written later
!
END SELECT

RETURN
END SUBROUTINE compute_alpha_anis
!
!-----------------------------------------------------------------------
SUBROUTINE compute_alpha_anis_p(celldm_ptt, celldm_p1, celldm_m1, &
            alpha_anis_ptt, press, npress, ibrav)
!-----------------------------------------------------------------------
!
!  This routine computes the thermal expansion tensor from temperature
!  derivatives of the crystal parameters. The input and output 
!  quantities are function of pressure. The input celldm are at temperature
!  T celldm_ptt, T+dt celldm_p1, T-dt celldm_m1.
!
USE kinds, ONLY : DP
USE temperature, ONLY : deltat

IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav, npress
REAL(DP), INTENT(IN) :: celldm_ptt(6,npress), celldm_p1(6,npress), &
                        celldm_m1(6,npress), press(npress)
REAL(DP), INTENT(INOUT) :: alpha_anis_ptt(6,npress)

INTEGER :: ipress
REAL(DP) :: fact1, fact2, deriv1, deriv2, dt

alpha_anis_ptt=0.0_DP
dt=2.0_DP * deltat
SELECT CASE (ibrav) 
   CASE(1,2,3) 
      DO ipress = 1, npress
         alpha_anis_ptt(1,ipress)=(celldm_p1(1,ipress)-celldm_m1(1,ipress))/ &
                         dt / celldm_ptt(1,ipress)
         alpha_anis_ptt(2,ipress) = alpha_anis_ptt(1,ipress)
         alpha_anis_ptt(3,ipress) = alpha_anis_ptt(1,ipress)
      END DO
   CASE(4,6,7)
      DO ipress = 1, npress
         alpha_anis_ptt(1,ipress) = (celldm_p1(1,ipress)-celldm_m1(1,ipress))/&
                         dt / celldm_ptt(1,ipress)
         alpha_anis_ptt(2,ipress) = alpha_anis_ptt(1,ipress)
         alpha_anis_ptt(3,ipress) = ( celldm_p1(3,ipress)*celldm_p1(1,ipress)-&
                               celldm_m1(3,ipress)*celldm_m1(1,ipress) )/ &
                               dt/(celldm_ptt(3,ipress)*celldm_ptt(1,ipress)) 
      END DO
   CASE(5)
      DO ipress = 1, npress
         fact1 = 2.0_DP * ( 1.0_DP - celldm_ptt(4,ipress) )
         fact2 = 1.0_DP + 2.0_DP * celldm_ptt(4,ipress)
         deriv1 = ( celldm_p1(1,ipress) - celldm_m1(1,ipress) ) / dt
         deriv2 = ( celldm_p1(4,ipress) - celldm_m1(4,ipress))  / dt 
         alpha_anis_ptt(1,ipress) = deriv1 / celldm_ptt(1, ipress) - &
                                                        deriv2 / fact1 
         alpha_anis_ptt(2,ipress) = alpha_anis_ptt(1,ipress)
         alpha_anis_ptt(3,ipress) = deriv1 / celldm_ptt(1, ipress) + &
                                    deriv2 / fact2
      END DO
   CASE (8,9,10,11)
      DO ipress = 1, npress
         alpha_anis_ptt(1,ipress) = (celldm_p1(1,ipress)-celldm_m1(1,ipress))/&
                        dt / celldm_ptt(1,ipress)
         alpha_anis_ptt(2,ipress) = ( celldm_p1(2,ipress)*celldm_p1(1,ipress)-&
                               celldm_m1(2,ipress)*celldm_m1(1,ipress) )/ &
                            dt/(celldm_ptt(2,ipress)*celldm_ptt(1,ipress)) 
         alpha_anis_ptt(3,ipress) = ( celldm_p1(3,ipress)*celldm_p1(1,ipress)-&
                                   celldm_m1(3,ipress)*celldm_m1(1,ipress) )/ &
                         dt/(celldm_ptt(3,ipress)*celldm_ptt(1,ipress)) 
      END DO
   CASE DEFAULT
! 
!   In this case do nothing. The thermal expansion tensor is not written later
!
END SELECT

RETURN
END SUBROUTINE compute_alpha_anis_p

!-----------------------------------------------------------------------
SUBROUTINE convert_ac_alpha(alpha_aux, alpha, cm, ibrav)
!-----------------------------------------------------------------------
!
!  this subroutine receives the thermal expansion calculated with the
!  gruneisen parameters which are derivatives of the frequencies with 
!  respect to the crystal parameters and transforms it into a thermal
!  expansion tensor.
!
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER :: ibrav
REAL(DP) :: alpha(6), alpha_aux(6), cm(6)

alpha(:)=0.0_DP
SELECT CASE(ibrav)
   CASE(1,2,3)
!
!  cubic 
!
      alpha(1)=alpha_aux(1) * cm(1) / 3.0_DP              
      alpha(2)=alpha_aux(1) * cm(1) / 3.0_DP            
      alpha(3)=alpha_aux(1) * cm(1) / 3.0_DP            
   CASE(4,6,7)
!
!  hexagonal or tetragonal
!
      alpha(1)=(alpha_aux(1) * cm(1) - alpha_aux(2) * cm(3) ) / 2.0_DP   
      alpha(2)=alpha(1)             
      alpha(3)=alpha_aux(2) * cm(3)             
   CASE(5)
!
!  trigonal 
!
      alpha(3)=( 1.0_DP + 2.0_DP * cm(4) ) * ( alpha_aux(1) * cm(1) + &
                 2.0_DP * ( 1.0_DP - cm(4) ) * alpha_aux(2) ) / 3.0_DP

      alpha(1)=(alpha_aux(1)*cm(1) - alpha(3) ) / 2.0_DP   
      alpha(2)=alpha(1)             
   CASE(8,9,10,11)
!
!   orthorhombic case
!
      alpha(1)=alpha_aux(1)*cm(1) - alpha_aux(2)*cm(2) - alpha_aux(3)*cm(3)
      alpha(2)=alpha_aux(2)*cm(2)             
      alpha(3)=alpha_aux(3)*cm(3)             
CASE DEFAULT
   CALL errore('convert_ac_alpha','ibrav not programmed',1)
END SELECT

RETURN
END SUBROUTINE convert_ac_alpha
!
!-----------------------------------------------------------------------
SUBROUTINE write_heat_anhar_anis(temp, cet, cvt, cpt, ntemp, filename)
!-----------------------------------------------------------------------
!
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), cet(ntemp), cvt(ntemp), cpt(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# All heat capacities in (Ry/cell/K)")')
   WRITE(iu_therm,'("# T (K)      (C_s-C_e)(T)              C_P(T) &
                              &                C_V-C_e(T) ")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,3e22.13)') temp(itemp), &
                 cpt(itemp)-cet(itemp), cpt(itemp), cvt(itemp)-cet(itemp)
   ENDDO
   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_heat_anhar_anis
!
!-----------------------------------------------------------------------
SUBROUTINE write_heat_anhar_anis_p(press, cet, cvt, cpt, npress, filename)
!-----------------------------------------------------------------------
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER,  INTENT(IN) :: npress
REAL(DP), INTENT(IN) :: press(npress), cet(npress), cvt(npress), cpt(npress)
CHARACTER(LEN=*) :: filename

INTEGER :: ipress, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# All heat capacities in (Ry/cell/K)")')
   WRITE(iu_therm,'("# T (K)      (C_s-C_e)(T)              C_P(T) &
                              &                C_V-C_e(T) ")')

   DO ipress = 1, npress
      WRITE(iu_therm, '(e12.5,3e22.13)') press(ipress), &
                 cpt(ipress)-cet(ipress), cpt(ipress), cvt(ipress)-cet(ipress)
   ENDDO
   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_heat_anhar_anis_p

!-----------------------------------------------------------------------
SUBROUTINE write_heat_anhar_small(temp, cet, ntemp, filename)
!-----------------------------------------------------------------------
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), cet(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# ")')
   WRITE(iu_therm,'("# T (K)        C_e(T) (Ry/cell/K) ")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,e22.13)') temp(itemp), cet(itemp)
   ENDDO
   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_heat_anhar_small
!
!-----------------------------------------------------------------------
SUBROUTINE write_heat_anhar_small_p(press, ce_p, npress, filename)
!-----------------------------------------------------------------------
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER,  INTENT(IN) :: npress
REAL(DP), INTENT(IN) :: press(npress), ce_p(npress)
CHARACTER(LEN=*) :: filename

INTEGER :: ipress, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# ")')
   WRITE(iu_therm,'("# p (kbar)     C_e(T) (Ry/cell/K) ")')

   DO ipress = 1, npress
      WRITE(iu_therm, '(e12.5,e22.13)') press(ipress), ce_p(ipress)
   ENDDO
   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_heat_anhar_small_p

!-----------------------------------------------------------------------
SUBROUTINE write_thermal_stress(temp, bths_t, ntemp, filename)
!-----------------------------------------------------------------------
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), bths_t(3,3,ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# Thermal stresses in (kbar)      ")')
   WRITE(iu_therm,'("# T (K)",5x," b_11", 9x," b_12", 9x," b_13", 9x,&
                      &" b_22", 9x, " b_23", 9x, " b_33")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,6e14.6)') temp(itemp), bths_t(1,1,itemp),   &
                  bths_t(1,2,itemp), bths_t(1,3,itemp), bths_t(2,2,itemp), &
                  bths_t(2,3,itemp), bths_t(3,3,itemp) 
   ENDDO
   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_thermal_stress
!
!-----------------------------------------------------------------------
SUBROUTINE write_thermal_stress_p(press, bths_p, npress, filename)
!-----------------------------------------------------------------------
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER,  INTENT(IN) :: npress
REAL(DP), INTENT(IN) :: press(npress), bths_p(3,3,npress)
CHARACTER(LEN=*) :: filename

INTEGER :: ipress, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# Thermal stresses in (kbar)      ")')
   WRITE(iu_therm,'("# p (kbar)",5x," b_11", 9x," b_12", 9x," b_13", 9x,&
                      &" b_22", 9x, " b_23", 9x, " b_33")')

   DO ipress = 1, npress
      WRITE(iu_therm, '(e12.5,6e14.6)') press(ipress), bths_p(1,1,ipress),   &
                  bths_p(1,2,ipress), bths_p(1,3,ipress), bths_p(2,2,ipress),&
                  bths_p(2,3,ipress), bths_p(3,3,ipress) 
   ENDDO
   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_thermal_stress_p

!-----------------------------------------------------------------------
SUBROUTINE write_generalized_gamma(temp, ggamma_t, ntemp, filename)
!-----------------------------------------------------------------------
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode

IMPLICIT NONE
INTEGER,  INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), ggamma_t(3,3,ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# Generalized Gruneisen parameter      ")')
   WRITE(iu_therm,'("# T (K)",5x," g_11", 9x," g_12", 9x," g_13", 9x,&
                       &" g_22", 9x, " g_23", 9x, " g_33")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,6e14.6)') temp(itemp),   ggamma_t(1,1,itemp), &
             ggamma_t(1,2,itemp), ggamma_t(1,3,itemp), ggamma_t(2,2,itemp), &
             ggamma_t(2,3,itemp), ggamma_t(3,3,itemp) 
   ENDDO
   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_generalized_gamma
!
!-----------------------------------------------------------------------
SUBROUTINE write_generalized_gamma_p(press, ggamma_p, npress, filename)
!-----------------------------------------------------------------------
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode

IMPLICIT NONE
INTEGER,  INTENT(IN) :: npress
REAL(DP), INTENT(IN) :: press(npress), ggamma_p(3,3,npress)
CHARACTER(LEN=*) :: filename

INTEGER :: ipress, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# Generalized Gruneisen parameter      ")')
   WRITE(iu_therm,'("# p(kbar)",5x," g_11", 9x," g_12", 9x," g_13", 9x,&
                       &" g_22", 9x, " g_23", 9x, " g_33")')

   DO ipress = 1, npress
      WRITE(iu_therm, '(e12.5,6e14.6)') press(ipress), ggamma_p(1,1,ipress),  &
             ggamma_p(1,2,ipress), ggamma_p(1,3,ipress), ggamma_p(2,2,ipress),&
             ggamma_p(2,3,ipress), ggamma_p(3,3,ipress) 
   ENDDO
   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_generalized_gamma_p
!
!-----------------------------------------------------------------------
SUBROUTINE set_elastic_grun()
!-----------------------------------------------------------------------
!
USE control_elastic_constants, ONLY : lelastic, lelasticf
USE grun_anharmonic,    ONLY : el_cons_grun_t, el_comp_grun_t, lelastic_grun
USE temperature,        ONLY : ntemp
USE anharmonic,         ONLY : el_cons_t, el_comp_t
USE ph_freq_anharmonic, ONLY : el_consf_t, el_compf_t

IMPLICIT NONE

INTEGER :: itemp

lelastic_grun=.FALSE.
IF (lelasticf) THEN
   DO itemp=1, ntemp
      el_cons_grun_t(:,:,itemp)=el_consf_t(:,:,itemp)
      el_comp_grun_t(:,:,itemp)=el_compf_t(:,:,itemp)
   ENDDO
   lelastic_grun=.TRUE.
ELSEIF(lelastic) THEN
   DO itemp=1, ntemp
      el_cons_grun_t(:,:,itemp)=el_cons_t(:,:,itemp)
      el_comp_grun_t(:,:,itemp)=el_comp_t(:,:,itemp)
   ENDDO
   lelastic_grun=.TRUE.
ENDIF

RETURN
END SUBROUTINE set_elastic_grun
!
!-----------------------------------------------------------------------
SUBROUTINE fit_free_energy_anis_t()
  !-----------------------------------------------------------------------
  !
  !   This routine fits the vibrational free energy (plus possibly the
  !   electronic one) with a polynomial of degree poly_degree_ph
  !   at all temperatures.
  !
  !   The output of this routine are the polynomial coefficients
  !   p2t_t, and one among p1t_t, p3t_t, or p4t_t depending on which is 
  !   poly_degree_ph
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE thermo_mod,  ONLY : celldm_geo, no_ph
  USE control_quartic_energy, ONLY :  lquartic, poly_degree_ph, lsolve
  USE temperature, ONLY : ntemp
  USE thermodynamics,    ONLY : ph_free_ener
  USE el_thermodynamics, ONLY : el_free_ener
  USE control_eldos,     ONLY : lel_free_energy
  USE anharmonic,  ONLY : p1t_t, p2t_t, p3t_t, p4t_t
  USE lattices,    ONLY : compress_celldm, crystal_parameters
  USE io_global,   ONLY : stdout

  USE linear_surfaces,    ONLY : fit_multi_linear
  USE quadratic_surfaces, ONLY : fit_multi_quadratic
  USE cubic_surfaces,     ONLY : fit_multi_cubic
  USE quartic_surfaces,   ONLY : fit_multi_quartic 

  USE polynomial,  ONLY : init_poly

  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: x(:,:), f(:)
  INTEGER  :: itemp, idata, ncoeff, nvar, ndata, ndatatot
  INTEGER  :: compute_nwork, compute_nwork_ph
  !
  nvar=crystal_parameters(ibrav)
  !
  ndatatot = compute_nwork()
  ndata = compute_nwork_ph(no_ph,ndatatot)

  ALLOCATE(x(nvar, ndata))
  ALLOCATE(f(ndata))
  DO itemp=1, ntemp
     CALL init_poly(nvar,p2t_t(itemp))
     ncoeff=1+nvar+p2t_t(itemp)%ncoeff2
     ndata=0
     DO idata=1,ndatatot
        IF (no_ph(idata)) CYCLE
        ndata=ndata+1
        CALL compress_celldm(celldm_geo(1,idata), x(1,ndata), nvar, ibrav)
        f(ndata)=ph_free_ener(itemp,idata)
        if (lel_free_energy) f(ndata)=f(ndata)+el_free_ener(itemp,idata)
     END DO
  !
  !    CALL summarize_fitting_data(nvar, ndata, x, f)
  !
     CALL fit_multi_quadratic(ndata, nvar, lsolve, x, f, p2t_t(itemp))

!     CALL print_quadratic_polynomial(nvar, p2t_t(itemp))

!  WRITE(stdout,'(/,7x,"Energy (1)      Fitted energy (2)   DeltaE (1)-(2)")') 
!     CALL print_chisq_quadratic(ndata, nvar, x, f, p2t_t(itemp))

     IF (lquartic) THEN
        IF (poly_degree_ph==4) THEN
           CALL init_poly(nvar,p4t_t(itemp))
           CALL fit_multi_quartic(ndata, nvar, lsolve, x, f, p4t_t(itemp))
        ELSEIF (poly_degree_ph==3) THEN
           CALL init_poly(nvar,p3t_t(itemp))
           CALL fit_multi_cubic(ndata, nvar, lsolve, x, f, p3t_t(itemp))
        ELSEIF (poly_degree_ph==1) THEN
           CALL init_poly(nvar,p1t_t(itemp))
           CALL fit_multi_linear(ndata, nvar, lsolve, x, f, p1t_t(itemp))
           WRITE(stdout,'(/,5x,"Extremum of the quartic+linear found at:")')
        ENDIF
     ENDIF
  ENDDO
  DEALLOCATE(f)
  DEALLOCATE(x)
  !
  RETURN
  !
END SUBROUTINE fit_free_energy_anis_t

