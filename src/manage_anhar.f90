!
! Copyright (C) 2016-2023 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE manage_anhar()
!---------------------------------------------------------------------

USE kinds,                 ONLY : DP
USE temperature,           ONLY : ntemp
USE thermo_mod,            ONLY : tot_ngeo
USE control_thermo,        ONLY : ltherm_dos, ltherm_freq, ltherm_glob
USE control_elastic_constants, ONLY : el_cons_qha_available, &
                                  el_consf_qha_available
USE control_eldos,         ONLY : lel_free_energy, hot_electrons
USE control_emp_free_ener, ONLY : add_empirical, emp_ener, emp_free_ener, &
                                  emp_entr, emp_ce
USE temperature,           ONLY : temp, ntemp, ntemp_plot, itemp_plot
USE internal_files_names,  ONLY : flfrq_thermo, flvec_thermo
USE el_thermodynamics,     ONLY : el_ener, el_free_ener, el_entr, el_ce
USE data_files,            ONLY : flanhar, fleltherm

USE io_global,             ONLY : stdout

IMPLICIT NONE

INTEGER :: itemp, igeom
CHARACTER(LEN=256) :: filedata, filerap, fileout, gnu_filename, filenameps, &
                      filepbs
LOGICAL :: all_geometry_done, all_el_free, ldummy

CALL check_all_geometries_done(all_geometry_done)
IF (.NOT.all_geometry_done) RETURN

IF (hot_electrons) THEN
   CALL manage_hot_electrons()
ELSEIF (lel_free_energy) THEN
   CALL check_all_el_free_ener_done(all_el_free)
   IF (.NOT.all_el_free) CALL errore('manage_anhar',&
                        'missing electron thermodynamics',1)
   DO igeom=1, tot_ngeo
      CALL set_el_files_names(igeom)
      filedata="therm_files/"//TRIM(fleltherm)
      CALL read_thermo(ntemp, temp, el_ener(:,igeom),           &
                       el_free_ener(:,igeom), el_entr(:,igeom), &
                       el_ce(:,igeom), ldummy, filedata)
   ENDDO
   CALL restore_el_file_names()
ELSE
   el_ener=0.0_DP
   el_free_ener=0.0_DP
   el_entr=0.0_DP
   el_ce=0.0_DP
ENDIF

IF (add_empirical) THEN
   CALL empirical_free_energy()
ELSE
   emp_ener=0.0_DP
   emp_free_ener=0.0_DP
   emp_entr=0.0_DP
   emp_ce=0.0_DP
ENDIF


IF (ltherm_dos) THEN
   WRITE(stdout,'(/,2x,76("-"))')
   WRITE(stdout,'(5x,"Computing the anharmonic properties within ")')
   WRITE(stdout,'(5x,"the QHA approximation using phonon dos.")') 
   WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(flanhar)
   WRITE(stdout,'(2x,76("-"),/)')
!
!  Fit the vibrational (and possibly electronic free energy) with a polynomial
!  fit_free_energy_noe fits only the vibrational part. The fit is made
!  also when ltherm_glob is true to calculate the thermal pressure below.
!
   CALL fit_free_energy()
   CALL fit_free_energy_noe()
!
!  first the crystal parameters as a function of temperature
!  at the input pressure. Then write the parameters on output, 
!  and interpolate the harmonic quantities at the temperature 
!  dependent crystal parameters
!
   IF (ltherm_glob) THEN
      CALL anhar_ev_glob_t()
      CALL anhar_ev_glob_noe_t()
   ELSE
      CALL anhar_ev_t()
      CALL anhar_ev_noe_t()
   ENDIF
   CALL compute_density_t()
   CALL compute_density_noe_t()
   CALL interpolate_harmonic()
   CALL interpolate_harmonic_noe_t()
!
!  then the crystal parameters as a function of temperature 
!  at several pressures
!
   IF (ltherm_glob) THEN
      CALL anhar_ev_glob_pt()
   ELSE
      CALL anhar_ev_pt()
   ENDIF
   CALL compute_density_pt()
   CALL interpolate_harmonic_pt()
!
!  then the crystal parameters as a function of pressure 
!  at several temperatures
!
   IF (ltherm_glob) THEN
      CALL anhar_ev_glob_ptt()
   ELSE
      CALL anhar_ev_ptt()
      CALL anhar_ev_ptt_pm()
   ENDIF
   CALL compute_density_ptt()
   CALL interpolate_harmonic_ptt()
!
!  some quantities as a function of temperature are needed 
!  at constant volume, they are computed here
!
   IF (ltherm_glob) THEN
      CALL anhar_ev_glob_vt()
   ELSE
      CALL anhar_ev_vt()
   ENDIF
!
!  calculate and writes several anharmonic quantities at the input pressure
!  (beta, b0, cp, gamma)
!
   CALL write_anhar()
!
!  writes anharmonic quantities obtained from equation of state
!  (b01, b02)
!
   CALL write_anhar_mur()
!
!  writes other anharmonic quantities (b_fact)
!
   CALL write_anhar_aux()
!
!  writes anharmonic quantities obtained without electronic contribution
!
   CALL write_anhar_noe()
!
!  writes electronic contribution to anharmonic quantities 
!
   CALL write_anhar_el_cont()
!
!   if requested in input writes on files the anharmonic quantities
!   at several pressures. No need here to use global routines, these
!   routines do not depend on the representation of the free energy.
!
    CALL write_anhar_pt()
    CALL write_anhar_mur_pt()
!
!  if requested in input writes on files anharmonic quantities 
!  at several temperatures
!
   IF (ltherm_glob) THEN
      CALL write_anhar_glob_ptt() 
   ELSE
      CALL write_anhar_ptt() 
      CALL write_anhar_mur_ptt() 
   ENDIF
   CALL write_tp_ptt()
!
!   if requested in input writes on files the anharmonic quantities
!   at several volumes
!
   CALL write_anhar_v()
!
!   write on output the electronic contributions if computed
!
   CALL write_anhar_el()
   CALL write_anhar_el_pt()
!
!  for diagnostic purposes write on file only the vibrational free energy
!  and the electronic one if available
!
   CALL write_free_energy()
!
!  write on output the Hugoniot
!   
   CALL write_hugoniot()
!
ENDIF

IF (ltherm_freq) THEN
   WRITE(stdout,'(/,2x,76("+"))')
   WRITE(stdout,'(5x,"Computing the anharmonic properties within ")')
   WRITE(stdout,'(5x,"the QHA approximation using phonon frequencies.")') 
   WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(flanhar)//'_ph'
   WRITE(stdout,'(2x,76("+"),/)')
!
!  Fit the vibrational (and possibly electronic free energy) with a polynomial
!
   CALL fit_free_energy_ph()
   CALL fit_free_energy_ph_noe()
!
!  the crystal parameters as a function of temperature
!
   IF (ltherm_glob) THEN
      CALL anhar_ev_glob_t_ph()
      CALL anhar_ev_glob_noe_t_ph()
   ELSE
      CALL anhar_ev_t_ph()
      CALL anhar_ev_noe_t_ph()
   ENDIF
   CALL compute_densityf_t()
   CALL compute_densityf_noe_t()
!  CALL summarize_anhar_param_ph()
   CALL interpolate_harmonic_ph()
   CALL interpolate_harmonic_noe_t_ph()

   IF (ltherm_glob) THEN
      CALL anhar_ev_glob_ph_pt()
   ELSE
      CALL anhar_ev_pt_ph()
   ENDIF
   CALL compute_densityf_pt()
   CALL interpolate_harmonicf_pt()
!
!  then the crystal parameters as a function of pressure 
!  at several temperatures
!
   IF (ltherm_glob) THEN
      CALL anhar_ev_glob_ph_ptt()
   ELSE
      CALL anhar_ev_ph_ptt()
      CALL anhar_ev_ph_ptt_pm()
   ENDIF
   CALL compute_densityf_ptt()
   CALL interpolate_harmonicf_ptt()
!
!  some quantities as a function of temperature are needed 
!  at constant volume, they are computed here
!
   IF (ltherm_glob) THEN
      CALL ph_freq_anhar_ev_glob_vt()
   ELSE
      CALL ph_freq_anhar_ev_vt()
   ENDIF
!
!  calculate several anharmonic quantities 
!
   CALL write_ph_freq_anhar()
!
!  writes anharmonic quantities obtained from equation of state
!
   CALL write_ph_freq_anhar_mur()
!
!  writes other anharmonic quantities (b_factf)
!
   CALL write_ph_freq_anhar_aux()
!
!  writes anharmonic quantities obtained without electronic contribution
!
   CALL write_ph_freq_anhar_noe()
!
!  writes electronic contribution to anharmonic quantities 
!
   CALL write_ph_freq_anhar_el_cont()
!
!   if requested in input writes on files the anharmonic quantities
!   at several pressures. No need here to use global routines, these
!   routines do not depend on the representation of the free energy.
!
    CALL write_ph_freq_anhar_pt()
    CALL write_ph_freq_anhar_mur_pt()
!
!  if requested in input writes on files anharmonic quantities 
!  at several temperatures
!
   IF (ltherm_glob) THEN
      CALL write_ph_freq_anhar_glob_ptt()
   ELSE
      CALL write_ph_freq_anhar_ptt()
      CALL write_ph_freq_anhar_mur_ptt()
   ENDIF
   CALL write_ph_freq_tp_ptt()
!
!   if requested in input writes on files the anharmonic quantities
!   at several volumes
!
   CALL write_ph_freq_anhar_v()
!
!   write on output the electronic contributions if computed to energy,
!   free_energy, entropy and cv.
!
   CALL write_ph_freq_anhar_el()
   CALL write_ph_freq_anhar_el_pt()
!
!  for diagnostic purposes write on file only the vibrational free energy
!  and the electronic one if available
!
   CALL write_free_energy_ph()
ENDIF
!
!  Check if the elastic constants are on file. 
!  First look for the quasi-harmonic ones
!
CALL check_el_cons_qha()
!
!  If not found search those at T=0 in the elastic_constants directory
!
IF (.NOT.(el_cons_qha_available.OR.el_consf_qha_available)) &
                                                CALL check_el_cons()
!
!  If the elastic constants are on file and the user allows it, the code 
!  computes the elastic constants as a function of temperature interpolating 
!  at the crystal parameters found in the quadratic/quartic fit
!
CALL set_elastic_constants_t()

IF (ltherm_dos) THEN
   CALL write_anhar_anis()
   CALL write_anhar_anis_pt()
   CALL write_anhar_anis_ptt()
ENDIF

IF (ltherm_freq) THEN
   CALL write_ph_freq_anhar_anis()
   CALL write_ph_freq_anhar_anis_pt()
   CALL write_ph_freq_anhar_anis_ptt()
ENDIF

WRITE(stdout,'(/,2x,76("-"))')
WRITE(stdout,'(5x,"Computing the anharmonic properties within ")')
WRITE(stdout,'(5x,"the QHA approximation using Gruneisen parameters.")') 
WRITE(stdout,'(2x,76("-"),/)')
!
!    calculate and plot the Gruneisen parameters along the given path
!
CALL write_gruneisen_band(flfrq_thermo,flvec_thermo)
CALL set_files_for_plot(3, flfrq_thermo, filedata, filerap, &
                                      fileout, gnu_filename, filenameps, filepbs)
CALL plotband_sub(3, filedata, filerap, fileout, gnu_filename, filenameps, filepbs)
CALL set_files_for_plot(4, flfrq_thermo, filedata, filerap,  &
                                       fileout, gnu_filename, filenameps, filepbs)
CALL plotband_sub(4, filedata, filerap, fileout, gnu_filename, filenameps, filepbs)
!
!    fit the frequencies of the dos mesh with a polynomial
!
CALL fit_frequencies()
!
!    calculate the Gruneisen parameters and the anharmonic quantities
!
CALL set_volume_b0_grun()
CALL write_grun_anharmonic()
!
!   now plot all the computed quantities
!
CALL manage_plot_anhar()

CALL manage_plot_elastic()
!
!   summarize on output the main anharmonic quantities 
!
IF (ltherm_dos) CALL summarize_anhar()
IF (ltherm_freq.AND.(.NOT.ltherm_dos)) CALL summarize_anhar_ph()
RETURN
END SUBROUTINE manage_anhar
!
!-------------------------------------------------------------------------
SUBROUTINE manage_anhar_anis()
!-------------------------------------------------------------------------
!
USE kinds,                 ONLY : DP
USE thermo_mod,            ONLY : reduced_grid, tot_ngeo
USE temperature,           ONLY : ntemp, temp, ntemp_plot, itemp_plot
USE control_thermo,        ONLY : ltherm_dos, ltherm_freq
USE control_elastic_constants, ONLY : el_cons_qha_available, &
                                  el_consf_qha_available
USE control_eldos,         ONLY : lel_free_energy
USE thermodynamics,        ONLY : ph_free_ener
USE ph_freq_thermodynamics, ONLY : phf_free_ener
USE el_thermodynamics,     ONLY : el_ener, el_free_ener, el_entr, &
                                  el_ce
USE data_files,            ONLY : fleltherm
USE internal_files_names,  ONLY : flfrq_thermo, flvec_thermo
USE io_global,             ONLY : stdout
USE mp,                    ONLY : mp_sum
USE mp_world,              ONLY : world_comm

IMPLICIT NONE
INTEGER :: itemp, itempp, igeom, startt, lastt, idata, ndata
CHARACTER(LEN=256) :: filedata, filerap, fileout, gnu_filename, filenameps, &
                      filepbs
REAL(DP), ALLOCATABLE :: phf(:)
LOGICAL :: all_geometry_done, all_el_free, ldummy
INTEGER :: compute_nwork

CALL check_all_geometries_done(all_geometry_done)
IF (.NOT.all_geometry_done) RETURN

IF (lel_free_energy) THEN
   CALL check_all_el_free_ener_done(all_el_free)
   IF (.NOT.all_el_free) CALL errore('manage_anhar_anis',&
                        'missing electron thermodynamics',1)
   DO igeom=1, tot_ngeo
      CALL set_el_files_names(igeom)
      filedata="therm_files/"//TRIM(fleltherm)
      CALL read_thermo(ntemp, temp, el_ener(:,igeom),           &
                       el_free_ener(:,igeom), el_entr(:,igeom), &
                       el_ce(:,igeom), ldummy, filedata)
   ENDDO
   CALL restore_el_file_names()
ELSE
   el_ener=0.0_DP
   el_free_ener=0.0_DP
   el_entr=0.0_DP
   el_ce=0.0_DP
ENDIF
!
!    Anisotropic solid. Compute the crystal parameters, the thermal expansion 
!    tensor and the Helmholtz (or Gibbs) free energy as a function 
!    of temperature.
!
ndata= compute_nwork()

ALLOCATE(phf(ndata))
CALL divide(world_comm, ntemp, startt, lastt)
IF (ltherm_dos) THEN
!
!  fit the free energy with a polynomial
!
   CALL fit_free_energy_anis_t()
   CALL fit_free_energy_noe_anis_t()
!
!  Use the polynomial to find the celldm_t and the energy at the minimum
!  of the free energy (Gibbs energy if pressure is given in input)
!  and the corresponding volumes
!
   CALL quadratic_fit_t_run()
   CALL quadratic_fit_t_noe_run()

   CALL compute_volume_t()
   CALL compute_density_t()
   CALL compute_volume_noe_t()
   CALL compute_density_noe_t()
!
!  Find celldm_t at pressure p+dp and p-dp
!
   CALL quadratic_fit_t_pm()
   CALL quadratic_fit_noe_t_pm()
!
!  Here compute the celldm_pt that minimize the Gibbs energy for 
!  selected pressures and the corresponding volumes
!
   CALL quadratic_fit_pt()
   CALL compute_volume_pt()
   CALL compute_density_pt()
!
!  Here compute the celldm_pt that minimize the Gibbs energy for 
!  selected pressures +- delta p for computing the bulk modulus
!  at selected pressures 
!
   CALL quadratic_fit_pt_pm()
!
!  Here compute the celldm_ptt that minimize the Gibbs energy for 
!  selected temperatures and the corresponding volumes
!
   CALL quadratic_fit_ptt()
   CALL compute_volume_ptt()
   CALL compute_density_ptt()
   CALL write_tp_anis_ptt()

!  Here compute the celldm_ptt_pm that minimize the Gibbs energy for 
!  selected temperatures +- deltaT and the corresponding volumes
!  needed to compute thermal expansion
!
   CALL compute_volume_ptt_pm()
!
!  Compute the bulk modulus as a funtion of T at the input pressure,
!  for selected pressures and for selected temperatures 
!
   CALL compute_bulk_modulus_t()
   CALL compute_bulk_modulus_noe_t()
   CALL compute_bulk_modulus_pt()
   CALL compute_bulk_modulus_ptt()
!
!  Intepolate the thermal energy, free_energy, entropy and c_e at the
!  minimum of the free energy, at the minimum of the gibbs energy for 
!  selected pressures and for selected temperatures
!
   CALL interpolate_harmonic()
   CALL interpolate_harmonic_noe_t()
   CALL interpolate_harmonic_pt()
   CALL interpolate_harmonic_ptt()
!
!  This is obsolete but presently the celldm as a function of pressure
!  are recalculated and written on file here. Only the writing on file
!  should remain.
!
   DO itempp=1,ntemp_plot
      itemp=itemp_plot(itempp)
      DO idata=1,ndata
         phf(idata)=ph_free_ener(itemp,idata)
         IF (lel_free_energy) phf(idata)=phf(idata)+ &
                                    el_free_ener(itemp, idata)
      ENDDO
      CALL write_e_omega_t(itemp, phf, ndata, '.mur_temp', '_mur_celldm.' )
   ENDDO
!
!  for diagnostic purposes write on file only the vibrational free energy
!  and the electronic one if available (only for cubic systems)
!
   CALL write_free_energy()
!
!  Here we plot the paths at constant temperature for which the
!  stress tensor is a uniform pressure
!
   CALL plot_multi_energy_t()

ENDIF

IF (ltherm_freq) THEN
!
!  fit the free energy with a polynomial
!
   CALL fit_free_energyf_anis_t()
   CALL fit_free_energyf_noe_anis_t()
!
!  Use the polynomial to find the celldm_t and the energy at the minimum
!  of the free energy (Gibbs energy if pressure is given in input)
!  and the corresponding volumes
!
   CALL quadratic_fitf_t_run()
   CALL quadratic_fitf_t_noe_run()

   CALL compute_volumef_t()
   CALL compute_densityf_t()
   CALL compute_volumef_noe_t()
   CALL compute_densityf_noe_t()
!
!  Find celldmf_t at pressure p+dp and p-dp
!
   CALL quadratic_fitf_t_pm()
   CALL quadratic_fitf_noe_t_pm()
!
!  Here compute the celldm_pt that minimize the Gibbs energy for 
!  selected pressures and the corresponding volumes
!
   CALL quadratic_fitf_pt()
   CALL compute_volumef_pt()
   CALL compute_densityf_pt()
!
!  Here compute the celldm_pt that minimize the Gibbs energy for 
!  selected pressures +- delta p for computing the bulk modulus
!  at selected pressures 
!
   CALL quadratic_fitf_pt_pm()
!
!  Here compute the celldm_ptt that minimize the Gibbs energy for 
!  selected temperatures and the corresponding volumes
!
   CALL quadratic_fitf_ptt()
   CALL compute_volumef_ptt()
   CALL compute_densityf_ptt()
   CALL write_ph_freq_tp_anis_ptt()
!
!  Here compute the celldm_ptt_pm that minimize the Gibbs energy for 
!  selected temperatures +- deltaT and the corresponding volumes
!  needed to compute thermal expansion
!
   CALL compute_volumef_ptt_pm()
!
!  Compute the bulk modulus as a funtion of T at the input pressure
!
   CALL compute_bulk_modulusf_t()
   CALL compute_bulk_modulusf_noe_t()
   CALL compute_bulk_modulusf_pt()
   CALL compute_bulk_modulusf_ptt()
!
!  Intepolate the thermal energy, free_energy, entropy and c_e at the
!  minimum of the free energy, at the minimum of the gibbs energy for 
!  selected pressures and for selected temperatures
!
   CALL interpolate_harmonic_ph()
   CALL interpolate_harmonicf_noe_t()
   CALL interpolate_harmonicf_pt()
   CALL interpolate_harmonicf_ptt()
!
!  This is obsolete but presently the celldm as a function of pressure
!  are recalculated and written on file here. Only the writing on file
!  should remain.
!
   DO itempp=1,ntemp_plot
      itemp=itemp_plot(itempp)
      DO idata=1,ndata
         phf(idata)=phf_free_ener(itemp,idata)
         IF (lel_free_energy) phf(idata)=phf(idata)+ &
                                    el_free_ener(itemp, idata)
      ENDDO
      CALL write_e_omega_t(itemp, phf, ndata, '.mur_ph_temp', '_mur_ph_celldm.' )
   ENDDO
!
!  for diagnostic purposes write on file only the vibrational free energy
!  and the electronic one if available (only for cubic systems)
!
   CALL write_free_energy_ph()
ENDIF
DEALLOCATE(phf)
!
!  Check if the elastic constants are on file. 
!  First look for the quasi-harmonic ones
!
CALL check_el_cons_qha()
!
!  If not found search those at T=0 in the elastic_constants directory
!
IF (.NOT.(el_cons_qha_available.OR.el_consf_qha_available)) &
                                                CALL check_el_cons()
!
!  If the elastic constants are on file and the user allows it, the code 
!  computes the elastic constants as a function of temperature interpolating 
!  at the crystal parameters found in the quadratic/quartic fit
!
CALL set_elastic_constants_t()

IF (ltherm_dos) THEN
!
!  calculates and writes several anharmonic quantities at the input pressure
!
   CALL write_anhar()
   CALL write_anhar_aux()
   CALL write_anhar_anis()

   CALL write_anhar_noe()
   CALL write_anhar_el_cont()

   CALL write_anhar_pt()
   CALL write_anhar_anis_pt()

   CALL write_anhar_ptt()
   CALL write_anhar_anis_ptt()
!
!   write on output the electronic contributions if computed, to
!   energy, free_energy, entropy and cv
!
   CALL write_anhar_el()
   CALL write_anhar_el_pt()
ENDIF
IF (ltherm_freq) THEN
   CALL write_ph_freq_anhar()
!
!  writes other anharmonic quantities (b_factf)
!
   CALL write_ph_freq_anhar_aux()
   CALL write_ph_freq_anhar_anis()

   CALL write_ph_freq_anhar_noe()
   CALL write_ph_freq_anhar_el_cont()

   CALL write_ph_freq_anhar_pt()
   CALL write_ph_freq_anhar_anis_pt()

   CALL write_ph_freq_anhar_ptt()
   CALL write_ph_freq_anhar_anis_ptt()
!
!   write on output the electronic contributions if computed
!
   CALL write_ph_freq_anhar_el()
   CALL write_ph_freq_anhar_el_pt()
ENDIF
!
!  Plot elastic constants and compliances when what='mur_lc_t'
!
CALL manage_plot_elastic()
!
!    calculate and plot the Gruneisen parameters along the given path.
!
WRITE(stdout,'(/,2x,76("-"))')
WRITE(stdout,'(5x,"Computing the anharmonic properties within ")')
WRITE(stdout,'(5x,"the QHA approximation using Gruneisen parameters.")')
WRITE(stdout,'(2x,76("-"),/)')

CALL write_gruneisen_band_anis(flfrq_thermo,flvec_thermo)
CALL set_files_for_plot(4, flfrq_thermo, filedata, filerap, &
                                          fileout, gnu_filename, filenameps, filepbs)
CALL plotband_sub(4, filedata, filerap, fileout, gnu_filename, filenameps, filepbs)
CALL plot_gruneisen_band_anis(flfrq_thermo)
!
!    fit the frequencies of the dos mesh with a polynomial
!
CALL set_volume_b0_grun()
CALL set_elastic_grun()
IF (reduced_grid) THEN
   CALL fit_frequencies_anis_reduced()
ELSE
   CALL fit_frequencies_anis()
ENDIF
!
!    calculate the Gruneisen parameters and the anharmonic quantities
!
CALL write_grun_anhar_anis()
!
!    and plot them
!
CALL manage_plot_anhar_anis()
!
!   summarize on output the main anharmonic quantities 
!
IF (ltherm_dos) CALL summarize_anhar()
IF (ltherm_freq.AND.(.NOT.ltherm_dos)) CALL summarize_anhar_ph()
!
RETURN
END SUBROUTINE manage_anhar_anis
!
!-------------------------------------------------------------------
SUBROUTINE manage_el_anhar()
!-------------------------------------------------------------------
!
USE kinds,                 ONLY : DP
USE temperature,           ONLY : ntemp
USE data_files,            ONLY : flelanhar
USE el_anharmonic,         ONLY : vmine_t, b0e_t, b01e_t, b02e_t, &
                                  free_e_mine_t
USE io_global,             ONLY : stdout
USE mp_images,             ONLY : inter_image_comm
USE mp,                    ONLY : mp_sum

IMPLICIT NONE

INTEGER :: itemp

WRITE(stdout,'(/,2x,76("-"))')
WRITE(stdout,'(5x,"Computing the crystal parameters adding ")')
WRITE(stdout,'(5x,"the electronic energy.")') 
WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(flelanhar)
WRITE(stdout,'(2x,76("-"),/)')
!
!  first the crystal parameters as a function of temperature
!
vmine_t=0.0_DP
b0e_t=0.0_DP
b01e_t=0.0_DP
b02e_t=0.0_DP
free_e_mine_t=0.0_DP
DO itemp=1, ntemp
   CALL do_ev_t_el(itemp)
ENDDO
CALL mp_sum(vmine_t, inter_image_comm)
CALL mp_sum(b0e_t, inter_image_comm)
CALL mp_sum(b01e_t, inter_image_comm)
CALL mp_sum(b02e_t, inter_image_comm)
CALL mp_sum(free_e_mine_t, inter_image_comm)
!
!    fit some electronic harmonic quantities 
!
CALL write_el_fit_harmonic()

RETURN
END SUBROUTINE manage_el_anhar
