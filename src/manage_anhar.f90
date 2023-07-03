!
! Copyright (C) 2016 Andrea Dal Corso
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
USE control_eldos,         ONLY : lel_free_energy, hot_electrons
USE control_quartic_energy, ONLY : poly_degree_ph
USE control_emp_free_ener, ONLY : add_empirical, emp_ener, emp_free_ener, &
                                  emp_entr, emp_ce
USE temperature,           ONLY : temp, ntemp, ntemp_plot, itemp_plot
USE internal_files_names,  ONLY : flfrq_thermo, flvec_thermo
USE anharmonic,            ONLY : a_t
USE el_thermodynamics,     ONLY : el_ener, el_free_ener, el_entr, el_ce
USE data_files,            ONLY : flanhar, fleltherm

USE control_mur,           ONLY : vmin, b0, b01, b02, emin
USE io_global,             ONLY : stdout
USE mp_images,             ONLY : inter_image_comm
USE mp,                    ONLY : mp_sum

IMPLICIT NONE

INTEGER :: itemp, itempp, igeom, ivol, ivolp, m1
CHARACTER(LEN=256) :: filedata, filerap, fileout, gnu_filename, filenameps
LOGICAL :: all_geometry_done, all_el_free, ldummy

CALL check_all_geometries_done(all_geometry_done)
IF (.NOT.all_geometry_done) RETURN

IF (lel_free_energy.AND..NOT.hot_electrons) THEN
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

IF (hot_electrons) CALL manage_hot_electrons()

IF (ltherm_dos) THEN
   WRITE(stdout,'(/,2x,76("-"))')
   WRITE(stdout,'(5x,"Computing the anharmonic properties within ")')
   WRITE(stdout,'(5x,"the QHA approximation using phonon dos.")') 
   WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(flanhar)
   WRITE(stdout,'(2x,76("-"),/)')
!
!  Fit the vibrational (and possibly electronic free energy) with a polynomial
!  fit_free_energy_noe fits only the vibrational part.
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
   CALL summarize_anhar_param()
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
   CALL interpolate_harmonic_pt()
!
!  then the crystal parameters as a function of pressure 
!  at several temperatures
!
   IF (ltherm_glob) THEN
      CALL anhar_ev_glob_ptt()
   ELSE
      CALL anhar_ev_ptt()
   ENDIF
   CALL interpolate_harmonic_ptt()
!
!  some quantities as a function of temperature are needed 
!  at constant volume, they are computed here
!
   CALL anhar_ev_vt()
!
!  calculate several anharmonic quantities at the input pressure
!
   CALL write_anhar()
!
!   if requested in input writes on files the anharmonic quantities
!   at several pressures
!
   CALL write_anhar_pt()
!
!  if requested in input writes on files anharmonic quantities 
!  at several temperatures
!
   IF (ltherm_glob) THEN
      CALL write_anhar_glob_t() 
   ELSE
      CALL write_anhar_ptt() 
   ENDIF
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
!
!  the crystal parameters as a function of temperature
!
   IF (ltherm_glob) THEN
      CALL anhar_ev_glob_t_ph()
   ELSE
      CALL anhar_ev_t_ph()
   ENDIF
   CALL summarize_anhar_param_ph()
   CALL interpolate_harmonic_ph()
!
!  calculate several anharmonic quantities 
!
   CALL write_ph_freq_anhar()
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
                                      fileout, gnu_filename, filenameps)
CALL plotband_sub(3, filedata, filerap, fileout, gnu_filename, filenameps)
CALL set_files_for_plot(4, flfrq_thermo, filedata, filerap,  &
                                       fileout, gnu_filename, filenameps)
CALL plotband_sub(4, filedata, filerap, fileout, gnu_filename, filenameps)
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
!
!   summarize the main anharmonic quantities of a solid
!
CALL summarize_anhar()
RETURN
END SUBROUTINE manage_anhar
!
!-------------------------------------------------------------------------
SUBROUTINE manage_anhar_anis()
!-------------------------------------------------------------------------

USE kinds,                 ONLY : DP
USE thermo_mod,            ONLY : reduced_grid, what, tot_ngeo
USE temperature,           ONLY : ntemp, temp, ntemp_plot, itemp_plot
USE control_pressure,      ONLY : pressure_kb
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
USE anharmonic,            ONLY : celldm_t, free_e_min_t
USE anharmonic_ptt,        ONLY : celldm_ptt, celldm_ptt_p1, &
                                  celldm_ptt_m1, emin_ptt, emin_ptt_p1, &
                                  emin_ptt_m1
USE ph_freq_anharmonic,    ONLY : celldmf_t, free_e_minf_t
USE io_global,             ONLY : ionode, stdout
USE mp,                    ONLY : mp_sum
USE mp_world,              ONLY : world_comm

USE control_elastic_constants, ONLY : lelastic, lelasticf

IMPLICIT NONE
INTEGER :: itemp, itempp, igeom, startt, lastt, idata, ndata
CHARACTER(LEN=256) :: filedata, filerap, fileout, gnu_filename, filenameps
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
   celldm_t=0.0_DP
   free_e_min_t=0.0_DP
   DO itemp = startt, lastt
      WRITE(stdout,'(/,5x,70("-"))')
      IF (pressure_kb > 0.0_DP) THEN
         WRITE(stdout,'(5x, "Gibbs energy from phdos, at T= ", f12.6)') &
                                                                  temp(itemp)
         WRITE(stdout,'(5x, "Pressure is :",f12.6)') pressure_kb
      ELSE
         WRITE(stdout,'(5x, "Helmholtz free energy from phdos, at &
                                                 &T= ", f12.6)') temp(itemp)
      ENDIF
      DO idata=1,ndata
         phf(idata)=ph_free_ener(itemp,idata)
         IF (lel_free_energy) phf(idata)=phf(idata)+el_free_ener(itemp, idata)
      ENDDO
      CALL quadratic_fit_t(itemp, celldm_t(:,itemp), free_e_min_t(itemp), phf,&
                                                                        ndata)
   ENDDO
   CALL mp_sum(celldm_t, world_comm)
   CALL mp_sum(free_e_min_t, world_comm)

   CALL fit_free_energy_anis_t()

   CALL quadratic_fit_pt()
   CALL interpolate_harmonic_pt()

   DO itempp=1,ntemp_plot
!
!  to obtain the thermal expansion we need the equilibrium celldm at
!  temperature itemp, itemp+1, and itemp-1
!
      itemp=itemp_plot(itempp)
      CALL quadratic_fit_ptt(celldm_ptt(:,:,itempp), emin_ptt(:,itempp), itemp)
      CALL quadratic_fit_ptt(celldm_ptt_p1(:,:,itempp), &
                             emin_ptt_p1(:,itempp), itemp+1)
      CALL quadratic_fit_ptt(celldm_ptt_m1(:,:,itempp), &
                             emin_ptt_m1(:,itempp), itemp-1)
   ENDDO
   CALL interpolate_harmonic_ptt()

   DO itempp=1,ntemp_plot
      DO idata=1,ndata
         phf(idata)=ph_free_ener(itemp_plot(itempp),idata)
         IF (lel_free_energy) phf(idata)=phf(idata)+ &
                                    el_free_ener(itemp_plot(itempp), idata)
      ENDDO
      CALL write_e_omega_t(itemp_plot(itempp), phf, ndata)
   ENDDO
!
!  for diagnostic purposes write on file only the vibrational free energy
!  and the electronic one if available
!
   CALL write_free_energy()

ENDIF

IF (ltherm_freq) THEN
   celldmf_t=0.0_DP
   free_e_minf_t=0.0_DP
   DO itemp = startt, lastt
      WRITE(stdout,'(/,5x,70("+"))')
      IF (pressure_kb > 0.0_DP) THEN
         WRITE(stdout,'(5x, "Gibbs energy from integration, at T= ", f12.6)') &
                                                                   temp(itemp)
         WRITE(stdout,'(5x, "Pressure is :",f12.6)') pressure_kb
      ELSE
         WRITE(stdout,'(5x, "Helmholtz Free energy from integration, at T= ", &
                                                      &f12.6)') temp(itemp)
      ENDIF
      DO idata=1,ndata
         phf(idata)=phf_free_ener(itemp,idata)
         IF (lel_free_energy) phf(idata)=phf(idata)+el_free_ener(itemp, idata)
      ENDDO
      CALL quadratic_fit_t(itemp, celldmf_t(:,itemp), free_e_minf_t(itemp), &
                                                      phf, ndata)
   ENDDO
   CALL mp_sum(celldmf_t, world_comm)
   CALL mp_sum(free_e_minf_t, world_comm)
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
   CALL write_anhar_anis()
   CALL write_anhar_anis_pt()
   CALL write_anhar_anis_ptt()
!
!   write on output the electronic contributions if computed
!
   CALL write_anhar_el()
   CALL write_anhar_el_pt()
ENDIF
IF (ltherm_freq) CALL write_ph_freq_anhar_anis()
!
!  Plot elastic constants and compliances
!
IF (what=='mur_lc_t') THEN
   CALL plot_elastic_t(0,.TRUE.)
   CALL plot_elastic_t(1,.TRUE.)
   CALL plot_elastic_pt(0,.TRUE.)
   CALL plot_elastic_pt(1,.TRUE.)
   CALL plot_elastic_ptt(0,.TRUE.)
   CALL plot_elastic_ptt(1,.TRUE.)

   CALL plot_macro_el_t()
!   CALL plot_macro_el_pt()
!   CALL plot_macro_el_ptt()
ENDIF
!
!    calculate and plot the Gruneisen parameters along the given path.
!
WRITE(stdout,'(/,2x,76("-"))')
WRITE(stdout,'(5x,"Computing the anharmonic properties within ")')
WRITE(stdout,'(5x,"the QHA approximation using Gruneisen parameters.")')
WRITE(stdout,'(2x,76("-"),/)')

CALL write_gruneisen_band_anis(flfrq_thermo,flvec_thermo)
CALL set_files_for_plot(4, flfrq_thermo, filedata, filerap, &
                                          fileout, gnu_filename, filenameps)
CALL plotband_sub(4, filedata, filerap, fileout, gnu_filename, filenameps)
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
!    and plot them
!
CALL write_grun_anhar_anis()

CALL manage_plot_anhar_anis()
!
!   summarize the main anharmonic quantities of a solid
!
CALL summarize_anhar()
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

!CALL plot_el_anhar() 

RETURN
END SUBROUTINE manage_el_anhar
