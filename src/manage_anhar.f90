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
USE control_thermo,        ONLY : ltherm_dos, ltherm_freq
USE control_eldos,         ONLY : lel_free_energy
USE temperature,           ONLY : temp, ntemp
USE internal_files_names,  ONLY : flfrq_thermo, flvec_thermo
USE el_thermodynamics,     ONLY : el_ener, el_free_ener, el_entr, &
                                  el_ce
USE data_files,            ONLY : flanhar, fleltherm
USE anharmonic,            ONLY : vmin_t, b0_t, b01_t, free_e_min_t
USE ph_freq_anharmonic,    ONLY : vminf_t, b0f_t, b01f_t, free_e_minf_t
USE io_global,             ONLY : stdout
USE mp_images,             ONLY : inter_image_comm
USE mp,                    ONLY : mp_sum

IMPLICIT NONE

INTEGER :: itemp, igeom
CHARACTER(LEN=256) :: filedata, filerap, fileout, gnu_filename, filenameps
LOGICAL :: all_geometry_done, all_el_free, ldummy

CALL check_all_geometries_done(all_geometry_done)
IF (.NOT.all_geometry_done) RETURN

IF (lel_free_energy) THEN
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
ENDIF

IF (ltherm_dos) THEN
   WRITE(stdout,'(/,2x,76("-"))')
   WRITE(stdout,'(5x,"Computing the anharmonic properties within ")')
   WRITE(stdout,'(5x,"the QHA approximation using phonon dos.")') 
   WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(flanhar)
   WRITE(stdout,'(2x,76("-"),/)')
   !
   !  first the crystal parameters as a function of temperature
   !
   vmin_t=0.0_DP
   b0_t=0.0_DP
   b01_t=0.0_DP
   free_e_min_t=0.0_DP
   DO itemp=1, ntemp
      CALL do_ev_t(itemp)
   ENDDO
   CALL mp_sum(vmin_t, inter_image_comm)
   CALL mp_sum(b0_t, inter_image_comm)
   CALL mp_sum(b01_t, inter_image_comm)
   CALL mp_sum(free_e_min_t, inter_image_comm)
!
!    calculate several anharmonic quantities 
!
   CALL write_anharmonic()
ENDIF

IF (ltherm_freq) THEN
   WRITE(stdout,'(/,2x,76("+"))')
   WRITE(stdout,'(5x,"Computing the anharmonic properties within ")')
   WRITE(stdout,'(5x,"the QHA approximation using phonon frequencies.")') 
   WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(flanhar)//'_ph'
   WRITE(stdout,'(2x,76("+"),/)')
   !
   !  first the crystal parameters as a function of temperature
   !
   vminf_t=0.0_DP
   b0f_t=0.0_DP
   b01f_t=0.0_DP
   free_e_minf_t=0.0_DP
   DO itemp = 1, ntemp
      CALL do_ev_t_ph(itemp)
   ENDDO
   CALL mp_sum(vminf_t, inter_image_comm)
   CALL mp_sum(b0f_t, inter_image_comm)
   CALL mp_sum(b01f_t, inter_image_comm)
   CALL mp_sum(free_e_minf_t, inter_image_comm)
!
!    calculate several anharmonic quantities 
!
   CALL write_ph_freq_anharmonic()
ENDIF

WRITE(stdout,'(/,2x,76("-"))')
WRITE(stdout,'(5x,"Computing the anharmonic properties within ")')
WRITE(stdout,'(5x,"the QHA approximation using Gruneisen parameters.")') 
WRITE(stdout,'(2x,76("-"),/)')
   !
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
CALL plot_anhar() 

RETURN
END SUBROUTINE manage_anhar

!-------------------------------------------------------------------------
SUBROUTINE manage_anhar_anis()
!-------------------------------------------------------------------------

USE kinds,                 ONLY : DP
USE thermo_mod,            ONLY : reduced_grid, what, tot_ngeo
USE temperature,           ONLY : ntemp, temp
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
USE ph_freq_anharmonic,    ONLY : celldmf_t, free_e_minf_t
USE io_global,             ONLY : ionode, stdout
USE mp,                    ONLY : mp_sum
USE mp_world,              ONLY : world_comm

IMPLICIT NONE
INTEGER :: itemp, igeom, startt, lastt, idata, ndata
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

IF (ltherm_dos) CALL write_anhar_anis()
IF (ltherm_freq) CALL write_ph_freq_anhar_anis()
!
!  Plot elastic constants and compliances
!
IF (what=='mur_lc_t') THEN
   CALL plot_elastic_t(0,.TRUE.)
   CALL plot_elastic_t(1,.TRUE.)

   CALL plot_macro_el_t()
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
CALL plot_anhar_anis()

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
USE el_anharmonic,         ONLY : vmine_t, b0e_t, b01e_t, free_e_mine_t
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
free_e_mine_t=0.0_DP
DO itemp=1, ntemp
   CALL do_ev_t_el(itemp)
ENDDO
CALL mp_sum(vmine_t, inter_image_comm)
CALL mp_sum(b0e_t, inter_image_comm)
CALL mp_sum(b01e_t, inter_image_comm)
CALL mp_sum(free_e_mine_t, inter_image_comm)
!
!    calculate several anharmonic quantities 
!
CALL write_el_anharmonic()

!CALL plot_el_anhar() 

RETURN
END SUBROUTINE manage_el_anhar
