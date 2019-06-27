!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE manage_anhar()

USE kinds,                 ONLY : DP
USE temperature,           ONLY : ntemp
USE control_thermo,        ONLY : ltherm_dos, ltherm_freq
USE internal_files_names,  ONLY : flfrq_thermo, flvec_thermo
USE anharmonic,            ONLY : vmin_t, b0_t, b01_t, free_e_min_t
USE ph_freq_anharmonic,    ONLY : vminf_t, b0f_t, b01f_t, free_e_minf_t
USE mp_images,             ONLY : inter_image_comm
USE mp,                    ONLY : mp_sum

IMPLICIT NONE

INTEGER :: itemp
CHARACTER(LEN=256) :: filedata, filerap, fileout, gnu_filename, filenameps
LOGICAL :: all_geometry_done

CALL check_all_geometry_done(all_geometry_done)
IF (.NOT.all_geometry_done) RETURN

IF (ltherm_dos) THEN
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
ENDIF

IF (ltherm_freq) THEN
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
ENDIF
!
!    calculate several anharmonic quantities 
!
IF (ltherm_dos) CALL write_anharmonic()
IF (ltherm_freq) CALL write_ph_freq_anharmonic()
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
!
SUBROUTINE manage_anhar_anis()

USE kinds,                 ONLY : DP
USE thermo_mod,            ONLY : reduced_grid, what
USE temperature,           ONLY : ntemp, temp
USE control_pressure,      ONLY : pressure_kb
USE control_thermo,        ONLY : ltherm_dos, ltherm_freq
USE thermodynamics,        ONLY : ph_free_ener
USE ph_freq_thermodynamics, ONLY : phf_free_ener
USE internal_files_names,  ONLY : flfrq_thermo, flvec_thermo
USE anharmonic,            ONLY : celldm_t, free_e_min_t
USE ph_freq_anharmonic,    ONLY : celldmf_t, free_e_minf_t
USE io_global,             ONLY : ionode, stdout
USE mp,                    ONLY : mp_sum
USE mp_world,              ONLY : world_comm

IMPLICIT NONE
INTEGER :: itemp, startt, lastt, idata, ndata
CHARACTER(LEN=256) :: filedata, filerap, fileout, gnu_filename, filenameps
REAL(DP), ALLOCATABLE :: phf(:)
LOGICAL :: all_geometry_done
INTEGER :: compute_nwork

CALL check_all_geometry_done(all_geometry_done)
IF (.NOT.all_geometry_done) RETURN
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
      ENDDO
      CALL quadratic_fit_t(itemp, celldmf_t(:,itemp), free_e_minf_t(itemp), &
                                                      phf, ndata)
   ENDDO
   CALL mp_sum(celldmf_t, world_comm)
   CALL mp_sum(free_e_minf_t, world_comm)
ENDIF
DEALLOCATE(phf)
!
!  Check if the elastic constants are on file. If they are, the code
!  computes the elastic constants as a function of temperature interpolating
!  at the crystal parameters found in the quadratic/quartic fit
!
CALL check_el_cons()
CALL write_elastic_t()

IF (ltherm_dos) CALL write_anhar_anis()
IF (ltherm_freq) CALL write_ph_freq_anhar_anis()
!
!  Plot elastic constants and compliances
!
IF (what=='mur_lc_t') THEN
   CALL plot_elastic_t(0,.TRUE.)
   CALL plot_elastic_t(1,.TRUE.)
ENDIF
!
!    calculate and plot the Gruneisen parameters along the given path.
!
CALL write_gruneisen_band_anis(flfrq_thermo,flvec_thermo)
CALL set_files_for_plot(4, flfrq_thermo, filedata, filerap, &
                                          fileout, gnu_filename, filenameps)
CALL plotband_sub(4, filedata, filerap, fileout, gnu_filename, filenameps)
CALL plot_gruneisen_band_anis(flfrq_thermo)
!
!    fit the frequencies of the dos mesh with a polynomial
!
CALL set_volume_b0_grun()
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
