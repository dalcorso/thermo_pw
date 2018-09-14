!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE manage_anhar()

USE kinds,                 ONLY : DP
USE temperature,           ONLY : ntemp, temp
USE control_thermo,        ONLY : ltherm_dos, ltherm_freq, set_internal_path
USE internal_files_names,  ONLY : flfrq_thermo, flvec_thermo
USE anharmonic,            ONLY : vmin_t, b0_t, b01_t, free_e_min_t
USE ph_freq_anharmonic,    ONLY : vminf_t, b0f_t, b01f_t, free_e_minf_t
USE mp_images,             ONLY : inter_image_comm
USE mp,                    ONLY : mp_sum
USE io_global,             ONLY : stdout

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
!    here we calculate several anharmonic quantities 
!
IF (ltherm_dos) CALL write_anharmonic()
IF (ltherm_freq) CALL write_ph_freq_anharmonic()
!
!    here we calculate and plot the Gruneisen parameters along the given path
!
CALL write_gruneisen_band(flfrq_thermo,flvec_thermo)
CALL set_files_for_plot(3, flfrq_thermo, filedata, filerap, &
                                      fileout, gnu_filename, filenameps)
CALL plotband_sub(3, filedata, filerap, fileout, gnu_filename, filenameps)
CALL set_files_for_plot(4, flfrq_thermo, filedata, filerap,  &
                                       fileout, gnu_filename, filenameps)
CALL plotband_sub(4, filedata, filerap, fileout, gnu_filename, filenameps)
!
!    here we compute the Gruneisen parameters on the dos mesh
!
CALL compute_gruneisen()
!
!    here we calculate several anharmonic quantities and plot them
!
CALL write_grun_anharmonic()
CALL plot_anhar() 

RETURN
END SUBROUTINE manage_anhar

!-------------------------------------------------------------------------
!
SUBROUTINE manage_anhar_anis()

USE kinds,                 ONLY : DP
USE temperature,           ONLY : ntemp
USE control_thermo,        ONLY : ltherm_dos, ltherm_freq
USE internal_files_names,  ONLY : flfrq_thermo, flvec_thermo
USE anharmonic,            ONLY : celldm_t
USE ph_freq_anharmonic,    ONLY : celldmf_t
USE io_global,             ONLY : ionode
USE mp,                    ONLY : mp_sum
USE mp_world,              ONLY : world_comm

IMPLICIT NONE
INTEGER :: itemp
CHARACTER(LEN=256) :: filedata, filerap, fileout, gnu_filename, filenameps
LOGICAL :: all_geometry_done

CALL check_all_geometry_done(all_geometry_done)
IF (.NOT.all_geometry_done) RETURN
!
!    Anisotropic solid. Compute only the crystal parameters as a function
!    of temperature and the thermal expansion tensor
!
IF (ltherm_dos) THEN
   celldm_t=0.0_DP
   DO itemp = 1, ntemp
      CALL quadratic_fit_t(itemp)
   ENDDO
   IF (.NOT.ionode) celldm_t=0.0_DP
   CALL mp_sum(celldm_t, world_comm)
ENDIF

IF (ltherm_freq) THEN
   celldmf_t=0.0_DP
   DO itemp = 1, ntemp
      CALL quadratic_fit_t_ph(itemp)
   ENDDO
   IF (.NOT.ionode) celldmf_t=0.0_DP
   CALL mp_sum(celldmf_t, world_comm)
ENDIF
!
!  Check if the elastic constants are on file. If they are, the code
!  computes the elastic constants as a function of temperature interpolating
!  at the crystal parameters found in the quadratic/quartic fit
!
CALL check_el_cons()
CALL write_elastic_t()
CALL plot_elastic_t()

IF (ltherm_dos) CALL write_anhar_anis()
IF (ltherm_freq) CALL write_ph_freq_anhar_anis()
!
!    here we calculate and plot the gruneisen parameters along the given path.
!
CALL write_gruneisen_band_anis(flfrq_thermo,flvec_thermo)
!
!    plotband_sub writes the interpolated phonon at the requested temperature
!
CALL set_files_for_plot(4, flfrq_thermo, filedata, filerap, &
                                          fileout, gnu_filename, filenameps)
CALL plotband_sub(4, filedata, filerap, fileout, gnu_filename, filenameps)
!
!   and the next driver plots all the gruneisen parameters depending on the
!   lattice.
!
CALL plot_gruneisen_band_anis(flfrq_thermo)
!
!   here we calculate the gruneisen parameters on a mesh and then recompute
!   the thermal expansion from the gruneisen parameters
!
CALL fit_frequencies_anis()

CALL write_grun_anharmonic_anis()
CALL plot_anhar_anis()

RETURN
END SUBROUTINE manage_anhar_anis
