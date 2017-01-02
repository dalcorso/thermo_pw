!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE manage_anhar()

USE temperature,           ONLY : ntemp
USE control_thermo,        ONLY : ltherm_dos, ltherm_freq
USE internal_files_names,  ONLY : flfrq_thermo, flvec_thermo

IMPLICIT NONE

INTEGER :: itemp
CHARACTER(LEN=256) :: filedata, filerap, fileout, gnu_filename, filenameps

DO itemp = 1, ntemp
   IF (ltherm_dos) CALL do_ev_t(itemp)
   IF (ltherm_freq) CALL do_ev_t_ph(itemp)
ENDDO
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

USE temperature,           ONLY : ntemp
USE control_thermo,        ONLY : ltherm_dos, ltherm_freq
USE internal_files_names,  ONLY : flfrq_thermo, flvec_thermo

IMPLICIT NONE
INTEGER :: itemp
CHARACTER(LEN=256) :: filedata, filerap, fileout, gnu_filename, filenameps
!
!    Anisotropic solid. Compute only the crystal parameters as a function
!    of temperature and the thermal expansion tensor
!
DO itemp = 1, ntemp
   IF (ltherm_dos) CALL quadratic_fit_t(itemp)
   IF (ltherm_freq) CALL quadratic_fit_t_ph(itemp)
ENDDO
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
