!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE manage_bands()
!
!   this is a driver that controls the band structure / electronic dos
!   calculation
!
USE control_thermo,   ONLY : ldos_syn_1, spin_component
USE control_bands,    ONLY : nbnd_bands
USE control_2d_bands, ONLY : only_bands_plot

USE wvfct,            ONLY : nbnd
USE lsda_mod,         ONLY : nspin
USE control_flags,    ONLY : lbands

USE io_global,        ONLY : stdout

IMPLICIT NONE

LOGICAL :: exit_status
INTEGER :: nspin0, ierr
CHARACTER(LEN=256) :: filedata, filerap, fileout, gnu_filename, filenameps
 
IF (.NOT.only_bands_plot) THEN
   IF (ldos_syn_1) THEN
      CALL set_dos_kpoints()
   ELSE
      CALL set_paths_disp()
      CALL set_k_points()
   ENDIF
!
!   by default in a band structure calculation we double the number of
!   computed bands
!
   IF (nbnd_bands == 0) nbnd_bands = 2*nbnd
   IF (nbnd_bands > nbnd) nbnd = nbnd_bands
   WRITE(stdout,'(/,2x,76("+"))')
   WRITE(stdout,'(5x,"Doing a non self-consistent calculation", i5)') 
   WRITE(stdout,'(2x,76("+"),/)')
   IF (ldos_syn_1) THEN
      lbands=.FALSE.
   ELSE
      lbands=.TRUE.
   ENDIF
   CALL set_fft_mesh()
   CALL do_pwscf(exit_status, .FALSE.)
   IF (ldos_syn_1) THEN
      CALL dos_sub()
      CALL plot_dos()
      CALL write_el_thermo()
      CALL plot_el_thermo()
   ELSE
      nspin0=nspin
      IF (nspin==4) nspin0=1
      DO spin_component = 1, nspin0
         CALL bands_sub()
         CALL read_minimal_info(.FALSE.,ierr)
         CALL set_files_for_plot(1, ' ', filedata, filerap, &
                                           fileout, gnu_filename, filenameps)
         CALL plotband_sub(1,filedata, filerap, fileout, &
                                           gnu_filename, filenameps)
      ENDDO
   ENDIF
ELSE
   CALL read_minimal_info(.TRUE., ierr)
   IF (ierr /= 0) THEN
!
!    The code might have only the punch files but not the bands files
!
      CALL set_paths_disp()
      nspin0=nspin
      IF (nspin==4) nspin0=1
      DO spin_component = 1, nspin0
         CALL bands_sub()
      END DO
   END IF
   CALL set_files_for_plot(1, ' ', filedata, filerap, fileout,  &
                                                gnu_filename, filenameps)
   CALL plotband_sub(1, filedata, filerap, fileout, &
                                                gnu_filename, filenameps)
ENDIF
RETURN
END SUBROUTINE manage_bands
