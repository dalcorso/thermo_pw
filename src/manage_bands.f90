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
!   calculation. It manage also the recovering of the bands from file.
!
USE control_thermo,   ONLY : ldos_syn_1, spin_component
USE control_bands,    ONLY : nbnd_bands
USE control_2d_bands, ONLY : only_bands_plot
USE control_paths,    ONLY : q2d, is_a_path

USE wvfct,            ONLY : nbnd
USE lsda_mod,         ONLY : nspin
USE control_flags,    ONLY : lbands

USE io_global,        ONLY : stdout

IMPLICIT NONE

LOGICAL :: exit_status
INTEGER :: nspin0, ierr
CHARACTER(LEN=256) :: filedata, filerap, fileout, gnu_filename, filenameps
 
ierr=0
IF (only_bands_plot) CALL read_minimal_info(.TRUE., ierr)

IF (ldos_syn_1) THEN
   IF (.NOT.only_bands_plot) CALL set_dos_kpoints()
   lbands=.FALSE.
ELSE
   CALL set_paths_disp()
   CALL set_k_points()
   lbands=.TRUE.
ENDIF
!
!   by default in a band structure calculation we double the number of
!   computed bands
!
IF (nbnd_bands == 0) nbnd_bands = 2*nbnd
IF (nbnd_bands > nbnd) nbnd = nbnd_bands

IF (.NOT.only_bands_plot) THEN
   WRITE(stdout,'(/,2x,76("+"))')
   WRITE(stdout,'(5x,"Doing a non self-consistent calculation", i5)') 
   WRITE(stdout,'(2x,76("+"),/)')
   CALL set_fft_mesh()
   CALL do_pwscf(exit_status, .FALSE.)
ENDIF

IF (ldos_syn_1) THEN
   IF (.NOT.only_bands_plot.OR.ierr/=0) THEN
      CALL dos_sub()
      CALL read_minimal_info(.FALSE.,ierr)
   ENDIF
   CALL plot_dos()
   CALL write_el_thermo()
   CALL plot_el_thermo()
ELSE
   nspin0=nspin
   IF (nspin==4) nspin0=1
!
!  If ierr/=0 the only possibility to make the calculation is to have 
!  the outdir directory and reanalyze the bands
!
   IF (.NOT.only_bands_plot.OR.ierr/=0) THEN
      DO spin_component = 1, nspin0
         CALL bands_sub()
         CALL read_minimal_info(.FALSE.,ierr)
      END DO
   ENDIF
!
!  If the code arrives here it could read the bands and possibly analyze
!  their symmetry. Plot them on output if they are on a path
!
   IF (is_a_path) THEN
      DO spin_component = 1, nspin0
         CALL set_files_for_plot(1, ' ', filedata, filerap, &
                                           fileout, gnu_filename, filenameps)
         CALL plotband_sub(1,filedata, filerap, fileout, &
                                           gnu_filename, filenameps)
      ENDDO
   ELSEIF (q2d) THEN
      spin_component=1
      CALL set_files_for_plot(1, ' ', filedata, filerap, &
                                           fileout, gnu_filename, filenameps)
      CALL plot_ef(filedata, gnu_filename, filenameps)
   ENDIF
ENDIF

RETURN
END SUBROUTINE manage_bands
