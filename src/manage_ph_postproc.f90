!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE manage_ph_postproc(igeom)
!
!  This driver computes the interatomic force constants and
!  interpolates the phonon frequencies to make a phonon dispersion
!  and then computes them on a thick mesh to compute the harmonic
!  thermodynamic quantities.
!
USE control_thermo,   ONLY : ltherm, ltherm_dos, ltherm_freq, set_internal_path
USE control_paths,    ONLY : disp_nqs
USE control_phrun,    ONLY : auxdyn

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom

CHARACTER(LEN=256) :: filedata, filerap, fileout, gnu_filename, filenameps
!
!   Compute the interatomic force constants from the dynamical matrices
!   written on file
!
CALL q2r_sub(auxdyn) 
!
!    compute interpolated dispersions
!
IF (set_internal_path) CALL set_bz_path()
CALL write_ph_dispersions()
CALL set_files_for_plot(2, ' ', filedata, filerap, fileout, &
                                                    gnu_filename, filenameps)
IF (disp_nqs>0) CALL plotband_sub(2, filedata, filerap, fileout, &
                                                    gnu_filename, filenameps)
!
!   Compute the harmonic thermodynamic quantities
!
IF (ltherm) THEN
!
!    the frequencies on a uniform mesh are interpolated or read from disk 
!    if available and saved in ph_freq_save 
!
   CALL write_ph_freq(igeom)
   IF (ltherm_freq.OR..NOT.ltherm_dos) CALL write_thermo_ph(igeom)
 
   IF (ltherm_dos) THEN
      CALL write_phdos(igeom)
      CALL plot_phdos()
      CALL write_thermo(igeom)
   ENDIF
   CALL plot_thermo()
ENDIF
CALL clean_ifc_variables()

RETURN
END SUBROUTINE manage_ph_postproc

