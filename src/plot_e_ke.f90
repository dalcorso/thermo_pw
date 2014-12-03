!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_e_ke()
!
!  This is a driver to plot the quantities written inside flkeconv
!  
!
USE kinds,           ONLY : DP
USE control_gnuplot, ONLY : flgnuplot, flpskeconv, lgnuplot, gnuplot_command
USE gnuplot,         ONLY : gnuplot_start, gnuplot_end,  &
                            gnuplot_write_header,        &
                            gnuplot_set_gfact,           &
                            gnuplot_set_eref,            &
                            gnuplot_ylabel,              &
                            gnuplot_xlabel,              &
                            gnuplot_write_file_mul_data
USE control_thermo,  ONLY : flkeconv
USE control_conv,    ONLY : ke, nke, nkeden
USE thermo_mod,      ONLY : energy_geo
USE mp_images,       ONLY : my_image_id, root_image
USE io_global,       ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename
CHARACTER(LEN=6), EXTERNAL :: int_to_char
INTEGER :: system
INTEGER :: iden, ierr

IF ( my_image_id /= root_image ) RETURN

gnu_filename=TRIM(flgnuplot)//'_keconv'
CALL gnuplot_start(gnu_filename)

filename=TRIM(flpskeconv)
CALL gnuplot_write_header(filename, ke(1), ke(nke), 0.0_DP, 0.0_DP, 1.0_DP ) 

CALL gnuplot_xlabel('Kinetic energy (Ry)',.FALSE.) 
CALL gnuplot_ylabel('Total energy (mRy)',.FALSE.) 
CALL gnuplot_set_eref(energy_geo(nke*nkeden),.FALSE.) 
CALL gnuplot_set_gfact(1000._DP,.FALSE.) 

DO iden=1,nkeden
   IF (nkeden > 1) THEN
      filename=TRIM(flkeconv)//TRIM(int_to_char(iden))//'/'//TRIM(flkeconv)
      IF (iden==1) THEN
         CALL gnuplot_write_file_mul_data(filename,1,2,'red',.TRUE.,&
                                          .FALSE.,.FALSE.)
      ELSEIF (iden==nkeden) THEN
         CALL gnuplot_write_file_mul_data(filename,1,2,'blue',.FALSE.,&
                                          .TRUE.,.FALSE.)
      ELSE
         CALL gnuplot_write_file_mul_data(filename,1,2,'green',.FALSE.,&
                                          .FALSE.,.FALSE.)
      ENDIF
   ELSE
      filename=TRIM(flkeconv)
      CALL gnuplot_write_file_mul_data(filename,1,2,'red',.TRUE.,.TRUE.,&
                                                 .FALSE.)
   END IF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_e_ke
