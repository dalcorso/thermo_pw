!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_e_nk()
!
!  This is a driver to plot the quantities written inside flnkconv
!
USE kinds,           ONLY : DP
USE control_gnuplot, ONLY : flgnuplot, lgnuplot, gnuplot_command
USE postscript_files, ONLY : flpsnkconv
USE gnuplot,         ONLY : gnuplot_start, gnuplot_end,  &
                            gnuplot_write_header,        &
                            gnuplot_ylabel,              &
                            gnuplot_xlabel,              &
                            gnuplot_set_eref,            &
                            gnuplot_set_gfact,           &
                            gnuplot_write_file_mul_data
USE thermo_mod,      ONLY : energy_geo
USE data_files,      ONLY : flnkconv
USE control_conv,    ONLY : nk_test, nnk, nsigma
USE mp_images,       ONLY : my_image_id, root_image
USE io_global,       ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename
CHARACTER(LEN=6), EXTERNAL :: int_to_char
INTEGER :: system
INTEGER :: isigma, ierr
REAL(DP) :: xmin, xmax

IF ( my_image_id /= root_image ) RETURN

gnu_filename=TRIM(flgnuplot)//'_nkconv'
CALL gnuplot_start(gnu_filename)

filename=TRIM(flpsnkconv)
xmin = nk_test(1)
xmax = nk_test(nnk)
CALL gnuplot_write_header(filename, xmin, xmax, 0.0_DP, 0.0_DP, 1.0_DP ) 

CALL gnuplot_xlabel(' nk ',.FALSE.) 
CALL gnuplot_ylabel('Total energy error (mRy)',.FALSE.) 
CALL gnuplot_set_eref(energy_geo(nnk),.FALSE.) 
CALL gnuplot_set_gfact(1000._DP,.FALSE.) 

DO isigma=1,nsigma
   IF (nsigma > 1) THEN
      filename=TRIM(flnkconv)//TRIM(int_to_char(isigma))//'/'//TRIM(flnkconv)
      IF (isigma==1) THEN
         CALL gnuplot_write_file_mul_data(filename,1,2,'color_red',.TRUE.,.FALSE.,&
                                                             .FALSE.)
      ELSEIF (isigma==nsigma) THEN
         CALL gnuplot_write_file_mul_data(filename,1,2,'color_blue',.FALSE.,.TRUE.,&
                                          .FALSE.)
      ELSE
         CALL gnuplot_write_file_mul_data(filename,1,2,'color_green',.FALSE.,.FALSE.,&
                                  .FALSE.)
      ENDIF
   ELSE
      filename=TRIM(flnkconv)
      CALL gnuplot_write_file_mul_data(filename,1,2,'color_red',.TRUE.,.TRUE.,.FALSE.)
   END IF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '// TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_e_nk
