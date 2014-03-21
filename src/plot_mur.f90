!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_mur()
!
!  This is a driver to plot the energy and pressure as a function of volume
!  
!
USE kinds,           ONLY : DP
USE control_gnuplot, ONLY : flgnuplot, flpsmur, lgnuplot, gnuplot_command
USE gnuplot,         ONLY : gnuplot_start, gnuplot_end,  &
                            gnuplot_write_header,        &
                            gnuplot_ylabel,              &
                            gnuplot_xlabel,              &
                            gnuplot_write_horizontal_line, &
                            gnuplot_write_file_mul_data, &
                            gnuplot_write_file_mul_point
USE control_thermo,  ONLY : flevdat
USE control_mur,     ONLY : vmin_input, vmax_input
USE mp_images,       ONLY : my_image_id, root_image
USE io_global,       ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: filename, filename1, filename2, gnu_filename
INTEGER :: system
INTEGER :: ierr

IF ( my_image_id /= root_image ) RETURN

gnu_filename=TRIM(flgnuplot)//'_mur'
CALL gnuplot_start(gnu_filename)

filename=TRIM(flpsmur)
CALL gnuplot_write_header(filename, vmin_input, vmax_input, 0.0_DP, 0.0_DP ) 

filename1=TRIM(flevdat)//'_mur'
filename2=TRIM(flevdat)//'_mur1'
CALL gnuplot_xlabel('Volume ((a.u.)^3)',.FALSE.) 
CALL gnuplot_ylabel('Energy (Ry)',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename1,1,2,'red',.TRUE.,.FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_point(filename2,1,2,'red',.FALSE.,.TRUE.,.FALSE.)

CALL gnuplot_ylabel('Pressure (kbar)',.FALSE.) 
CALL gnuplot_write_horizontal_line(0.0_DP, 2, 'front', 'black',.FALSE.)
CALL gnuplot_write_file_mul_data(filename1,1,3,'red',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_mur

