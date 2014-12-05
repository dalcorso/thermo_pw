!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_anhar()
!
!  This is a driver to plot the quantities written inside flanhar
!  
!
USE kinds,           ONLY : DP
USE control_gnuplot, ONLY : flgnuplot, flpsanhar, gnuplot_command, lgnuplot
USE gnuplot,         ONLY : gnuplot_start, gnuplot_end,  &
                            gnuplot_write_header,        &
                            gnuplot_ylabel,              &
                            gnuplot_xlabel,              &
                            gnuplot_write_file_mul_data, &
                            gnuplot_write_file_mul_data_sum, &
                            gnuplot_write_file_mul_point,  &
                            gnuplot_write_horizontal_line, &
                            gnuplot_write_command,         &
                            gnuplot_set_fact
USE control_thermo,  ONLY : flanhar
USE temperature,     ONLY : tmin, tmax
USE constants,       ONLY : ry_kbar
USE mp_images,       ONLY : my_image_id, root_image
USE io_global,       ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename1, filename2, filename3
CHARACTER(LEN=6), EXTERNAL :: int_to_char
INTEGER :: system
INTEGER :: ierr

IF ( my_image_id /= root_image ) RETURN

gnu_filename=TRIM(flgnuplot)//'_anhar'
CALL gnuplot_start(gnu_filename)

filename=TRIM(flpsanhar)
IF (tmin /= 1.0_DP) THEN
   CALL gnuplot_write_header(filename, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ELSE
   CALL gnuplot_write_header(filename, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ENDIF

CALL gnuplot_write_command('color_red="red"',.FALSE.)
CALL gnuplot_write_command('color_blue="blue"',.FALSE.)
CALL gnuplot_write_command('color_green="green"',.FALSE.)

filename=TRIM(flanhar)//'_ph'
filename1=TRIM(flanhar)//'.aux'
filename2=TRIM(flanhar)//'.aux_ph'
filename3=TRIM(flanhar)//'.aux_grun'

CALL gnuplot_xlabel('T (K)',.FALSE.) 
CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('Volume ((a.u.)^3)',.FALSE.) 
CALL gnuplot_write_file_mul_data(flanhar,1,2,'color_red',.TRUE.,.FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filename,1,2,'color_blue',.FALSE.,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('Bulk modulus (kbar)',.FALSE.) 
CALL gnuplot_write_file_mul_data(flanhar,1,3,'color_red',.TRUE.,.FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filename,1,3,'color_blue',.FALSE.,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('d B / d p',.FALSE.) 
CALL gnuplot_write_file_mul_data(flanhar,1,4,'color_red',.TRUE.,.FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',.FALSE.,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('Thermal expansion ({/Symbol b} x 10^{6})',.FALSE.) 
CALL gnuplot_write_file_mul_data(flanhar,1,5,'color_red',.TRUE.,.FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filename,1,5,'color_blue',.FALSE.,.FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filename3,1,2,'color_green',.FALSE.,.TRUE.,.FALSE.)
!
!  put as a comment the possibility to plot also the experimental data
!
CALL gnuplot_write_file_mul_data(filename3,1,2,'color_green',.FALSE.,.FALSE.,.TRUE.)
CALL gnuplot_write_file_mul_point('anhar.exp',1,2,'color_red',.FALSE.,.TRUE.,.TRUE.)

CALL gnuplot_set_fact(1313313.0_DP,.FALSE.)
CALL gnuplot_ylabel('Heat capacity C_v (J / K / N / mol)',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename1,1,3,'color_red',.TRUE.,.FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filename2,1,3,'color_blue',.FALSE.,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1313313.0_DP,.FALSE.)
CALL gnuplot_ylabel('Heat capacity C_p (J / K / N / mol)',.FALSE.) 
CALL gnuplot_write_file_mul_data_sum(filename1,1,3,4,'color_red',.TRUE.,.FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data_sum(filename2,1,3,4,'color_blue',.FALSE.,.TRUE.,.FALSE.)
!
!  put as a comment the possibility to plot also the experimental data
!
CALL gnuplot_write_file_mul_data_sum(filename2,1,3,4,'color_blue',.FALSE.,.FALSE.,.TRUE.)
CALL gnuplot_write_file_mul_point('cv.exp',1,2,'color_red',.FALSE.,.TRUE.,.TRUE.)

CALL gnuplot_set_fact(1313313.0_DP,.FALSE.)
CALL gnuplot_ylabel('C_p - C_v (J / K / N / mol)',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename1,1,4,'color_red',.TRUE.,.FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filename2,1,4,'color_blue',.FALSE.,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1._DP,.FALSE.)
CALL gnuplot_ylabel('B_S - B_T (kbar)',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename1,1,5,'color_red',.TRUE.,.FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filename2,1,5,'color_blue',.FALSE.,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('{/Symbol g}',.FALSE.) 
CALL gnuplot_write_horizontal_line(0.0_DP, 2, 'front', 'black', .FALSE.)
CALL gnuplot_write_file_mul_data(filename1,1,2,'color_red',.TRUE.,.FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filename2,1,2,'color_blue',.FALSE.,.FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filename3,1,3,'color_green',.FALSE.,.TRUE.,.FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_anhar

