!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_thermo(igeom)
!
!  This is a driver to plot the quantities written inside fltherm
!  
!
USE kinds,           ONLY : DP
USE control_gnuplot, ONLY : flgnuplot, flpstherm, gnuplot_command, lgnuplot
USE gnuplot,         ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                            gnuplot_ylabel, gnuplot_write_command, &
                            gnuplot_xlabel, gnuplot_write_file_mul_data, &
                            gnuplot_set_fact
USE control_thermo,  ONLY : fltherm
USE temperature,     ONLY : tmin, tmax
USE mp_images,       ONLY : root_image, my_image_id
USE io_global,       ONLY : ionode

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom
CHARACTER(LEN=256) :: gnu_filename, filename
CHARACTER(LEN=6), EXTERNAL :: int_to_char
INTEGER :: system
INTEGER :: ierr

IF ( my_image_id /= root_image ) RETURN

gnu_filename=TRIM(flgnuplot)//'_therm'
CALL gnuplot_start(gnu_filename)

IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(flpstherm, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ELSE
   CALL gnuplot_write_header(flpstherm, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ENDIF
filename=TRIM(fltherm)//'_ph'
CALL gnuplot_xlabel('T (K)', .FALSE.) 
CALL gnuplot_ylabel('Vibrational energy (kJ / N / mol)',.FALSE.) 
CALL gnuplot_set_fact(1313.3130_DP, .FALSE.) 

CALL gnuplot_write_command('color_red="red"',.FALSE.)
CALL gnuplot_write_command('color_blue="blue"',.FALSE.)

CALL gnuplot_write_file_mul_data(fltherm,1,2,'color_red',.TRUE.,.FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filename,1,2,'color_blue',.FALSE.,.TRUE.,.FALSE.)

CALL gnuplot_ylabel('Vibrational free energy (kJ / N / mol)', .FALSE.) 
CALL gnuplot_write_file_mul_data(fltherm,1,3,'color_red',.TRUE.,.FALSE., .FALSE.)
CALL gnuplot_write_file_mul_data(filename,1,3,'color_blue',.FALSE.,.TRUE., .FALSE.)

CALL gnuplot_set_fact(1313313.0_DP, .FALSE.) 
CALL gnuplot_ylabel('Entropy (J / K / N / mol))',.FALSE.) 
CALL gnuplot_write_file_mul_data(fltherm,1,4,'color_red',.TRUE.,.FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',.FALSE.,.TRUE.,.FALSE.)

CALL gnuplot_ylabel('Heat capacity C_v (J / K / N / mol)',.FALSE.) 
CALL gnuplot_write_file_mul_data(fltherm,1,5,'color_red',.TRUE.,.FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filename,1,5,'color_blue',.FALSE.,.TRUE.,.FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_thermo
