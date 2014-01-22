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
USE control_gnuplot, ONLY : flgnuplot, flpstherm
USE gnuplot,         ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                            gnuplot_ylabel, &
                            gnuplot_xlabel, gnuplot_write_file_mul_data, &
                            gnuplot_set_fact
USE control_thermo,  ONLY : fltherm
USE thermodynamics,  ONLY : tmin, tmax
USE mp_images,       ONLY : root_image, my_image_id

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=6), EXTERNAL :: int_to_char

IF ( my_image_id /= root_image ) RETURN

filename=TRIM(flgnuplot)//'_therm'
CALL gnuplot_start(filename)

IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(flpstherm, 0.0_DP, tmax, 0.0_DP, 0.0_DP ) 
ELSE
   CALL gnuplot_write_header(flpstherm, tmin, tmax, 0.0_DP, 0.0_DP ) 
ENDIF
CALL gnuplot_xlabel('T (K)', .FALSE.) 
CALL gnuplot_ylabel('Vibrational energy (kJ / N / mol)',.FALSE.) 
CALL gnuplot_set_fact(1313.3130_DP, .FALSE.) 

CALL gnuplot_write_file_mul_data(fltherm,1,2,'red',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_ylabel('Vibrational free energy (kJ / N / mol)', .FALSE.) 
CALL gnuplot_write_file_mul_data(fltherm,1,3,'red',.TRUE.,.TRUE., .FALSE.)

CALL gnuplot_set_fact(1313313.0_DP, .FALSE.) 
CALL gnuplot_ylabel('Entropy (J / K / N / mol))',.FALSE.) 
CALL gnuplot_write_file_mul_data(fltherm,1,4,'blue',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_ylabel('Heat capacity C_v (J / K / N / mol)',.FALSE.) 
CALL gnuplot_write_file_mul_data(fltherm,1,5,'blue',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_end()

RETURN
END SUBROUTINE plot_thermo

