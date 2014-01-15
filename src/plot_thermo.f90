!
! Copyright (C) 2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_thermo()
!
!  This is a driver to plot the quantities written inside fltherm
!  
!
USE kinds,           ONLY : DP
USE control_gnuplot, ONLY : flgnuplot, ncount, flpstherm
USE gnuplot,         ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                            gnuplot_ylabel, &
                            gnuplot_xlabel, gnuplot_write_file_mul_data, &
                            gnuplot_write_vertical_line, &
                            gnuplot_set_fact
USE control_thermo,  ONLY : fltherm
USE thermodynamics,  ONLY : tmin, tmax
IMPLICIT NONE

CHARACTER(LEN=256) :: filename
CHARACTER(LEN=6), EXTERNAL :: int_to_char

ncount=ncount+1

filename=TRIM(flgnuplot)//TRIM(int_to_char(ncount))
CALL gnuplot_start(filename)

filename=TRIM(flpstherm)
CALL gnuplot_write_header(filename, tmin, tmax, 0.0_DP, 0.0_DP ) 

CALL gnuplot_xlabel('T (K)') 
CALL gnuplot_ylabel('Vibrational energy (J / mol)') 
CALL gnuplot_set_fact(1313313.0_DP) 

CALL gnuplot_write_file_mul_data(fltherm,1,2,'red',.TRUE.,.TRUE.)

CALL gnuplot_ylabel('Vibrational free energy (J / mol)') 
CALL gnuplot_write_file_mul_data(fltherm,1,3,'red',.TRUE.,.TRUE.)

CALL gnuplot_ylabel('Entropy (J / K / mol)') 
CALL gnuplot_write_file_mul_data(fltherm,1,4,'blue',.TRUE.,.TRUE.)

CALL gnuplot_ylabel('Specific heat C_v (J / K / mol)') 
CALL gnuplot_write_file_mul_data(fltherm,1,5,'blue',.TRUE.,.TRUE.)

CALL gnuplot_end()

RETURN
END SUBROUTINE plot_thermo

