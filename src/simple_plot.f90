!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE simple_plot(data_filename, psfilename, xlabel, ylabel, &
                       colore, xmin, xmax, ymin, ymax)
!
!  This is a simple routine which write a gnuplot script in flgnuplot,
!  The script create a graph of the data contained in data_filename,
!  in a postscript file called psfilename.
!
USE kinds,           ONLY : DP
USE control_gnuplot, ONLY : flgnuplot, ncount
USE gnuplot,       ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                          gnuplot_write_file_data, gnuplot_ylabel, &
                          gnuplot_xlabel, &
                          gnuplot_write_vertical_line, gnuplot_write_label
IMPLICIT NONE
REAL(DP), INTENT(IN) :: xmin, xmax, ymin, ymax
CHARACTER(LEN=*), INTENT(IN) :: data_filename, psfilename, xlabel, ylabel

CHARACTER(LEN=256) :: filename
CHARACTER(LEN=*) :: colore
CHARACTER(LEN=6), EXTERNAL :: int_to_char

ncount=ncount+1
filename=TRIM(flgnuplot)//TRIM(int_to_char(ncount))
CALL gnuplot_start(filename)

filename=TRIM(psfilename)
CALL gnuplot_write_header(filename, xmin, xmax, ymin, ymax ) 

CALL gnuplot_ylabel(TRIM(ylabel)) 
CALL gnuplot_xlabel(TRIM(xlabel)) 

CALL gnuplot_write_file_data(data_filename,colore,.TRUE.,.TRUE.)

CALL gnuplot_end()

RETURN
END SUBROUTINE simple_plot

