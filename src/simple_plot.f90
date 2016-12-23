!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE simple_plot(ext, data_filename, psfilename, xlabel, ylabel, &
                       colore, xmin, xmax, ymin, ymax)
!
!  This is a simple routine which write a gnuplot script in flgnuplot,
!  The script create a graph of the data contained in data_filename,
!  in a postscript file called psfilename.
!
USE kinds,           ONLY : DP
USE control_gnuplot, ONLY : flgnuplot, gnuplot_command, lgnuplot
USE gnuplot,         ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                            gnuplot_write_file_data, gnuplot_ylabel, &
                            gnuplot_xlabel, gnuplot_write_command
USE mp_images,       ONLY : root_image, my_image_id
USE io_global,       ONLY : ionode

IMPLICIT NONE
REAL(DP), INTENT(IN) :: xmin, xmax, ymin, ymax
CHARACTER(LEN=*), INTENT(IN) :: data_filename, psfilename, xlabel, ylabel

INTEGER :: ierr
CHARACTER(LEN=256) :: gnu_filename, filename
CHARACTER(LEN=*) :: colore, ext
CHARACTER(LEN=6), EXTERNAL :: int_to_char
INTEGER :: system

IF ( my_image_id /= root_image ) RETURN

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//TRIM(ext)
CALL gnuplot_start(gnu_filename)

filename=TRIM(psfilename)
CALL gnuplot_write_header(filename, xmin, xmax, ymin, ymax, 1.0_DP ) 

CALL gnuplot_ylabel(TRIM(ylabel), .FALSE.) 
CALL gnuplot_xlabel(TRIM(xlabel), .FALSE.) 
CALL gnuplot_write_command('plot_width=2',.FALSE.)

CALL gnuplot_write_file_data(data_filename,'plot_width',colore,.TRUE.,.TRUE.,&
                                                                      .FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE simple_plot

