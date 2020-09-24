!
! Copyright (C) 2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE plot_epsilon_omega_q()
!-------------------------------------------------------------------------
!
!  This is a driver to plot the quantities inside the epsilon files. It makes 
!  four plots. The first two contain the real and the imaginary part of the 
!  inverse of the dielectric constant (q w), the other two the real and 
!  imaginary part of the dielectric constant of (q w).
!
USE kinds,            ONLY : DP
USE constants,        ONLY : rytoev
USE control_gnuplot,  ONLY : flgnuplot, lgnuplot, gnuplot_command, flext
USE postscript_files, ONLY : flpsepsilon
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,   &
                            gnuplot_write_header,          &
                            gnuplot_ylabel,                &
                            gnuplot_xlabel,                &
                            gnuplot_write_command,         &
                            gnuplot_write_file_mul_data_minus, &
                            gnuplot_write_file_mul_data
USE data_files,       ONLY : flepsilon
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, string
INTEGER :: im
INTEGER :: ierr, system

IF ( my_image_id /= root_image ) RETURN

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_epsilon'
CALL gnuplot_start(gnu_filename)

filename=TRIM(flpsepsilon)//TRIM(flext)
CALL gnuplot_write_header(filename, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
                                                               rytoev, flext ) 

CALL gnuplot_xlabel('{/Symbol w}  (eV)',.FALSE.) 

CALL gnuplot_ylabel('Re 1 / {/Symbol e} (q, {/Symbol w})',.FALSE.) 
filename='dynamical_matrices/'//TRIM(flepsilon)
CALL gnuplot_write_file_mul_data(filename,2,4,'color_red',.TRUE.,&
                                                .TRUE.,.FALSE.)

CALL gnuplot_ylabel('- Im 1 / {/Symbol e} (q, {/Symbol w})',.FALSE.) 
CALL gnuplot_write_file_mul_data_minus(filename,2,5,'color_red',.TRUE.,&
                                                .TRUE.,.FALSE.)

CALL gnuplot_ylabel('{/Symbol e}_1 (q, {/Symbol w})',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename,2,6,'color_red',.TRUE.,&
                                                .TRUE.,.FALSE.)

CALL gnuplot_ylabel('{/Symbol e}_2 (q, {/Symbol w})',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename,2,7,'color_red',.TRUE.,&
                                                .TRUE.,.FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_epsilon_omega_q
