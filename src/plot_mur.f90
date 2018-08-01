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
USE kinds,            ONLY : DP
USE control_gnuplot,  ONLY : flgnuplot, lgnuplot, gnuplot_command
USE postscript_files, ONLY : flpsmur
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_write_horizontal_line, &
                             gnuplot_write_file_mul_data, &
                             gnuplot_write_file_mul_point
USE data_files,       ONLY : flevdat
USE control_mur,      ONLY : lmurn, vmin_input, vmax_input
USE control_pressure, ONLY : pressure_kb
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: filename, filename1, filename2, gnu_filename, label
CHARACTER(LEN=8) :: float_to_char
INTEGER :: ierr, system

IF ( my_image_id /= root_image ) RETURN

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_mur'
filename=TRIM(flpsmur)//'.ps'
IF (pressure_kb /= 0.0_DP) THEN
   gnu_filename=TRIM(gnu_filename)//'.'//TRIM(float_to_char(pressure_kb,1))
   filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))
END IF

CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filename, vmin_input, vmax_input, 0.0_DP, 0.0_DP, &
                          1.0_DP ) 

CALL gnuplot_xlabel('Volume ((a.u.)^3)',.FALSE.) 

filename1="energy_files/"//TRIM(flevdat)//'_mur'
filename2="energy_files/"//TRIM(flevdat)

IF (pressure_kb /= 0.0_DP) THEN
   filename1=TRIM(filename1)//'.'//TRIM(float_to_char(pressure_kb,1))
   filename2=TRIM(filename2)//'.'//TRIM(float_to_char(pressure_kb,1))
   label='Enthalpy (Ry)    p= '//TRIM(float_to_char(pressure_kb,1))//' kbar'
   CALL gnuplot_ylabel(TRIM(label),.FALSE.) 
ELSE
   CALL gnuplot_ylabel('Energy (Ry)',.FALSE.) 
END IF
CALL gnuplot_write_file_mul_data(filename1,1,2,'color_red',.TRUE.,.NOT.lmurn, &
                                                                        .FALSE.)
IF (lmurn) &
   CALL gnuplot_write_file_mul_point(filename2,1,2,'color_red',.FALSE.,.TRUE.,&
                                                                        .FALSE.)

CALL gnuplot_ylabel('Pressure (kbar)',.FALSE.) 
IF (pressure_kb /= 0.0_DP) &
   CALL gnuplot_write_horizontal_line(pressure_kb, 2, 'front', 'color_green',&
                                                                     .FALSE.)
CALL gnuplot_write_horizontal_line(0.0_DP, 2, 'front', 'color_black',.FALSE.)
CALL gnuplot_write_file_mul_data(filename1,1,3,'color_red',.TRUE.,.TRUE.,&
                                                                        .FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_mur

