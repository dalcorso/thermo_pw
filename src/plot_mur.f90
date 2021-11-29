!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------
SUBROUTINE plot_mur()
!-------------------------------------------------------------------
!
!  This is a driver to plot the energy and pressure as a function of volume
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE control_gnuplot,  ONLY : flgnuplot, lgnuplot, gnuplot_command, flext
USE postscript_files, ONLY : flpsmur
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_write_command,       &
                             gnuplot_write_horizontal_line, &
                             gnuplot_write_file_mul_data, &
                             gnuplot_write_file_mul_point
USE data_files,       ONLY : flevdat
USE control_mur,      ONLY : lmurn
USE control_vol,      ONLY : vmin_input, vmax_input
USE control_pressure, ONLY : pressure_kb, pmax, pmin
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: filename, filename1, filename2, gnu_filename, label
CHARACTER(LEN=8) :: float_to_char
INTEGER :: ierr, system

IF ( my_image_id /= root_image ) RETURN

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_mur'
CALL add_pressure(gnu_filename)
filename=TRIM(flpsmur)
CALL add_pressure(filename)
filename=TRIM(filename)//TRIM(flext)

CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filename, vmin_input, vmax_input, 0.0_DP, 0.0_DP, &
                          1.0_DP, flext ) 

CALL gnuplot_xlabel('Volume ((a.u.)^3)',.FALSE.) 
!
!  Energy or enthalpy as a function of the volume.
!
filename1="energy_files/"//TRIM(flevdat)//'_mur'
filename2="energy_files/"//TRIM(flevdat)
CALL add_pressure(filename1)
CALL add_pressure(filename2)
IF (pressure_kb /= 0.0_DP) THEN
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
!
!  Pressure as a function of the volume
!
CALL gnuplot_ylabel('Pressure (kbar)',.FALSE.) 
IF (pressure_kb /= 0.0_DP) &
   CALL gnuplot_write_horizontal_line(pressure_kb, 2, 'front', 'color_green',&
                                                                     .FALSE.)
CALL gnuplot_write_horizontal_line(0.0_DP, 2, 'front', 'color_black',.FALSE.)
CALL gnuplot_write_file_mul_data(filename1,1,4,'color_red',.TRUE.,.TRUE.,&
                                                                     .FALSE.)
!
!  Enthalpy as a function of pressure
!
CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 
WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                    pmax*ry_kbar
CALL gnuplot_write_command(TRIM(label),.FALSE.)
CALL gnuplot_ylabel('Enthalpy (Ry)',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename1,4,3,'color_red',.TRUE.,.TRUE.,&
                                                                   .FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_mur

!-------------------------------------------------------------------
SUBROUTINE plot_mur_p()
!-------------------------------------------------------------------
!
!  This is a driver to plot the energy and pressure as a function of volume
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE control_gnuplot,  ONLY : flgnuplot, lgnuplot, gnuplot_command, flext
USE postscript_files, ONLY : flpsmur
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_write_command,       &
                             gnuplot_write_file_mul_data
USE control_ev,       ONLY : ieos
USE data_files,       ONLY : flevdat
USE temperature,      ONLY : ntemp_plot
USE control_pressure, ONLY : press, npress, npress_plot
USE control_pressure, ONLY : pressure_kb, pmax, pmin
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: filename, filename1, gnu_filename, label
CHARACTER(LEN=8) :: float_to_char
INTEGER :: ierr, system

IF ( my_image_id /= root_image ) RETURN
IF (npress_plot==0.AND.ntemp_plot==0) RETURN

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_mur_press'
CALL add_pressure(gnu_filename)
filename=TRIM(flpsmur)//'_press'
CALL add_pressure(filename)
filename=TRIM(filename)//TRIM(flext)

CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filename, press(1), press(npress), 0.0_DP, 0.0_DP, &
                          1.0_DP, flext ) 

CALL gnuplot_xlabel('Pressure (kbar)',.FALSE.) 
!
!  Volume as a function of pressure
!
filename1="energy_files/"//TRIM(flevdat)//'_mur_press'
CALL add_pressure(filename1)
CALL gnuplot_ylabel('Volume (a.u.)^3',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename1,1,2,'color_red',.TRUE.,.TRUE., &
                                                                  .FALSE.)
!
!  Bulk modulus as a function of pressure
!
CALL gnuplot_ylabel('Bulk modulus B_0 (kbar)',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename1,1,3,'color_red',.TRUE.,.TRUE., &
                                                                  .FALSE.)
!
!  Pressure derivative of the bulk modulus as a function of pressure
!
CALL gnuplot_ylabel('Pressure derivative of B_0',.FALSE.) 
WRITE(label,'("set yrange [3.5:5.0]")') 
CALL gnuplot_write_command(TRIM(label),.FALSE.)
CALL gnuplot_write_file_mul_data(filename1,1,4,'color_red',.TRUE.,.TRUE., &
                                                                  .FALSE.)
!
!  Second pressure derivative of the bulk modulus as a function of pressure
!
IF (ieos==2) THEN
   CALL gnuplot_ylabel('Second pressure derivative of B_0',.FALSE.) 
   WRITE(label,'("set yrange [-0.1:0.1]")') 
   CALL gnuplot_write_command(TRIM(label),.FALSE.)
   CALL gnuplot_write_file_mul_data(filename1,1,5,'color_red',.TRUE.,.TRUE., &
                                                                  .FALSE.)
ENDIF

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_mur_p
