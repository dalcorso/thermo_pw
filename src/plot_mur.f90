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
USE control_mur,      ONLY : lmurn, vmin_input, vmax_input, press_max, &
                             press_min
USE control_pressure, ONLY : pressure_kb
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
WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') press_min*ry_kbar, &
                                                    press_max*ry_kbar
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
SUBROUTINE plot_mur_t()
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
USE control_mur,      ONLY : lmurn, vmin_input, vmax_input, press_max, &
                             press_min
USE temperature,      ONLY : ntemp, temp_nstep
USE control_pressure, ONLY : pressure_kb
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: filename, filename1, filename2, gnu_filename, label
CHARACTER(LEN=12) :: color(8)
CHARACTER(LEN=8) :: float_to_char
CHARACTER(LEN=6) :: int_to_char
INTEGER :: itemp, istep
INTEGER :: ierr, system
LOGICAL :: first_step, last_step

IF ( my_image_id /= root_image ) RETURN

IF (temp_nstep>ntemp) RETURN

color(1)='color_red'
color(2)='color_green'
color(3)='color_blue'
color(4)='color_yellow'
color(5)='color_pink'
color(6)='color_cyan'
color(7)='color_orange'
color(8)='color_black'

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_mur_t'
CALL add_pressure(gnu_filename)
filename=TRIM(flpsmur)//'_t'
CALL add_pressure(filename)
filename=TRIM(filename)//TRIM(flext)

CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filename, vmin_input, vmax_input, &
                                0.0_DP, 0.0_DP, 1.0_DP, flext ) 

CALL gnuplot_xlabel('Volume ((a.u.)^3)',.FALSE.) 
!
!  Energy or enthalpy as a function of the volume.
!
istep=0
DO itemp=1,ntemp,temp_nstep
   first_step=(itemp==1)
   last_step=((itemp+temp_nstep)>ntemp)
   istep=MOD(istep,8)+1
   filename1="therm_files/"//TRIM(flevdat)//'_mur'//TRIM(int_to_char(itemp))
   CALL add_pressure(filename1)
   IF (first_step) THEN
      IF (pressure_kb /= 0.0_DP) THEN
         label='Gibbs free-energy (Ry)    p= '//&
                    &TRIM(float_to_char(pressure_kb,1))//' kbar'
         CALL gnuplot_ylabel(TRIM(label),.FALSE.) 
      ELSE
         CALL gnuplot_ylabel('Helmholtz Free Energy (Ry)',.FALSE.) 
      END IF
   ENDIF
   CALL gnuplot_write_file_mul_data(filename1,1,2,color(istep),first_step, &
                                                        last_step, .FALSE.)
ENDDO
!
!  Pressure as a function of the volume
!
istep=0
DO itemp=1,ntemp,temp_nstep
   first_step=(itemp==1)
   last_step=((itemp+temp_nstep)>ntemp)
   istep=MOD(istep,8)+1

   filename1="therm_files/"//TRIM(flevdat)//'_mur'//TRIM(int_to_char(itemp))
   CALL add_pressure(filename1)
   IF (first_step) THEN
      CALL gnuplot_ylabel('Pressure (kbar)',.FALSE.) 
      CALL gnuplot_write_horizontal_line(0.0_DP, 2, 'front', 'color_black', &
                                          .FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename1,1,4,color(istep),first_step, &
                                                        last_step, .FALSE.)
ENDDO
!
!  Enthalpy as a function of pressure
!
istep=0
DO itemp=1,ntemp,temp_nstep
   first_step=(itemp==1)
   last_step=((itemp+temp_nstep)>ntemp)
   istep=MOD(istep,8)+1
   filename1="therm_files/"//TRIM(flevdat)//'_mur'//TRIM(int_to_char(itemp))
   CALL add_pressure(filename1)

   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') press_min*ry_kbar, &
                                                    press_max*ry_kbar
      CALL gnuplot_write_command(TRIM(label),.FALSE.)

      CALL gnuplot_ylabel('Gibbs free-energy (Ry)',.FALSE.) 
   ENDIF
   CALL gnuplot_write_file_mul_data(filename1,4,3,color(istep),first_step, &
                                                  last_step, .FALSE.)
ENDDO
!
!   close the file and make the plot
!
CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_mur_t

