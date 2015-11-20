!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_anhar_anis()
!
!  This is a driver to plot the quantities written inside flanhar for the
!  case of anisotropic systems. Presently plots celldm as a function of
!  temperature (one parameter per plot), the thermal expansion tensor
!  components (all in the same plot), the volume as a function of temperature,
!  the volume thermal expansion as a function of temperature.
!  
!
USE kinds,           ONLY : DP
USE control_gnuplot, ONLY : flgnuplot, gnuplot_command, lgnuplot
USE postscript_files, ONLY : flpsanhar
USE gnuplot,         ONLY : gnuplot_start, gnuplot_end,  &
                            gnuplot_write_header,        &
                            gnuplot_ylabel,              &
                            gnuplot_xlabel,              &
                            gnuplot_write_file_mul_data, &
                            gnuplot_write_file_mul_point, &
                            gnuplot_set_fact
USE data_files,  ONLY : flanhar
USE grun_anharmonic, ONLY : done_grun
USE control_pwrun, ONLY : ibrav_save
USE temperature,     ONLY : tmin, tmax
USE control_pressure, ONLY : pressure, pressure_kb
USE mp_images,       ONLY : my_image_id, root_image
USE io_global,       ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename1, filename2, &
                      filename3, filename4, filenameps
CHARACTER(LEN=8) :: float_to_char
INTEGER :: system
INTEGER :: ierr

IF ( my_image_id /= root_image ) RETURN

gnu_filename=TRIM(flgnuplot)//'_anhar'
IF (pressure /= 0.0_DP) &
   gnu_filename=TRIM(gnu_filename)//'.'//float_to_char(pressure_kb,1)

CALL gnuplot_start(gnu_filename)

filenameps=TRIM(flpsanhar)
IF (pressure /= 0.0_DP) &
   filenameps=TRIM(filenameps)//'.'//float_to_char(pressure_kb,1)

IF (tmin /= 1.0_DP) THEN
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ELSE
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ENDIF

filename=TRIM(flanhar)
filename1=TRIM(flanhar)//'_ph'
filename2=TRIM(flanhar)//'.celldm'
filename3=TRIM(flanhar)//'.celldm_ph'
filename4=TRIM(flanhar)//'.aux_grun'

IF (pressure /= 0.0_DP) THEN
   filename=TRIM(filename)//'.'//float_to_char(pressure_kb,1)
   filename1=TRIM(filename1)//'.'//float_to_char(pressure_kb,1)
   filename2=TRIM(filename2)//'.'//float_to_char(pressure_kb,1)
   filename3=TRIM(filename3)//'.'//float_to_char(pressure_kb,1)
   filename4=TRIM(filename4)//'.'//float_to_char(pressure_kb,1)
END IF

CALL gnuplot_xlabel('T (K)',.FALSE.) 
CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('a (a.u.)',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename2,1,2,'color_red',.TRUE., &
                                                      .FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filename3,1,2,'color_blue',.FALSE., &
                                                      .TRUE.,.FALSE.)

IF (ibrav_save==4.OR.ibrav_save==5.OR.ibrav_save==6.OR.ibrav_save==7) THEN
   CALL gnuplot_ylabel('c/a ',.FALSE.) 
   CALL gnuplot_write_file_mul_data(filename2,1,3,'color_red',.TRUE., &
                                                              .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filename3,1,3,'color_blue',.FALSE., &
                                                              .TRUE.,.FALSE.)
ELSEIF (ibrav_save==8.OR.ibrav_save==9.OR.ibrav_save==10.OR.ibrav_save==11) THEN
   CALL gnuplot_ylabel('b/a ',.FALSE.) 
   CALL gnuplot_write_file_mul_data(filename2,1,3,'color_red',.TRUE., &
                                                              .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filename3,1,3,'color_blue',.FALSE., &
                                                              .TRUE.,.FALSE.)
   CALL gnuplot_ylabel('c/a ',.FALSE.) 
   CALL gnuplot_write_file_mul_data(filename2,1,4,'color_red',.TRUE., &
                                                              .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filename3,1,4,'color_blue',.FALSE., &
                                                              .TRUE.,.FALSE.)
ENDIF
CALL gnuplot_ylabel('Volume ((a.u.)^3)',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename,1,2,'color_red',.TRUE.,.FALSE.,&
                                                               .FALSE.)
CALL gnuplot_write_file_mul_data(filename1,1,2,'color_blue',.FALSE.,.TRUE.,&
                                                               .FALSE.)
CALL gnuplot_ylabel('Linear thermal expansion {/Symbol a} x 10^6 (K^{-1})',.FALSE.) 
IF (ibrav_save==1.OR.ibrav_save==2.OR.ibrav_save==3) THEN
   CALL gnuplot_write_file_mul_data(filename2,1,3,'color_red',.TRUE.,.FALSE.,&
                                                                .FALSE.)
   IF (done_grun) &
      CALL gnuplot_write_file_mul_data(filename4,1,3,'color_green',.FALSE., &
                                                              .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filename3,1,3,'color_blue',.FALSE., &
                                                              .TRUE.,.FALSE.)
!
!  put as a comment the possibility to plot also the experimental data
!
   CALL gnuplot_write_file_mul_data(filename3,1,3,'color_blue',.FALSE.,&
                                                                .FALSE.,.TRUE.)
   CALL gnuplot_write_file_mul_point('anhar.exp',1,2,'color_red',.FALSE.,&
                                                                 .TRUE.,.TRUE.)

ELSEIF (ibrav_save==4.OR.ibrav_save==5.OR.ibrav_save==6.OR.ibrav_save==7) THEN
   CALL gnuplot_write_file_mul_data(filename2,1,4,'color_red',.TRUE.,.FALSE.,&
                                                                .FALSE.)
   CALL gnuplot_write_file_mul_data(filename3,1,4,'color_blue',.FALSE.,.FALSE.,&
                                                                .FALSE.)
   IF (done_grun) &
      CALL gnuplot_write_file_mul_data(filename4,1,4,'color_cyan',.FALSE.,&
                                              .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filename2,1,5,'color_pink',.FALSE.,.FALSE.,&
                                                                .FALSE.)
   IF (done_grun) THEN
      CALL gnuplot_write_file_mul_data(filename3,1,5,'color_green',.FALSE.,&
                                                             .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename4,1,5,'color_orange',.FALSE.,&
                                                              .TRUE.,.FALSE.)
   ELSE
      CALL gnuplot_write_file_mul_data(filename3,1,5,'color_green',.FALSE.,&
                                                         .TRUE.,.FALSE.)
   END IF
ELSEIF (ibrav_save==8.OR.ibrav_save==9.OR.ibrav_save==10.OR.ibrav_save==11) THEN
   CALL gnuplot_write_file_mul_data(filename2,1,5,'color_red',.TRUE.,.FALSE.,&
                                                                .FALSE.)
   IF (done_grun) &
      CALL gnuplot_write_file_mul_data(filename4,1,5,'color_cyan',.FALSE.,&
                                              .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filename2,1,6,'color_green',.FALSE.,&
                                                    .FALSE.,.FALSE.)
   IF (done_grun) &
      CALL gnuplot_write_file_mul_data(filename4,1,6,'color_orange',.FALSE.,&
                                              .FALSE.,.FALSE.)
   IF (done_grun) THEN
      CALL gnuplot_write_file_mul_data(filename2,1,7,'color_blue',.FALSE., &
                                                            .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename4,1,7,'color_ligth_blue',&
                                                     .FALSE.,.TRUE.,.FALSE.)
   ELSE
      CALL gnuplot_write_file_mul_data(filename2,1,7,'color_blue',.FALSE.,&
                                                      .TRUE.,.FALSE.)
   ENDIF
END IF

CALL gnuplot_ylabel('Volume thermal expansion {/Symbol b} x 10^6 (K^{-1})',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename,1,3,'color_red',.TRUE.,.FALSE.,&
                                                               .FALSE.)
CALL gnuplot_write_file_mul_data(filename1,1,3,'color_blue',.FALSE.,.TRUE.,&
                                                               .FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_anhar_anis
