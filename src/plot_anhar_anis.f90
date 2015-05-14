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
!  temperature (one parameter per plot), the thermal expansions (all in the
!  same plot).
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
                            gnuplot_write_file_mul_data_sum, &
                            gnuplot_write_file_mul_point,  &
                            gnuplot_write_horizontal_line, &
                            gnuplot_set_fact
USE data_files,      ONLY : flanhar
USE control_pwrun,   ONLY : ibrav_save
USE temperature,     ONLY : tmin, tmax
USE constants,       ONLY : ry_kbar
USE mp_images,       ONLY : my_image_id, root_image
USE io_global,       ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename1, filenameps
                      
CHARACTER(LEN=6), EXTERNAL :: int_to_char
CHARACTER(LEN=8), EXTERNAL :: float_to_char
INTEGER :: system
INTEGER :: ierr

IF ( my_image_id /= root_image ) RETURN

gnu_filename=TRIM(flgnuplot)//'_anhar'

CALL gnuplot_start(gnu_filename)

filenameps=TRIM(flpsanhar)
IF (tmin /= 1.0_DP) THEN
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ELSE
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ENDIF

filename=TRIM(flanhar)
filename1=TRIM(flanhar)//'.celldm'

CALL gnuplot_xlabel('T (K)',.FALSE.) 
CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('a (a.u.)',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename1,1,2,'color_red',.TRUE., &
                                                      .TRUE.,.FALSE.)

IF (ibrav_save==4.OR.ibrav_save==5.OR.ibrav_save==6.OR.ibrav_save==7) THEN
   CALL gnuplot_ylabel('c/a ',.FALSE.) 
   CALL gnuplot_write_file_mul_data(filename1,1,3,'color_blue',.TRUE., &
                                                              .TRUE.,.FALSE.)
ELSEIF (ibrav_save==8.OR.ibrav_save==9.OR.ibrav_save==10.OR.ibrav_save==11) THEN
   CALL gnuplot_ylabel('b/a ',.FALSE.) 
   CALL gnuplot_write_file_mul_data(filename1,1,3,'color_green',.TRUE., &
                                                              .FALSE.,.FALSE.)
   CALL gnuplot_ylabel('c/a ',.FALSE.) 
   CALL gnuplot_write_file_mul_data(filename1,1,4,'color_green',.FALSE., &
                                                              .TRUE.,.FALSE.)

ENDIF
CALL gnuplot_ylabel('Linear thermal expansion {/Symbol a} x 10^6 (K^{-1})',.FALSE.) 
IF (ibrav_save==4.OR.ibrav_save==5.OR.ibrav_save==6.OR.ibrav_save==7) THEN
   CALL gnuplot_write_file_mul_data(filename1,1,4,'color_red',.TRUE.,.FALSE.,&
                                                                .FALSE.)
   CALL gnuplot_write_file_mul_data(filename1,1,5,'color_green',.FALSE.,.TRUE.,&
                                                                .FALSE.)
ELSEIF (ibrav_save==8.OR.ibrav_save==9.OR.ibrav_save==10.OR.ibrav_save==11) THEN
   CALL gnuplot_write_file_mul_data(filename1,1,5,'color_red',.TRUE.,.FALSE.,&
                                                                .FALSE.)
   CALL gnuplot_write_file_mul_data(filename1,1,6,'color_green',.FALSE.,.FALSE.,&
                                                                .FALSE.)
   CALL gnuplot_write_file_mul_data(filename1,1,7,'color_blue',.FALSE.,.TRUE.,&
                                                                .FALSE.)
END IF

CALL gnuplot_ylabel('Volume ((a.u.)^3)',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename,1,2,'color_red',.TRUE.,.TRUE.,&
                                                               .FALSE.)
CALL gnuplot_ylabel('Volume thermal expansion {/Symbol b} x 10^6 (K^{-1})',&
                                                               .FALSE.) 
CALL gnuplot_write_file_mul_data(filename,1,3,'color_red',.TRUE.,.TRUE.,&
                                                               .FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_anhar_anis
