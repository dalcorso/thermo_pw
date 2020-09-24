!
! Copyright (C) 2020 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
SUBROUTINE plot_optical_constants()
!---------------------------------------------------------------
!
!  This is a driver to plot the quantities inside the optical constant 
!  file. For the general crystal it produces four plots.
!  It plots the real and imaginary parts of the dielectric constant
!  the real and imaginary part of the refractive index
!  ibrav = 1,2,3  it plots only epsilon_xx 
!  ibrav = 4,5,6,7  it plots only epsilon_xx (red) and epsilon_zz (green)
!  ibrav = 8,9,10,11 it plots epsilon_xx (red), epsilon_yy (green), 
!                   epsilon_zz (blue).
!  other ibrav are presently not programmed
!  For the cubic crystals it produces two additional plots that contain
!  the reflectivity for normal incidence and 
!  the logarithm with basis 10 or the absorption coefficient in cm^-1.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : rytoev
USE control_gnuplot,  ONLY : flgnuplot, lgnuplot, gnuplot_command, flext
USE postscript_files, ONLY : flpsoptical
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,   &
                             gnuplot_write_header,         &
                             gnuplot_ylabel,               &
                             gnuplot_xlabel,               &
                             gnuplot_write_command,        &
                             gnuplot_write_file_mul_data,  &
                             gnuplot_write_file_mul_data_log10
USE data_files,       ONLY : floptical
USE cell_base,        ONLY : ibrav
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename1, filename2
INTEGER :: ierr, system, col

IF ( my_image_id /= root_image ) RETURN

DO col=1, 2
   IF (col==1) THEN
      gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_optical'
      filename=TRIM(flpsoptical)//TRIM(flext)
   ELSE
      gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_optical_lambda'
      filename=TRIM(flpsoptical)//'_lambda'//TRIM(flext)
   ENDIF
   CALL gnuplot_start(gnu_filename)

   IF (col==1) THEN
      CALL gnuplot_write_header(filename, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
                                                       rytoev, flext ) 
      CALL gnuplot_xlabel('{/Symbol w}  (eV)',.FALSE.) 
   ELSE
      CALL gnuplot_write_header(filename, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext ) 
      CALL gnuplot_xlabel('{/Symbol l}  (nm)',.FALSE.) 
   ENDIF

   IF (ibrav==1.OR.ibrav==2.OR.ibrav==3) THEN
      CALL gnuplot_ylabel('{/Symbol e}_1 ({/Symbol w})',.FALSE.) 
      filename='dynamical_matrices/'//TRIM(floptical)//'_xx'
      CALL gnuplot_write_file_mul_data(filename,col,3,'color_red',.TRUE.,&
                                                .TRUE.,.FALSE.)

      CALL gnuplot_ylabel('{/Symbol e}_2 ({/Symbol w})',.FALSE.) 
      CALL gnuplot_write_file_mul_data(filename,col,4,'color_red',.TRUE.,&
                                                .TRUE.,.FALSE.)

      CALL gnuplot_ylabel('n ({/Symbol w})',.FALSE.) 
      CALL gnuplot_write_file_mul_data(filename,col,5,'color_red',.TRUE.,&
                                                .TRUE.,.FALSE.)

      CALL gnuplot_ylabel('k ({/Symbol w})',.FALSE.) 
      CALL gnuplot_write_file_mul_data(filename,col,6,'color_red',.TRUE.,&
                                                .TRUE.,.FALSE.)

      CALL gnuplot_ylabel('Reflectivity ({/Symbol w})',.FALSE.) 
      CALL gnuplot_write_command('set yrange [0:1]',.FALSE.)
      CALL gnuplot_write_file_mul_data(filename,col,7,'color_red',.TRUE.,&
                                                .TRUE.,.FALSE.)

      CALL gnuplot_write_command('unset yrange',.FALSE.)
      CALL gnuplot_ylabel('log_{10} ({/Symbol a} ({/Symbol w}) (cm^{-1}))', &
                                                                 .FALSE.) 
      filename='dynamical_matrices/'//TRIM(floptical)//'_xx'
      CALL gnuplot_write_file_mul_data_log10(filename,col,8,'color_red',    &
                                                   .TRUE.,.TRUE.,.FALSE.)
   ELSEIF (ibrav==4.OR.ibrav==5.OR.ibrav==6.OR.ibrav==7) THEN
      CALL gnuplot_ylabel('{/Symbol e}_1 ({/Symbol w})',.FALSE.) 
      filename='dynamical_matrices/'//TRIM(floptical)//'_xx'
      filename1='dynamical_matrices/'//TRIM(floptical)//'_zz'
      CALL gnuplot_write_file_mul_data(filename,col,3,'color_red',.TRUE.,&
                                                   .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename1,col,3,'color_green',.FALSE.,&
                                                   .TRUE.,.FALSE.)

      CALL gnuplot_ylabel('{/Symbol e}_2 ({/Symbol w})',.FALSE.) 
      CALL gnuplot_write_file_mul_data(filename,col,4,'color_red',.TRUE.,&
                                                .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename1,col,4,'color_green',.FALSE.,&
                                                .TRUE.,.FALSE.)

      CALL gnuplot_ylabel('n ({/Symbol w})',.FALSE.) 
      CALL gnuplot_write_file_mul_data(filename,col,5,'color_red',.TRUE.,&
                                                .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename1,col,5,'color_green',.FALSE.,&
                                                .TRUE.,.FALSE.)

      CALL gnuplot_ylabel('k ({/Symbol w})',.FALSE.) 
      CALL gnuplot_write_file_mul_data(filename,col,6,'color_red',.TRUE.,&
                                                .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename1,col,6,'color_green',.FALSE.,&
                                                .TRUE.,.FALSE.)
   ELSEIF (ibrav==8.OR.ibrav==9.OR.ibrav==10.OR.ibrav==11) THEN
      CALL gnuplot_ylabel('{/Symbol e}_1 ({/Symbol w})',.FALSE.) 
      filename='dynamical_matrices/'//TRIM(floptical)//'_xx'
      filename1='dynamical_matrices/'//TRIM(floptical)//'_yy'
      filename2='dynamical_matrices/'//TRIM(floptical)//'_zz'
      CALL gnuplot_write_file_mul_data(filename,col,3,'color_red',.TRUE.,&
                                                   .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename1,col,3,'color_green',.FALSE.,&
                                                   .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename2,col,3,'color_blue',.FALSE.,&
                                                    .TRUE.,.FALSE.)
      CALL gnuplot_ylabel('{/Symbol e}_2 ({/Symbol w})',.FALSE.) 
      CALL gnuplot_write_file_mul_data(filename,col,4,'color_red',.TRUE.,&
                                                   .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename1,col,4,'color_green',.FALSE.,&
                                                   .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename2,col,4,'color_blue',.FALSE.,&
                                                   .TRUE.,.FALSE.)

      CALL gnuplot_ylabel('n ({/Symbol w})',.FALSE.) 
      CALL gnuplot_write_file_mul_data(filename,col,5,'color_red',.TRUE.,&
                                                  .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename1,col,5,'color_green',.FALSE.,&
                                                  .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename2,col,5,'color_blue',.FALSE.,&
                                                  .TRUE.,.FALSE.)

      CALL gnuplot_ylabel('k ({/Symbol w})',.FALSE.) 
      CALL gnuplot_write_file_mul_data(filename,col,6,'color_red',.TRUE.,&
                                                .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename1,col,6,'color_green',.FALSE.,&
                                                .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename2,col,6,'color_blue',.FALSE.,&
                                                .TRUE.,.FALSE.)
   ENDIF

   CALL gnuplot_end()

   IF (lgnuplot.AND.ionode) &
      ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

   !IF (lgnuplot.AND.ionode) &
   !   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
   !                                       //TRIM(gnu_filename), WAIT=.FALSE.)
ENDDO

RETURN
END SUBROUTINE plot_optical_constants
