!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_anhar()
!
!  This is a driver to plot the quantities written inside flanhar,
!  flanhar//'_ph', flanhar//'.aux', flanhar//'.aux_ph' and 
!  flanhar//'.aux_grun'
!  
!
USE kinds,            ONLY : DP
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE control_thermo,   ONLY : ltherm, ltherm_dos, ltherm_freq, with_eigen
USE postscript_files, ONLY : flpsanhar
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_write_file_mul_data, &
                             gnuplot_write_file_mul_point,  &
                             gnuplot_write_horizontal_line, &
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar
USE temperature,      ONLY : tmin, tmax
USE control_pressure, ONLY : pressure_kb
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename0, filename3,  &
                      filename_bulk, filename_bulk_ph,               &
                      filename_dbulk, filename_dbulk_ph,             &
                      filename_heat, filename_heat_ph 
CHARACTER(LEN=8) :: float_to_char
INTEGER :: ierr, system

IF ( my_image_id /= root_image ) RETURN

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_anhar'
CALL add_pressure(gnu_filename)

CALL gnuplot_start(gnu_filename)

filename=TRIM(flpsanhar)
CALL add_pressure(filename)
filename=TRIM(filename)//TRIM(flext)
IF (tmin /= 1.0_DP) THEN
   CALL gnuplot_write_header(filename, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext ) 
ELSE
   CALL gnuplot_write_header(filename, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext ) 
ENDIF

filename="anhar_files/"//TRIM(flanhar)//'_ph'
CALL add_pressure(filename)
filename0="anhar_files/"//TRIM(flanhar)
CALL add_pressure(filename0)
filename3="anhar_files/"//TRIM(flanhar)//'.aux_grun'
CALL add_pressure(filename3)
filename_bulk="anhar_files/"//TRIM(flanhar)//'.bulk_mod'
CALL add_pressure(filename_bulk)
filename_bulk_ph="anhar_files/"//TRIM(flanhar)//'.bulk_mod_ph'
CALL add_pressure(filename_bulk_ph)
filename_dbulk="anhar_files/"//TRIM(flanhar)//'.dbulk_mod'
CALL add_pressure(filename_dbulk)
filename_dbulk_ph="anhar_files/"//TRIM(flanhar)//'.dbulk_mod_ph'
CALL add_pressure(filename_dbulk_ph)
filename_heat="anhar_files/"//TRIM(flanhar)//'.heat'
CALL add_pressure(filename_heat)
filename_heat_ph="anhar_files/"//TRIM(flanhar)//'.heat_ph'
CALL add_pressure(filename_heat_ph)

CALL gnuplot_xlabel('T (K)',.FALSE.) 
CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('Volume ((a.u.)^3)',.FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename0,1,2,'color_red',.TRUE.,&
                                              .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,2,'color_blue',&
                                        .NOT.ltherm_dos,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1.0_DP,.FALSE.)

IF (pressure_kb /= 0.0_DP) THEN
   CALL gnuplot_ylabel('Gibbs Free Energy (Ry)',.FALSE.) 
ELSE
   CALL gnuplot_ylabel('Helmholtz Free Energy (Ry)',.FALSE.) 
ENDIF

IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename0,1,3,'color_red',.TRUE.,&
                                              .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,3,'color_blue',&
                                        .NOT.ltherm_dos,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('Bulk modulus (kbar)',.FALSE.) 

IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename_bulk,1,3,'color_red',.TRUE., &
                                                .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename_bulk_ph,1,3,'color_blue',&
                                                .NOT.ltherm_dos,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('d B / d p',.FALSE.) 

IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename_dbulk,1,2,'color_red',.TRUE., &
                                                .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename_dbulk_ph,1,2,'color_blue',&
                                               .NOT.ltherm_dos,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('Thermal expansion ({/Symbol b} x 10^{6}) (K^{-1})',.FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename0,1,4,'color_red',.TRUE.,&
                                                            .FALSE.,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',&
                                            .NOT.ltherm_dos,.FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filename3,1,2,'color_green',&
                                .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)
!
!  put as a comment the possibility to plot also the experimental data
!
CALL gnuplot_write_file_mul_data(filename3,1,2,'color_green',&
                                               .NOT.ltherm,.FALSE.,.TRUE.)
CALL gnuplot_write_file_mul_point('anhar.exp',1,2,'color_red',.FALSE.,&
                                                            .TRUE.,.TRUE.)

CALL gnuplot_set_fact(1313313.0_DP,.FALSE.)
CALL gnuplot_ylabel('Heat capacity C_v (J / K / N / mol)',.FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename_heat,1,2,'color_red',.TRUE.,&
                                            .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename_heat_ph,1,2,'color_blue',&
                                            .NOT.ltherm_dos,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1313313.0_DP,.FALSE.)
CALL gnuplot_ylabel('Heat capacity C_p (J / K / N / mol)',.FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename_heat,1,3,'color_red',.TRUE.,&
                                          .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) THEN
   CALL gnuplot_write_file_mul_data(filename_heat_ph,1,3,'color_blue', &
                                    .NOT.ltherm_dos,.TRUE.,.FALSE.)
!
!  put as a comment the possibility to plot also the experimental data
!
   CALL gnuplot_write_file_mul_data(filename_heat,1,3,'color_blue',&
                                       .NOT.ltherm_dos,.FALSE.,.TRUE.)
   CALL gnuplot_write_file_mul_point('cv.exp',1,2,'color_red',.FALSE.,&
                                                              .TRUE.,.TRUE.)
ENDIF

CALL gnuplot_set_fact(1313313.0_DP,.FALSE.)
CALL gnuplot_ylabel('C_p - C_v (J / K / N / mol)',.FALSE.) 

IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename_heat,1,4,'color_red',.TRUE., &
                                                   .FALSE.,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename_heat_ph,1,4,'color_blue',&
                                               .NOT.ltherm_dos,.FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filename3,1,4,'color_green', &
                              .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1._DP,.FALSE.)
CALL gnuplot_ylabel('B_S - B_T (kbar)',.FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename_bulk,1,4,'color_red',.TRUE.,&
                                             .FALSE.,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename_bulk_ph,1,4,'color_blue',&
                                           .NOT.ltherm_dos,.FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filename3,1,5,'color_green', &
                              .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('Gr\374neisen parameter ({/Symbol g})',.FALSE.) 
CALL gnuplot_write_horizontal_line(0.0_DP, 2, 'front', 'color_black', .FALSE.)
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename_bulk,1,2,'color_red',.TRUE.,&
                                                     .FALSE.,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename_bulk_ph,1,2,'color_blue',&
                         .NOT.ltherm_dos,.FALSE.,.FALSE.)

   CALL gnuplot_write_file_mul_data(filename3,1,3,'color_green', &
                     .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)
IF (with_eigen) CALL plot_dw_anhar()

RETURN
END SUBROUTINE plot_anhar

! Copyright (C) 2018 Cristiano Malica

SUBROUTINE plot_dw_anhar()
!
!  This is a driver to plot the quantities written inside 
!
USE kinds,            ONLY : DP
USE ions_base,        ONLY : nat
USE cell_base,        ONLY : ibrav
USE control_thermo,   ONLY : ltherm_dos, ltherm_freq
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE postscript_files, ONLY : flpsanhar
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                             gnuplot_ylabel, &
                             gnuplot_xlabel, &
                             gnuplot_write_file_mul_data
USE gnuplot_color,    ONLY : gnuplot_set_greens
USE data_files,       ONLY : flanhar
USE temperature,      ONLY : tmin, tmax
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: gnu_filename, filename, psfilename, filetherm
INTEGER :: ierr, system, na
CHARACTER(LEN=6) :: int_to_char

IF ( my_image_id /= root_image ) RETURN

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_anhar_dw'
CALL gnuplot_start(gnu_filename)

psfilename=TRIM(flpsanhar)//'_dw'//TRIM(flext)
IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(psfilename, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ELSE
   CALL gnuplot_write_header(psfilename, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ENDIF
CALL gnuplot_set_greens()
CALL gnuplot_xlabel('T (K)', .FALSE.)
DO na=1,nat

   filename='anhar_files/'//TRIM(flanhar)//'_ph.'//TRIM(int_to_char(na))//'.dw'
   filetherm='anhar_files/'//TRIM(flanhar)//'.'//TRIM(int_to_char(na))//'.dw'

!
!   First the diagonal components
!
   CALL gnuplot_ylabel('B_{ii} ({\305}^2) (atom '// &
                                        TRIM(int_to_char(na))//')',.FALSE.)

   IF (ltherm_dos) THEN
      CALL gnuplot_write_file_mul_data(filetherm,1,2,'color_red',.TRUE., &
                                                   .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filetherm,1,5,'color_blue',.FALSE., &
                                                   .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filetherm,1,7,'color_dark_spring_green',&
                                   .FALSE., .NOT.ltherm_freq,.FALSE.)
   ENDIF

   IF (ltherm_freq) THEN
      CALL gnuplot_write_file_mul_data(filename,1,2,'color_pink', &
                                               .NOT.ltherm_dos,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filename,1,5,'color_light_blue', &
                                               .FALSE., .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filename,1,7,'color_green', &
                                               .FALSE.,.TRUE.,.FALSE.)
   ENDIF
!
!  And then the off diagonal, only for noncubic solids
!
   IF (ibrav/=1.AND.ibrav/=2.AND.ibrav/=3) THEN

      IF (ltherm_dos) THEN
         CALL gnuplot_write_file_mul_data(filetherm,1,3,'color_red',.TRUE., &
                                                   .FALSE.,.FALSE.)
         CALL gnuplot_write_file_mul_data(filetherm,1,4,'color_blue',.FALSE., &
                                                   .FALSE.,.FALSE.)
         CALL gnuplot_write_file_mul_data(filetherm,1,6,&
                      'color_dark_spring_green',.FALSE., &
                                                 .NOT.ltherm_freq,.FALSE.)
      ENDIF

      IF (ltherm_freq) THEN
         CALL gnuplot_write_file_mul_data(filename,1,3,'color_pink', &
                                               .NOT.ltherm_dos,.FALSE.,.FALSE.)

         CALL gnuplot_write_file_mul_data(filename,1,4,'color_light_blue', &
                                               .FALSE.,.FALSE.,.FALSE.)

         CALL gnuplot_write_file_mul_data(filename,1,6,'color_green', &
                                               .FALSE.,.TRUE.,.FALSE.)
      ENDIF

   ENDIF

ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_dw_anhar
