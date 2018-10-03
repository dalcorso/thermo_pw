!
! Copyright (C) 2015-2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_anhar_anis()
!
!  This is a driver to plot the quantities written inside flanhar for the
!  case of anisotropic systems. Presently plots celldm (one parameter &
!  per plot), the thermal expansion tensor !  components (all in the same 
!  plot), the volume, and the volume thermal expansion as a function 
!  of temperature.
!  
!
USE kinds,           ONLY : DP
USE control_gnuplot, ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE postscript_files, ONLY : flpsanhar
USE gnuplot,         ONLY : gnuplot_start, gnuplot_end,  &
                            gnuplot_write_header,        &
                            gnuplot_ylabel,              &
                            gnuplot_xlabel,              &
                            gnuplot_write_file_mul_data, &
                            gnuplot_write_file_mul_point, &
                            gnuplot_write_horizontal_line, &
                            gnuplot_set_fact
USE data_files,  ONLY : flanhar
USE grun_anharmonic, ONLY : done_grun
USE control_grun,  ONLY : lb0_t
USE initial_conf,  ONLY : ibrav_save
USE control_thermo,  ONLY : ltherm_dos, ltherm_freq, with_eigen
USE control_elastic_constants, ONLY : el_cons_t_available
USE anharmonic, ONLY : lelastic
USE ph_freq_anharmonic, ONLY : lelasticf
USE temperature,     ONLY : tmin, tmax
USE control_pressure, ONLY : pressure_kb
USE mp_images,       ONLY : my_image_id, root_image
USE io_global,       ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename1, filename2, &
                      filename3, filename4, filename7, filename8,   &
                      filename9, filename_bulk, filename_bulk_ph,   &
                      filename_heat, filename_heat_ph, filenameps
LOGICAL :: lgrun
INTEGER :: ierr, system

IF ( my_image_id /= root_image ) RETURN

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_anhar'
CALL add_pressure(gnu_filename)

CALL gnuplot_start(gnu_filename)

filenameps=TRIM(flpsanhar)
CALL add_pressure(filenameps)
filenameps=TRIM(filenameps)//TRIM(flext)
IF (tmin /= 1.0_DP) THEN
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                            flext ) 
ELSE
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                            flext ) 
ENDIF

filename='anhar_files/'//TRIM(flanhar)
CALL add_pressure(filename)
filename1='anhar_files/'//TRIM(flanhar)//'_ph'
CALL add_pressure(filename1)
filename2='anhar_files/'//TRIM(flanhar)//'.celldm'
CALL add_pressure(filename2)
filename3='anhar_files/'//TRIM(flanhar)//'.celldm_ph'
CALL add_pressure(filename3)
filename4='anhar_files/'//TRIM(flanhar)//'.celldm_grun'
CALL add_pressure(filename4)
filename_bulk='anhar_files/'//TRIM(flanhar)//'.bulk_mod'
CALL add_pressure(filename_bulk)
filename_bulk_ph='anhar_files/'//TRIM(flanhar)//'.bulk_mod_ph'
CALL add_pressure(filename_bulk_ph)
filename_heat='anhar_files/'//TRIM(flanhar)//'.heat'
CALL add_pressure(filename_heat)
filename_heat_ph='anhar_files/'//TRIM(flanhar)//'.heat_ph'
CALL add_pressure(filename_heat_ph)
filename7='anhar_files/'//TRIM(flanhar)//'.aux_grun'
CALL add_pressure(filename7)
filename8='anhar_files/'//TRIM(flanhar)//'.anis'
CALL add_pressure(filename8)
filename9='anhar_files/'//TRIM(flanhar)//'.anis_ph'
CALL add_pressure(filename9)

lgrun = lelastic .OR. lelasticf

CALL gnuplot_xlabel('T (K)',.FALSE.) 
CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('a (a.u.)',.FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename2,1,2,'color_red',.TRUE., &
                                                .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename3,1,2,'color_blue',&
                                             .NOT.ltherm_dos, .TRUE., .FALSE.)

IF (ibrav_save==4.OR.ibrav_save==5.OR.ibrav_save==6.OR.ibrav_save==7) THEN
   IF (ibrav_save==5) THEN
      CALL gnuplot_ylabel('cos({/Symbol a})',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('c/a ',.FALSE.) 
   ENDIF
   IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename2,1,3,'color_red',.TRUE., &
                                                     .NOT.ltherm_freq,.FALSE.)
   IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename3,1,3,'color_blue',&
                               .NOT.ltherm_dos, .TRUE.,.FALSE.)
ELSEIF (ibrav_save==8.OR.ibrav_save==9.OR.ibrav_save==10.OR.ibrav_save==11) THEN
   CALL gnuplot_ylabel('b/a ',.FALSE.) 
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename2,1,3,'color_red',.TRUE., &
                                                   .NOT.ltherm_freq,.FALSE.)
   IF (ltherm_freq) & 
      CALL gnuplot_write_file_mul_data(filename3,1,3,'color_blue',&
                                             .NOT.ltherm_dos,.TRUE.,.FALSE.)
   CALL gnuplot_ylabel('c/a ',.FALSE.) 
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename2,1,4,'color_red',.TRUE., &
                                             .NOT.ltherm_freq,.FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data(filename3,1,4,'color_blue', &
                                             .NOT.ltherm_dos, .TRUE.,.FALSE.)
ENDIF
CALL gnuplot_ylabel('Volume ((a.u.)^3)',.FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename,1,2,'color_red',.TRUE., &
                                                  .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename1,1,2,'color_blue',&
                                             .NOT.ltherm_dos, .TRUE., .FALSE.)

IF (pressure_kb /= 0.0_DP) THEN
   CALL gnuplot_ylabel('Gibbs free energy (Ry)',.FALSE.) 
ELSE
   CALL gnuplot_ylabel('Helmholtz free energy (Ry)',.FALSE.) 
ENDIF
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename,1,3,'color_red',.TRUE., &
                                                  .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename1,1,3,'color_blue',&
                                             .NOT.ltherm_dos, .TRUE., .FALSE.)


CALL gnuplot_ylabel('Linear thermal expansion {/Symbol a} x 10^6 (K^{-1})',.FALSE.) 
IF (ibrav_save==1.OR.ibrav_save==2.OR.ibrav_save==3) THEN
   IF (done_grun) &
      CALL gnuplot_write_file_mul_data(filename4,1,3,'color_green',.TRUE., &
                              .NOT.(ltherm_dos.OR.ltherm_freq),.FALSE.)
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename2,1,3,'color_red',&
                                .NOT.done_grun,.NOT.ltherm_freq,.FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data(filename3,1,3,'color_blue', &
                         .NOT.(ltherm_dos.OR.done_grun), .TRUE.,.FALSE.)
!
!  put as a comment the possibility to plot also the experimental data
!
   CALL gnuplot_write_file_mul_data(filename3,1,3,'color_blue', &
                         .NOT.(ltherm_dos.OR.done_grun), .FALSE.,.TRUE.)
   CALL gnuplot_write_file_mul_point('anhar.exp',1,2,'color_red',.FALSE.,&
                                                                 .TRUE.,.TRUE.)

ELSEIF (ibrav_save==4.OR.ibrav_save==5.OR.ibrav_save==6.OR.ibrav_save==7) THEN
   IF (done_grun) THEN
      CALL gnuplot_write_file_mul_data(filename4,1,4,'color_green',.TRUE.,&
                                              .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename4,1,5,'color_orange',.FALSE.,&
                                 .NOT.(ltherm_dos.OR.ltherm_freq),.FALSE.)
   ENDIF
   IF (ltherm_dos) THEN
      CALL gnuplot_write_file_mul_data(filename2,1,4,'color_red', &
                                     .NOT.done_grun,.FALSE., .FALSE.)
      CALL gnuplot_write_file_mul_data(filename2,1,5,'color_pink',.FALSE., &
                                     .NOT.ltherm_freq,.FALSE.)
   ENDIF
   IF (ltherm_freq) THEN
      CALL gnuplot_write_file_mul_data(filename3,1,4,'color_blue', &
             .NOT.(done_grun.OR.ltherm_dos),.FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename3,1,5,'color_cyan',.FALSE.,&
                                                         .TRUE.,.FALSE.)
   END IF
ELSEIF (ibrav_save==8.OR.ibrav_save==9.OR.ibrav_save==10.OR.ibrav_save==11) THEN
   IF (done_grun) THEN
      CALL gnuplot_write_file_mul_data(filename4,1,5,'color_green',.TRUE.,&
                                              .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename4,1,6,'color_light_blue',&
                                              .FALSE.,.FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename4,1,7,'color_orange',&
                        .FALSE.,.NOT.(ltherm_dos.OR.ltherm_freq),.FALSE.)
   ENDIF
   IF (ltherm_dos) THEN 
      CALL gnuplot_write_file_mul_data(filename2,1,5,'color_red', &
                                      .NOT.done_grun,.FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename2,1,6,'color_gold',.FALSE.,&
                                                    .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename2,1,7,'color_pink',.FALSE.,&
                                               .NOT.ltherm_freq,.FALSE.)
   ENDIF
   IF (ltherm_freq) THEN 
      CALL gnuplot_write_file_mul_data(filename3,1,5,'color_blue', &
                           .NOT.(done_grun.OR.ltherm_dos),.FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename3,1,6,'color_olive',.FALSE.,&
                                                    .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename3,1,7,'color_cyan',.FALSE.,&
                                               .TRUE.,.FALSE.)
   ENDIF
END IF

CALL gnuplot_ylabel('Volume thermal expansion {/Symbol b} x 10^6 (K^{-1})',.FALSE.) 

IF (ltherm_dos) THEN
   CALL gnuplot_write_file_mul_data(filename,1,4,'color_red',.TRUE., &
                           .NOT.(ltherm_freq.OR.lgrun), .FALSE.)
END IF
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename1,1,4,'color_blue', &
                                    .NOT.ltherm_dos,.NOT.lgrun, .FALSE.)

IF (lgrun) THEN
      CALL gnuplot_write_file_mul_data(filename7,1,2,'color_green', &
                  .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE., .FALSE.)

   CALL gnuplot_set_fact(1313313.0_DP,.FALSE.)
   CALL gnuplot_ylabel('Heat capacity C_v (J / K / N / mol)',.FALSE.)
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename_heat,1,2,'color_red',.TRUE.,&
                                            .NOT.ltherm_freq,.FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data(filename_heat_ph,1,2,'color_blue',&
                                            .NOT.ltherm_dos,.TRUE.,.FALSE.)

   CALL gnuplot_set_fact(1313313.0_DP,.FALSE.)
   CALL gnuplot_ylabel('C_p - C_v (J / K / N / mol)',.FALSE.)

   IF (ltherm_dos) THEN
      CALL gnuplot_write_file_mul_data(filename_heat,1,4,'color_red',.TRUE., &
                                  .FALSE.,.FALSE.)
      IF (lgrun) &
      CALL gnuplot_write_file_mul_data(filename8,1,2,'color_gold',.FALSE., &
                                  .FALSE.,.FALSE.)
   ENDIF   
   IF (ltherm_freq) THEN
      CALL gnuplot_write_file_mul_data(filename_heat_ph,1,4,'color_blue',&
                                  .NOT.ltherm_dos,.NOT.lgrun,.FALSE.)

      IF (lgrun) &
         CALL gnuplot_write_file_mul_data(filename9,1,2,'color_orange',&
                                                .FALSE., .FALSE.,.FALSE.)

   ENDIF
   IF (lgrun) &
      CALL gnuplot_write_file_mul_data(filename7,1,4,'color_green',&
                    .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)

   CALL gnuplot_set_fact(1.0_DP,.FALSE.)
   CALL gnuplot_ylabel('Gr\374neisen parameter ({/Symbol g})',.FALSE.)
   CALL gnuplot_write_horizontal_line(0.0_DP, 2, 'front', 'color_black', &
                                                                  .FALSE.)
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename_bulk,1,2,'color_red',.TRUE.,&
                                                        .FALSE.,.FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data(filename_bulk_ph,1,2,'color_blue',&
                                        .NOT.ltherm_dos,.NOT.lgrun,.FALSE.)

      CALL gnuplot_write_file_mul_data(filename7,1,3,'color_green',&
                          .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)
END IF

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!  CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_anhar_anis
