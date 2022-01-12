!
! Copyright (C) 2015-2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------
SUBROUTINE plot_anhar_anis()
!-----------------------------------------------------------------
!
!  This is a driver to plot the quantities written inside flanhar,
!  flanhar//'_ph', flanhar//'.celldm', flanhar//'.celldm_ph', 
!  flanhar//'.bulk_mod', flanhar//'.bulk_mod_ph', 
!  flanhar//'.heat', flanhar//'.heat_ph', and flanhar//'.aux_grun'.
!  flanhar//'.heat_anis', flanhar//'.heat_anis_ph', flanhar//'.heat_anis_grun'
!  flanhar//'.gamma', flanhar//'.gamma_ph', flanhar//'.gamma_grun'
!
USE kinds,           ONLY : DP
USE constants,       ONLY : rydberg_si, avogadro
USE control_gnuplot, ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE postscript_files, ONLY : flpsanhar
USE gnuplot,         ONLY : gnuplot_start, gnuplot_end,  &
                            gnuplot_write_header,        &
                            gnuplot_ylabel,              &
                            gnuplot_xlabel,              &
                            gnuplot_write_file_mul_data, &
                            gnuplot_write_file_mul_data_times, &
                            gnuplot_write_file_mul_point, &
                            gnuplot_write_horizontal_line, &
                            gnuplot_set_fact
USE data_files,  ONLY : flanhar
USE grun_anharmonic, ONLY : done_grun
USE initial_conf,  ONLY : ibrav_save
USE control_thermo,  ONLY : ltherm_dos, ltherm_freq
USE control_elastic_constants, ONLY : lelastic, lelasticf
USE temperature,     ONLY : tmin, tmax
USE control_pressure, ONLY : pressure_kb
USE mp_images,       ONLY : my_image_id, root_image
USE io_global,       ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename1, filename2, &
                      filename3, filename4, filename_aux_grun,    &
                      filename_bulk, filename_bulk_ph,   &
                      filename_heat, filename_heat_ph, filename_heat_anis,  &
                      filename_heat_anis_ph, filename_heat_anis_grun,       &
                      filename_gamma, filename_gamma_ph, filename_gamma_grun, &
                      filenameps
INTEGER :: ierr, system
REAL(DP) :: factor
LOGICAL :: isoent_avail, noncubic

IF ( my_image_id /= root_image ) RETURN

isoent_avail=lelastic.OR.lelasticf
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
filename_aux_grun='anhar_files/'//TRIM(flanhar)//'.aux_grun'
CALL add_pressure(filename_aux_grun)
filename_heat_anis='anhar_files/'//TRIM(flanhar)//'.heat_anis'
CALL add_pressure(filename_heat_anis)
filename_heat_anis_ph='anhar_files/'//TRIM(flanhar)//'.heat_anis_ph'
CALL add_pressure(filename_heat_anis_ph)
filename_heat_anis_grun='anhar_files/'//TRIM(flanhar)//'.heat_anis_grun'
CALL add_pressure(filename_heat_anis_grun)
filename_gamma='anhar_files/'//TRIM(flanhar)//'.gamma'
CALL add_pressure(filename_gamma)
filename_gamma_ph='anhar_files/'//TRIM(flanhar)//'.gamma_ph'
CALL add_pressure(filename_gamma_ph)
filename_gamma_grun='anhar_files/'//TRIM(flanhar)//'.gamma_grun'
CALL add_pressure(filename_gamma_grun)

CALL gnuplot_xlabel('T (K)',.FALSE.) 
!
!  Part 1: celldm parameters
!
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
   IF (ibrav_save==4.OR.ibrav_save==6.OR.ibrav_save==7) THEN
      CALL gnuplot_ylabel('c (a.u.) ',.FALSE.) 
      IF (ltherm_dos) &
         CALL gnuplot_write_file_mul_data_times(filename2,1,2,3,&
                               'color_red',.TRUE., .NOT.ltherm_freq,.FALSE.)
      IF (ltherm_freq) &
         CALL gnuplot_write_file_mul_data_times(filename3,1,2,3,'color_blue',&
                               .NOT.ltherm_dos, .TRUE.,.FALSE.)
   ENDIF
ELSEIF (ibrav_save==8.OR.ibrav_save==9.OR.ibrav_save==10.OR.ibrav_save==11&
        .OR.ibrav_save==12.OR.ibrav_save==-12.OR.ibrav_save==13.OR. &
            ibrav_save==-13.OR.ibrav_save==14) THEN
   CALL gnuplot_ylabel('b/a ',.FALSE.) 
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename2,1,3,'color_red',.TRUE., &
                                                   .NOT.ltherm_freq,.FALSE.)
   IF (ltherm_freq) & 
      CALL gnuplot_write_file_mul_data(filename3,1,3,'color_blue',&
                                             .NOT.ltherm_dos,.TRUE.,.FALSE.)
   CALL gnuplot_ylabel('b (a.u.) ',.FALSE.) 
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data_times(filename2,1,2,3,&
                               'color_red',.TRUE., .NOT.ltherm_freq,.FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data_times(filename3,1,2,3,'color_blue',&
                               .NOT.ltherm_dos, .TRUE.,.FALSE.)
   CALL gnuplot_ylabel('c/a ',.FALSE.) 
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename2,1,4,'color_red',.TRUE., &
                                             .NOT.ltherm_freq,.FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data(filename3,1,4,'color_blue', &
                                             .NOT.ltherm_dos, .TRUE.,.FALSE.)
   CALL gnuplot_ylabel('c (a.u.) ',.FALSE.) 
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data_times(filename2,1,2,4,&
                               'color_red',.TRUE., .NOT.ltherm_freq,.FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data_times(filename3,1,2,4,'color_blue',&
                               .NOT.ltherm_dos, .TRUE.,.FALSE.)

   IF (ibrav_save==12.OR.ibrav_save==13) THEN
      CALL gnuplot_ylabel('cos({/Symbol a})',.FALSE.) 
      IF (ltherm_dos) &
         CALL gnuplot_write_file_mul_data(filename2,1,5,'color_red',.TRUE., &
                                             .NOT.ltherm_freq,.FALSE.)
      IF (ltherm_freq) &
         CALL gnuplot_write_file_mul_data(filename3,1,5,'color_blue', &
                                             .NOT.ltherm_dos, .TRUE.,.FALSE.)
   ELSEIF (ibrav_save==-12.OR.ibrav_save==-13) THEN
      CALL gnuplot_ylabel('cos({/Symbol b})',.FALSE.) 
      IF (ltherm_dos) &
         CALL gnuplot_write_file_mul_data(filename2,1,5,'color_red',.TRUE., &
                                             .NOT.ltherm_freq,.FALSE.)
      IF (ltherm_freq) &
         CALL gnuplot_write_file_mul_data(filename3,1,5,'color_blue', &
                                             .NOT.ltherm_dos, .TRUE.,.FALSE.)
   ELSEIF (ibrav_save==14) THEN
      CALL gnuplot_ylabel('cos({/Symbol a})',.FALSE.) 
      IF (ltherm_dos) &
         CALL gnuplot_write_file_mul_data(filename2,1,5,'color_red',.TRUE., &
                                             .NOT.ltherm_freq,.FALSE.)
      IF (ltherm_freq) &
         CALL gnuplot_write_file_mul_data(filename3,1,5,'color_blue', &
                                             .NOT.ltherm_dos, .TRUE.,.FALSE.)
      CALL gnuplot_ylabel('cos({/Symbol b})',.FALSE.) 
      IF (ltherm_dos) &
         CALL gnuplot_write_file_mul_data(filename2,1,6,'color_red',.TRUE., &
                                             .NOT.ltherm_freq,.FALSE.)
      IF (ltherm_freq) &
         CALL gnuplot_write_file_mul_data(filename3,1,6,'color_blue', &
                                             .NOT.ltherm_dos, .TRUE.,.FALSE.)
      CALL gnuplot_ylabel('cos({/Symbol c})',.FALSE.) 
      IF (ltherm_dos) &
         CALL gnuplot_write_file_mul_data(filename2,1,7,'color_red',.TRUE., &
                                             .NOT.ltherm_freq,.FALSE.)
      IF (ltherm_freq) &
         CALL gnuplot_write_file_mul_data(filename3,1,7,'color_blue', &
                                             .NOT.ltherm_dos, .TRUE.,.FALSE.)
   ENDIF
ENDIF
!
!  Part 2: Volume
!
CALL gnuplot_ylabel('Volume ((a.u.)^3)',.FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename,1,2,'color_red',.TRUE., &
                                                  .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename1,1,2,'color_blue',&
                                             .NOT.ltherm_dos, .TRUE., .FALSE.)
!
!  Part 3: Helmholtz (or Gibbs) free energy
!
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
!
!  Part 4: Thermal expansion
!
CALL gnuplot_ylabel('Thermal expansion {/Symbol a} x 10^6 (K^{-1})',.FALSE.) 
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
!
!  Part 5: Volume thermal expansion
!
CALL gnuplot_ylabel('Volume thermal expansion {/Symbol b} x 10^6 (K^{-1})',.FALSE.) 

IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename,1,4,'color_red',.TRUE., &
                           .NOT.(ltherm_freq.OR.done_grun), .FALSE.)

IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename1,1,4,'color_blue', &
                                    .NOT.ltherm_dos,.NOT.done_grun, .FALSE.)

IF (done_grun) &
   CALL gnuplot_write_file_mul_data(filename_aux_grun,1,2,'color_green', &
                  .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE., .FALSE.)
!
!  Part 6: C_e heat capacity
!
factor = rydberg_si*avogadro 
CALL gnuplot_set_fact(factor,.FALSE.)
CALL gnuplot_ylabel('Heat capacity C_{/Symbol e} (J / K / N / mol)',.FALSE.)
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename_heat,1,2,'color_red',.TRUE.,&
                                         .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename_heat_ph,1,2,'color_blue',&
                                         .NOT.ltherm_dos,.TRUE.,.FALSE.)
!
!  Part 7: C_p Heat capacity
!
IF (isoent_avail) THEN
   CALL gnuplot_set_fact(factor,.FALSE.)
   CALL gnuplot_ylabel('Heat capacity C_p (J / K / N / mol)',.FALSE.)
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename_heat_anis,1,3,'color_red',  &
                            .TRUE.,.NOT.(ltherm_freq.OR.done_grun),.FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data(filename_heat_anis,1,3,'color_blue', &
                                 .NOT.ltherm_dos,.NOT.done_grun,.FALSE.)
   IF (done_grun) &
      CALL gnuplot_write_file_mul_data(filename_heat_anis_grun,1,3,         &
             'color_green',.NOT.(ltherm_dos.OR.ltherm_freq), .TRUE.,.FALSE.)
ENDIF
!
!  Part 8: C_p -C_e Heat capacity
!
IF (isoent_avail) THEN
   CALL gnuplot_set_fact(factor,.FALSE.)
   CALL gnuplot_ylabel('C_{/Symbol s} - C_{/Symbol e} (J / K / N / mol)',&
                                                                   .FALSE.)
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename_heat_anis,1,2,'color_red',&
                             .TRUE., .NOT.(ltherm_freq.OR.done_grun),.FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data(filename_heat_anis_ph,1,2,'color_blue',&
                                      .NOT.ltherm_dos, .NOT.done_grun,.FALSE.)
   IF (done_grun) &
      CALL gnuplot_write_file_mul_data(filename_heat_anis_grun,1,2,&
             'color_green',.NOT.(ltherm_dos.OR.ltherm_freq), .TRUE.,.FALSE.)
ENDIF
!
!  Part 10: Difference C_V-C_e heat capacity
!
noncubic=(ibrav_save/=1.AND.ibrav_save/=2.AND.ibrav_save/=3)
IF (isoent_avail.AND.noncubic) THEN
   CALL gnuplot_set_fact(factor,.FALSE.)
   CALL gnuplot_ylabel('C_{V} - C_{/Symbol e} (J / K / N / mol)',.FALSE.)
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename_heat_anis,1,4,'color_red',&
                              .TRUE., .NOT.(ltherm_freq.OR.done_grun),.FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data(filename_heat_anis_ph,1,4,'color_blue',&
                                      .NOT.ltherm_dos, .NOT.done_grun,.FALSE.)
   IF (done_grun) &
      CALL gnuplot_write_file_mul_data(filename_heat_anis_grun,1,4,&
             'color_green',.NOT.(ltherm_dos.OR.ltherm_freq), .TRUE.,.FALSE.)
ENDIF
!
!  Part 10: Difference between isoentropic and isothermal bulk modulus
!
IF (isoent_avail) THEN
   CALL gnuplot_set_fact(1._DP,.FALSE.)
   CALL gnuplot_ylabel('B_S - B_T (kbar)',.FALSE.)
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename_bulk,1,4,'color_red',.TRUE.,&
                               .NOT.(ltherm_freq.OR.done_grun),.FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data(filename_bulk_ph,1,4,'color_blue',&
                                  .NOT.ltherm_dos,.NOT.done_grun,.FALSE.)

   IF (done_grun) &
      CALL gnuplot_write_file_mul_data(filename_aux_grun,1,4,'color_green', &
                              .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)
ENDIF
!
!  Part 11: Average gruneisen parameter
!
IF (isoent_avail) THEN
   CALL gnuplot_set_fact(1.0_DP,.FALSE.)
   CALL gnuplot_ylabel('Gr\374neisen parameter ({/Symbol g})',.FALSE.)
   CALL gnuplot_write_horizontal_line(0.0_DP, 2, 'front', 'color_black', &
                                                                  .FALSE.)
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename_gamma,1,2,'color_red',.TRUE.,&
                                .NOT.(done_grun.OR.ltherm_freq),.FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data(filename_gamma_ph,1,2,'color_blue',&
                                     .NOT.ltherm_dos,.NOT.done_grun,.FALSE.)
   IF (done_grun) &
      CALL gnuplot_write_file_mul_data(filename_gamma_grun,1,2,'color_green',&
                          .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)
ENDIF
CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!  CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_anhar_anis

!-----------------------------------------------------------------
SUBROUTINE plot_thermal_stress()
!-----------------------------------------------------------------
!
!  This routine plot the thermal stress tensor
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
USE control_elastic_constants, ONLY : lelastic, lelasticf
USE gnuplot_color,    ONLY : gnuplot_set_greens
USE data_files,       ONLY : flanhar
USE temperature,      ONLY : tmin, tmax
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: gnu_filename, filename, psfilename, filenameph
INTEGER :: ierr, system, na

IF ( my_image_id /= root_image ) RETURN
IF (.NOT.(lelastic.OR.lelasticf)) RETURN

IF (.NOT.(ltherm_freq.OR.ltherm_dos)) RETURN

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'.tstress'
CALL gnuplot_start(gnu_filename)

psfilename=TRIM(flpsanhar)//'.tstress'//TRIM(flext)
IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(psfilename, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ELSE
   CALL gnuplot_write_header(psfilename, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ENDIF
CALL gnuplot_set_greens()
CALL gnuplot_xlabel('T (K)', .FALSE.)

filename='anhar_files/'//TRIM(flanhar)//'.tstress'
filenameph='anhar_files/'//TRIM(flanhar)//'.tstress_ph'
!
!   First the diagonal components
!
CALL gnuplot_ylabel('Thermal stress b_{ii} (kbar)',.FALSE.)

IF (ltherm_dos) THEN
   CALL gnuplot_write_file_mul_data(filename,1,2,'color_red',.TRUE., &
                                                   .FALSE.,.FALSE.)

   CALL gnuplot_write_file_mul_data(filename,1,5,'color_blue',.FALSE., &
                                                   .FALSE.,.FALSE.)

   CALL gnuplot_write_file_mul_data(filename,1,7,'color_dark_spring_green',&
                                .FALSE., .NOT.ltherm_freq,.FALSE.)
ENDIF

IF (ltherm_freq) THEN
   CALL gnuplot_write_file_mul_data(filenameph,1,2,'color_pink', &
                                               .NOT.ltherm_dos,.FALSE.,.FALSE.)

   CALL gnuplot_write_file_mul_data(filenameph,1,5,'color_light_blue', &
                                               .FALSE., .FALSE.,.FALSE.)

   CALL gnuplot_write_file_mul_data(filenameph,1,7,'color_green', &
                                               .FALSE.,.TRUE.,.FALSE.)
ENDIF
!
!  And then the off diagonal, only for noncubic solids
!
IF (ibrav/=1.AND.ibrav/=2.AND.ibrav/=3) THEN
   CALL gnuplot_ylabel('Thermal stress b_{ij} (kbar)',.FALSE.)

   IF (ltherm_dos) THEN
      CALL gnuplot_write_file_mul_data(filename,1,3,'color_red',.TRUE., &
                                               .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',.FALSE., &
                                                .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename,1,6,&
                   'color_dark_spring_green',.FALSE.,.NOT.ltherm_freq,.FALSE.)
   ENDIF

   IF (ltherm_freq) THEN
      CALL gnuplot_write_file_mul_data(filenameph,1,3,'color_pink', &
                                               .NOT.ltherm_dos,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filenameph,1,4,'color_light_blue', &
                                               .FALSE.,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filenameph,1,6,'color_green', &
                                               .FALSE.,.TRUE.,.FALSE.)
   ENDIF
ENDIF

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_thermal_stress

!-----------------------------------------------------------------
SUBROUTINE plot_generalized_gruneisen()
!-----------------------------------------------------------------
!
!  This routine plot the thermal stress tensor
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
USE control_elastic_constants, ONLY : lelastic, lelasticf
USE gnuplot_color,    ONLY : gnuplot_set_greens
USE data_files,       ONLY : flanhar
USE temperature,      ONLY : tmin, tmax
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: gnu_filename, filename, psfilename, filenameph
INTEGER :: ierr, system, na

IF ( my_image_id /= root_image ) RETURN
IF (.NOT.(ltherm_freq.OR.ltherm_dos)) RETURN
IF (.NOT.(lelastic.OR.lelasticf)) RETURN

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'.ggamma'
CALL gnuplot_start(gnu_filename)

psfilename=TRIM(flpsanhar)//'.ggamma'//TRIM(flext)
IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(psfilename, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ELSE
   CALL gnuplot_write_header(psfilename, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ENDIF
CALL gnuplot_set_greens()
CALL gnuplot_xlabel('T (K)', .FALSE.)

filename='anhar_files/'//TRIM(flanhar)//'.ggamma'
filenameph='anhar_files/'//TRIM(flanhar)//'.ggamma_ph'
!
!   First the diagonal components
!
CALL gnuplot_ylabel('Gr\374neisen parameters ({/Symbol g}_{ii})',.FALSE.)

IF (ltherm_dos) THEN
   CALL gnuplot_write_file_mul_data(filename,1,2,'color_red',.TRUE., &
                                                   .FALSE.,.FALSE.)

   CALL gnuplot_write_file_mul_data(filename,1,5,'color_blue',.FALSE., &
                                                   .FALSE.,.FALSE.)

   CALL gnuplot_write_file_mul_data(filename,1,7,'color_dark_spring_green',&
                                .FALSE., .NOT.ltherm_freq,.FALSE.)
ENDIF

IF (ltherm_freq) THEN
   CALL gnuplot_write_file_mul_data(filenameph,1,2,'color_pink', &
                                               .NOT.ltherm_dos,.FALSE.,.FALSE.)

   CALL gnuplot_write_file_mul_data(filenameph,1,5,'color_light_blue', &
                                               .FALSE., .FALSE.,.FALSE.)

   CALL gnuplot_write_file_mul_data(filenameph,1,7,'color_green', &
                                               .FALSE.,.TRUE.,.FALSE.)
ENDIF
!
!  And then the off diagonal, only for noncubic solids
!
IF (ibrav/=1.AND.ibrav/=2.AND.ibrav/=3) THEN

   CALL gnuplot_ylabel('Gr\374neisen parameters ({/Symbol g}_{ij})',.FALSE.)

   IF (ltherm_dos) THEN
      CALL gnuplot_write_file_mul_data(filename,1,3,'color_red',.TRUE., &
                                               .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',.FALSE., &
                                                .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename,1,6,&
                   'color_dark_spring_green',.FALSE.,.NOT.ltherm_freq,.FALSE.)
   ENDIF

   IF (ltherm_freq) THEN
      CALL gnuplot_write_file_mul_data(filenameph,1,3,'color_pink', &
                                               .NOT.ltherm_dos,.FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filenameph,1,4,'color_light_blue', &
                                               .FALSE.,.FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filenameph,1,6,'color_green', &
                                               .FALSE.,.TRUE.,.FALSE.)
   ENDIF
ENDIF

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_generalized_gruneisen

! Copyright (C) 2018 Cristiano Malica

!-----------------------------------------------------------------
SUBROUTINE plot_dw_anhar_anis()
!-----------------------------------------------------------------
!
!  This is a driver to plot the quantities written inside 
!
USE kinds,            ONLY : DP
USE ions_base,        ONLY : nat
USE cell_base,        ONLY : ibrav
USE control_thermo,   ONLY : ltherm_dos, ltherm_freq, with_eigen
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
IF (.NOT.with_eigen) RETURN

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'.anhar_anis_dw'
CALL gnuplot_start(gnu_filename)

psfilename=TRIM(flpsanhar)//'.anis_dw'//TRIM(flext)
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

   filename='anhar_files/'//TRIM(flanhar)//'.anis_ph.'//TRIM(int_to_char(na))//'.dw'
   filetherm='anhar_files/'//TRIM(flanhar)//'.anis.'//TRIM(int_to_char(na))//'.dw'
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
!
!   Then the off diagonal components
!
      CALL gnuplot_ylabel('B_{ij} ({\305}^2) (atom '// &
                                        TRIM(int_to_char(na))//')',.FALSE.)

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
END SUBROUTINE plot_dw_anhar_anis
