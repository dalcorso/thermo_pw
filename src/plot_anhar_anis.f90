!
! Copyright (C) 2015-2022 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------
SUBROUTINE plot_anhar_anis_celldm()
!-----------------------------------------------------------------
!
!  This is a driver to plot the quantities written inside 
!  flanhar//'.celldm', flanhar//'.celldm_press'. In the same postscript
!  file it plots celldm as a function of temperature at the input pressure
!  and for several pressures.
!
USE kinds,           ONLY : DP
USE control_gnuplot, ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE postscript_files, ONLY : flpsanhar
USE gnuplot,         ONLY : gnuplot_start, gnuplot_end,  &
                            gnuplot_write_header,        &
                            gnuplot_ylabel,              &
                            gnuplot_xlabel,              &
                            gnuplot_write_command,       &
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
USE lattices,        ONLY : celldm_gnuplot_name, needed_celldm
USE nye,             ONLY : thermal_gnuplot_name, needed_tensor2
USE temperature,     ONLY : tmin, tmax, temp, ntemp_plot, itemp_plot
USE control_pressure, ONLY : pressure_kb, press, npress_plot, ipress_plot, &
                             npress, press
USE color_mod,       ONLY : color
USE mp_images,       ONLY : my_image_id, root_image
USE io_global,       ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename_ph, filename1, &
                      filename2, filename3, filenameps, label
INTEGER :: i, istep, ipressp, ipress, itempp, itemp, iusing, ierr, system, &
           last_iusing
REAL(DP) :: factor
LOGICAL :: isoent_avail, noncubic, celldm_in_use(6), tensor2_in_use(6), &
           first_step, last_step, first, last

IF ( my_image_id /= root_image ) RETURN

isoent_avail=lelastic.OR.lelasticf
gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_anhar_celldm'
CALL add_pressure(gnu_filename)

CALL gnuplot_start(gnu_filename)

filenameps=TRIM(flpsanhar)//'.celldm'
CALL add_pressure(filenameps)
filenameps=TRIM(filenameps)//TRIM(flext)
IF (tmin /= 1.0_DP) THEN
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                            flext ) 
ELSE
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP,&
                                                            flext ) 
ENDIF

filename1='anhar_files/'//TRIM(flanhar)//'.celldm'
CALL add_pressure(filename1)
filename2='anhar_files/'//TRIM(flanhar)//'.celldm_ph'
CALL add_pressure(filename2)
filename3='anhar_files/'//TRIM(flanhar)//'.celldm_grun'
CALL add_pressure(filename3)

CALL gnuplot_xlabel('T (K)',.FALSE.) 
!
!  Part 1: celldm parameters
!
CALL gnuplot_set_fact(1.0_DP,.FALSE.)

CALL needed_celldm(ibrav_save, celldm_in_use)

iusing=1
DO i=1,6
   IF (celldm_in_use(i)) THEN
      iusing=iusing+1
      last_iusing=iusing
   ENDIF
ENDDO

iusing=1
DO i=1,6
   IF (celldm_in_use(i)) THEN
      iusing=iusing+1
      CALL gnuplot_ylabel(TRIM(celldm_gnuplot_name(i)),.FALSE.) 
      IF (ltherm_dos) &
         CALL gnuplot_write_file_mul_data(filename1,1,iusing,'color_red',   &
                        .TRUE., .NOT.ltherm_freq,.FALSE.)
      IF (ltherm_freq) &
         CALL gnuplot_write_file_mul_data(filename2,1,iusing,'color_blue',  &
                                 (.NOT.ltherm_dos), .TRUE., .FALSE.)
!
!  plot also b and c
!
      IF (i==2.OR.i==3) THEN
         IF (i==2) CALL gnuplot_ylabel("b (a.u.)",.FALSE.) 
         IF (i==3) CALL gnuplot_ylabel("c (a.u.)",.FALSE.) 
         IF (ltherm_dos) &
            CALL gnuplot_write_file_mul_data_times(filename1,1,2,iusing, &
                         'color_red', .TRUE., .NOT.ltherm_freq,.FALSE.)
         IF (ltherm_freq) &
            CALL gnuplot_write_file_mul_data_times(filename2,1,2,iusing, &
                        'color_blue', (.NOT.ltherm_dos), .TRUE., .FALSE.)
      ENDIF
   ENDIF
ENDDO
!
!   celldm as a function of T for several pressures
!
iusing=1
DO i=1,6
   IF (celldm_in_use(i)) THEN
      iusing=iusing+1
      istep=0
      DO ipressp=1,npress_plot
         first_step=(ipressp==1)
         last_step=(ipressp==npress_plot)
         ipress=ipress_plot(ipressp)
         istep=MOD(istep,8)+1
         filename="anhar_files/"//TRIM(flanhar)//'.celldm_press'
         CALL add_value(filename,press(ipress))
         filename_ph="anhar_files/"//TRIM(flanhar)//'.celldm_ph_press'
         CALL add_value(filename_ph,press(ipress))
         IF (first_step) THEN
            CALL gnuplot_xlabel('T (K)',.FALSE.)
            CALL gnuplot_set_fact(1.0_DP,.FALSE.)
            CALL gnuplot_ylabel(TRIM(celldm_gnuplot_name(i)),.FALSE.)
         ENDIF
         IF (ltherm_dos) &
            CALL gnuplot_write_file_mul_data(filename,1,iusing,color(istep), &
                     first_step, (last_step.AND..NOT.ltherm_freq), .FALSE.)
         IF (ltherm_freq) &
            CALL gnuplot_write_file_mul_data(filename_ph,1,iusing,   &
                    color(istep), (first_step.AND..NOT.ltherm_dos),  &
                                                     last_step, .FALSE.)
      ENDDO
      IF (i==2.OR.i==3) THEN
         istep=0
         DO ipressp=1,npress_plot
            first_step=(ipressp==1)
            last_step=(ipressp==npress_plot)
            ipress=ipress_plot(ipressp)
            istep=MOD(istep,8)+1
            filename="anhar_files/"//TRIM(flanhar)//'.celldm_press'
            CALL add_value(filename,press(ipress))
            filename_ph="anhar_files/"//TRIM(flanhar)//'.celldm_ph_press'
            CALL add_value(filename_ph,press(ipress))
            IF (first_step) THEN
               CALL gnuplot_xlabel('T (K)',.FALSE.)
               CALL gnuplot_set_fact(1.0_DP,.FALSE.)
               IF (i==2) CALL gnuplot_ylabel("b (a.u.)",.FALSE.) 
               IF (i==3) CALL gnuplot_ylabel("c (a.u.)",.FALSE.) 
            ENDIF
            IF (ltherm_dos) &
               CALL gnuplot_write_file_mul_data_times(filename,1,2,iusing,&
                color(istep), first_step, (last_step.AND..NOT.ltherm_freq),&
                                                                     .FALSE.)
            IF (ltherm_freq) &
               CALL gnuplot_write_file_mul_data_times(filename_ph,1,2,iusing,&
                    color(istep), (first_step.AND..NOT.ltherm_dos),  &
                                                     last_step, .FALSE.)
         ENDDO
      ENDIF
   ENDIF
ENDDO
!
!   celldm as a function of pressure for several temperatures
!
iusing=1
DO i=1,6
   IF (celldm_in_use(i)) THEN
      iusing=iusing+1
      istep=0
      DO itempp=1,ntemp_plot
         first_step=(itempp==1)
         last_step=(itempp==ntemp_plot)
         itemp=itemp_plot(itempp)
         istep=MOD(istep,8)+1
         filename="anhar_files/"//TRIM(flanhar)//'.celldm_temp'
         CALL add_value(filename,temp(itemp))
         filename_ph="anhar_files/"//TRIM(flanhar)//'.celldm_ph_temp'
         CALL add_value(filename_ph,temp(itemp))
         IF (first_step) THEN
            WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') &
                                  MAX(0.0_DP,press(1)), press(npress)
            CALL gnuplot_write_command(TRIM(label),.FALSE.)
            CALL gnuplot_xlabel('p (kbar)',.FALSE.)
            CALL gnuplot_set_fact(1.0_DP,.FALSE.)
            CALL gnuplot_ylabel(TRIM(celldm_gnuplot_name(i)),.FALSE.)
         ENDIF
         IF (ltherm_dos) &
            CALL gnuplot_write_file_mul_data(filename,1,iusing,color(istep), &
                    first_step, (last_step.AND..NOT.ltherm_freq), .FALSE.)
         IF (ltherm_freq) &
            CALL gnuplot_write_file_mul_data(filename_ph,1,iusing, &
             color(istep), (first_step.AND..NOT.ltherm_dos), last_step, &
                                                                 .FALSE.)
      ENDDO
      IF (i==2.OR.i==3) THEN
         istep=0
         DO itempp=1,ntemp_plot
            first_step=(itempp==1)
            last_step=(itempp==ntemp_plot)
            itemp=itemp_plot(itempp)
            istep=MOD(istep,8)+1
            filename="anhar_files/"//TRIM(flanhar)//'.celldm_temp'
            CALL add_value(filename,temp(itemp))
            filename_ph="anhar_files/"//TRIM(flanhar)//'.celldm_ph_temp'
            CALL add_value(filename_ph,temp(itemp))
            IF (first_step) THEN
               WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') &
                                  MAX(0.0_DP,press(1)), press(npress)
               CALL gnuplot_write_command(TRIM(label),.FALSE.)
               CALL gnuplot_xlabel('p (kbar)',.FALSE.)
               CALL gnuplot_set_fact(1.0_DP,.FALSE.)
               IF (i==2) CALL gnuplot_ylabel("b (a.u.)",.FALSE.) 
               IF (i==3) CALL gnuplot_ylabel("c (a.u.)",.FALSE.) 
            ENDIF
            IF (ltherm_dos) &
               CALL gnuplot_write_file_mul_data_times(filename,1,2,iusing, &
              color(istep), first_step, (last_step.AND..NOT.ltherm_freq),  &
                                                                   .FALSE.)
            IF (ltherm_freq) &
               CALL gnuplot_write_file_mul_data_times(filename_ph,1,2,iusing,&
                color(istep), (first_step.AND..NOT.ltherm_dos), last_step, &
                                                                 .FALSE.)
         ENDDO
      ENDIF
   ENDIF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!  CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)
RETURN
END SUBROUTINE plot_anhar_anis_celldm
!
!-----------------------------------------------------------------
SUBROUTINE plot_anhar_anis_alpha()
!-----------------------------------------------------------------
!
!  This is a driver to plot the quantities written inside 
!  flanhar//'.celldm', flanhar//'.celldm_press'. In the same postscript
!  file it plots the thermal expansion tensor as a function of temperature 
!  at the input pressure and for several pressures. It writes also
!  the thermal expansion as a function of pressure at several temperatures.
!
USE kinds,           ONLY : DP
USE control_gnuplot, ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE postscript_files, ONLY : flpsanhar
USE gnuplot,         ONLY : gnuplot_start, gnuplot_end,  &
                            gnuplot_write_header,        &
                            gnuplot_ylabel,              &
                            gnuplot_xlabel,              &
                            gnuplot_write_command,       &
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
USE control_pressure, ONLY : press
USE lattices,        ONLY : celldm_gnuplot_name, needed_celldm
USE nye,             ONLY : thermal_gnuplot_name, needed_tensor2
USE temperature,     ONLY : tmin, tmax, temp, ntemp_plot, itemp_plot
USE control_pressure, ONLY : pressure_kb, npress_plot, ipress_plot, &
                             npress, press
USE color_mod,       ONLY : color
USE mp_images,       ONLY : my_image_id, root_image
USE io_global,       ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename_ph, filename1, &
                      filename2, filename3, filenameps, label
INTEGER :: i, istep, ipressp, ipress, itempp, itemp, iusing, ierr, system, &
           last_iusing
REAL(DP) :: factor
LOGICAL :: isoent_avail, noncubic, celldm_in_use(6), tensor2_in_use(6), &
           first_step, last_step, first, last

IF ( my_image_id /= root_image ) RETURN

isoent_avail=lelastic.OR.lelasticf
gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_anhar_alpha'
CALL add_pressure(gnu_filename)

CALL gnuplot_start(gnu_filename)

filenameps=TRIM(flpsanhar)//'.alpha'
CALL add_pressure(filenameps)
filenameps=TRIM(filenameps)//TRIM(flext)
IF (tmin /= 1.0_DP) THEN
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                            flext ) 
ELSE
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP,&
                                                            flext ) 
ENDIF

filename1='anhar_files/'//TRIM(flanhar)//'.celldm'
CALL add_pressure(filename1)
filename2='anhar_files/'//TRIM(flanhar)//'.celldm_ph'
CALL add_pressure(filename2)
filename3='anhar_files/'//TRIM(flanhar)//'.celldm_grun'
CALL add_pressure(filename3)

CALL gnuplot_xlabel('T (K)',.FALSE.) 
!
!  Part 1: celldm parameters
!
CALL gnuplot_set_fact(1.0_DP,.FALSE.)

CALL needed_celldm(ibrav_save, celldm_in_use)

iusing=1
DO i=1,6
   IF (celldm_in_use(i)) THEN
      iusing=iusing+1
      last_iusing=iusing
   ENDIF
ENDDO
!
!  Thermal expansion at input pressure as a function of T
!
CALL needed_tensor2(ibrav_save, tensor2_in_use)

iusing=last_iusing
DO i=1,6
   IF (tensor2_in_use(i)) THEN 
      WRITE(label,'("Thermal expansion ", a," x 10^6 (K^{-1})")') &
                                              TRIM(thermal_gnuplot_name(i))
      CALL gnuplot_ylabel(TRIM(label),.FALSE.) 
      iusing=iusing+1
      IF (done_grun) &
         CALL gnuplot_write_file_mul_data(filename3,1,iusing,'color_green',&
                   .TRUE., (.NOT.(ltherm_dos.OR.ltherm_freq)),.FALSE.)
      IF (ltherm_dos) &
         CALL gnuplot_write_file_mul_data(filename1,1,iusing,'color_red',&
                       .NOT.done_grun,(.NOT.ltherm_freq),.FALSE.)
      IF (ltherm_freq) &
         CALL gnuplot_write_file_mul_data(filename2,1,iusing,'color_blue', &
                         .NOT.(ltherm_dos.OR.done_grun), .TRUE., .FALSE.)
   ENDIF
ENDDO
!
!  Temperature dependent thermal expansion at several pressures
!
iusing=last_iusing
DO i=1,6
   IF (tensor2_in_use(i)) THEN 
      iusing=iusing+1
      istep=0
      DO ipressp=1,npress_plot
         first_step=(ipressp==1)
         last_step=(ipressp==npress_plot)
         ipress=ipress_plot(ipressp)
         istep=MOD(istep,8)+1
         filename="anhar_files/"//TRIM(flanhar)//'.celldm_press'
         CALL add_value(filename,press(ipress))
         filename_ph="anhar_files/"//TRIM(flanhar)//'.celldm_ph_press'
         CALL add_value(filename_ph,press(ipress))
         IF (first_step) THEN
            CALL gnuplot_xlabel('T (K)',.FALSE.)
            CALL gnuplot_set_fact(1.0_DP,.FALSE.)
            WRITE(label,'("Thermal expansion ", a," x 10^6 (K^{-1})")') &
                                              TRIM(thermal_gnuplot_name(i))
            CALL gnuplot_ylabel(TRIM(label),.FALSE.) 
         ENDIF
         IF (ltherm_dos) &
            CALL gnuplot_write_file_mul_data(filename,1,iusing,color(istep), &
                   first_step,(last_step.AND..NOT.ltherm_freq),.FALSE.)
         IF (ltherm_freq) &
            CALL gnuplot_write_file_mul_data(filename_ph,1,iusing,&
               color(istep),(first_step.AND..NOT.ltherm_dos),last_step,.FALSE.)
      ENDDO
   ENDIF
ENDDO
!
!  Pressure dependent thermal expansion at several temperatures
!
iusing=last_iusing
DO i=1,6
   IF (tensor2_in_use(i)) THEN 
      iusing=iusing+1
      istep=0
      DO itempp=1,ntemp_plot
         first_step=(itempp==1)
         last_step=(itempp==ntemp_plot)
         itemp=itemp_plot(itempp)
         istep=MOD(istep,8)+1
         filename="anhar_files/"//TRIM(flanhar)//'.celldm_temp'
         CALL add_value(filename,temp(itemp))
         filename_ph="anhar_files/"//TRIM(flanhar)//'.celldm_ph_temp'
         CALL add_value(filename_ph,temp(itemp))
         IF (first_step) THEN
            WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') &
                                   MAX(0.0_DP,press(1)), press(npress)
            CALL gnuplot_write_command(TRIM(label),.FALSE.) 
            CALL gnuplot_xlabel('p (kbar)',.FALSE.)
            CALL gnuplot_set_fact(1.0_DP,.FALSE.)
            WRITE(label,'("Thermal expansion ", a," x 10^6 (K^{-1})")') &
                                              TRIM(thermal_gnuplot_name(i))
            CALL gnuplot_ylabel(TRIM(label),.FALSE.) 
         ENDIF
         IF (ltherm_dos) &
            CALL gnuplot_write_file_mul_data(filename,1,iusing,color(istep), &
                         first_step,(last_step.AND..NOT.ltherm_freq),.FALSE.)
         IF (ltherm_freq) &
            CALL gnuplot_write_file_mul_data(filename_ph,1,iusing,&
             color(istep), (first_step.AND..NOT.ltherm_dos),last_step,.FALSE.)
      ENDDO
   ENDIF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!  CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_anhar_anis_alpha
!
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

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_anhar_tstress'
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

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_anhar_ggamma'
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
SUBROUTINE plot_anhar_anis_dw()
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
END SUBROUTINE plot_anhar_anis_dw
