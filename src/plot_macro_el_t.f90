!
! Copyright (C) 2020 Cristiano Malica 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE plot_sound_t()
!-------------------------------------------------------------------------

USE kinds,               ONLY : DP
USE temperature,         ONLY : tmin, tmax
USE data_files,          ONLY : flanhar
USE postscript_files,    ONLY : flpsanhar
USE control_elastic_constants, ONLY : lelastic, lelasticf
USE control_gnuplot,     ONLY : gnuplot_command, flgnuplot, lgnuplot, flext
USE gnuplot,             ONLY : gnuplot_start, gnuplot_end,           &
                                gnuplot_write_header, gnuplot_xlabel, &
                                gnuplot_ylabel, gnuplot_write_file_mul_data
USE io_global,        ONLY : ionode
USE mp_images,        ONLY : root_image, my_image_id

IMPLICIT NONE

CHARACTER(LEN=256) :: filelastic, filelastic_s, &
                      filelastic_ph, filelastic_s_ph, gnu_filename, &
                      filenameps, filedebye_ph, filedebye
INTEGER :: ierr, system

IF (.NOT.(lelastic.OR.lelasticf)) RETURN

filelastic=TRIM(filelastic)
gnu_filename=TRIM(gnu_filename)
filenameps=TRIM(filenameps)

IF (lelasticf) THEN
   filelastic_ph="anhar_files/"//TRIM(flanhar)//".sound_vel_ph"
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//".sound_vel_s_ph"
   filedebye_ph="anhar_files/"//TRIM(flanhar)//".macro_el_debye_ph"
END IF

IF (lelastic) THEN
   filelastic="anhar_files/"//TRIM(flanhar)//".sound_vel"
   filelastic_s="anhar_files/"//TRIM(flanhar)//".sound_vel_s"
   filedebye="anhar_files/"//TRIM(flanhar)//".macro_el_debye"
END IF

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_sound_vel"
filenameps=TRIM(flpsanhar)//".sound_vel"//TRIM(flext)

CALL gnuplot_start(gnu_filename)

IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ELSE
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ENDIF

CALL gnuplot_xlabel('T (K)', .FALSE.)

CALL gnuplot_ylabel('V_{P} (m/s)',.FALSE.)

IF (lelastic) THEN
   CALL gnuplot_write_file_mul_data(filelastic,1,2,'color_red',.TRUE.,  &
                                                             .FALSE.,.FALSE.)

   CALL gnuplot_write_file_mul_data(filelastic_s,1,2,'color_green',.FALSE., &
                                                       .NOT.lelasticf,.FALSE.)
ENDIF
IF (lelasticf) THEN
   CALL gnuplot_write_file_mul_data(filelastic_ph,1,2,'color_blue',&
                                                .NOT.lelastic,.FALSE.,.FALSE.)

   CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,2,'color_orange',&
                                             .FALSE.,.TRUE.,.FALSE.)
ENDIF
CALL gnuplot_ylabel('V_{B} (m/s)',.FALSE.)

IF (lelastic) THEN
   CALL gnuplot_write_file_mul_data(filelastic,1,3,'color_red',.TRUE.,  &
                                                             .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filelastic_s,1,3,'color_green',.FALSE., &
                                                      .NOT.lelasticf,.FALSE.)
ENDIF
IF (lelasticf) THEN
   CALL gnuplot_write_file_mul_data(filelastic_ph,1,3,'color_blue', &
                                   .NOT.lelastic, .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,3,'color_orange',&
                            .FALSE.,.TRUE.,.FALSE.)
ENDIF
CALL gnuplot_ylabel('V_{G} (m/s)',.FALSE.)

IF (lelastic) THEN
   CALL gnuplot_write_file_mul_data(filelastic,1,4,'color_red',.TRUE., &
                                                             .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filelastic_s,1,4,'color_green',.FALSE., &
                                                     .NOT.lelasticf,.FALSE.)
ENDIF
IF (lelasticf) THEN
   CALL gnuplot_write_file_mul_data(filelastic_ph,1,4,'color_blue', &
                     .NOT.lelastic,.FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,4,'color_orange',&
                                 .FALSE., .TRUE.,.FALSE.)
ENDIF

CALL gnuplot_ylabel('{/Symbol Q}_D (K)',.FALSE.)

IF (lelastic) THEN
   CALL gnuplot_write_file_mul_data(filedebye,1,2,'color_red',.TRUE., &
                                                             .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filedebye,1,3,'color_green',.FALSE., &
                                                     .NOT.lelasticf,.FALSE.)
ENDIF
IF (lelasticf) THEN
   CALL gnuplot_write_file_mul_data(filedebye_ph,1,2,'color_blue', &
                     .NOT.lelastic,.FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filedebye_ph,1,3,'color_orange',&
                                 .FALSE., .TRUE.,.FALSE.)
ENDIF

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
           ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN

END SUBROUTINE plot_sound_t

! Copyright (C) 2023 Andrea Dal Corso
!
!  This routines generalize the previous one by plotting quantities
!  as a function of temperature at several pressures or as a function
!  of pressure at several temperatures
!

!-------------------------------------------------------------------------
SUBROUTINE plot_macro_el_new_t()
!-------------------------------------------------------------------------

USE kinds,               ONLY : DP
USE temperature,         ONLY : tmin, tmax
USE data_files,          ONLY : flanhar
USE postscript_files,    ONLY : flpsanhar
USE control_elastic_constants, ONLY : lelastic, lelasticf
USE control_gnuplot,     ONLY : gnuplot_command, flgnuplot, lgnuplot, flext
USE gnuplot,             ONLY : gnuplot_start, gnuplot_end,           &
                                gnuplot_write_header, gnuplot_xlabel, &
                                gnuplot_ylabel, gnuplot_write_file_mul_data, &
                                gnuplot_write_file_mul_data_div,      &
                                gnuplot_write_file_mul_data_linear
USE io_global,        ONLY : ionode
USE mp_images,        ONLY : root_image, my_image_id

IMPLICIT NONE

CHARACTER(LEN=256) :: filelastic, filelastic_s, &
                      filelastic_ph, filelastic_s_ph, gnu_filename, filenameps
INTEGER :: ierr, system

IF (.NOT.(lelastic.OR.lelasticf)) RETURN

filelastic=TRIM(filelastic)
gnu_filename=TRIM(gnu_filename)
filenameps=TRIM(filenameps)

IF (lelasticf) THEN
   filelastic_ph="anhar_files/"//TRIM(flanhar)//".macro_el_ph_aver"
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//".macro_el_s_ph_aver"
END IF

IF (lelastic) THEN
   filelastic="anhar_files/"//TRIM(flanhar)//".macro_el_aver"
   filelastic_s="anhar_files/"//TRIM(flanhar)//".macro_el_s_aver"
END IF

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_macro_el"
filenameps=TRIM(flpsanhar)//".macro_el"//TRIM(flext)

CALL gnuplot_start(gnu_filename)

IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ELSE
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ENDIF

CALL gnuplot_xlabel('T (K)', .FALSE.)

CALL gnuplot_ylabel('Bulk Modulus (kbar)',.FALSE.)

IF (lelastic) THEN
   CALL gnuplot_write_file_mul_data(filelastic,1,2,'color_red',.TRUE.,  &
                                                             .FALSE.,.FALSE.)

   CALL gnuplot_write_file_mul_data(filelastic_s,1,2,'color_green',.FALSE., &
                                                       .NOT.lelasticf,.FALSE.)
ENDIF
IF (lelasticf) THEN
   CALL gnuplot_write_file_mul_data(filelastic_ph,1,2,'color_blue',&
                                                .NOT.lelastic,.FALSE.,.FALSE.)

   CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,2,'color_orange',&
                                             .FALSE.,.TRUE.,.FALSE.)
ENDIF
CALL gnuplot_ylabel('Young modulus (kbar)',.FALSE.)

IF (lelastic) THEN
   CALL gnuplot_write_file_mul_data(filelastic,1,3,'color_red',.TRUE.,  &
                                                             .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filelastic_s,1,3,'color_green',.FALSE., &
                                                      .NOT.lelasticf,.FALSE.)
ENDIF
IF (lelasticf) THEN
   CALL gnuplot_write_file_mul_data(filelastic_ph,1,3,'color_blue', &
                                   .NOT.lelastic, .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,3,'color_orange',&
                            .FALSE.,.TRUE.,.FALSE.)
ENDIF
CALL gnuplot_ylabel('Shear modulus (kbar)',.FALSE.)

IF (lelastic) THEN
   CALL gnuplot_write_file_mul_data(filelastic,1,4,'color_red',.TRUE., &
                                                             .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filelastic_s,1,4,'color_green',.FALSE., &
                                                     .NOT.lelasticf,.FALSE.)
ENDIF
IF (lelasticf) THEN
   CALL gnuplot_write_file_mul_data(filelastic_ph,1,4,'color_blue', &
                     .NOT.lelastic,.FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,4,'color_orange',&
                                 .FALSE., .TRUE.,.FALSE.)
ENDIF

CALL gnuplot_ylabel('Longitudinal modulus (kbar)',.FALSE.)

IF (lelastic) THEN
   CALL gnuplot_write_file_mul_data_linear(filelastic,1,2,4,1.0_DP,&
                 4.0_DP/3.0_DP, 'color_red',.TRUE.,.FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data_linear(filelastic_s,1,2,4,1.0_DP,&
                 4.0_DP/3.0_DP, 'color_green',.FALSE.,.NOT.lelasticf,.FALSE.)
ENDIF
IF (lelasticf) THEN
   CALL gnuplot_write_file_mul_data_linear(filelastic_ph,1,2,4,1.0_DP,&
       4.0_DP/3.0_DP, 'color_blue', .NOT.lelastic, .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data_linear(filelastic_s_ph,1,2,4,1.0_DP,&
       4.0_DP/3.0_DP,'color_orange',.FALSE.,.TRUE.,.FALSE.)
ENDIF

CALL gnuplot_ylabel('Poisson ratio',.FALSE.)

IF (lelastic) THEN
   CALL gnuplot_write_file_mul_data(filelastic,1,5,'color_red',.TRUE., &
                                                             .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filelastic_s,1,5,'color_green',.FALSE., &
                                                     .NOT.lelasticf,.FALSE.)
ENDIF
IF (lelasticf) THEN
   CALL gnuplot_write_file_mul_data(filelastic_ph,1,5,'color_blue', &
                     .NOT.lelastic,.FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,5,'color_orange',&
                                 .FALSE., .TRUE.,.FALSE.)
ENDIF

CALL gnuplot_ylabel('Pugh ratio',.FALSE.)

IF (lelastic) THEN
   CALL gnuplot_write_file_mul_data_div(filelastic,1,4,2,'color_red',&
                                        .TRUE.,.FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data_div(filelastic_s,1,4,2,'color_green',&
                                 .FALSE.,.NOT.lelasticf,.FALSE.)
ENDIF
IF (lelasticf) THEN
   CALL gnuplot_write_file_mul_data_div(filelastic_ph,1,4,2,'color_blue', &
                     .NOT.lelastic,.FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data_div(filelastic_s_ph,1,4,2,'color_orange',&
                                 .FALSE., .TRUE.,.FALSE.)
ENDIF

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
           ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN

END SUBROUTINE plot_macro_el_new_t
!
!-------------------------------------------------------------------------
SUBROUTINE plot_macro_el_new_pt()
!-------------------------------------------------------------------------
!
USE kinds,               ONLY : DP
USE temperature,         ONLY : tmin, tmax
USE data_files,          ONLY : flanhar
USE postscript_files,    ONLY : flpsanhar
USE control_elastic_constants, ONLY : lelastic, lelasticf
USE control_pressure,    ONLY : npress_plot, ipress_plot, press
USE control_gnuplot,     ONLY : gnuplot_command, flgnuplot, lgnuplot, flext
USE gnuplot,             ONLY : gnuplot_start, gnuplot_end,           &
                                gnuplot_write_header, gnuplot_xlabel, &
                                gnuplot_ylabel, gnuplot_write_file_mul_data, &
                                gnuplot_write_file_mul_data_div, &
                                gnuplot_set_fact, &
                                gnuplot_write_file_mul_data_linear
USE color_mod,           ONLY : color
USE io_global,           ONLY : ionode
USE mp_images,           ONLY : root_image, my_image_id

IMPLICIT NONE

CHARACTER(LEN=256) :: filelastic, filelastic_s, &
                      filelastic_ph, filelastic_s_ph, gnu_filename, filenameps
LOGICAL :: first_step, last_step
INTEGER :: istep, ipressp, ipress, ierr, system

IF (npress_plot==0) RETURN
IF (.NOT.(lelastic.OR.lelasticf)) RETURN

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_macro_el_p"
filenameps=TRIM(flpsanhar)//".macro_el_p"//TRIM(flext)

CALL gnuplot_start(gnu_filename)

IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ELSE
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ENDIF

CALL gnuplot_xlabel('T (K)', .FALSE.)

istep=0
DO ipressp=1, npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.macro_el_press'
   CALL add_value(filelastic,press(ipress))
   filelastic=TRIM(filelastic)//'_aver'
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.macro_el_s_press'
   CALL add_value(filelastic_s,press(ipress))
   filelastic_s=TRIM(filelastic_s)//'_aver'
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_ph_press'
   CALL add_value(filelastic_ph,press(ipress))
   filelastic_ph=TRIM(filelastic_ph)//'_aver'
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_s_ph_press'
   CALL add_value(filelastic_s_ph,press(ipress))
   filelastic_s_ph=TRIM(filelastic_s_ph)//'_aver'
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('Bulk Modulus (kbar)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,2,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s,1,2,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filelastic_ph,1,2,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,2,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

istep=0
DO ipressp=1, npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.macro_el_press'
   CALL add_value(filelastic,press(ipress))
   filelastic=TRIM(filelastic)//'_aver'
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.macro_el_s_press'
   CALL add_value(filelastic_s,press(ipress))
   filelastic_s=TRIM(filelastic_s)//'_aver'
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_ph_press'
   CALL add_value(filelastic_ph,press(ipress))
   filelastic_ph=TRIM(filelastic_ph)//'_aver'
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_s_ph_press'
   CALL add_value(filelastic_s_ph,press(ipress))
   filelastic_s_ph=TRIM(filelastic_s_ph)//'_aver'
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('Young modulus (kbar)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,3,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s,1,3,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filelastic_ph,1,3,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,3,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

istep=0
DO ipressp=1, npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.macro_el_press'
   CALL add_value(filelastic,press(ipress))
   filelastic=TRIM(filelastic)//'_aver'
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.macro_el_s_press'
   CALL add_value(filelastic_s,press(ipress))
   filelastic_s=TRIM(filelastic_s)//'_aver'
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_ph_press'
   CALL add_value(filelastic_ph,press(ipress))
   filelastic_ph=TRIM(filelastic_ph)//'_aver'
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_s_ph_press'
   CALL add_value(filelastic_s_ph,press(ipress))
   filelastic_s_ph=TRIM(filelastic_s_ph)//'_aver'
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('Shear modulus (kbar)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,4,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s,1,4,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filelastic_ph,1,4,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,4,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

istep=0
DO ipressp=1, npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.macro_el_press'
   CALL add_value(filelastic,press(ipress))
   filelastic=TRIM(filelastic)//'_aver'
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.macro_el_s_press'
   CALL add_value(filelastic_s,press(ipress))
   filelastic_s=TRIM(filelastic_s)//'_aver'
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_ph_press'
   CALL add_value(filelastic_ph,press(ipress))
   filelastic_ph=TRIM(filelastic_ph)//'_aver'
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_s_ph_press'
   CALL add_value(filelastic_s_ph,press(ipress))
   filelastic_s_ph=TRIM(filelastic_s_ph)//'_aver'
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('Longitudinal modulus (kbar)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data_linear(filelastic,1,2,4, &
            1.0_DP, 4.0_DP/3.0_DP, color(istep),first_step,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data_linear(filelastic_s,1,2,4,  &
                      1.0_DP, 4.0_DP/3.0_DP, color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data_linear(filelastic_ph,1,2,4, &
                   1.0_DP, 4.0_DP/3.0_DP, color(istep),            &
                   first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data_linear(filelastic_s_ph,1,2,4,&
                    1.0_DP, 4.0_DP/3.0_DP,color(istep),.FALSE.,last_step,&
                    .FALSE.)
   ENDIF
ENDDO

istep=0
DO ipressp=1, npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.macro_el_press'
   CALL add_value(filelastic,press(ipress))
   filelastic=TRIM(filelastic)//'_aver'
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.macro_el_s_press'
   CALL add_value(filelastic_s,press(ipress))
   filelastic_s=TRIM(filelastic_s)//'_aver'
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_ph_press'
   CALL add_value(filelastic_ph,press(ipress))
   filelastic_ph=TRIM(filelastic_ph)//'_aver'
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_s_ph_press'
   CALL add_value(filelastic_s_ph,press(ipress))
   filelastic_s_ph=TRIM(filelastic_s_ph)//'_aver'
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('Poisson ratio',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,5,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s,1,5,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filelastic_ph,1,5,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,5,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

istep=0
DO ipressp=1, npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.macro_el_press'
   CALL add_value(filelastic,press(ipress))
   filelastic=TRIM(filelastic)//'_aver'
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.macro_el_s_press'
   CALL add_value(filelastic_s,press(ipress))
   filelastic_s=TRIM(filelastic_s)//'_aver'
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_ph_press'
   CALL add_value(filelastic_ph,press(ipress))
   filelastic_ph=TRIM(filelastic_ph)//'_aver'
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_s_ph_press'
   CALL add_value(filelastic_s_ph,press(ipress))
   filelastic_s_ph=TRIM(filelastic_s_ph)//'_aver'
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('Pugh ratio',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data_div(filelastic,1,4,2,color(istep),&
                           first_step,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data_div(filelastic_s,1,4,2,color(istep),&
                          .FALSE., last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data_div(filelastic_ph,1,4,2,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data_div(filelastic_s_ph,1,4,2,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
           ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN

END SUBROUTINE plot_macro_el_new_pt

!-------------------------------------------------------------------------
SUBROUTINE plot_sound_pt()
!-------------------------------------------------------------------------

USE kinds,               ONLY : DP
USE temperature,         ONLY : tmin, tmax
USE data_files,          ONLY : flanhar
USE postscript_files,    ONLY : flpsanhar
USE control_elastic_constants, ONLY : lelastic, lelasticf
USE control_pressure,    ONLY : npress_plot, ipress_plot, press
USE control_gnuplot,     ONLY : gnuplot_command, flgnuplot, lgnuplot, flext
USE gnuplot,             ONLY : gnuplot_start, gnuplot_end,           &
                                gnuplot_write_header, gnuplot_xlabel, &
                                gnuplot_ylabel, gnuplot_write_file_mul_data, &
                                gnuplot_set_fact
USE color_mod,           ONLY : color
USE io_global,           ONLY : ionode
USE mp_images,           ONLY : root_image, my_image_id

IMPLICIT NONE

CHARACTER(LEN=256) :: filelastic, filelastic_s, filelastic_ph, &
                      filelastic_s_ph, filedebye, filedebye_ph, &
                      gnu_filename, filenameps
LOGICAL :: first_step, last_step
INTEGER :: istep, ipressp, ipress, ierr, system

IF (npress_plot==0) RETURN
IF (.NOT.(lelastic.OR.lelasticf)) RETURN

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_sound_vel_p"
filenameps=TRIM(flpsanhar)//".sound_vel_p"//TRIM(flext)

CALL gnuplot_start(gnu_filename)

IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ELSE
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ENDIF

CALL gnuplot_xlabel('T (K)', .FALSE.)

istep=0
DO ipressp=1, npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.sound_vel_press'
   CALL add_value(filelastic,press(ipress))
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.sound_vel_s_press'
   CALL add_value(filelastic_s,press(ipress))
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.sound_vel_ph_press'
   CALL add_value(filelastic_ph,press(ipress))
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.sound_vel_s_ph_press'
   CALL add_value(filelastic_s_ph,press(ipress))
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('V_{P} (m/s)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,2,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s,1,2,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filelastic_ph,1,2,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,2,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

istep=0
DO ipressp=1, npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.sound_vel_press'
   CALL add_value(filelastic,press(ipress))
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.sound_vel_s_press'
   CALL add_value(filelastic_s,press(ipress))
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.sound_vel_ph_press'
   CALL add_value(filelastic_ph,press(ipress))
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.sound_vel_s_ph_press'
   CALL add_value(filelastic_s_ph,press(ipress))
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('V_{B} (m/s)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,3,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s,1,3,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filelastic_ph,1,3,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,3,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

istep=0
DO ipressp=1, npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.sound_vel_press'
   CALL add_value(filelastic,press(ipress))
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.sound_vel_s_press'
   CALL add_value(filelastic_s,press(ipress))
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.sound_vel_ph_press'
   CALL add_value(filelastic_ph,press(ipress))
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.sound_vel_s_ph_press'
   CALL add_value(filelastic_s_ph,press(ipress))
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('V_{G} (m/s)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,4,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s,1,4,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filelastic_ph,1,4,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,4,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

istep=0
DO ipressp=1, npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filedebye="anhar_files/"//TRIM(flanhar)//'.macro_el_debye_press'
   CALL add_value(filedebye,press(ipress))
   filedebye_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_debye_ph_press'
   CALL add_value(filedebye_ph,press(ipress))
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('{/Symbol Q}_D (K)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filedebye,1,2,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filedebye,1,3,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filedebye_ph,1,2,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filedebye_ph,1,3,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
           ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN

END SUBROUTINE plot_sound_pt

!-------------------------------------------------------------------------
SUBROUTINE plot_macro_el_new_ptt()
!-------------------------------------------------------------------------

USE kinds,               ONLY : DP
USE constants,           ONLY : ry_kbar
USE temperature,         ONLY : ntemp_plot, itemp_plot, temp
USE data_files,          ONLY : flanhar
USE postscript_files,    ONLY : flpsanhar
USE control_elastic_constants, ONLY : lelastic, lelasticf
USE control_pressure,    ONLY : pmin, pmax
USE control_gnuplot,     ONLY : gnuplot_command, flgnuplot, lgnuplot, flext
USE gnuplot,             ONLY : gnuplot_start, gnuplot_end,           &
                                gnuplot_write_header, gnuplot_xlabel, &
                                gnuplot_ylabel, gnuplot_write_file_mul_data, &
                                gnuplot_write_file_mul_data_div,    &
                                gnuplot_write_file_mul_data_linear, &
                                gnuplot_set_fact
USE color_mod,           ONLY : color
USE io_global,           ONLY : ionode
USE mp_images,           ONLY : root_image, my_image_id

IMPLICIT NONE

CHARACTER(LEN=256) :: filelastic, filelastic_s, &
                      filelastic_ph, filelastic_s_ph, gnu_filename, filenameps
LOGICAL :: first_step, last_step
INTEGER :: istep, itempp, itemp, ierr, system

IF ( my_image_id /= root_image ) RETURN
IF (ntemp_plot==0) RETURN
IF (.NOT.(lelastic.OR.lelasticf)) RETURN

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_macro_el_t"
filenameps=TRIM(flpsanhar)//".macro_el_t"//TRIM(flext)

CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filenameps, pmin*ry_kbar, pmax*ry_kbar, 0.0_DP, &
                              0.0_DP, 1.0_DP, flext )
CALL gnuplot_xlabel('p (kbar)', .FALSE.)

istep=0
DO itempp=1, ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.macro_el_temp'
   CALL add_value(filelastic,temp(itemp))
   filelastic=TRIM(filelastic)//'_aver'
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.macro_el_s_temp'
   CALL add_value(filelastic_s,temp(itemp))
   filelastic_s=TRIM(filelastic_s)//'_aver'
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_ph_temp'
   CALL add_value(filelastic_ph,temp(itemp))
   filelastic_ph=TRIM(filelastic_ph)//'_aver'
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_s_ph_temp'
   CALL add_value(filelastic_s_ph,temp(itemp))
   filelastic_s_ph=TRIM(filelastic_s_ph)//'_aver'
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('Bulk Modulus (kbar)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,2,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s,1,2,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filelastic_ph,1,2,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,2,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO
istep=0
DO itempp=1, ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.macro_el_temp'
   CALL add_value(filelastic,temp(itemp))
   filelastic=TRIM(filelastic)//'_aver'
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.macro_el_s_temp'
   CALL add_value(filelastic_s,temp(itemp))
   filelastic_s=TRIM(filelastic_s)//'_aver'
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_ph_temp'
   CALL add_value(filelastic_ph,temp(itemp))
   filelastic_ph=TRIM(filelastic_ph)//'_aver'
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_s_ph_temp'
   CALL add_value(filelastic_s_ph,temp(itemp))
   filelastic_s_ph=TRIM(filelastic_s_ph)//'_aver'
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('Young modulus (kbar)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,3,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s,1,3,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filelastic_ph,1,3,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,3,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

istep=0
DO itempp=1, ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.macro_el_temp'
   CALL add_value(filelastic,temp(itemp))
   filelastic=TRIM(filelastic)//'_aver'
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.macro_el_s_temp'
   CALL add_value(filelastic_s,temp(itemp))
   filelastic_s=TRIM(filelastic_s)//'_aver'
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_ph_temp'
   CALL add_value(filelastic_ph,temp(itemp))
   filelastic_ph=TRIM(filelastic_ph)//'_aver'
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_s_ph_temp'
   CALL add_value(filelastic_s_ph,temp(itemp))
   filelastic_s_ph=TRIM(filelastic_s_ph)//'_aver'
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('Shear modulus (kbar)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,4,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s,1,4,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filelastic_ph,1,4,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,4,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

istep=0
DO itempp=1, ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.macro_el_temp'
   CALL add_value(filelastic,temp(itemp))
   filelastic=TRIM(filelastic)//'_aver'
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.macro_el_s_temp'
   CALL add_value(filelastic_s,temp(itemp))
   filelastic_s=TRIM(filelastic_s)//'_aver'
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_ph_temp'
   CALL add_value(filelastic_ph,temp(itemp))
   filelastic_ph=TRIM(filelastic_ph)//'_aver'
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_s_ph_temp'
   CALL add_value(filelastic_s_ph,temp(itemp))
   filelastic_s_ph=TRIM(filelastic_s_ph)//'_aver'
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('Longitudinal modulus (kbar)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data_linear(filelastic,1,2,4,       &
          1.0_DP, 4.0_DP/3.0_DP, color(istep),first_step,.FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data_linear(filelastic_s,1,2,4,     &
          1.0_DP, 4.0_DP/3.0_DP,color(istep),.FALSE.,                 &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data_linear(filelastic_ph,1,2,4,    &
         1.0_DP, 4.0_DP/3.0_DP, color(istep),                         &
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data_linear(filelastic_s_ph,1,2,4,  &
         1.0_DP, 4.0_DP/3.0_DP,color(istep),.FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

istep=0
DO itempp=1, ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.macro_el_temp'
   CALL add_value(filelastic,temp(itemp))
   filelastic=TRIM(filelastic)//'_aver'
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.macro_el_s_temp'
   CALL add_value(filelastic_s,temp(itemp))
   filelastic_s=TRIM(filelastic_s)//'_aver'
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_ph_temp'
   CALL add_value(filelastic_ph,temp(itemp))
   filelastic_ph=TRIM(filelastic_ph)//'_aver'
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_s_ph_temp'
   CALL add_value(filelastic_s_ph,temp(itemp))
   filelastic_s_ph=TRIM(filelastic_s_ph)//'_aver'
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('Poisson ratio',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,5,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s,1,5,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filelastic_ph,1,5,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,5,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

istep=0
DO itempp=1, ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.macro_el_temp'
   CALL add_value(filelastic,temp(itemp))
   filelastic=TRIM(filelastic)//'_aver'
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.macro_el_s_temp'
   CALL add_value(filelastic_s,temp(itemp))
   filelastic_s=TRIM(filelastic_s)//'_aver'
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_ph_temp'
   CALL add_value(filelastic_ph,temp(itemp))
   filelastic_ph=TRIM(filelastic_ph)//'_aver'
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_s_ph_temp'
   CALL add_value(filelastic_s_ph,temp(itemp))
   filelastic_s_ph=TRIM(filelastic_s_ph)//'_aver'
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('Pugh ratio',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data_div(filelastic,1,4,2,color(istep), & 
                               first_step,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data_div(filelastic_s,1,4,2,color(istep), &
                              .FALSE.,last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data_div(filelastic_ph,1,4,2,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data_div(filelastic_s_ph,1,4,2,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
           ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN

END SUBROUTINE plot_macro_el_new_ptt

!-------------------------------------------------------------------------
SUBROUTINE plot_sound_ptt()
!-------------------------------------------------------------------------

USE kinds,               ONLY : DP
USE constants,           ONLY : ry_kbar
USE temperature,         ONLY : ntemp_plot, itemp_plot, temp
USE data_files,          ONLY : flanhar
USE postscript_files,    ONLY : flpsanhar
USE control_elastic_constants, ONLY : lelastic, lelasticf
USE control_pressure,    ONLY : pmin, pmax
USE control_gnuplot,     ONLY : gnuplot_command, flgnuplot, lgnuplot, flext
USE gnuplot,             ONLY : gnuplot_start, gnuplot_end,           &
                                gnuplot_write_header, gnuplot_xlabel, &
                                gnuplot_ylabel, gnuplot_write_file_mul_data, &
                                gnuplot_set_fact
USE color_mod,           ONLY : color
USE io_global,           ONLY : ionode
USE mp_images,           ONLY : root_image, my_image_id

IMPLICIT NONE

CHARACTER(LEN=256) :: filelastic, filelastic_s, filelastic_ph, &
                      filelastic_s_ph, filedebye, filedebye_ph, &
                      gnu_filename, filenameps
LOGICAL :: first_step, last_step
INTEGER :: istep, itempp, itemp, ierr, system

IF ( my_image_id /= root_image ) RETURN
IF (ntemp_plot==0) RETURN
IF (.NOT.(lelastic.OR.lelasticf)) RETURN

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_sound_vel_t"
filenameps=TRIM(flpsanhar)//".sound_vel_t"//TRIM(flext)

CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filenameps, pmin*ry_kbar, pmax*ry_kbar, 0.0_DP, &
                              0.0_DP, 1.0_DP, flext )
CALL gnuplot_xlabel('p (kbar)', .FALSE.)

istep=0
DO itempp=1, ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.sound_vel_temp'
   CALL add_value(filelastic,temp(itemp))
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.sound_vel_s_temp'
   CALL add_value(filelastic_s,temp(itemp))
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.sound_vel_ph_temp'
   CALL add_value(filelastic_ph,temp(itemp))
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.sound_vel_s_ph_temp'
   CALL add_value(filelastic_s_ph,temp(itemp))
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('V_{P} (m/s)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,2,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s,1,2,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filelastic_ph,1,2,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,2,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO
istep=0
DO itempp=1, ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.sound_vel_temp'
   CALL add_value(filelastic,temp(itemp))
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.sound_vel_s_temp'
   CALL add_value(filelastic_s,temp(itemp))
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.sound_vel_ph_temp'
   CALL add_value(filelastic_ph,temp(itemp))
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.sound_vel_s_ph_temp'
   CALL add_value(filelastic_s_ph,temp(itemp))
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('V_{B} (m/s)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,3,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s,1,3,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filelastic_ph,1,3,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,3,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

istep=0
DO itempp=1, ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filelastic="anhar_files/"//TRIM(flanhar)//'.sound_vel_temp'
   CALL add_value(filelastic,temp(itemp))
   filelastic_s="anhar_files/"//TRIM(flanhar)//'.sound_vel_s_temp'
   CALL add_value(filelastic_s,temp(itemp))
   filelastic_ph="anhar_files/"//TRIM(flanhar)//'.sound_vel_ph_temp'
   CALL add_value(filelastic_ph,temp(itemp))
   filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.sound_vel_s_ph_temp'
   CALL add_value(filelastic_s_ph,temp(itemp))
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('V_{G} (m/s)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,4,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s,1,4,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filelastic_ph,1,4,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,4,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

istep=0
DO itempp=1, ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filedebye="anhar_files/"//TRIM(flanhar)//'.macro_el_debye_temp'
   CALL add_value(filedebye,temp(itemp))
   filedebye_ph="anhar_files/"//TRIM(flanhar)//'.macro_el_debye_ph_temp'
   CALL add_value(filedebye_ph,temp(itemp))
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP, .FALSE.)
      CALL gnuplot_ylabel('{/Symbol Q}_{D} (K)',.FALSE.)
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filedebye,1,2,color(istep),first_step,&
                                                             .FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filedebye,1,3,color(istep),.FALSE., &
                                       last_step.AND..NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filedebye_ph,1,2,color(istep),&
                              first_step.AND..NOT.lelastic,.FALSE.,.FALSE.)

      CALL gnuplot_write_file_mul_data(filedebye_ph,1,3,color(istep),&
                                               .FALSE.,last_step,.FALSE.)
   ENDIF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
           ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN

END SUBROUTINE plot_sound_ptt
