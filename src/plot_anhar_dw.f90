!
! Copyright (C) 2018 Cristiano Malica
!               2022 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE plot_anhar_thermo()
!-----------------------------------------------------------------------
!
!  This is a driver to plot the quantities written inside flanhar.therm
!  and flanhar.therm_ph in the directory anhar_files
!  
USE kinds,            ONLY : DP
USE constants,        ONLY : rydberg_si, avogadro
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE control_thermo,   ONLY : ltherm_dos, ltherm_freq
USE control_pressure, ONLY : npress_plot, ipress_plot, press
USE postscript_files, ONLY : flpsanhar
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                             gnuplot_ylabel, &
                             gnuplot_xlabel, &
                             gnuplot_write_file_mul_data, &
                             gnuplot_set_fact
USE data_files,      ONLY : flanhar 
USE color_mod,       ONLY : color
USE temperature,     ONLY : tmin, tmax
USE mp_images,       ONLY : root_image, my_image_id
USE io_global,       ONLY : ionode

IMPLICIT NONE
INTEGER :: ierr, system
CHARACTER(LEN=256) :: gnu_filename, filename, filename1, filename1_ph, &
                      filetherm, filepstherm
LOGICAL :: first_step, last_step
INTEGER :: istep, ipressp, ipress
REAL(DP) :: factor

IF ( my_image_id /= root_image ) RETURN

IF (.NOT.(ltherm_freq.OR.ltherm_dos)) RETURN

filetherm="anhar_files/"//TRIM(flanhar)//'.therm'
CALL add_pressure(filetherm)
filename="anhar_files/"//TRIM(flanhar)//'.therm_ph'
CALL add_pressure(filename)

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_anhar_therm'
CALL add_pressure(gnu_filename)

CALL gnuplot_start(gnu_filename)

filepstherm=TRIM(flpsanhar)//'.therm'//TRIM(flext)
CALL add_pressure(filepstherm)

IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(filepstherm, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                        1.0_DP, flext )
ELSE
   CALL gnuplot_write_header(filepstherm, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                        1.0_DP, flext )
ENDIF

CALL gnuplot_xlabel('T (K)', .FALSE.)
CALL gnuplot_ylabel('Thermal Energy (kJ / (N mol))',.FALSE.)
factor = rydberg_si*avogadro / 1.D3
CALL gnuplot_set_fact(factor, .FALSE.)

IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filetherm,1,2,'color_red',.TRUE.,&
                                                     .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,2,'color_blue', &
                                                .NOT.ltherm_dos,.TRUE.,.FALSE.)
!
! ADC addition
!
istep=0
factor = rydberg_si*avogadro/1.D3 
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename1="anhar_files/"//TRIM(flanhar)//'.therm_press'
   CALL add_value(filename1,press(ipress))
   filename1_ph="anhar_files/"//TRIM(flanhar)//'.therm_ph_press'
   CALL add_value(filename1_ph,press(ipress))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      CALL gnuplot_set_fact(factor,.FALSE.)
      CALL gnuplot_ylabel('Thermal Energy (kJ / (N mol))',.FALSE.)
   ENDIF
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename1,1,2,color(istep),first_step, &
                            (last_step.AND..NOT.ltherm_freq), .FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data(filename1_ph,1,2,color(istep), &
                  (first_step.AND..NOT.ltherm_dos), last_step, .FALSE.)
ENDDO
!
! End ADC addition
!

CALL gnuplot_ylabel('Vibrational free energy (kJ / (N mol))', .FALSE.)
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filetherm,1,3,'color_red',.TRUE.,&
                                                     .NOT.ltherm_freq, .FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,3,'color_blue',&
                                                .NOT.ltherm_dos,.TRUE.,.FALSE.)
!
! ADC addition
!
istep=0
factor = rydberg_si*avogadro/1.D3 
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename1="anhar_files/"//TRIM(flanhar)//'.therm_press'
   CALL add_value(filename1,press(ipress))
   filename1_ph="anhar_files/"//TRIM(flanhar)//'.therm_ph_press'
   CALL add_value(filename1_ph,press(ipress))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      CALL gnuplot_set_fact(factor,.FALSE.)
      CALL gnuplot_ylabel('Thermal Free Energy (kJ / (N mol))',.FALSE.)
   ENDIF
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename1,1,3,color(istep),first_step, &
                          (last_step.AND..NOT.ltherm_freq), .FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data(filename1_ph,1,3,color(istep),  &
                          (first_step.AND..NOT.ltherm_dos), last_step, .FALSE.)
ENDDO
!
! End ADC addition
!
 
CALL gnuplot_set_fact(factor*1.D3, .FALSE.)
CALL gnuplot_ylabel('Entropy (J / K / (N mol))',.FALSE.)
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filetherm,1,4,'color_red',.TRUE., &
                                                   .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',.NOT.ltherm_dos,&
                                                             .TRUE.,.FALSE.)
!
! ADC addition
!
istep=0
factor = rydberg_si*avogadro 
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename1="anhar_files/"//TRIM(flanhar)//'.therm_press'
   CALL add_value(filename1,press(ipress))
   filename1_ph="anhar_files/"//TRIM(flanhar)//'.therm_ph_press'
   CALL add_value(filename1_ph,press(ipress))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      CALL gnuplot_set_fact(factor,.FALSE.)
      CALL gnuplot_ylabel('Entropy (J / K / (N mol))',.FALSE.)
   ENDIF
   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename1,1,4,color(istep),first_step, &
                  (last_step.AND..NOT.ltherm_freq), .FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data(filename1_ph,1,4,color(istep), &
                  (first_step.AND..NOT.ltherm_dos), last_step, .FALSE.)
ENDDO
!
! End ADC addition
!
CALL gnuplot_ylabel('Heat capacity C_v (J / K / (N mol))',.FALSE.)
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filetherm,1,5,'color_red',.TRUE.,&
                 .NOT.ltherm_freq,.FALSE.)

IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,5,'color_blue',.NOT.ltherm_dos,&
                                                           .TRUE.,.FALSE.)
CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_anhar_thermo

!-----------------------------------------------------------------------
SUBROUTINE plot_anhar_dw()
!-----------------------------------------------------------------------
!
!  This is a driver to plot the quantities written inside 
!  
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

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_anhar_dw'
CALL gnuplot_start(gnu_filename)

psfilename=TRIM(flpsanhar)//'.dw'//TRIM(flext)
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
END SUBROUTINE plot_anhar_dw
