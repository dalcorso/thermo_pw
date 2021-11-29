!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE plot_anhar()
!-----------------------------------------------------------------------
!
!  This is a driver to plot the quantities written inside flanhar,
!  flanhar//'_ph', flanhar//'.bulk_mod', flanhar//'.bulk_mod_ph', 
!  flanhar//'.dbulk_mod', flanhar//'.dbulk_mod_ph', flanhar//'.heat', 
!  flanhar//'.heat_ph', flanhar//'.aux_grun', flanhar//'.gamma', 
!  flanhar//'.gamma_ph', flanhar//'.gamma_grun'.
!  
!
USE kinds,            ONLY : DP
USE constants,        ONLY : rydberg_si, avogadro
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

CHARACTER(LEN=256) :: gnu_filename, filename, filename0, filename_aux_grun,  &
                      filename_bulk, filename_bulk_ph,               &
                      filename_dbulk, filename_dbulk_ph,             &
                      filename_heat, filename_heat_ph,               &
                      filename_gamma, filename_gamma_ph, filename_gamma_grun 
INTEGER :: ierr, system
REAL(DP) :: factor

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

filename0="anhar_files/"//TRIM(flanhar)
CALL add_pressure(filename0)
filename="anhar_files/"//TRIM(flanhar)//'_ph'
CALL add_pressure(filename)
filename_aux_grun="anhar_files/"//TRIM(flanhar)//'.aux_grun'
CALL add_pressure(filename_aux_grun)
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

filename_gamma="anhar_files/"//TRIM(flanhar)//'.gamma'
CALL add_pressure(filename_gamma)
filename_gamma_ph="anhar_files/"//TRIM(flanhar)//'.gamma_ph'
CALL add_pressure(filename_gamma_ph)
filename_gamma_grun="anhar_files/"//TRIM(flanhar)//'.gamma_grun'
CALL add_pressure(filename_gamma_grun)
!
!  Part 1: Volume
!
CALL gnuplot_xlabel('T (K)',.FALSE.) 
CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('Volume ((a.u.)^3)',.FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename0,1,2,'color_red',.TRUE.,&
                                              .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,2,'color_blue',&
                                        .NOT.ltherm_dos,.TRUE.,.FALSE.)

!
!  Part 2: Helmholtz (or Gibbs) free energy
!
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
!
!  Part 3: bulk modulus
!
CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('Bulk modulus (kbar)',.FALSE.) 

IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename_bulk,1,2,'color_red',.TRUE., &
                                                .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename_bulk_ph,1,2,'color_blue',&
                                                .NOT.ltherm_dos,.TRUE.,.FALSE.)
!
!  Part 4: pressure derivative of the bulk modulus
!
CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('d B / d p',.FALSE.) 

IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename_dbulk,1,2,'color_red',.TRUE., &
                                                .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename_dbulk_ph,1,2,'color_blue',&
                                               .NOT.ltherm_dos,.TRUE.,.FALSE.)
!
!  Part 5: Volume thermal expansion
!
CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('Thermal expansion ({/Symbol b} x 10^{6}) (K^{-1})',.FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename0,1,4,'color_red',.TRUE.,&
                                                            .FALSE.,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',&
                                            .NOT.ltherm_dos,.FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filename_aux_grun,1,2,'color_green',&
                                .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)
!
!  put as a comment the possibility to plot also the experimental data
!
CALL gnuplot_write_file_mul_data(filename_aux_grun,1,2,'color_green',&
                                               .NOT.ltherm,.FALSE.,.TRUE.)
CALL gnuplot_write_file_mul_point('anhar.exp',1,2,'color_red',.FALSE.,&
                                                            .TRUE.,.TRUE.)
!
!  Part 6: isochoric heat capacity
!
factor = rydberg_si*avogadro 
CALL gnuplot_set_fact(factor,.FALSE.)
CALL gnuplot_ylabel('Heat capacity C_v (J / K / N / mol)',.FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename_heat,1,2,'color_red',.TRUE.,&
                                            .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename_heat_ph,1,2,'color_blue',&
                                            .NOT.ltherm_dos,.TRUE.,.FALSE.)
!
!  Part 7: isobaric heat capacity
!
CALL gnuplot_set_fact(factor,.FALSE.)
CALL gnuplot_ylabel('Heat capacity C_p (J / K / N / mol)',.FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename_heat,1,4,'color_red',.TRUE.,&
                                          .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename_heat_ph,1,4,'color_blue', &
                                    .NOT.ltherm_dos,.TRUE.,.FALSE.)
!
!  put as a comment the possibility to plot also the experimental data
!
   CALL gnuplot_write_file_mul_data(filename_heat,1,4,'color_blue',&
                                       .NOT.ltherm_dos,.FALSE.,.TRUE.)
   CALL gnuplot_write_file_mul_point('cv.exp',1,2,'color_red',.FALSE.,&
                                                              .TRUE.,.TRUE.)
!
!  Part 8: isobaric -isochoric heat capacity
!
CALL gnuplot_set_fact(factor,.FALSE.)
CALL gnuplot_ylabel('C_p - C_v (J / K / N / mol)',.FALSE.) 

IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename_heat,1,3,'color_red',.TRUE., &
                                                   .FALSE.,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename_heat_ph,1,3,'color_blue',&
                                               .NOT.ltherm_dos,.FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filename_aux_grun,1,3,'color_green', &
                              .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)
!
!  Part 8: isoentropic-isothermal bulk modulus
!
CALL gnuplot_set_fact(1._DP,.FALSE.)
CALL gnuplot_ylabel('B_S - B_T (kbar)',.FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename_bulk,1,4,'color_red',.TRUE.,&
                                             .FALSE.,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename_bulk_ph,1,4,'color_blue',&
                                           .NOT.ltherm_dos,.FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filename_aux_grun,1,4,'color_green', &
                              .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)

!
!  Part 10: average gruneisen parameter
!
CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('Gr\374neisen parameter ({/Symbol g})',.FALSE.) 
CALL gnuplot_write_horizontal_line(0.0_DP, 2, 'front', 'color_black', .FALSE.)
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename_gamma,1,2,'color_red',.TRUE.,&
                                                     .FALSE.,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename_gamma_ph,1,2,'color_blue',&
                         .NOT.ltherm_dos,.FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filename_gamma_grun,1,2,'color_green', &
                     .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

CALL plot_thermo_anhar()

IF (with_eigen) CALL plot_dw_anhar()

RETURN
END SUBROUTINE plot_anhar
!
!-----------------------------------------------------------------------
SUBROUTINE plot_anhar_p()
!-----------------------------------------------------------------------
!
! This routine plots a few quantities dependent on the temperature
! for several pressures
!
USE kinds,            ONLY : DP
USE constants,        ONLY : rydberg_si, avogadro
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_write_file_mul_data, &
                             gnuplot_write_command,       &
                             gnuplot_set_fact
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE control_ev,       ONLY : ieos
USE data_files,       ONLY : flanhar
USE postscript_files, ONLY : flpsanhar
USE temperature,      ONLY : tmin, tmax
USE control_pressure, ONLY : press, npress_plot, ipress_plot
USE color_mod,        ONLY : color
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: gnu_filename, filename
CHARACTER(LEN=8) :: float_to_char
REAL(DP) :: factor
INTEGER  :: istep, ipress, ipressp
INTEGER  :: ierr, system
LOGICAL  :: first_step, last_step
!
IF ( my_image_id /= root_image ) RETURN
!
gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_anhar_press'
CALL gnuplot_start(gnu_filename)

filename=TRIM(flpsanhar)//'.press'
filename=TRIM(filename)//TRIM(flext)
IF (tmin /= 1.0_DP) THEN
   CALL gnuplot_write_header(filename, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext )
ELSE
   CALL gnuplot_write_header(filename, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext )
ENDIF
!
!  V(T)
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.press.'//&
                            TRIM(float_to_char(press(ipress),1))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_ylabel('V(T) ((a.u.)^3)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step, &
                                                        last_step, .FALSE.)
ENDDO
!
!  V(T)/V(T=300 K)
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.press.'//&
                            TRIM(float_to_char(press(ipress),1))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_ylabel('V(T)/V(T=300 K)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,3,color(istep),first_step, &
                                                        last_step, .FALSE.)
ENDDO
!
!   The bulk modulus
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.bulk_press.'//&
                            TRIM(float_to_char(press(ipress),1))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_ylabel('B_T(T) (kbar))',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step, &
                                                        last_step, .FALSE.)
ENDDO
!
!   The derivative of the bulk modulus with respect to pressure
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.dbulk_press.'//&
                            TRIM(float_to_char(press(ipress),1))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_ylabel('d B_T/dp (T)',.FALSE.)
      CALL gnuplot_write_command('set yrange [3.0:5.0]',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step, &
                                                        last_step, .FALSE.)
ENDDO
CALL gnuplot_write_command('unset yrange',.FALSE.)
!
!   The second derivative of the bulk modulus with respect to pressure
!
IF (ieos==2) THEN
   istep=0
   DO ipressp=1,npress_plot
      first_step=(ipressp==1)
      last_step=(ipressp==npress_plot)
      ipress=ipress_plot(ipressp)
      istep=MOD(istep,8)+1
      filename="anhar_files/"//TRIM(flanhar)//'.dbulk_press.'//&
                               TRIM(float_to_char(press(ipress),1))
      IF (first_step) THEN
         CALL gnuplot_xlabel('T (K)',.FALSE.)
         CALL gnuplot_set_fact(1.0_DP,.FALSE.)
         CALL gnuplot_ylabel('d^2 B_T/dp^2 (T) (1/kbar)',.FALSE.)
         CALL gnuplot_write_command('set yrange [-0.1:0.1]',.FALSE.)
      ENDIF
      CALL gnuplot_write_file_mul_data(filename,1,3,color(istep),first_step, &
                                                        last_step, .FALSE.)
   ENDDO
   CALL gnuplot_write_command('unset yrange',.FALSE.)
ENDIF
!
!  Thermal expansion
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.press.'//& 
                   TRIM(float_to_char(press(ipress),1))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_ylabel('Thermal expansion ({/Symbol b} x 10^{6}) &
                                                       &(K^{-1})',.FALSE.)

   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,5,color(istep),first_step, &
                                                        last_step,.FALSE.)
ENDDO
!
!  C_V
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.heat_press.'//&
                            TRIM(float_to_char(press(ipress),1))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      factor = rydberg_si*avogadro
      CALL gnuplot_set_fact(factor,.FALSE.)
      CALL gnuplot_ylabel('Heat capacity C_v (J / K / N / mol)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step, &
                                                        last_step,.FALSE.)
ENDDO
!
!  C_P
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.heat_press.'//&
                      TRIM(float_to_char(press(ipress),1))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      factor = rydberg_si*avogadro
      CALL gnuplot_set_fact(factor,.FALSE.)
      CALL gnuplot_ylabel('Heat capacity C_p (J / K / N / mol)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,4,color(istep),first_step, &
                                                        last_step,.FALSE.)
ENDDO
!
!  C_P - C_V
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.heat_press.'//&
                      TRIM(float_to_char(press(ipress),1))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      factor = rydberg_si*avogadro
      CALL gnuplot_set_fact(factor,.FALSE.)
      CALL gnuplot_ylabel('Heat capacity C_p -C_V(J / K / N / mol)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,3,color(istep),first_step, &
                                                        last_step,.FALSE.)
ENDDO
!
!  B_S - B_T
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.bulk_press.'//&
                      TRIM(float_to_char(press(ipress),1))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      factor = 1.0_DP
      CALL gnuplot_set_fact(factor,.FALSE.)
      CALL gnuplot_ylabel('B_S - B_T (kbar)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,4,color(istep),first_step, &
                                                        last_step,.FALSE.)
ENDDO
!
!   The average Gruneisen parameter 
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.gamma_press.'//&
                      TRIM(float_to_char(press(ipress),1))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      factor = 1.0_DP
      CALL gnuplot_set_fact(factor,.FALSE.)
      CALL gnuplot_ylabel('Gr\374neisen parameter ({/Symbol g})',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step, &
                                                        last_step,.FALSE.)
ENDDO
!
CALL gnuplot_end()
!
IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))
!
RETURN
END SUBROUTINE plot_anhar_p

!-------------------------------------------------------------------
SUBROUTINE plot_anhar_t()
!-------------------------------------------------------------------
!
!  This is a driver to plot the energy, the pressure, the thermal pressure,
!  the Gibbs free energy, the vibrational contribution to the free energy 
!  as a function of volume for several temperatures.
!  It writes also the volume thermal expansion, the average Gruneisen 
!  parameter and the product of the volume thermal expansion and the 
!  bulk modulus as a function of pressure for several temperatures.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar, rytoev
USE control_gnuplot,  ONLY : flgnuplot, lgnuplot, gnuplot_command, flext
USE postscript_files, ONLY : flpsanhar
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_set_fact,            &
                             gnuplot_write_command,       &
                             gnuplot_write_horizontal_line,     &
                             gnuplot_write_file_mul_data,       &
                             gnuplot_write_file_mul_line_point, &
                             gnuplot_write_file_mul_point_sum,  &
                             gnuplot_write_file_mul_data_diff,  &
                             gnuplot_write_file_mul_point
USE control_ev,       ONLY : ieos
USE data_files,       ONLY : flevdat, flanhar
USE control_mur,      ONLY : lmurn
USE control_vol,      ONLY : vmin_input, vmax_input
USE control_eldos,    ONLY : lel_free_energy
USE temperature,      ONLY : temp, ntemp, ntemp_plot, itemp_plot
USE control_pressure, ONLY : pressure_kb, pmax, pmin
USE color_mod,        ONLY : color
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: filename, filename1, gnu_filename, label
CHARACTER(LEN=8) :: float_to_char
INTEGER :: itemp, itempp, istep
INTEGER :: ierr, system
LOGICAL :: first_step, last_step

IF ( my_image_id /= root_image ) RETURN

IF (ntemp_plot==0) RETURN

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_anhar_temp'
filename=TRIM(flpsanhar)//'.temp'
filename=TRIM(filename)//TRIM(flext)

CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filename, vmin_input, vmax_input, &
                                0.0_DP, 0.0_DP, 1.0_DP, flext ) 
!
!  Energy or enthalpy as a function of the volume.
!
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flevdat)//'_mur.'//&
                                     TRIM(float_to_char(temp(itemp),1))
   CALL add_pressure(filename)
   IF (first_step) THEN
      CALL gnuplot_xlabel('Volume ((a.u.)^3)',.FALSE.) 
      IF (pressure_kb /= 0.0_DP) THEN
         label='Gibbs free-energy (Ry)    p= '//&
                    &TRIM(float_to_char(pressure_kb,1))//' kbar'
         CALL gnuplot_ylabel(TRIM(label),.FALSE.) 
      ELSE
         CALL gnuplot_ylabel('Helmholtz Free Energy (Ry)',.FALSE.) 
      END IF
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step, &
                                                        last_step, .FALSE.)
ENDDO
!
!  Pressure as a function of the volume
!
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1

   filename="anhar_files/"//TRIM(flevdat)//'_mur.'//&
                                    TRIM(float_to_char(temp(itemp),1))
   CALL add_pressure(filename)
   IF (first_step) THEN
      CALL gnuplot_xlabel('Volume ((a.u.)^3)',.FALSE.) 
      CALL gnuplot_ylabel('Pressure (kbar)',.FALSE.) 
      CALL gnuplot_write_horizontal_line(0.0_DP, 2, 'front', 'color_black', &
                                          .FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,4,color(istep),first_step, &
                                                        last_step, .FALSE.)
ENDDO
!
!  The thermal pressure as a function of the volume
!
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1

   filename="anhar_files/"//TRIM(flevdat)//'_mur.'//&
                                   TRIM(float_to_char(temp(itemp),1))
   CALL add_pressure(filename)
   IF (first_step) THEN
      CALL gnuplot_xlabel('Volume ((a.u.)^3)',.FALSE.) 
      CALL gnuplot_ylabel('Thermal pressure (kbar)',.FALSE.) 
   ENDIF
   CALL gnuplot_write_file_mul_data_diff(filename,1,4,5,color(istep),&
                             first_step, last_step, .FALSE.)
ENDDO

!
!   Vibrational free energy (+ electronic if computed) as a function of the 
!   volume
!
istep=0
CALL gnuplot_set_fact(rytoev, .FALSE.)
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1

   filename="anhar_files/"//TRIM(flevdat)//'_free.'//&
                                          TRIM(float_to_char(temp(itemp),1))
   CALL add_pressure(filename)
   filename1="anhar_files/"//TRIM(flevdat)//'_poly_free.'//&
                                          TRIM(float_to_char(temp(itemp),1))
   IF (first_step) THEN
      CALL gnuplot_xlabel('Volume ((a.u.)^3)',.FALSE.) 
      IF (lel_free_energy) THEN
         CALL gnuplot_ylabel('Vibrat. + elec. free energy (eV)',.FALSE.) 
      ELSE
         CALL gnuplot_ylabel('Vibrational free energy (eV)',.FALSE.) 
      ENDIF
   ENDIF
   IF (lel_free_energy) THEN
      CALL gnuplot_write_file_mul_point_sum(filename,1,3,4,color(istep), &
                                        first_step, .FALSE., .FALSE.)
   ELSE
      CALL gnuplot_write_file_mul_point(filename,1,3,color(istep), &
                                        first_step, .FALSE., .FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename1,1,2,color(istep), &
                                        .FALSE., last_step, .FALSE.)
ENDDO
!
!   Electronic free energy as a function of the volume when computed
!
IF (lel_free_energy) THEN
   istep=0
   CALL gnuplot_set_fact(rytoev, .FALSE.)
   DO itempp=1,ntemp_plot
      first_step=(itempp==1)
      last_step=(itempp==ntemp_plot)
      itemp=itemp_plot(itempp)
      istep=MOD(istep,8)+1

      filename="anhar_files/"//TRIM(flevdat)//'_free.'// &
                                           TRIM(float_to_char(temp(itemp),1))
      CALL add_pressure(filename)
      IF (first_step) THEN
         CALL gnuplot_xlabel('Volume ((a.u.)^3)',.FALSE.) 
         CALL gnuplot_ylabel('Electronic free energy (eV)',.FALSE.) 
      ENDIF
      CALL gnuplot_write_file_mul_line_point(filename,1,4,color(istep), &
                                        first_step, last_step, .FALSE.)
   ENDDO
ENDIF
!
!  Gibbs free energy as a function of pressure
!
istep=0
CALL gnuplot_set_fact(1.0_DP, .FALSE.)
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flevdat)//'_mur.'//&
                                    TRIM(float_to_char(temp(itemp),1))
   CALL add_pressure(filename)

   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_write_command(TRIM(label),.FALSE.)

      CALL gnuplot_ylabel('Gibbs free-energy (Ry)',.FALSE.) 
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,4,3,color(istep),first_step, &
                                                  last_step, .FALSE.)
ENDDO
!
!  Bulk modulus as a function of pressure
!
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.temp.'// &
                                 TRIM(float_to_char(temp(itemp),1))
   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_write_command(TRIM(label),.FALSE.)

      CALL gnuplot_ylabel('Bulk modulus (kbar)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,3,color(istep),first_step,&
                                                            last_step,.FALSE.)
ENDDO
!
!  Derivative of the bulk modulus as a function of pressure
!
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.dbulk_temp.'// &
                                 TRIM(float_to_char(temp(itemp),1))
   IF (first_step) THEN
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_write_command(TRIM(label),.FALSE.)
      WRITE(label,'("set yrange [3.0:5.0]")') 
      CALL gnuplot_write_command(TRIM(label),.FALSE.)
      CALL gnuplot_ylabel('dB/dp',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step,&
                                                            last_step,.FALSE.)
ENDDO
CALL gnuplot_write_command('unset yrange',.FALSE.)
!
!  Second derivative of the bulk modulus with respect to pressure, if available
!
IF (ieos==2) THEN
   istep=0
   DO itempp=1,ntemp_plot
      first_step=(itempp==1)
      last_step=(itempp==ntemp_plot)
      itemp=itemp_plot(itempp)
      istep=MOD(istep,8)+1
      filename="anhar_files/"//TRIM(flanhar)//'.dbulk_temp.'// &
                                 TRIM(float_to_char(temp(itemp),1))
      IF (first_step) THEN
         CALL gnuplot_set_fact(1.0_DP,.FALSE.)
         CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 
         WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                             pmax*ry_kbar
         CALL gnuplot_write_command(TRIM(label),.FALSE.)
         WRITE(label,'("set yrange [-0.1:0.1]")') 
         CALL gnuplot_write_command(TRIM(label),.FALSE.)
         CALL gnuplot_ylabel('d^2B/dp^2     (1/kbar)',.FALSE.)
      ENDIF
      CALL gnuplot_write_file_mul_data(filename,1,3,color(istep),first_step,&
                                                            last_step,.FALSE.)
   ENDDO
   CALL gnuplot_write_command('unset yrange',.FALSE.)
ENDIF

!
!  Thermal expansion as a function of pressure
!
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.temp.'// &
                                 TRIM(float_to_char(temp(itemp),1))
   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_write_command(TRIM(label),.FALSE.)

      CALL gnuplot_ylabel('Thermal expansion ({/Symbol b} x 10^{6}) &
                                                    (K^{-1})',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,4,color(istep),first_step,&
                                                            last_step,.FALSE.)
ENDDO
!
!  Average Gruneisen parameter as a function of pressure
!
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.gamma_temp.'//&
                                 TRIM(float_to_char(temp(itemp),1))
   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_write_command(TRIM(label),.FALSE.)

      CALL gnuplot_ylabel('Average Gr\374neisen parameter',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step,&
                                                            last_step,.FALSE.)
ENDDO
!
!  Product beta B_T as a function of pressure
!
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.gamma_temp.'//&
                                 TRIM(float_to_char(temp(itemp),1))
   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_write_command(TRIM(label),.FALSE.)

      CALL gnuplot_ylabel('{/Symbol b} B_T (kbar/K)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,4,color(istep),first_step,&
                                                            last_step,.FALSE.)
ENDDO

CALL gnuplot_end()
!
!   close the file and make the plot
!
IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_anhar_t

!-------------------------------------------------------------------
SUBROUTINE plot_anhar_v()
!-------------------------------------------------------------------
!
!  This is a driver to plot the the thermal pressure,
!  as a function of temperature for several volumes.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE thermo_mod,       ONLY : omega_geo
USE control_gnuplot,  ONLY : flgnuplot, lgnuplot, gnuplot_command, flext
USE postscript_files, ONLY : flpsanhar
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_set_fact,            &
                             gnuplot_write_command,       &
                             gnuplot_write_file_mul_data_diff
USE data_files,       ONLY : flanhar
USE control_vol,      ONLY : nvol_plot, ivol_plot
USE temperature,      ONLY : temp, ntemp, tmin, tmax
USE color_mod,        ONLY : color
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: filename, gnu_filename, label
CHARACTER(LEN=8) :: float_to_char
INTEGER :: ivol, ivolp, istep
INTEGER :: ierr, system
REAL(DP) :: omega
LOGICAL :: first_step, last_step

IF ( my_image_id /= root_image ) RETURN

IF (nvol_plot==0) RETURN

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_anhar_vol'
filename=TRIM(flpsanhar)//'.vol'
filename=TRIM(filename)//TRIM(flext)

CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filename, tmin, tmax, &
                                0.0_DP, 0.0_DP, 1.0_DP, flext ) 
!
!  The thermal pressure as a function of the volume
!
istep=0
DO ivolp=1,nvol_plot
   first_step=(ivolp==1)
   last_step=(ivolp==nvol_plot)
   ivol=ivol_plot(ivolp)
   istep=MOD(istep,8)+1

   omega=omega_geo(ivol)
   filename="anhar_files/"//TRIM(flanhar)//'.vol.'//&
                                      TRIM(float_to_char(omega,2))
   CALL add_pressure(filename)
   IF (first_step) THEN
      CALL gnuplot_xlabel('Temperature (K)',.FALSE.) 
      CALL gnuplot_ylabel('Thermal pressure (kbar)',.FALSE.) 
   ENDIF
   CALL gnuplot_write_file_mul_data_diff(filename,1,2,3,color(istep),&
                                     first_step, last_step, .FALSE.)
ENDDO

CALL gnuplot_end()
!
!   close the file and make the plot
!
IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_anhar_v

! Copyright (C) 2018 Cristiano Malica

!-----------------------------------------------------------------------
SUBROUTINE plot_thermo_anhar()
!-----------------------------------------------------------------------
!
!  This is a driver to plot the quantities written inside flanhar.therm
!  and flanhar.therm_ph in the directory anhar_files
!  
USE kinds,            ONLY : DP
USE constants,        ONLY : rydberg_si, avogadro
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE control_thermo,   ONLY : ltherm_dos, ltherm_freq, with_eigen
USE postscript_files, ONLY : flpsanhar
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                             gnuplot_ylabel, &
                             gnuplot_xlabel, &
                             gnuplot_write_file_mul_data, &
                             gnuplot_set_fact
USE data_files,      ONLY : flanhar 
USE temperature,     ONLY : tmin, tmax
USE mp_images,       ONLY : root_image, my_image_id
USE io_global,       ONLY : ionode

IMPLICIT NONE
INTEGER :: ierr, system
CHARACTER(LEN=256) :: gnu_filename, filename, filetherm, filepstherm
REAL(DP) :: factor

IF ( my_image_id /= root_image ) RETURN

IF (.NOT.(ltherm_freq.OR.ltherm_dos)) RETURN

filetherm="anhar_files/"//TRIM(flanhar)//'.therm'
CALL add_pressure(filetherm)
filename="anhar_files/"//TRIM(flanhar)//'.therm_ph'
CALL add_pressure(filename)

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_therm_anhar'
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
CALL gnuplot_ylabel('Vibrational energy (kJ / (N mol))',.FALSE.)
factor = rydberg_si*avogadro / 1.D3
CALL gnuplot_set_fact(factor, .FALSE.)

IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filetherm,1,2,'color_red',.TRUE.,&
                                                     .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,2,'color_blue', &
                                                .NOT.ltherm_dos,.TRUE.,.FALSE.)

CALL gnuplot_ylabel('Vibrational free energy (kJ / (N mol))', .FALSE.)
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filetherm,1,3,'color_red',.TRUE.,&
                                                     .NOT.ltherm_freq, .FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,3,'color_blue',&
                                                .NOT.ltherm_dos,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(factor*1.D3, .FALSE.)
CALL gnuplot_ylabel('Entropy (J / K / (N mol))',.FALSE.)
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filetherm,1,4,'color_red',.TRUE., &
                                                   .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',.NOT.ltherm_dos,&
                                                             .TRUE.,.FALSE.)

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
END SUBROUTINE plot_thermo_anhar

!-----------------------------------------------------------------------
SUBROUTINE plot_dw_anhar()
!-----------------------------------------------------------------------
!
!  This is a driver to plot the quantities written inside 
!  
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
END SUBROUTINE plot_dw_anhar
