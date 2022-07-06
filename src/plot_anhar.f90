!
! Copyright (C) 2014-2022 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE plot_anhar_energy()
!-----------------------------------------------------------------------
!
!  This routine plots the energies in the anharmonic case  
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar, rytoev
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE control_thermo,   ONLY : ltherm_dos, ltherm_freq
USE postscript_files, ONLY : flpsanhar
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_write_file_mul_data, &
                             gnuplot_write_command,       &
                             gnuplot_write_file_mul_point,&
                             gnuplot_write_file_mul_point_sum,&
                             gnuplot_write_file_mul_line_point,&
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar
USE temperature,      ONLY : temp, tmin, tmax, ntemp_plot, itemp_plot
USE control_pressure, ONLY : pmin, pmax, pressure_kb
USE control_vol,      ONLY : vmin_input, vmax_input
USE control_eldos,    ONLY : lel_free_energy
USE color_mod,        ONLY : color
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename1, filename2, label
CHARACTER(LEN=8) :: float_to_char

INTEGER :: ierr, system, istep, itemp, itempp
LOGICAL :: first_step, last_step

IF ( my_image_id /= root_image ) RETURN
!
!   gnuplot script
!
gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_anhar_energy'
CALL add_pressure(gnu_filename)
!
!  name of the postcript file
!
filename=TRIM(flpsanhar)//'.energy'
CALL add_pressure(filename)
filename=TRIM(filename)//TRIM(flext)
!
!  Files with the data
!
filename1="anhar_files/"//TRIM(flanhar)
CALL add_pressure(filename1)
filename2="anhar_files/"//TRIM(flanhar)//'_ph'
CALL add_pressure(filename2)
!
!  open the script
!
CALL gnuplot_start(gnu_filename)
!
!  set the ranges and axis
!
IF (tmin /= 1.0_DP) THEN
   CALL gnuplot_write_header(filename, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext ) 
ELSE
   CALL gnuplot_write_header(filename, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext ) 
ENDIF
!
CALL gnuplot_xlabel('T (K)',.FALSE.) 
!
!  Helmholtz (or Gibbs) free energy
!
CALL gnuplot_set_fact(1.0_DP,.FALSE.)

IF (pressure_kb /= 0.0_DP) THEN
   CALL gnuplot_ylabel('Gibbs Free Energy (Ry)',.FALSE.) 
ELSE
   CALL gnuplot_ylabel('Helmholtz Free Energy (Ry)',.FALSE.) 
ENDIF

IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename1,1,3,'color_red',.TRUE.,&
                                              .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename2,1,3,'color_blue',&
                                        .NOT.ltherm_dos,.TRUE.,.FALSE.)
!
!  Helmholtz or Gibbs free energy as a function of the volume for several 
!  temperatures
!
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.free_temp'
   CALL add_value(filename,temp(itemp))
   CALL add_pressure(filename)
   filename1="anhar_files/"//TRIM(flanhar)//'.mur_temp'
   CALL add_value(filename1, temp(itemp))
   CALL add_pressure(filename1)
   IF (first_step) THEN
      CALL gnuplot_xlabel('Volume ((a.u.)^3)',.FALSE.) 
      IF (pressure_kb /= 0.0_DP) THEN
         label='Gibbs free-energy (Ry)    p= '//&
                    &TRIM(float_to_char(pressure_kb,1))//' kbar'
         CALL gnuplot_ylabel(TRIM(label),.FALSE.) 
      ELSE
         CALL gnuplot_ylabel('Helmholtz Free Energy (Ry)',.FALSE.) 
      END IF
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') vmin_input, &
                                                          vmax_input
      CALL gnuplot_write_command(TRIM(label),.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename1,1,2,color(istep),first_step, &
                                                        .FALSE., .FALSE.)
   CALL gnuplot_write_file_mul_point_sum(filename,1,2,3,color(istep), &
                                              .FALSE., last_step, .FALSE.)
ENDDO
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
   filename="anhar_files/"//TRIM(flanhar)//'.mur_temp'
   CALL add_value(filename,temp(itemp))
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
!   Vibrational free energy (+ electronic if computed) as a function of the 
!   volume
!
WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') vmin_input, vmax_input
CALL gnuplot_write_command(TRIM(label),.FALSE.)

CALL gnuplot_set_fact(rytoev, .FALSE.)
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1

   filename="anhar_files/"//TRIM(flanhar)//'.free_temp'
   CALL add_value(filename,temp(itemp))
   CALL add_pressure(filename)
   filename1="anhar_files/"//TRIM(flanhar)//'.poly_free_temp'
   CALL add_value(filename1,temp(itemp))
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
      filename="anhar_files/"//TRIM(flanhar)//'.free_temp'
      CALL add_value(filename,temp(itemp))
      CALL add_pressure(filename)
      IF (first_step) THEN
         CALL gnuplot_xlabel('Volume ((a.u.)^3)',.FALSE.) 
         CALL gnuplot_ylabel('Electronic free energy (eV)',.FALSE.) 
      ENDIF
      CALL gnuplot_write_file_mul_line_point(filename,1,4,color(istep), &
                                        first_step, last_step, .FALSE.)
   ENDDO
ENDIF

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_anhar_energy
!
!-----------------------------------------------------------------------
SUBROUTINE plot_anhar_volume()
!-----------------------------------------------------------------------
!
!  This routine plots the volume in the anharmonic case in
!  several forms 
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE control_thermo,   ONLY : ltherm_dos, ltherm_freq
USE postscript_files, ONLY : flpsanhar
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_write_file_mul_data, &
                             gnuplot_write_command,       &
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar
USE temperature,      ONLY : temp, tmin, tmax, ntemp_plot, itemp_plot
USE control_pressure, ONLY : press, pmin, pmax, pressure_kb, &
                             ipress_plot, npress_plot
USE control_vol,      ONLY : vmin_input, vmax_input
USE color_mod,        ONLY : color
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename1, filename2, label

INTEGER :: ierr, system, istep, ipress, ipressp, itemp, itempp
REAL(DP) :: factor
LOGICAL :: first_step, last_step

IF ( my_image_id /= root_image ) RETURN
!
!   gnuplot script
!
gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_anhar_volume'
CALL add_pressure(gnu_filename)
!
!  name of the postcript file
!
filename=TRIM(flpsanhar)//'.volume'
CALL add_pressure(filename)
filename=TRIM(filename)//TRIM(flext)
!
!  Files with the data
!
filename1="anhar_files/"//TRIM(flanhar)
CALL add_pressure(filename1)
filename2="anhar_files/"//TRIM(flanhar)//'_ph'
CALL add_pressure(filename2)
!
!  open the script
!
CALL gnuplot_start(gnu_filename)
!
!  set the ranges and axis
!
IF (tmin /= 1.0_DP) THEN
   CALL gnuplot_write_header(filename, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext ) 
ELSE
   CALL gnuplot_write_header(filename, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext ) 
ENDIF
!
CALL gnuplot_xlabel('T (K)',.FALSE.) 

CALL gnuplot_ylabel('Volume ((a.u.)^3)',.FALSE.)
!
!  Volume as a function of temperature at zero pressure
!
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename1,1,2,'color_red',.TRUE.,&
                                              .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename2,1,2,'color_blue',&
                                        .NOT.ltherm_dos,.TRUE.,.FALSE.)
!
!  Volume as a function of temperature at several pressures
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.press'
   CALL add_value(filename,press(ipress))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_ylabel('V(T) ((a.u.)^3)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step, &
                                                        last_step, .FALSE.)
ENDDO
!
!  V / V(300 K) as a function of temperature at several pressures
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.press'
   CALL add_value(filename,press(ipress))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_ylabel('V(T)/V(T=300 K)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,3,color(istep),first_step, &
                                                        last_step, .FALSE.)
ENDDO
!
!  Volume as a function of pressure at several temperatures
!
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.mur_temp'
   CALL add_value(filename,temp(itemp))
   CALL add_pressure(filename)
   IF (first_step) THEN
      CALL gnuplot_xlabel('Pressure (kbar)',.FALSE.)
      CALL gnuplot_ylabel('Volume ((a.u.)^3)',.FALSE.)
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_write_command(TRIM(label),.FALSE.)

   ENDIF
   CALL gnuplot_write_file_mul_data(filename,4,1,color(istep),first_step, &
                                                        last_step, .FALSE.)
ENDDO

istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.temp'
   CALL add_value(filename,temp(itemp))
   CALL add_pressure(filename)
   IF (first_step) THEN
      CALL gnuplot_xlabel('Pressure (kbar)',.FALSE.)
      CALL gnuplot_ylabel('V/V(T = 300 K) ',.FALSE.)
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_write_command(TRIM(label),.FALSE.)

   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,3,color(istep),first_step, &
                                                        last_step, .FALSE.)
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_anhar_volume
!
!-----------------------------------------------------------------------
SUBROUTINE plot_anhar_bulk()
!-----------------------------------------------------------------------
!
!  This routine plots the bulk modulus in the anharmonic case in
!  several forms 
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE control_thermo,   ONLY : ltherm_dos, ltherm_freq
USE postscript_files, ONLY : flpsanhar
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_write_file_mul_data, &
                             gnuplot_write_command,       &
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar
USE temperature,      ONLY : temp, tmin, tmax, ntemp_plot, itemp_plot
USE control_pressure, ONLY : press, pmin, pmax, pressure_kb, &
                             ipress_plot, npress_plot
USE control_vol,      ONLY : vmin_input, vmax_input
USE control_eldos,    ONLY : lel_free_energy
USE color_mod,        ONLY : color
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename1, filename2, &
                      filename3, filename_aux_grun, label

INTEGER :: ierr, system, istep, ipress, ipressp, itemp, itempp
LOGICAL :: first_step, last_step

IF ( my_image_id /= root_image ) RETURN
!
!   gnuplot script
!
gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_anhar_bulk'
CALL add_pressure(gnu_filename)
!
!  name of the postcript file
!
filename=TRIM(flpsanhar)//'.bulk'
CALL add_pressure(filename)
filename=TRIM(filename)//TRIM(flext)
!
!  Files with the data
!
filename1="anhar_files/"//TRIM(flanhar)//'.bulk'
CALL add_pressure(filename1)
filename2="anhar_files/"//TRIM(flanhar)//'.bulk_ph'
CALL add_pressure(filename2)
filename3="anhar_files/"//TRIM(flanhar)//'.el_anhar'
CALL add_pressure(filename3)
filename_aux_grun="anhar_files/"//TRIM(flanhar)//'.aux_grun'
CALL add_pressure(filename_aux_grun)
!
!  open the script
!
CALL gnuplot_start(gnu_filename)
!
!  set the ranges and axis
!
IF (tmin /= 1.0_DP) THEN
   CALL gnuplot_write_header(filename, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext ) 
ELSE
   CALL gnuplot_write_header(filename, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext ) 
ENDIF
!
CALL gnuplot_xlabel('T (K)',.FALSE.) 

CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('Bulk modulus (B_T) (kbar)',.FALSE.)

IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename1,1,2,'color_red',.TRUE., &
                                                .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename2,1,2,'color_blue',&
                              .NOT.ltherm_dos,.NOT.lel_free_energy,.FALSE.)
IF (lel_free_energy) &
   CALL gnuplot_write_file_mul_data(filename3,1,5,'color_orange',&
                             .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)
!
!   Bulk modulus as a function of temperature for several pressures
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1

   filename="anhar_files/"//TRIM(flanhar)//'.bulk_press'
   CALL add_value(filename,press(ipress))

   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_ylabel('B_T(T) (kbar))',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step, &
                                                        last_step, .FALSE.)
ENDDO
!
!   Bulk modulus as a function of pressure for several temperatures
!
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.bulk_temp'
   CALL add_value(filename,temp(itemp))
   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.)
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_write_command(TRIM(label),.FALSE.)

      CALL gnuplot_ylabel('B_T (kbar)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step,&
                                                            last_step,.FALSE.)
ENDDO

WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') tmin, tmax
CALL gnuplot_write_command(TRIM(label),.FALSE.)
CALL gnuplot_xlabel('T (K)',.FALSE.)

CALL gnuplot_set_fact(1._DP,.FALSE.)
CALL gnuplot_ylabel('B_S (kbar)',.FALSE.)
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename1,1,3,'color_red',.TRUE.,&
                                             .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename2,1,3,'color_blue',&
                                           .NOT.ltherm_dos,.TRUE.,.FALSE.)
!
!  B_S as a function of temperature, at several pressures
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.bulk_press'
   CALL add_value(filename,press(ipress))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_ylabel('B_S (kbar)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,3,color(istep),first_step, &
                                                        last_step,.FALSE.)
ENDDO
!
!   B_S as a function of pressure for several temperatures
!
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.bulk_temp'
   CALL add_value(filename,temp(itemp))
   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.)
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_write_command(TRIM(label),.FALSE.)

      CALL gnuplot_ylabel('B_S (kbar)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,3,color(istep),first_step,&
                                                            last_step,.FALSE.)
ENDDO

WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') tmin, tmax
CALL gnuplot_write_command(TRIM(label),.FALSE.)
CALL gnuplot_xlabel('T (K)',.FALSE.)

CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('B_S - B_T (kbar)',.FALSE.)
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename1,1,4,'color_red',.TRUE.,&
                                             .FALSE.,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename2,1,4,'color_blue',&
                                           .NOT.ltherm_dos,.FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filename_aux_grun,1,4,'color_green', &
                              .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)
!
!  B_S - B_T as a function of temperature, at several pressures
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.bulk_press'
   CALL add_value(filename,press(ipress))
   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_ylabel('B_S - B_T (kbar)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,4,color(istep),first_step, &
                                                        last_step,.FALSE.)
ENDDO
!
!   B_S - B_T as a function of pressure for several temperatures
!
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.bulk_temp'
   CALL add_value(filename,temp(itemp))
   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.)
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_write_command(TRIM(label),.FALSE.)

      CALL gnuplot_ylabel('B_S - B_T (kbar)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,4,color(istep),first_step,&
                                                            last_step,.FALSE.)
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_anhar_bulk
!
!-----------------------------------------------------------------------
SUBROUTINE plot_anhar_beta()
!-----------------------------------------------------------------------
!
!  This routine plots the volume thermal expansion in the anharmonic case in
!  several forms 
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE control_thermo,   ONLY : ltherm_dos, ltherm_freq
USE postscript_files, ONLY : flpsanhar
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_write_file_mul_data, &
                             gnuplot_write_command,       &
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar
USE temperature,      ONLY : temp, tmin, tmax, ntemp_plot, itemp_plot
USE control_pressure, ONLY : press, pmin, pmax, pressure_kb, &
                             ipress_plot, npress_plot
USE control_eldos,    ONLY : lel_free_energy
USE control_vol,      ONLY : vmin_input, vmax_input
USE color_mod,        ONLY : color
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename1, filename2, &
                      filename3, filename_aux_grun, label

INTEGER :: ierr, system, istep, ipress, ipressp, itemp, itempp
LOGICAL :: first_step, last_step

IF ( my_image_id /= root_image ) RETURN
!
!   gnuplot script
!
gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_anhar_beta'
CALL add_pressure(gnu_filename)
!
!  name of the postcript file
!
filename=TRIM(flpsanhar)//'.beta'
CALL add_pressure(filename)
filename=TRIM(filename)//TRIM(flext)
!
!  Files with the data
!
filename1="anhar_files/"//TRIM(flanhar)
CALL add_pressure(filename1)
filename2="anhar_files/"//TRIM(flanhar)//'_ph'
CALL add_pressure(filename2)
filename3="anhar_files/"//TRIM(flanhar)//'.el_anhar'
CALL add_pressure(filename3)
filename_aux_grun="anhar_files/"//TRIM(flanhar)//'.aux_grun'
CALL add_pressure(filename_aux_grun)
!
!  open the script
!
CALL gnuplot_start(gnu_filename)
!
!  set the ranges and axis
!
IF (tmin /= 1.0_DP) THEN
   CALL gnuplot_write_header(filename, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext ) 
ELSE
   CALL gnuplot_write_header(filename, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext ) 
ENDIF
!
CALL gnuplot_xlabel('T (K)',.FALSE.) 

CALL gnuplot_set_fact(1.0_DP,.FALSE.)
!
!  Volume thermal expansion
!
CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('Thermal expansion ({/Symbol b} x 10^{6}) (K^{-1})',.FALSE.)
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename1,1,4,'color_red',.TRUE.,&
                                                            .FALSE.,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename2,1,4,'color_blue',&
                                            .NOT.ltherm_dos,.FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filename_aux_grun,1,2,'color_green',&
                .NOT.(ltherm_dos.OR.ltherm_freq),.NOT.lel_free_energy,.FALSE.)

IF (lel_free_energy) &
   CALL gnuplot_write_file_mul_data(filename3,1,4,'color_orange',&
                             .FALSE.,.TRUE.,.FALSE.)
!
!  Thermal expansion as a function of temperature, at several pressures.
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1

   filename="anhar_files/"//TRIM(flanhar)//'.press'
   CALL add_value(filename,press(ipress))

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
!  Thermal expansion as a function of pressure, at several temperatures.
!
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.temp'
   CALL add_value(filename,temp(itemp))
   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_write_command(TRIM(label),.FALSE.)

      CALL gnuplot_ylabel('Thermal expansion ({/Symbol b} x 10^{6}) &
                                                    &(K^{-1})',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,5,color(istep),first_step,&
                                                            last_step,.FALSE.)
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_anhar_beta
!
!-----------------------------------------------------------------------
SUBROUTINE plot_anhar_heat()
!-----------------------------------------------------------------------
!
!  This routine plots the volume thermal expansion in the anharmonic case in
!  several forms 
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar, rydberg_si, avogadro
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE control_thermo,   ONLY : ltherm_dos, ltherm_freq
USE postscript_files, ONLY : flpsanhar
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,   &
                             gnuplot_write_header,         &
                             gnuplot_ylabel,               &
                             gnuplot_xlabel,               &
                             gnuplot_write_file_mul_data,  &
                             gnuplot_write_file_mul_point, &
                             gnuplot_write_command,        &
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar
USE temperature,      ONLY : temp, tmin, tmax, ntemp_plot, itemp_plot

USE control_pressure, ONLY : press, pmin, pmax, pressure_kb, &
                             ipress_plot, npress_plot
USE control_vol,      ONLY : vmin_input, vmax_input
USE control_eldos,    ONLY : lel_free_energy
USE color_mod,        ONLY : color
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename1, filename2,          &
                      filename3, filename4, filename_aux, filename_aux_grun, &
                      label

INTEGER :: ierr, system, istep, ipress, ipressp, itemp, itempp
REAL(DP) :: factor
LOGICAL :: first_step, last_step

IF ( my_image_id /= root_image ) RETURN
!
!   gnuplot script
!
gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_anhar_heat'
CALL add_pressure(gnu_filename)
!
!  name of the postcript file
!
filename=TRIM(flpsanhar)//'.heat'
CALL add_pressure(filename)
filename=TRIM(filename)//TRIM(flext)
!
!  Files with the data
!
filename1="anhar_files/"//TRIM(flanhar)//'.heat'
CALL add_pressure(filename1)
filename2="anhar_files/"//TRIM(flanhar)//'.heat_ph'
CALL add_pressure(filename2)
filename3="anhar_files/"//TRIM(flanhar)//'.el_anhar'
CALL add_pressure(filename3)
filename4="anhar_files/"//TRIM(flanhar)//'.el_press'
CALL add_pressure(filename4)
filename_aux_grun="anhar_files/"//TRIM(flanhar)//'.aux_grun'
CALL add_pressure(filename_aux_grun)
!
!  open the script
!
CALL gnuplot_start(gnu_filename)
!
!  set the ranges and axis
!
IF (tmin /= 1.0_DP) THEN
   CALL gnuplot_write_header(filename, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext ) 
ELSE
   CALL gnuplot_write_header(filename, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext ) 
ENDIF
!
CALL gnuplot_xlabel('T (K)',.FALSE.) 
CALL gnuplot_set_fact(1.0_DP,.FALSE.)
!
! isochoric heat capacity
!
factor = rydberg_si*avogadro
CALL gnuplot_set_fact(factor,.FALSE.)
CALL gnuplot_ylabel('Heat capacity C_v (J / K / N / mol)',.FALSE.)
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename1,1,2,'color_red',.TRUE.,&
                                            .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename2,1,2,'color_blue',&
                              .NOT.ltherm_dos,.NOT.lel_free_energy,.FALSE.)
IF (lel_free_energy) &
   CALL gnuplot_write_file_mul_data(filename4,1,2,'color_orange',&
                            .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)
!
! isochoric heat capacity as a function of temperature, at several pressures
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1

   filename="anhar_files/"//TRIM(flanhar)//'.heat_press'
   CALL add_value(filename,press(ipress))
   filename_aux="anhar_files/"//TRIM(flanhar)//'.el_press'
   CALL add_value(filename_aux,press(ipress))

   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      factor = rydberg_si*avogadro
      CALL gnuplot_set_fact(factor,.FALSE.)
      CALL gnuplot_ylabel('Heat capacity C_v (J / K / N / mol)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step, &
                           (last_step.AND.(.NOT.lel_free_energy)),.FALSE.)
   IF (lel_free_energy) &
      CALL gnuplot_write_file_mul_data(filename_aux,1,2,color(istep),.FALSE., &
                                                        last_step,.FALSE.)
ENDDO
!
! isochoric heat capacity as a function of pressure, at several temperatures
! Only if the required temperatures are larger than T=50 K
!
istep=0
first_step=.TRUE.
DO itempp=1,ntemp_plot
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   IF (temp(itemp)<49.9_DP) CYCLE
   filename="anhar_files/"//TRIM(flanhar)//'.heat_temp'
   CALL add_value(filename,temp(itemp))
   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_write_command(TRIM(label),.FALSE.)
      factor = rydberg_si*avogadro
      CALL gnuplot_set_fact(factor,.FALSE.)
      CALL gnuplot_ylabel('Heat capacity C_v (J / K / N / mol)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step,&
                                                            last_step,.FALSE.)
   first_step=.FALSE.
ENDDO

WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') tmin, tmax
CALL gnuplot_write_command(TRIM(label),.FALSE.)

CALL gnuplot_ylabel('Heat capacity C_p (J / K / N / mol)',.FALSE.)
factor = rydberg_si*avogadro
CALL gnuplot_set_fact(factor,.FALSE.)
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename1,1,4,'color_red',.TRUE.,    &
                                          .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename2,1,4,'color_blue',          &
                            .NOT.ltherm_dos,.NOT.lel_free_energy,.FALSE.)

IF (lel_free_energy) &
   CALL gnuplot_write_file_mul_data(filename3,1,3,'color_orange',&
                   .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)

!
!  put as a comment the possibility to plot also the experimental data
!
   CALL gnuplot_write_file_mul_data(filename1,1,4,'color_blue',&
                                       .NOT.ltherm_dos,.FALSE.,.TRUE.)
   CALL gnuplot_write_file_mul_point('cv.exp',1,2,'color_red',.FALSE.,&
                                                              .TRUE.,.TRUE.)
!
! isobaric heat capacity as a function of temperature, at several pressures
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1

   filename="anhar_files/"//TRIM(flanhar)//'.heat_press'
   CALL add_value(filename,press(ipress))

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
! isobaric heat capacity as a function of pressure, at several temperatures
!
istep=0
first_step=.TRUE.
DO itempp=1,ntemp_plot
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   IF (temp(itemp)<49.9_DP) CYCLE
   filename="anhar_files/"//TRIM(flanhar)//'.heat_temp'
   CALL add_value(filename,temp(itemp))
   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_write_command(TRIM(label),.FALSE.)
      factor = rydberg_si*avogadro
      CALL gnuplot_set_fact(factor,.FALSE.)
      CALL gnuplot_ylabel('Heat capacity C_p (J / K / N / mol)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,4,color(istep),first_step,&
                                                            last_step,.FALSE.)
   first_step=.FALSE.
ENDDO
CALL gnuplot_write_command('unset yrange',.FALSE.)
!
!  Now the difference C_p-C_V as a function of temperature, at zero pressure.
!
WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') tmin, tmax
CALL gnuplot_write_command(TRIM(label),.FALSE.)
!
!  isobaric -isochoric heat capacity
!
factor = rydberg_si*avogadro
CALL gnuplot_set_fact(factor,.FALSE.)
CALL gnuplot_ylabel('C_p - C_v (J / K / N / mol)',.FALSE.) 

IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename1,1,3,'color_red',.TRUE., &
                                                   .FALSE.,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename2,1,3,'color_blue',&
                                               .NOT.ltherm_dos,.FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filename_aux_grun,1,3,'color_green', &
                              .NOT.(ltherm_dos.OR.ltherm_freq),.TRUE.,.FALSE.)
!
!  isobaric - isochoric heat capacity as a function of temperature at several
!  pressures
!
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.heat_press'
   CALL add_value(filename,press(ipress))
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
! isobaric -isochoric heat capacity as a function of pressure, 
! at several temperatures
!
istep=0
first_step=.TRUE.
DO itempp=1,ntemp_plot
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   IF (temp(itemp)<49.9_DP) CYCLE
   filename="anhar_files/"//TRIM(flanhar)//'.heat_temp'
   CALL add_value(filename,temp(itemp))
   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_write_command(TRIM(label),.FALSE.)
      factor = rydberg_si*avogadro
      CALL gnuplot_set_fact(factor,.FALSE.)
      CALL gnuplot_ylabel('Heat capacity C_p -C_V(J / K / N / mol)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,3,color(istep),first_step,&
                                                            last_step,.FALSE.)
   first_step=.FALSE.
ENDDO
CALL gnuplot_write_command('unset yrange',.FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_anhar_heat
!
!-----------------------------------------------------------------------
SUBROUTINE plot_anhar_gamma()
!-----------------------------------------------------------------------
!
!  This routine plots the volume thermal expansion in the anharmonic case in
!  several forms 
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE control_thermo,   ONLY : ltherm_dos, ltherm_freq
USE postscript_files, ONLY : flpsanhar
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_write_file_mul_data, &
                             gnuplot_write_command,       &
                             gnuplot_write_horizontal_line, &
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar
USE temperature,      ONLY : temp, tmin, tmax, ntemp_plot, itemp_plot

USE control_pressure, ONLY : press, pmin, pmax, pressure_kb, &
                             ipress_plot, npress_plot
USE control_vol,      ONLY : vmin_input, vmax_input
USE control_eldos,    ONLY : lel_free_energy
USE color_mod,        ONLY : color
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename1, filename2, &
                      filename3, filename_gamma_grun, label

INTEGER  :: ierr, system, istep, ipress, ipressp, itemp, itempp
REAL(DP) :: factor
LOGICAL  :: first_step, last_step

IF ( my_image_id /= root_image ) RETURN
!
!   gnuplot script
!
gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_anhar_gamma'
CALL add_pressure(gnu_filename)
!
!  name of the postcript file
!
filename=TRIM(flpsanhar)//'.gamma'
CALL add_pressure(filename)
filename=TRIM(filename)//TRIM(flext)
!
!  Files with the data
!
filename1="anhar_files/"//TRIM(flanhar)//'.gamma'
CALL add_pressure(filename1)
filename2="anhar_files/"//TRIM(flanhar)//'.gamma_ph'
CALL add_pressure(filename2)
filename3="anhar_files/"//TRIM(flanhar)//'.el_anhar'
CALL add_pressure(filename3)
filename_gamma_grun="anhar_files/"//TRIM(flanhar)//'.gamma_grun'
CALL add_pressure(filename_gamma_grun)
!
!  open the script
!
CALL gnuplot_start(gnu_filename)
!
!  set the ranges of temperatures (below 50 K the Average Gruneisen parameter
!  is inaccurate)
!
IF (tmin /= 1.0_DP) THEN
   CALL gnuplot_write_header(filename, MAX(tmin, 50.0_DP), tmax, 0.0_DP, &
                                                 0.0_DP, 1.0_DP, flext ) 
ELSE
   CALL gnuplot_write_header(filename, 50.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                        1.0_DP, flext ) 
ENDIF
!
CALL gnuplot_xlabel('T (K)',.FALSE.) 

CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('Gr\374neisen parameter ({/Symbol g})',.FALSE.)
CALL gnuplot_write_horizontal_line(0.0_DP, 2, 'front', 'color_black', .FALSE.)
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename1,1,2,'color_red',.TRUE.,&
                                                     .FALSE.,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename2,1,2,'color_blue',&
                         .NOT.ltherm_dos,.FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filename_gamma_grun,1,2,'color_green', &
              .NOT.(ltherm_dos.OR.ltherm_freq),.NOT.lel_free_energy,.FALSE.)
!
!   Average Gruneisen parameter as a function of temperature, for 
!   several pressures
!
IF (lel_free_energy) &
   CALL gnuplot_write_file_mul_data(filename3,1,2,'color_orange',&
                             .FALSE.,.TRUE.,.FALSE.)

istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1

   filename="anhar_files/"//TRIM(flanhar)//'.gamma_press'
   CALL add_value(filename,press(ipress))

   IF (first_step) THEN
      CALL gnuplot_xlabel('T (K)',.FALSE.)
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_ylabel('Gr\374neisen parameter ({/Symbol g})',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step, &
                                                        last_step,.FALSE.)
ENDDO
!
!   Average Gruneisen parameter as a function of pressure, for 
!   several temperatures
!   Do not plot the average Gruneisen parameters for temperatures lower than
!   T = 50 K
!
istep=0
first_step=.TRUE.
DO itempp=1,ntemp_plot
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   IF (temp(itemp)<49.9_DP) CYCLE
   filename="anhar_files/"//TRIM(flanhar)//'.gamma_temp'
   CALL add_value(filename,temp(itemp))
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
   first_step=.FALSE.
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

   filename="anhar_files/"//TRIM(flanhar)//'.gamma_temp'
   CALL add_value(filename,temp(itemp))

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

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_anhar_gamma

!-----------------------------------------------------------------------
SUBROUTINE plot_anhar_press()
!-----------------------------------------------------------------------
!
!  This routine plots the pressure as a function of volume or of 
!  temperature in several forms 
!
USE kinds,            ONLY : DP
USE thermo_mod,       ONLY : omega_geo
USE constants,        ONLY : ry_kbar
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE control_thermo,   ONLY : ltherm_dos, ltherm_freq
USE postscript_files, ONLY : flpsanhar
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_write_file_mul_data, &
                             gnuplot_write_file_mul_data_diff, &
                             gnuplot_write_command,       &
                             gnuplot_write_horizontal_line, &
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar
USE temperature,      ONLY : temp, tmin, tmax, ntemp_plot, itemp_plot

USE control_pressure, ONLY : press, pmin, pmax, pressure_kb, &
                             ipress_plot, npress_plot
USE control_vol,      ONLY : vmin_input, vmax_input, nvol_plot, ivol_plot
USE color_mod,        ONLY : color
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename1, filename2, &
                      filename_aux, filename_gamma_grun, label

INTEGER :: ierr, system, istep, ipress, ipressp, itemp, itempp, ivol, ivolp
REAL(DP) :: factor, omega
LOGICAL :: first_step, last_step
CHARACTER(LEN=8) :: float_to_char

IF ( my_image_id /= root_image ) RETURN
IF ( ntemp_plot==0 .AND. nvol_plot==0 ) RETURN
!
!   gnuplot script
!
gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_pressure'
CALL add_pressure(gnu_filename)
!
!  name of the postcript file
!
filename=TRIM(flpsanhar)//'.pressure'
CALL add_pressure(filename)
filename=TRIM(filename)//TRIM(flext)
!
!  open the script
!
CALL gnuplot_start(gnu_filename)
!
!  set the ranges of volumes
!
CALL gnuplot_write_header(filename, vmin_input, vmax_input, 0.0_DP, &
                                                0.0_DP, 1.0_DP, flext ) 
!
CALL gnuplot_xlabel('T (K)',.FALSE.) 
!
!  Pressure as a function of the volume, at several temperatures
!
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.mur_temp'
   CALL add_value(filename, temp(itemp))
   CALL add_pressure(filename)
   IF (first_step) THEN
      CALL gnuplot_xlabel('Volume ((a.u.)^3)',.FALSE.) 
      CALL gnuplot_ylabel('Pressure (kbar)',.FALSE.) 
      CALL gnuplot_write_horizontal_line(0.0_DP, 2, 'front', &
                                              'color_black', .FALSE.)
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

   filename="anhar_files/"//TRIM(flanhar)//'.mur_temp'
   CALL add_value(filename,temp(itemp))
   CALL add_pressure(filename)
   filename_aux="anhar_files/"//TRIM(flanhar)//'.poly_free_temp'
   CALL add_value(filename_aux,temp(itemp))
   CALL add_pressure(filename_aux)
   IF (first_step) THEN
      CALL gnuplot_xlabel('Volume ((a.u.)^3)',.FALSE.) 
      CALL gnuplot_ylabel('Thermal pressure (kbar)',.FALSE.) 
   ENDIF
   CALL gnuplot_write_file_mul_data_diff(filename,1,4,5,color(istep),&
                             first_step, .FALSE., .FALSE.)
   CALL gnuplot_write_file_mul_data(filename_aux,1,3,color(istep),&
                             .FALSE., last_step, .FALSE.)
ENDDO
!
!   Thermal pressure as a function of temperature
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
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') tmin, tmax
      CALL gnuplot_write_command(TRIM(label),.FALSE.)
      CALL gnuplot_xlabel('Temperature (K)',.FALSE.)
      CALL gnuplot_ylabel('Thermal pressure (kbar)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data_diff(filename,1,2,3,color(istep),&
                                     first_step, last_step, .FALSE.)
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_anhar_press
!
!-----------------------------------------------------------------------
SUBROUTINE plot_anhar_dbulk()
!-----------------------------------------------------------------------
!
!  This routine plots the derivative of the bulk modulus with 
!  respect to pressure (and also the second derivative if ieos=2)
!  as a function of temperature
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE control_thermo,   ONLY : ltherm_dos, ltherm_freq
USE postscript_files, ONLY : flpsanhar
USE control_ev,       ONLY : ieos
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_write_file_mul_data, &
                             gnuplot_write_file_mul_data_diff, &
                             gnuplot_write_command,       &
                             gnuplot_write_horizontal_line, &
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar
USE temperature,      ONLY : temp, tmin, tmax, ntemp_plot, itemp_plot
USE control_pressure, ONLY : press, pmin, pmax, pressure_kb, &
                             ipress_plot, npress_plot
USE control_vol,      ONLY : vmin_input, vmax_input
USE color_mod,        ONLY : color
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, filename1, filename2, &
                      filename_gamma_grun, label

INTEGER :: ierr, system, istep, ipress, ipressp, itemp, itempp
REAL(DP) :: factor
LOGICAL :: first_step, last_step

IF ( my_image_id /= root_image ) RETURN
!
!   gnuplot script
!
gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'.dbulk'
CALL add_pressure(gnu_filename)
!
!  name of the postcript file
!
filename=TRIM(flpsanhar)//'.dbulk'
CALL add_pressure(filename)
filename=TRIM(filename)//TRIM(flext)
!
!  Names of the files with data
!
filename1="anhar_files/"//TRIM(flanhar)//'.dbulk'
CALL add_pressure(filename1)
filename2="anhar_files/"//TRIM(flanhar)//'.dbulk_ph'
CALL add_pressure(filename2)

!
!  open the script
!
CALL gnuplot_start(gnu_filename)
!
!  set the ranges of temperatures
!
IF (tmin /= 1.0_DP) THEN
   CALL gnuplot_write_header(filename, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext )
ELSE
   CALL gnuplot_write_header(filename, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP, &
                                                                   flext )
ENDIF
!
!
CALL gnuplot_xlabel('T (K)',.FALSE.) 

!
!  pressure derivative of the bulk modulus
!
CALL gnuplot_set_fact(1.0_DP,.FALSE.)
CALL gnuplot_ylabel('d B / d p',.FALSE.)

IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filename1,1,2,'color_red',.TRUE., &
                                                .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename2,1,2,'color_blue',&
                                               .NOT.ltherm_dos,.TRUE.,.FALSE.)
istep=0
DO ipressp=1,npress_plot
   first_step=(ipressp==1)
   last_step=(ipressp==npress_plot)
   ipress=ipress_plot(ipressp)
   istep=MOD(istep,8)+1

   filename="anhar_files/"//TRIM(flanhar)//'.dbulk_press'
   CALL add_value(filename,press(ipress))

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

istep=0
first_step=.TRUE.
DO itempp=1,ntemp_plot
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.dbulk_temp'
   CALL add_value(filename,temp(itemp))
   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.)
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_write_command(TRIM(label),.FALSE.)
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_ylabel('d B_T/dp (p)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step, &
                                                        last_step, .FALSE.)
   first_step=.FALSE.
ENDDO

!
!   The second derivative of the bulk modulus with respect to pressure
!
IF (ieos==2) THEN
   CALL gnuplot_set_fact(1.0_DP,.FALSE.)
   CALL gnuplot_ylabel('d^2 B / d p^2 (1/kbar)',.FALSE.)

   IF (ltherm_dos) &
      CALL gnuplot_write_file_mul_data(filename1,1,3,'color_red',.TRUE., &
                                                .NOT.ltherm_freq,.FALSE.)
   IF (ltherm_freq) &
      CALL gnuplot_write_file_mul_data(filename2,1,3,'color_blue',&
                                               .NOT.ltherm_dos,.TRUE.,.FALSE.)
!
!   The derivative of the bulk modulus with respect to pressure as a
!   function of temperature for several pressures
!
   istep=0
   DO ipressp=1,npress_plot
      first_step=(ipressp==1)
      last_step=(ipressp==npress_plot)
      ipress=ipress_plot(ipressp)
      istep=MOD(istep,8)+1

      filename="anhar_files/"//TRIM(flanhar)//'.dbulk_press'
      CALL add_value(filename,press(ipress))

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

   istep=0
   first_step=.TRUE.
   DO itempp=1,ntemp_plot
      last_step=(itempp==ntemp_plot)
      itemp=itemp_plot(itempp)
      istep=MOD(istep,8)+1
      filename="anhar_files/"//TRIM(flanhar)//'.dbulk_temp'
      CALL add_value(filename,temp(itemp))
      IF (first_step) THEN
         CALL gnuplot_xlabel('pressure (kbar)',.FALSE.)
         WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                             pmax*ry_kbar
         CALL gnuplot_write_command(TRIM(label),.FALSE.)
         CALL gnuplot_set_fact(1.0_DP,.FALSE.)
         CALL gnuplot_ylabel('d^2 B_T/dp^2 (p)',.FALSE.)
      ENDIF
      CALL gnuplot_write_file_mul_data(filename,1,3,color(istep),first_step, &
                                                        last_step, .FALSE.)
      first_step=.FALSE.
   ENDDO
ENDIF

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_anhar_dbulk

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
USE data_files,       ONLY : flanhar
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
   filename="anhar_files/"//TRIM(flanhar)//'_mur'
   CALL add_value(filename,temp(itemp))
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

   filename="anhar_files/"//TRIM(flanhar)//'_mur'
   CALL add_value(filename,temp(itemp))
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

   filename="anhar_files/"//TRIM(flanhar)//'_mur'
   CALL add_value(filename,temp(itemp))
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

   filename="anhar_files/"//TRIM(flanhar)//'.free_temp'
   CALL add_value(filename,temp(itemp))
   CALL add_pressure(filename)
   filename1="anhar_files/"//TRIM(flanhar)//'.poly_free'
   CALL add_value(filename1,temp(itemp))
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

      filename="anhar_files/"//TRIM(flanhar)//'.free_temp'
      CALL add_value(filename,temp(itemp))
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
   filename="anhar_files/"//TRIM(flanhar)//'_mur'
   CALL add_value(filename,temp(itemp))
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
   filename="anhar_files/"//TRIM(flanhar)//'.bulk_temp'
   CALL add_value(filename,temp(itemp))
   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_write_command(TRIM(label),.FALSE.)

      CALL gnuplot_ylabel('Bulk modulus (kbar)',.FALSE.)
   ENDIF
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step,&
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
   filename="anhar_files/"//TRIM(flanhar)//'.dbulk_temp'
   CALL add_value(filename,temp(itemp))
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
      filename="anhar_files/"//TRIM(flanhar)//'.dbulk_temp'
      CALL add_value(filename,temp(itemp))
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
   filename="anhar_files/"//TRIM(flanhar)//'.temp'
   CALL add_value(filename,temp(itemp))
   IF (first_step) THEN
      CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 
      WRITE(label,'("set xrange [",f12.5,":",f12.5,"]")') pmin*ry_kbar, &
                                                          pmax*ry_kbar
      CALL gnuplot_set_fact(1.0_DP,.FALSE.)
      CALL gnuplot_write_command(TRIM(label),.FALSE.)

      CALL gnuplot_ylabel('Thermal expansion ({/Symbol b} x 10^{6}) &
                                                    &(K^{-1})',.FALSE.)
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
   filename="anhar_files/"//TRIM(flanhar)//'.gamma_temp'
   CALL add_value(filename,temp(itemp))
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

   filename="anhar_files/"//TRIM(flanhar)//'.gamma_temp'
   CALL add_value(filename,temp(itemp))

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

