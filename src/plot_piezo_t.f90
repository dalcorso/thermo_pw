!
! Copyright (C) 2025 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE plot_piezo_t(iflag)
!---------------------------------------------------------------------
!
!  This is a driver to plot the piezoelectric tensor as a function of
!  temperature or pressure
!  iflag    0   piezoelectric tensor as a function of temperature
!           2   piezoelectric tensor as a function of pressure
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar, electron_si, bohr_radius_si
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar, fl_piezo
USE postscript_files, ONLY : flpsanhar
USE control_piezoelectric_tensor,  ONLY : lpiezo, lpiezof
USE piezoelectric_tensor, ONLY : get_pt_type, pt_present, pt_names, pt_types
USE control_pressure, ONLY : pmin, pmax
USE temperature,      ONLY : tmin, tmax
USE rap_point_group,  ONLY : code_group
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
INTEGER :: iflag
CHARACTER(LEN=256) :: gnu_filename, filenameps, filepiezo
CHARACTER(LEN=10) :: unit
INTEGER :: ierr, system, i, j, ipt, piezo_type
REAL(DP) :: fact

IF ( my_image_id /= root_image ) RETURN

IF (.NOT.(lpiezo.OR.lpiezof)) RETURN

IF (iflag==0) THEN
   gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_e_piezo"
   filenameps=TRIM(flpsanhar)//".e_piezo"//TRIM(flext)
   filepiezo="anhar_files/"//TRIM(flanhar)//".e_piezo"
ELSE
!
!  In this case plot the elastic constants as a function of pressure
!
   gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_e_piezo_p"
   filenameps=TRIM(flpsanhar)//".e_piezo_p"//TRIM(flext)
   filepiezo="elastic_constants/"//TRIM(fl_piezo)//".e_piezo_p"
ENDIF

CALL gnuplot_start(gnu_filename)

IF (iflag<2) THEN
   IF (tmin ==1._DP) THEN
      CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext ) 
   ELSE
      CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext ) 
   ENDIF
   CALL gnuplot_xlabel('T (K)', .FALSE.) 
ELSE
   CALL gnuplot_write_header(filenameps, pmin*ry_kbar, pmax*ry_kbar, &
                                              0.0_DP, 0.0_DP, 1.0_DP, flext ) 
   CALL gnuplot_xlabel('p (kbar)', .FALSE.) 
ENDIF

fact= electron_si / (bohr_radius_si)**2
piezo_type=get_pt_type(code_group)
unit=' (C/m^2)'
DO ipt=1,pt_types
   IF (pt_present(ipt, piezo_type)>0) THEN 
      CALL gnuplot_set_fact(fact, .FALSE.)
      CALL plot_one_piezo_element(1, pt_present(ipt, piezo_type)+1, &
              pt_names(ipt) // TRIM(unit), lpiezo, lpiezof, filepiezo)
   ENDIF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_piezo_t

!---------------------------------------------------------------------
SUBROUTINE plot_one_piezo_element(i, j, label, lpiezo, lpiezof, filepiezo)
!---------------------------------------------------------------------
USE thermo_mod, ONLY : what
USE control_thermo,   ONLY : ltherm_freq, ltherm_dos
USE gnuplot, ONLY : gnuplot_ylabel, &
                    gnuplot_write_file_mul_data
USE color_mod, ONLY : color
IMPLICIT NONE
INTEGER, INTENT(IN) :: i, j
LOGICAL, INTENT(IN) :: lpiezo, lpiezof
CHARACTER(LEN=*) :: label
CHARACTER(LEN=256), INTENT(IN) :: filepiezo

CHARACTER(LEN=256) :: filename, filename_ph, filename_s_ph

filename=TRIM(filepiezo)
filename_ph=TRIM(filepiezo)//'_ph'

CALL gnuplot_ylabel(TRIM(label),.FALSE.) 

IF (ltherm_dos) THEN
   CALL gnuplot_write_file_mul_data(filename,i,j,'color_red',.TRUE., &
                                                  .NOT.ltherm_freq,.FALSE.)
ENDIF
IF (ltherm_freq) THEN
   CALL gnuplot_write_file_mul_data(filename_ph,i,j,'color_red',     &
                                            .NOT.ltherm_dos,.TRUE.,.FALSE.)
ENDIF

RETURN
END SUBROUTINE plot_one_piezo_element
!
!---------------------------------------------------------------------
SUBROUTINE plot_piezo_pt()
!---------------------------------------------------------------------
!
!  This is a driver to plot the piezoelectric tensor 
!  as a function of temperature at a set of pressures.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar, electron_si, bohr_radius_si
USE control_thermo,   ONLY : ltherm_freq, ltherm_dos
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact, gnuplot_ylabel,     &
                             gnuplot_write_file_mul_data
USE color_mod,        ONLY : color
USE data_files,       ONLY : flanhar, fl_el_cons
USE piezoelectric_tensor, ONLY : get_pt_type, pt_present, pt_names, &
                                 pt_types
USE postscript_files, ONLY : flpsanhar, flps_el_cons
USE control_piezoelectric_tensor,  ONLY : lpiezo_pt, lpiezof_pt
USE control_pressure, ONLY : pmin, pmax, npress_plot, ipress_plot, press
USE temperature,      ONLY : tmin, tmax
USE color_mod,        ONLY : color
USE thermo_sym,       ONLY : laue
USE rap_point_group,  ONLY : code_group
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
INTEGER :: iflag
LOGICAL :: with_s
CHARACTER(LEN=256) :: gnu_filename, filenameps, filepiezo, filepiezo_ph, label
INTEGER :: piezo_type, ipt, istep, ipressp, ipress
INTEGER :: ierr, system
LOGICAL :: first, last, first_step, last_step
REAL(DP) :: fact

IF ( my_image_id /= root_image ) RETURN
IF (.NOT.(lpiezo_pt.OR.lpiezof_pt)) RETURN
IF (npress_plot==0) RETURN


gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_e_piezo_p"
filenameps=TRIM(flpsanhar)//".e_piezo_p"//TRIM(flext)

CALL gnuplot_start(gnu_filename)

IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                    1.0_DP, flext ) 
ELSE
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                    1.0_DP, flext ) 
ENDIF
CALL gnuplot_xlabel('T (K)', .FALSE.) 

piezo_type=get_pt_type(code_group)

fact= electron_si / (bohr_radius_si)**2
DO ipt=1,pt_types
   IF (pt_present(ipt, piezo_type)>0) THEN 
      istep=0
      DO ipressp=1, npress_plot
         first_step=(ipressp==1)
         last_step=(ipressp==npress_plot)
         ipress=ipress_plot(ipressp)
         istep=MOD(istep,8)+1
         filepiezo="anhar_files/"//TRIM(flanhar)//'.e_piezo_press'
         CALL add_value(filepiezo,press(ipress))
         filepiezo_ph="anhar_files/"//TRIM(flanhar)//'.e_piezo_ph_press'
         CALL add_value(filepiezo_ph,press(ipress))
         IF (first_step) THEN
            CALL gnuplot_set_fact(fact, .FALSE.) 
            WRITE(label,'(a," (C/m^2)")') TRIM(pt_names(ipt)) 
            CALL gnuplot_ylabel(TRIM(label),.FALSE.)
         ENDIF
         IF (ltherm_dos) THEN
            CALL gnuplot_write_file_mul_data(filepiezo,1,                &
                          pt_present(ipt, piezo_type)+1,color(istep), &
                          first_step,(last_step.AND..NOT.ltherm_freq),.FALSE.)
         ENDIF
         IF (ltherm_freq) THEN
            CALL gnuplot_write_file_mul_data(filepiezo_ph,1,&
                          pt_present(ipt, piezo_type)+1,color(istep), &
                          (first_step.AND..NOT.ltherm_dos),last_step,.FALSE.)
         ENDIF
      ENDDO
   ENDIF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_piezo_pt
!
!---------------------------------------------------------------------
SUBROUTINE plot_piezo_ptt()
!---------------------------------------------------------------------
!
!  This is a driver to plot the piezoelectric tensor
!  as a function of pressure at a set of temperatures.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar, electron_si, bohr_radius_si
USE control_thermo,   ONLY : ltherm_freq, ltherm_dos
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact, gnuplot_ylabel,     &
                             gnuplot_write_file_mul_data
USE color_mod,        ONLY : color
USE data_files,       ONLY : flanhar, fl_piezo
USE piezoelectric_tensor, ONLY : get_pt_type, pt_present, pt_names, pt_types
USE postscript_files, ONLY : flpsanhar
USE control_piezoelectric_tensor,  ONLY : lpiezo_ptt, lpiezof_ptt
USE control_pressure, ONLY : pmin, pmax
USE temperature,      ONLY : temp, ntemp_plot, itemp_plot
USE color_mod,        ONLY : color
USE rap_point_group,  ONLY : code_group
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: gnu_filename, filenameps, filepiezo, filepiezo_ph, label
INTEGER :: pt_type, ipt, istep, itempp, itemp
INTEGER :: ierr, system
LOGICAL :: first, last, first_step, last_step
REAL(DP) :: fact

IF ( my_image_id /= root_image ) RETURN
IF (.NOT.(lpiezo_ptt.OR.lpiezof_ptt)) RETURN
IF (ntemp_plot==0) RETURN


gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_e_piezo_t"
filenameps=TRIM(flpsanhar)//".e_piezo_t"//TRIM(flext)

CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filenameps, pmin*ry_kbar, pmax*ry_kbar, 0.0_DP, &
                              0.0_DP, 1.0_DP, flext ) 
CALL gnuplot_xlabel('p (kbar)', .FALSE.) 

pt_type=get_pt_type(code_group)
fact= electron_si / (bohr_radius_si)**2
DO ipt=1,pt_types
   IF (pt_present(ipt, pt_type)>0) THEN 
      istep=0
      DO itempp=1, ntemp_plot
         first_step=(itempp==1)
         last_step=(itempp==ntemp_plot)
         itemp=itemp_plot(itempp)
         istep=MOD(istep,8)+1
         filepiezo="anhar_files/"//TRIM(flanhar)//'.e_piezo_temp'
         CALL add_value(filepiezo,temp(itemp))
         filepiezo_ph="anhar_files/"//TRIM(flanhar)//'.e_piezo_ph_temp'
         CALL add_value(filepiezo_ph,temp(itemp))
         IF (first_step) THEN
            CALL gnuplot_set_fact(fact, .FALSE.) 
            WRITE(label,'(a," (C/m^2)")') TRIM(pt_names(ipt)) 
            CALL gnuplot_ylabel(TRIM(label),.FALSE.)
         ENDIF
         IF (ltherm_dos) THEN
            CALL gnuplot_write_file_mul_data(filepiezo,1,          &
                          pt_present(ipt, pt_type)+1,color(istep), &
                          first_step,(last_step.AND.               &
                           .NOT.ltherm_freq),.FALSE.)
         ENDIF
         IF (ltherm_freq) THEN
            CALL gnuplot_write_file_mul_data(filepiezo_ph,1,&
                          pt_present(ipt, pt_type)+1,color(istep), &
                          (first_step.AND..NOT.ltherm_dos),&
                          last_step,.FALSE.)
         ENDIF
      ENDDO
   ENDIF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_piezo_ptt

!---------------------------------------------------------------------
SUBROUTINE plot_piezo_d_t(iflag)
!---------------------------------------------------------------------
!
!  This is a driver to plot the d piezoelectric tensor as a function of
!  temperature or pressure
!  iflag    0   piezoelectric tensor as a function of temperature
!           2   piezoelectric tensor as a function of pressure
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar, electron_si, bohr_radius_si
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar, fl_piezo
USE postscript_files, ONLY : flpsanhar
USE control_piezoelectric_tensor,  ONLY : lpiezo_d, lpiezof_d
USE piezoelectric_tensor, ONLY : get_pt_type, pt_present, ptd_names, pt_types
USE control_pressure, ONLY : pmin, pmax
USE temperature,      ONLY : tmin, tmax
USE rap_point_group,  ONLY : code_group
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
INTEGER :: iflag
CHARACTER(LEN=256) :: gnu_filename, filenameps, filepiezo
CHARACTER(LEN=10) :: unit
INTEGER :: ierr, system, i, j, ipt, piezo_type
REAL(DP) :: fact

IF ( my_image_id /= root_image ) RETURN

IF (.NOT.(lpiezo_d.OR.lpiezof_d)) RETURN

IF (iflag==0) THEN
   gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_d_piezo"
   filenameps=TRIM(flpsanhar)//".d_piezo"//TRIM(flext)
   filepiezo="anhar_files/"//TRIM(flanhar)//".d_piezo"
ELSE
!
!  In this case plot the elastic constants as a function of pressure
!
   gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_d_piezo_p"
   filenameps=TRIM(flpsanhar)//".d_piezo_p"//TRIM(flext)
   filepiezo="elastic_constants/"//TRIM(fl_piezo)//".d_piezo_p"
ENDIF

CALL gnuplot_start(gnu_filename)

IF (iflag<2) THEN
   IF (tmin ==1._DP) THEN
      CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext ) 
   ELSE
      CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext ) 
   ENDIF
   CALL gnuplot_xlabel('T (K)', .FALSE.) 
ELSE
   CALL gnuplot_write_header(filenameps, pmin*ry_kbar, pmax*ry_kbar, &
                                              0.0_DP, 0.0_DP, 1.0_DP, flext ) 
   CALL gnuplot_xlabel('p (kbar)', .FALSE.) 
ENDIF

fact= electron_si / (bohr_radius_si)**2 * 1.D4
piezo_type=get_pt_type(code_group)
unit=' (pC/N)'
DO ipt=1,pt_types
   IF (pt_present(ipt, piezo_type)>0) THEN 
      CALL gnuplot_set_fact(fact, .FALSE.)
      CALL plot_one_piezo_element(1, pt_present(ipt, piezo_type)+1, &
              ptd_names(ipt) // TRIM(unit), lpiezo_d, lpiezof_d, filepiezo)
   ENDIF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_piezo_d_t

!---------------------------------------------------------------------
SUBROUTINE plot_piezo_d_pt()
!---------------------------------------------------------------------
!
!  This is a driver to plot the strain piezoelectric tensor 
!  as a function of temperature at a set of pressures.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar, electron_si, bohr_radius_si
USE control_thermo,   ONLY : ltherm_freq, ltherm_dos
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact, gnuplot_ylabel,     &
                             gnuplot_write_file_mul_data
USE color_mod,        ONLY : color
USE data_files,       ONLY : flanhar, fl_el_cons
USE piezoelectric_tensor, ONLY : get_pt_type, pt_present, ptd_names, &
                                 pt_types
USE postscript_files, ONLY : flpsanhar, flps_el_cons
USE control_piezoelectric_tensor,  ONLY : lpiezo_d_pt, lpiezof_d_pt
USE control_pressure, ONLY : pmin, pmax, npress_plot, ipress_plot, press
USE temperature,      ONLY : tmin, tmax
USE color_mod,        ONLY : color
USE thermo_sym,       ONLY : laue
USE rap_point_group,  ONLY : code_group
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
INTEGER :: iflag
LOGICAL :: with_s
CHARACTER(LEN=256) :: gnu_filename, filenameps, filepiezo, filepiezo_ph, label
INTEGER :: piezo_type, ipt, istep, ipressp, ipress
INTEGER :: ierr, system
LOGICAL :: first, last, first_step, last_step
REAL(DP) :: fact

IF ( my_image_id /= root_image ) RETURN
IF (.NOT.(lpiezo_d_pt.OR.lpiezof_d_pt)) RETURN
IF (npress_plot==0) RETURN


gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_d_piezo_p"
filenameps=TRIM(flpsanhar)//".d_piezo_p"//TRIM(flext)

CALL gnuplot_start(gnu_filename)

IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                    1.0_DP, flext ) 
ELSE
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                    1.0_DP, flext ) 
ENDIF
CALL gnuplot_xlabel('T (K)', .FALSE.) 

piezo_type=get_pt_type(code_group)

fact= electron_si / (bohr_radius_si)**2 * 1.D4
DO ipt=1,pt_types
   IF (pt_present(ipt, piezo_type)>0) THEN 
      istep=0
      DO ipressp=1, npress_plot
         first_step=(ipressp==1)
         last_step=(ipressp==npress_plot)
         ipress=ipress_plot(ipressp)
         istep=MOD(istep,8)+1
         filepiezo="anhar_files/"//TRIM(flanhar)//'.d_piezo_press'
         CALL add_value(filepiezo,press(ipress))
         filepiezo_ph="anhar_files/"//TRIM(flanhar)//'.d_piezo_ph_press'
         CALL add_value(filepiezo_ph,press(ipress))
         IF (first_step) THEN
            CALL gnuplot_set_fact(fact, .FALSE.) 
            WRITE(label,'(a," (pC/N)")') TRIM(ptd_names(ipt)) 
            CALL gnuplot_ylabel(TRIM(label),.FALSE.)
         ENDIF
         IF (ltherm_dos) THEN
            CALL gnuplot_write_file_mul_data(filepiezo,1,                &
                          pt_present(ipt, piezo_type)+1,color(istep), &
                          first_step,(last_step.AND..NOT.ltherm_freq),.FALSE.)
         ENDIF
         IF (ltherm_freq) THEN
            CALL gnuplot_write_file_mul_data(filepiezo_ph,1,&
                          pt_present(ipt, piezo_type)+1,color(istep), &
                          (first_step.AND..NOT.ltherm_dos),last_step,.FALSE.)
         ENDIF
      ENDDO
   ENDIF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_piezo_d_pt
!
!---------------------------------------------------------------------
SUBROUTINE plot_piezo_d_ptt()
!---------------------------------------------------------------------
!
!  This is a driver to plot the strain piezoelectric tensor
!  as a function of pressure at a set of temperatures.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar, electron_si, bohr_radius_si
USE control_thermo,   ONLY : ltherm_freq, ltherm_dos
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact, gnuplot_ylabel,     &
                             gnuplot_write_file_mul_data
USE color_mod,        ONLY : color
USE data_files,       ONLY : flanhar, fl_piezo
USE piezoelectric_tensor, ONLY : get_pt_type, pt_present, ptd_names, pt_types
USE postscript_files, ONLY : flpsanhar
USE control_piezoelectric_tensor,  ONLY : lpiezo_d_ptt, lpiezof_d_ptt
USE control_pressure, ONLY : pmin, pmax
USE temperature,      ONLY : temp, ntemp_plot, itemp_plot
USE color_mod,        ONLY : color
USE rap_point_group,  ONLY : code_group
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: gnu_filename, filenameps, filepiezo, filepiezo_ph, label
INTEGER :: pt_type, ipt, istep, itempp, itemp
INTEGER :: ierr, system
LOGICAL :: first, last, first_step, last_step
REAL(DP) :: fact

IF ( my_image_id /= root_image ) RETURN
IF (.NOT.(lpiezo_d_ptt.OR.lpiezof_d_ptt)) RETURN
IF (ntemp_plot==0) RETURN


gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_d_piezo_t"
filenameps=TRIM(flpsanhar)//".d_piezo_t"//TRIM(flext)

CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filenameps, pmin*ry_kbar, pmax*ry_kbar, 0.0_DP, &
                              0.0_DP, 1.0_DP, flext ) 
CALL gnuplot_xlabel('p (kbar)', .FALSE.) 

pt_type=get_pt_type(code_group)
fact= electron_si / (bohr_radius_si)**2 * 1.D4
DO ipt=1,pt_types
   IF (pt_present(ipt, pt_type)>0) THEN 
      istep=0
      DO itempp=1, ntemp_plot
         first_step=(itempp==1)
         last_step=(itempp==ntemp_plot)
         itemp=itemp_plot(itempp)
         istep=MOD(istep,8)+1
         filepiezo="anhar_files/"//TRIM(flanhar)//'.d_piezo_temp'
         CALL add_value(filepiezo,temp(itemp))
         filepiezo_ph="anhar_files/"//TRIM(flanhar)//'.d_piezo_ph_temp'
         CALL add_value(filepiezo_ph,temp(itemp))
         IF (first_step) THEN
            CALL gnuplot_set_fact(fact, .FALSE.) 
            WRITE(label,'(a," (pC/N)")') TRIM(ptd_names(ipt)) 
            CALL gnuplot_ylabel(TRIM(label),.FALSE.)
         ENDIF
         IF (ltherm_dos) THEN
            CALL gnuplot_write_file_mul_data(filepiezo,1,          &
                          pt_present(ipt, pt_type)+1,color(istep), &
                          first_step,(last_step.AND.               &
                           .NOT.ltherm_freq),.FALSE.)
         ENDIF
         IF (ltherm_freq) THEN
            CALL gnuplot_write_file_mul_data(filepiezo_ph,1,&
                          pt_present(ipt, pt_type)+1,color(istep), &
                          (first_step.AND..NOT.ltherm_dos),&
                          last_step,.FALSE.)
         ENDIF
      ENDDO
   ENDIF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_piezo_d_ptt

