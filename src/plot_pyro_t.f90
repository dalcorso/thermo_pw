!
! Copyright (C) 2025 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE plot_pyro_t(iflag)
!---------------------------------------------------------------------
!
!  This is a driver to plot the pyroelectric tensor as a function of
!  temperature or pressure
!  iflag    0   pyroelectric tensor as a function of temperature
!           2   pyroelectric tensor as a function of pressure
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar, electron_si, bohr_radius_si
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar, fl_piezo
USE postscript_files, ONLY : flpsanhar
USE control_pyroelectric_tensor,  ONLY : lpyro, lpyrof
USE polarization_vector, ONLY : get_py_type, py_present, py_names, py_types, &
                                py_elements
USE control_pressure, ONLY : pmin, pmax
USE temperature,      ONLY : tmin, tmax
USE rap_point_group,  ONLY : code_group
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
INTEGER :: iflag
CHARACTER(LEN=256) :: gnu_filename, filenameps, filepyro
CHARACTER(LEN=15) :: unit
INTEGER :: ierr, system, i, j, ipt, pyro_type
REAL(DP) :: fact

IF ( my_image_id /= root_image ) RETURN

IF (.NOT.(lpyro.OR.lpyrof)) RETURN

IF (iflag==0) THEN
   gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_pyro"
   filenameps=TRIM(flpsanhar)//".pyro"//TRIM(flext)
   filepyro="anhar_files/"//TRIM(flanhar)//".pyro"
ELSE
!
!  In this case plot the elastic constants as a function of pressure
!
   gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_pyro_p"
   filenameps=TRIM(flpsanhar)//".pyro_p"//TRIM(flext)
   filepyro="elastic_constants/"//TRIM(fl_piezo)//".pyro_p"
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

fact= electron_si * 1.D6 / (bohr_radius_si)**2
pyro_type=get_py_type(code_group)
unit='x10^6 (C/m^2/K)'
DO ipt=1,py_elements
   IF (py_present(ipt, pyro_type)>0) THEN 
      CALL gnuplot_set_fact(fact, .FALSE.)
      CALL plot_one_pyro_element(1, py_present(ipt, pyro_type)+1,   &
              'Pyroelectric tensor '//TRIM(py_names(ipt)) //        &
                               TRIM(unit), lpyro, lpyrof, filepyro)
   ENDIF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_pyro_t
!
!---------------------------------------------------------------------
SUBROUTINE plot_pyro_pt(iflag)
!---------------------------------------------------------------------
!
!  This is a driver to plot the pyroelectric tensor as a function of
!  temperature or pressure
!  iflag    0   pyroelectric tensor as a function of temperature
!           2   pyroelectric tensor as a function of pressure
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar, bohr_radius_si, electron_si
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact, gnuplot_ylabel, &
                             gnuplot_write_file_mul_data, &
                             gnuplot_write_file_mul_data_sum
USE data_files,       ONLY : flanhar, fl_piezo
USE postscript_files, ONLY : flpsanhar
USE control_pyroelectric_tensor,  ONLY : lpyro_pt, lpyrof_pt
USE polarization_vector, ONLY : get_py_type, py_present, py_names, py_types, &
                                py_elements
USE control_pressure, ONLY : pmin, pmax, press, npress_plot, ipress_plot
USE control_thermo,   ONLY : ltherm_freq, ltherm_dos
USE temperature,      ONLY : tmin, tmax
USE rap_point_group,  ONLY : code_group
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
INTEGER :: iflag
CHARACTER(LEN=256) :: gnu_filename, filenameps, filepyro, filepyro_ph, label
CHARACTER(LEN=15) :: aunit
INTEGER :: ierr, system, i, j, ipt, pyro_type, ipressp, ipress, istep
LOGICAL :: first_step, last_step
REAL(DP) :: fact

IF ( my_image_id /= root_image ) RETURN

IF (.NOT.(lpyro_pt.OR.lpyrof_pt)) RETURN

IF (iflag==0) THEN
   gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_pyro_p"
   filenameps=TRIM(flpsanhar)//".pyro_p"//TRIM(flext)
ELSE
!
!  In this case plot the elastic constants as a function of pressure
!
   gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_pyro_p"
   filenameps=TRIM(flpsanhar)//".pyro_p"//TRIM(flext)
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

fact= electron_si * 1.D6 / (bohr_radius_si)**2
pyro_type=get_py_type(code_group)
aunit='x10^6 (C/m^2/K)'
DO ipt=1,py_elements
   IF (py_present(ipt, pyro_type)>0) THEN 
      i=1
      j=py_present(ipt, pyro_type)+1 
      istep=0
      DO ipressp=1, npress_plot
         first_step=(ipressp==1)
         last_step=(ipressp==npress_plot)
         ipress=ipress_plot(ipressp)
         filepyro="anhar_files/"//TRIM(flanhar)//'.pyro_press'
         CALL add_value(filepyro,press(ipress))
         filepyro_ph="anhar_files/"//TRIM(flanhar)//'.pyro_ph_press'
         CALL add_value(filepyro_ph,press(ipress))
         IF (first_step) THEN
            CALL gnuplot_set_fact(fact, .FALSE.)
            WRITE(label,'(a)') "Pyroelectric tensor "//TRIM(py_names(ipt)) &
                                                     // TRIM(aunit)
            CALL gnuplot_ylabel(TRIM(label),.FALSE.)
         ENDIF
         IF (ltherm_dos) THEN
            CALL gnuplot_write_file_mul_data(filepyro,i,j,'color_green',     &
                                        first_step, .FALSE.,.FALSE.)
            CALL gnuplot_write_file_mul_data(filepyro,i,j+py_elements,       &
                               'color_orange',.FALSE.,.FALSE.,.FALSE.)
            CALL gnuplot_write_file_mul_data_sum(filepyro,i,j,j+py_elements, &
                 'color_blue',.FALSE.,last_step.AND..NOT.ltherm_freq,.FALSE.)
         ENDIF
         IF (ltherm_freq) THEN
            CALL gnuplot_write_file_mul_data(filepyro_ph,i,j,'color_green',   &
                              first_step.AND..NOT.ltherm_dos, .FALSE.,.FALSE.)
            CALL gnuplot_write_file_mul_data(filepyro_ph,i,j+py_elements,     &
                               'color_orange',.FALSE.,.FALSE.,.FALSE.)
            CALL gnuplot_write_file_mul_data_sum(filepyro_ph,i,j,j+py_elements,&
                 'color_blue',.FALSE.,last_step,.FALSE.)
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
END SUBROUTINE plot_pyro_pt
!
!---------------------------------------------------------------------
SUBROUTINE plot_pyro_ptt(iflag)
!---------------------------------------------------------------------
!
!  This is a driver to plot the pyroelectric tensor as a function of
!  temperature or pressure
!  iflag    0   pyroelectric tensor as a function of temperature
!           2   pyroelectric tensor as a function of pressure
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar, bohr_radius_si, electron_si
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact, gnuplot_ylabel, &
                             gnuplot_write_file_mul_data, &
                             gnuplot_write_file_mul_data_sum
USE data_files,       ONLY : flanhar, fl_piezo
USE postscript_files, ONLY : flpsanhar
USE control_pyroelectric_tensor,  ONLY : lpyro_ptt, lpyrof_ptt
USE polarization_vector, ONLY : get_py_type, py_present, py_names, py_types, &
                                py_elements
USE control_pressure, ONLY : pmin, pmax, deltap
USE control_thermo,   ONLY : ltherm_freq, ltherm_dos
USE temperature,      ONLY : tmin, tmax, temp, ntemp_plot, itemp_plot
USE rap_point_group,  ONLY : code_group
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
INTEGER :: iflag
CHARACTER(LEN=256) :: gnu_filename, filenameps, filepyro, filepyro_ph, label
CHARACTER(LEN=15) :: aunit
INTEGER :: ierr, system, i, j, ipt, pyro_type, itempp, itemp, istep
LOGICAL :: first_step, last_step
REAL(DP) :: fact

IF ( my_image_id /= root_image ) RETURN

IF (.NOT.(lpyro_ptt.OR.lpyrof_ptt)) RETURN

IF (iflag==0) THEN
   gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_pyro_t"
   filenameps=TRIM(flpsanhar)//".pyro_t"//TRIM(flext)
ELSE
!
!  In this case plot the pyroelectric tensor as a function of pressure
!
   gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_pyro_t"
   filenameps=TRIM(flpsanhar)//".pyro_t"//TRIM(flext)
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
   CALL gnuplot_write_header(filenameps, (pmin+deltap)*ry_kbar, &
                   (pmax-deltap)*ry_kbar, 0.0_DP, 0.0_DP, 1.0_DP, flext ) 
   CALL gnuplot_xlabel('p (kbar)', .FALSE.) 
ENDIF

fact= electron_si * 1.D6 / (bohr_radius_si)**2
pyro_type=get_py_type(code_group)
aunit='x10^6 (C/m^2/K)'
DO ipt=1,py_elements
   IF (py_present(ipt, pyro_type)>0) THEN 
      i=1
      j=py_present(ipt, pyro_type)+1 
      istep=0
      DO itempp=1, ntemp_plot
         first_step=(itempp==1)
         last_step=(itempp==ntemp_plot)
         itemp=itemp_plot(itempp)
         filepyro="anhar_files/"//TRIM(flanhar)//'.pyro_temp'
         CALL add_value(filepyro,temp(itemp))
         filepyro_ph="anhar_files/"//TRIM(flanhar)//'.pyro_ph_temp'
         CALL add_value(filepyro_ph,temp(itemp))
         IF (first_step) THEN
            CALL gnuplot_set_fact(fact, .FALSE.)
            WRITE(label,'(a)') "Pyroelectric tensor "//TRIM(py_names(ipt)) &
                                                     // TRIM(aunit)
            CALL gnuplot_ylabel(TRIM(label),.FALSE.)
         ENDIF
         IF (ltherm_dos) THEN
            CALL gnuplot_write_file_mul_data(filepyro,i,j,'color_green',     &
                                        .TRUE., .FALSE.,.FALSE.)
            CALL gnuplot_write_file_mul_data(filepyro,i,j+py_elements,       &
                               'color_orange',.FALSE.,.FALSE.,.FALSE.)
            CALL gnuplot_write_file_mul_data_sum(filepyro,i,j,j+py_elements, &
                 'color_blue',.FALSE.,.NOT.ltherm_freq,.FALSE.)
         ENDIF
         IF (ltherm_freq) THEN
            CALL gnuplot_write_file_mul_data(filepyro_ph,i,j,'color_green',   &
                              .NOT.ltherm_dos, .FALSE.,.FALSE.)
            CALL gnuplot_write_file_mul_data(filepyro_ph,i,j+py_elements,     &
                               'color_orange',.FALSE.,.FALSE.,.FALSE.)
            CALL gnuplot_write_file_mul_data_sum(filepyro_ph,i,j, &
                  j+py_elements,'color_blue',.FALSE.,.TRUE.,.FALSE.)
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
END SUBROUTINE plot_pyro_ptt

!---------------------------------------------------------------------
SUBROUTINE plot_one_pyro_element(i, j, label, lpyro, lpyrof, filepyro)
!---------------------------------------------------------------------
USE thermo_mod, ONLY : what
USE control_thermo,   ONLY : ltherm_freq, ltherm_dos
USE gnuplot, ONLY : gnuplot_ylabel, &
                    gnuplot_write_file_mul_data, &
                    gnuplot_write_file_mul_data_sum
USE color_mod, ONLY : color
USE control_elastic_constants, ONLY : ngeom, all_geometry_done_geo
USE polarization_vector, ONLY : py_elements
IMPLICIT NONE
INTEGER, INTENT(IN) :: i, j
LOGICAL, INTENT(IN) :: lpyro, lpyrof
CHARACTER(LEN=*) :: label
CHARACTER(LEN=256), INTENT(IN) :: filepyro

CHARACTER(LEN=256) :: filename, filename_ph
CHARACTER(LEN=6) :: int_to_char
INTEGER :: last_ngeom, first_ngeom, igeom, ic

filename=TRIM(filepyro)
filename_ph=TRIM(filepyro)//'_ph'

CALL gnuplot_ylabel(TRIM(label),.FALSE.) 

IF (lpyro.OR.lpyrof) THEN
   last_ngeom=1
   first_ngeom=0
   DO igeom=1,ngeom
      IF (all_geometry_done_geo(igeom).AND.first_ngeom==0) first_ngeom=igeom
      IF (all_geometry_done_geo(igeom)) last_ngeom=igeom
   ENDDO
ENDIF

IF (ltherm_dos) THEN
   IF (what=='polarization_geo') THEN
      DO igeom=1, ngeom
         IF (.NOT.all_geometry_done_geo(igeom)) CYCLE
         filename=TRIM(filepyro)//".g"//TRIM(int_to_char(igeom))
         ic=MOD(igeom-1,8)+1
         CALL gnuplot_write_file_mul_data(filename,i,j,color(ic),&
               (igeom==first_ngeom), .FALSE.,.FALSE.)
         CALL gnuplot_write_file_mul_data(filename,i,j+py_elements,color(ic),&
               .FALSE., (igeom==last_ngeom).AND. &
                                          .NOT.ltherm_freq,.FALSE.)
      ENDDO
   ELSE
      CALL gnuplot_write_file_mul_data(filename,i,j,'color_red',.TRUE., &
                                                  .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename,i,j+py_elements,&
                               'color_red',.FALSE.,.NOT.ltherm_freq,.FALSE.)
      CALL gnuplot_write_file_mul_data_sum(filename,i,j,j+py_elements, &
                                  'color_red',.FALSE.,.TRUE.,.FALSE.)
   ENDIF
ENDIF
IF (ltherm_freq) THEN
   IF (what=='polarization_geo') THEN
      DO igeom=1, ngeom
         IF (.NOT.all_geometry_done_geo(igeom)) CYCLE
         filename_ph=TRIM(filepyro)//".g"//TRIM(int_to_char(igeom))//"_ph"
         ic=MOD(igeom-1,8)+1
         CALL gnuplot_write_file_mul_data(filename_ph,i,j,color(ic), &
             .NOT.ltherm_dos.AND.(igeom==first_ngeom), .FALSE.,&
             .FALSE.)
         CALL gnuplot_write_file_mul_data(filename_ph,i,j+py_elements, &
              color(ic), .FALSE., (igeom==last_ngeom), .FALSE.)
      ENDDO
   ELSE
      CALL gnuplot_write_file_mul_data(filename_ph,i,j,'color_green',     &
                                            .NOT.ltherm_dos,.FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename_ph,i,j+py_elements, &
                                  'color_red',.FALSE.,.FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data_sum(filename_ph,i,j,j+py_elements, &
                                  'color_blue',.FALSE.,.TRUE.,.FALSE.)
   ENDIF
ENDIF

RETURN
END SUBROUTINE plot_one_pyro_element
