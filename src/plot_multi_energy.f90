!
! Copyright (C) 2014-2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE plot_multi_energy()
  !-----------------------------------------------------------------------
  !
  !  this routine writes a gnuplot script to make one or more 
  !  two dimensional countour plots of the energy as a function
  !  of the structural parameters.
  !
  !
  USE kinds,                ONLY : DP
  USE thermo_mod,           ONLY : ngeo, celldm_geo, energy_geo, omega_geo
  USE initial_conf,         ONLY : ibrav_save
  USE control_gnuplot,      ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
  USE data_files,           ONLY : flenergy, flevdat, flgeom
  USE postscript_files,     ONLY : flpsenergy
  USE control_thermo,       ONLY : lgeo_to_file
  USE control_energy_plot,  ONLY : ncontours, ene_levels, color_levels
  USE control_quadratic_energy, ONLY : x_pos_min, hessian_v, nvar, show_fit
  USE control_quartic_energy, ONLY : x_min_4, hessian4_v, lquartic
  USE control_vol,          ONLY : nvol
  USE control_pressure,     ONLY : pressure, pressure_kb
  USE mp_images,            ONLY : my_image_id, root_image
  USE io_global,            ONLY : ionode, stdout
  USE color_mod,            ONLY : color
  USE gnuplot,              ONLY : gnuplot_start, gnuplot_end,             &
                                   gnuplot_start_2dplot,                   &
                                   gnuplot_set_contour, gnuplot_do_2dplot, &
                                   gnuplot_xlabel,     &
                                   gnuplot_write_header,                   &
                                   gnuplot_write_file_mul_data,            &
                                   gnuplot_write_file_mul_point,           &
                                   gnuplot_write_horizontal_line,          &
                                   gnuplot_write_file_mul_line_point,      &
                                   gnuplot_ylabel, gnuplot_close_2dplot_prep, &
                                   gnuplot_line_v

  IMPLICIT NONE
  CHARACTER(LEN=256) :: gnu_filename, filename, filename1, label, filenameps, &
                        tablefile, fileout, filename2
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=8) :: float_to_char
  CHARACTER(LEN=50) :: xlabel, ylabel
  REAL(DP) :: emax, emin, deltae, ene_levels_int(ncontours)
  REAL(DP) :: xmin, xmax, ymin, ymax, x2
  INTEGER :: nx, ny, icont, ifile, tot_n, iwork
  INTEGER :: compute_nwork
  INTEGER :: ierr, system

  IF ( my_image_id /= root_image ) RETURN

  gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_energy'
  CALL add_pressure(gnu_filename)
  IF (show_fit) THEN
     filename='energy_files/'//TRIM(flevdat)//'_quadratic'
  ELSE
     filename='energy_files/'//TRIM(flenergy)//int_to_char(1)
  ENDIF
  CALL add_pressure(filename)
  filenameps=TRIM(flpsenergy)
  CALL add_pressure(filenameps)
  filenameps=TRIM(filenameps)//TRIM(flext)

  tablefile='gnuplot_files/table'
  CALL gnuplot_start(gnu_filename)

  SELECT CASE (ibrav_save) 
     CASE(1,2,3)
        nx=ngeo(1)
        xmax=0.0_DP
        xmin=1000000.0_DP
        DO iwork = 1, nx
           IF (celldm_geo(1,iwork) > xmax) xmax=celldm_geo(1,iwork)
           IF (celldm_geo(1,iwork) < xmin) xmin=celldm_geo(1,iwork)
        ENDDO
        xmin=xmin*0.99_DP
        xmax=xmax*1.01_DP
        CALL gnuplot_write_header(filenameps, xmin, xmax, 0.0_DP, 0.0_DP, &
                                                          1.0_DP, flext)
        CALL gnuplot_xlabel('a (a.u.)',.FALSE.)

        IF (pressure /= 0.0_DP) THEN
           label='Enthalpy (Ry)    p= '&
                             &//TRIM(float_to_char(pressure_kb,1))//' kbar'
        ELSE
           label='Energy (Ry)'
        END IF
        CALL gnuplot_ylabel(TRIM(label),.FALSE.)
        filename1='energy_files/'//TRIM(flevdat)//'_quadratic'
        CALL add_pressure(filename1)
        CALL gnuplot_write_file_mul_data(filename1,1,2,'color_red',.TRUE., &
                                                          .FALSE.,.FALSE.)
        CALL gnuplot_write_file_mul_point(filename,1,2,'color_red',      &
                                                   .FALSE.,.TRUE.,.FALSE.)

        CALL gnuplot_ylabel('Pressure (kbar)',.FALSE.)
        IF (pressure /= 0.0_DP) &
            CALL gnuplot_write_horizontal_line(pressure_kb, 2, 'front', &
                                                  'color_green',.FALSE.)
        CALL gnuplot_write_horizontal_line(0.0_DP, 2, 'front', 'color_black',&
                                                                     .FALSE.)
        CALL gnuplot_write_file_mul_data(filename1,1,3,'color_red',.TRUE.,&
                                                              .TRUE.,.FALSE.)

     CASE(4,5,6,7)
        IF (ncontours==0) RETURN
        ene_levels_int(:)=ene_levels(:)
        IF (show_fit) THEN
           nx=nvol
           ny=nvol
        ELSE
           nx=ngeo(1)
           IF (ibrav_save==5) THEN
              ny=ngeo(4)
           ELSE
              ny=ngeo(3)
           ENDIF
        ENDIF
        tot_n=compute_nwork()
        xmin=1.d10
        xmax=-1.d10
        ymin=1.d10
        ymax=-1.d10
        DO iwork = 1, tot_n
           IF (celldm_geo(1,iwork) > xmax) xmax=celldm_geo(1,iwork)  
           IF (celldm_geo(1,iwork) < xmin) xmin=celldm_geo(1,iwork)  
           IF (ibrav_save==5) THEN
              IF (celldm_geo(4,iwork) > ymax) ymax=celldm_geo(4,iwork)  
              IF (celldm_geo(4,iwork) < ymin) ymin=celldm_geo(4,iwork)  
           ELSE
              IF (celldm_geo(3,iwork) > ymax) ymax=celldm_geo(3,iwork)  
              IF (celldm_geo(3,iwork) < ymin) ymin=celldm_geo(3,iwork)  
           ENDIF
        END DO
        IF (ene_levels(1)==-1000._DP) THEN
           emin=1.d10
           emax=-1.d10
           DO iwork = 1, tot_n
              IF ( energy_geo(iwork) + pressure * omega_geo(iwork) > emax ) &
                 emax = energy_geo(iwork) + pressure * omega_geo(iwork)
              IF ( energy_geo(iwork) + pressure * omega_geo(iwork) < emin ) &
                 emin = energy_geo(iwork) + pressure * omega_geo(iwork)
           ENDDO
!
!    emax and emin are not used as contours 
!
           deltae = (emax-emin) / (ncontours+1)
           DO icont = 1, ncontours
              ene_levels_int(icont) = emin + icont * deltae
              color_levels(icont) = color(MOD((icont-1)/3,8)+1)
           END DO
        END IF

        WRITE(stdout,'(/,5x,"The plot will have ",i5," levels")') ncontours
        DO icont=1, ncontours
           WRITE(stdout,'(5x,"Level ",i5," Energy= ",f15.8," ", a)') icont, &
                             ene_levels_int(icont), TRIM(color_levels(icont))
        ENDDO

        CALL gnuplot_start_2dplot(ncontours, nx, ny)
        DO icont=1,ncontours
           CALL gnuplot_set_contour(filename,ene_levels_int(icont), &
                                            color_levels(icont),tablefile)
        ENDDO
        CALL gnuplot_close_2dplot_prep()
        xlabel='a (a.u.)'
        IF (ibrav_save==5) THEN
           ylabel='cos({/Symbol a})'
        ELSE
           ylabel='c/a '
        ENDIF

        CALL gnuplot_do_2dplot(filenameps, xmin, xmax, ymin, ymax, xlabel, &
                                                  ylabel, tablefile, flext)

        IF (lquartic) THEN
           x2=x_min_4(2)
           IF (ibrav_save==5) x2=COS(x_min_4(2))
           CALL gnuplot_line_v(hessian4_v(1,1), hessian4_v(2,1),              &
                                   x_min_4(1), x2, .FALSE., 'color_blue')
           CALL gnuplot_line_v(hessian4_v(1,2), hessian4_v(2,2), x_min_4(1),  &
                                       x2, .FALSE., 'color_blue')
        ELSE
           x2=x_pos_min(2)
           IF (ibrav_save==5) x2=COS(x_pos_min(2))
           CALL gnuplot_line_v(hessian_v(1,1), hessian_v(2,1), x_pos_min(1),  &
                                       x2,.FALSE.,'color_blue')
           CALL gnuplot_line_v(hessian_v(1,2), hessian_v(2,2), x_pos_min(1),  &
                                       x2,.FALSE.,'color_blue')
        ENDIF
        IF (lgeo_to_file) THEN
!
!   In this case a set of geometries is written on file. Here we plot
!   this set with points in the energy contour plot.
!
           filename2='./'//TRIM(flgeom)//'.dat'
           CALL gnuplot_write_file_mul_line_point(filename2, 3, 5, &
                           'color_orange', .FALSE., .TRUE., .TRUE., .FALSE.)

        ENDIF   
     CASE (8,9,91,10,11) 
        IF (ncontours==0) RETURN
        ene_levels_int(:)=ene_levels(:)
        IF (show_fit) THEN
           nx=nvol
           ny=nvol
        ELSE
           nx=ngeo(1)
           ny=ngeo(2)
        ENDIF
        tot_n=compute_nwork()
        xmin=1.d10
        xmax=-1.d10
        ymin=1.d10
        ymax=-1.d10
        DO iwork = 1, tot_n
           IF (celldm_geo(1,iwork) > xmax) xmax=celldm_geo(1,iwork)  
           IF (celldm_geo(1,iwork) < xmin) xmin=celldm_geo(1,iwork)  
           IF (celldm_geo(2,iwork) > ymax) ymax=celldm_geo(2,iwork)  
           IF (celldm_geo(2,iwork) < ymin) ymin=celldm_geo(2,iwork)  
        END DO
        IF (ene_levels(1)==-1000._DP) THEN
           emin=1.d10
           emax=-1.d10
           DO iwork = 1, tot_n
              IF ( energy_geo(iwork) + pressure * omega_geo(iwork) > emax ) &
                 emax = energy_geo(iwork) + pressure * omega_geo(iwork)
              IF ( energy_geo(iwork) + pressure * omega_geo(iwork) < emin ) &
                 emin = energy_geo(iwork) + pressure * omega_geo(iwork)
           ENDDO
!
!    emax and emin are not used as contours 
!
           deltae = (emax-emin) / (ncontours+1)
           DO icont = 1, ncontours
              ene_levels_int(icont) = emin + icont * deltae
              color_levels(icont) = color(MOD((icont-1)/3,8)+1)
           END DO
        END IF

        WRITE(stdout,'(/,5x,"The plot will have ",i5," levels")') ncontours
        DO icont=1, ncontours
           WRITE(stdout,'(5x,"Level ",i5," Energy= ",f15.8," ", a)') icont, &
                          ene_levels_int(icont), TRIM(color_levels(icont))
        ENDDO
        DO ifile=1,ngeo(3)
           filenameps=TRIM(flpsenergy)//int_to_char(ifile)
           CALL add_pressure(filenameps)
           fileout='energy_files/'//TRIM(flenergy)//int_to_char(ifile)
           CALL add_pressure(fileout)
           filenameps=TRIM(filenameps)//TRIM(flext)
           CALL gnuplot_start_2dplot(ncontours, nx, ny)
           DO icont=1,ncontours
              CALL gnuplot_set_contour(fileout,ene_levels_int(icont), &
                                            color_levels(icont),tablefile)
           ENDDO
           CALL gnuplot_close_2dplot_prep()
           xlabel='a (a.u.)'
           ylabel='b/a '
           CALL gnuplot_do_2dplot(filenameps, xmin, xmax, ymin, ymax, xlabel, &
                                                  ylabel, tablefile, flext)
        ENDDO
     CASE DEFAULT
        RETURN
  END SELECT

  CALL gnuplot_end()

  IF (lgnuplot.AND.ionode) &
     ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!  IF (lgnuplot.AND.ionode) &
!     CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

  RETURN
  END SUBROUTINE plot_multi_energy
!
!-----------------------------------------------------------------------
SUBROUTINE plot_multi_energy_t()
  !-----------------------------------------------------------------------
  !
  !  this routine writes a gnuplot script to make one or more 
  !  two dimensional countour plots of the free_energy as a function
  !  of the structural parameters. This plot is made at the ntemp_plot
  !  temperatures requested in input. The routine shows also the 
  !  pressure path at that temperature.
  !  At the end also a contour plot of the energy with the pressure paths
  !  at each of the ntemp_plot temperatures is shown.
  !
  !
  USE kinds,                ONLY : DP
  USE thermo_mod,           ONLY : ngeo, celldm_geo, energy_geo, omega_geo
  USE initial_conf,         ONLY : ibrav_save
  USE control_gnuplot,      ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
  USE data_files,           ONLY : flenergy, flevdat, flgeom, flanhar
  USE postscript_files,     ONLY : flpsenergy
  USE control_thermo,       ONLY : lgeo_to_file
  USE control_energy_plot,  ONLY : ncontours, ene_levels, color_levels
  USE control_quadratic_energy, ONLY : x_pos_min, hessian_v, nvar, show_fit
  USE control_quartic_energy, ONLY : x_min_4, hessian4_v, lquartic
  USE control_vol,          ONLY : nvol
  USE control_pressure,     ONLY : pressure, pressure_kb, press, npress_plot, &
                                   ipress_plot
  USE temperature,          ONLY : ntemp_plot, itemp_plot, temp
  USE mp_images,            ONLY : my_image_id, root_image
  USE io_global,            ONLY : ionode, stdout
  USE color_mod,            ONLY : color
  USE gnuplot,              ONLY : gnuplot_start, gnuplot_end,             &
                                   gnuplot_start_2dplot,                   &
                                   gnuplot_set_contour, gnuplot_do_2dplot, &
                                   gnuplot_xlabel,     &
                                   gnuplot_write_header,                   &
                                   gnuplot_write_file_mul_data,            &
                                   gnuplot_write_file_mul_point,           &
                                   gnuplot_write_horizontal_line,          &
                                   gnuplot_write_command,                  &
                                   gnuplot_ylabel, gnuplot_close_2dplot_prep, &
                                   gnuplot_line_v

  IMPLICIT NONE
  CHARACTER(LEN=256) :: gnu_filename, filename, filename1, label, filenameps, &
                        tablefile, fileout, filename2
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=8) :: float_to_char
  CHARACTER(LEN=50) :: xlabel, ylabel
  REAL(DP) :: emax, emin, deltae, ene_levels_int(ncontours)
  REAL(DP) :: xmin, xmax, ymin, ymax, x2
  LOGICAL :: first_step, last_step
  INTEGER :: nx, ny, icont, ifile, tot_n, iwork, itempp, itemp, istep
  INTEGER :: ipressp, ipress
  INTEGER :: compute_nwork
  INTEGER :: ierr, system

  IF ( my_image_id /= root_image ) RETURN

  gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_energy_t'
  CALL add_pressure(gnu_filename)
  IF (show_fit) THEN
     filename='energy_files/'//TRIM(flevdat)//'_quadratic'
  ELSE
     filename='energy_files/'//TRIM(flenergy)//int_to_char(1)
  ENDIF
  CALL add_pressure(filename)
  filenameps=TRIM(flpsenergy)//'_t'
  CALL add_pressure(filenameps)
  filenameps=TRIM(filenameps)//TRIM(flext)

  tablefile='gnuplot_files/table_t'
  CALL gnuplot_start(gnu_filename)

  SELECT CASE (ibrav_save) 
     CASE(1,2,3)
        nx=ngeo(1)
        xmax=0.0_DP
        xmin=1000000.0_DP
        DO iwork = 1, nx
           IF (celldm_geo(1,iwork) > xmax) xmax=celldm_geo(1,iwork)
           IF (celldm_geo(1,iwork) < xmin) xmin=celldm_geo(1,iwork)
        ENDDO
        xmin=xmin*0.99_DP
        xmax=xmax*1.01_DP
        CALL gnuplot_write_header(filenameps, xmin, xmax, 0.0_DP, 0.0_DP, &
                                                          1.0_DP, flext)
        CALL gnuplot_xlabel('a (a.u.)',.FALSE.)

        IF (pressure /= 0.0_DP) THEN
           label='Enthalpy (Ry)    p= '&
                             &//TRIM(float_to_char(pressure_kb,1))//' kbar'
        ELSE
           label='Energy (Ry)'
        END IF
        CALL gnuplot_ylabel(TRIM(label),.FALSE.)
        filename1='energy_files/'//TRIM(flevdat)//'_quadratic'
        CALL add_pressure(filename1)
        CALL gnuplot_write_file_mul_data(filename1,1,2,'color_red',.TRUE., &
                                                          .FALSE.,.FALSE.)
        CALL gnuplot_write_file_mul_point(filename,1,2,'color_red',      &
                                                   .FALSE.,.TRUE.,.FALSE.)

        CALL gnuplot_ylabel('Pressure (kbar)',.FALSE.)
        IF (pressure /= 0.0_DP) &
            CALL gnuplot_write_horizontal_line(pressure_kb, 2, 'front', &
                                                  'color_green',.FALSE.)
        CALL gnuplot_write_horizontal_line(0.0_DP, 2, 'front', 'color_black',&
                                                                     .FALSE.)
        CALL gnuplot_write_file_mul_data(filename1,1,3,'color_red',.TRUE.,&
                                                              .TRUE.,.FALSE.)

     CASE(4,5,6,7)
        IF (ncontours==0) RETURN
        ene_levels_int(:)=ene_levels(:)
        IF (show_fit) THEN
           nx=nvol
           ny=nvol
        ELSE
           nx=ngeo(1)
           IF (ibrav_save==5) THEN
              ny=ngeo(4)
           ELSE
              ny=ngeo(3)
           ENDIF
        ENDIF
        tot_n=compute_nwork()
        xmin=1.d10
        xmax=-1.d10
        ymin=1.d10
        ymax=-1.d10
        DO iwork = 1, tot_n
           IF (celldm_geo(1,iwork) > xmax) xmax=celldm_geo(1,iwork)  
           IF (celldm_geo(1,iwork) < xmin) xmin=celldm_geo(1,iwork)  
           IF (ibrav_save==5) THEN
              IF (celldm_geo(4,iwork) > ymax) ymax=celldm_geo(4,iwork)  
              IF (celldm_geo(4,iwork) < ymin) ymin=celldm_geo(4,iwork)  
           ELSE
              IF (celldm_geo(3,iwork) > ymax) ymax=celldm_geo(3,iwork)  
              IF (celldm_geo(3,iwork) < ymin) ymin=celldm_geo(3,iwork)  
           ENDIF
        END DO
        IF (ene_levels(1)==-1000._DP) THEN
           emin=1.d10
           emax=-1.d10
           DO iwork = 1, tot_n
              IF ( energy_geo(iwork) + pressure * omega_geo(iwork) > emax ) &
                 emax = energy_geo(iwork) + pressure * omega_geo(iwork)
              IF ( energy_geo(iwork) + pressure * omega_geo(iwork) < emin ) &
                 emin = energy_geo(iwork) + pressure * omega_geo(iwork)
           ENDDO
!
!    emax and emin are not used as contours 
!
           deltae = (emax-emin) / (ncontours+1)
           DO icont = 1, ncontours
              ene_levels_int(icont) = emin + icont * deltae
              color_levels(icont) = color(MOD((icont-1)/3,8)+1)
           END DO
        END IF

        WRITE(stdout,'(/,5x,"The plot will have ",i5," levels")') ncontours
        DO icont=1, ncontours
           WRITE(stdout,'(5x,"Level ",i5," Energy= ",f15.8," ", a)') icont, &
                             ene_levels_int(icont), TRIM(color_levels(icont))
        ENDDO

        CALL gnuplot_start_2dplot(ncontours, nx, ny)
        DO icont=1,ncontours
           CALL gnuplot_set_contour(filename,ene_levels_int(icont), &
                                            color_levels(icont),tablefile)
        ENDDO
        CALL gnuplot_close_2dplot_prep()
        xlabel='a (a.u.)'
        IF (ibrav_save==5) THEN
           ylabel='cos({/Symbol a})'
        ELSE
           ylabel='c/a '
        ENDIF

        CALL gnuplot_do_2dplot(filenameps, xmin, xmax, ymin, ymax, xlabel, &
                                                  ylabel, tablefile, flext)

        CALL gnuplot_write_command('set size ratio 1.0',.FALSE.)
        IF (lquartic) THEN
           x2=x_min_4(2)
           IF (ibrav_save==5) x2=COS(x_min_4(2))
           CALL gnuplot_line_v(hessian4_v(1,1), hessian4_v(2,1),              &
                                   x_min_4(1), x2, .FALSE., 'color_blue')
           CALL gnuplot_line_v(hessian4_v(1,2), hessian4_v(2,2), x_min_4(1),  &
                                       x2, .FALSE., 'color_blue')
        ELSE
           x2=x_pos_min(2)
           IF (ibrav_save==5) x2=COS(x_pos_min(2))
           CALL gnuplot_line_v(hessian_v(1,1), hessian_v(2,1), x_pos_min(1),  &
                                       x2,.FALSE.,'color_blue')
           CALL gnuplot_line_v(hessian_v(1,2), hessian_v(2,2), x_pos_min(1),  &
                                       x2,.FALSE.,'color_blue')
        ENDIF
        IF (ntemp_plot>0) THEN
!
!   In this case for each temp_plot there is file with the celldm. Here we plot
!   this set with lines in the energy contour plot.
!
           istep=0
           DO itempp=1, ntemp_plot
              first_step=(itempp==1)
              last_step=(itempp==ntemp_plot)
              itemp=itemp_plot(itempp)
              istep=MOD(istep,8)+1
              filename2='anhar_files/'//TRIM(flanhar)//'.celldm_temp'
              CALL add_value(filename2,temp(itemp))
              CALL gnuplot_write_file_mul_data(filename2, 2, 3, &
                        color(istep),first_step, last_step, .FALSE., .TRUE.)
           ENDDO
        ENDIF   
        IF (npress_plot>0) THEN
!
!   In this case for each press_plot there is file with the celldm. 
!   Here we plot this set with lines in the energy contour plot.
!
           istep=0
           DO ipressp=1, npress_plot
              first_step=(ipressp==1)
              last_step=(ipressp==npress_plot)
              ipress=ipress_plot(ipressp)
              istep=MOD(istep,8)+1
              filename2='anhar_files/'//TRIM(flanhar)//'.celldm_press'
              CALL add_value(filename2,press(ipress))
              CALL gnuplot_write_file_mul_data(filename2, 2, 3, &
                        color(istep),first_step, last_step, .FALSE., .TRUE.)
           ENDDO
        ENDIF   
     CASE (8,9,91,10,11) 
        IF (ncontours==0) RETURN
        ene_levels_int(:)=ene_levels(:)
        IF (show_fit) THEN
           nx=nvol
           ny=nvol
        ELSE
           nx=ngeo(1)
           ny=ngeo(2)
        ENDIF
        tot_n=compute_nwork()
        xmin=1.d10
        xmax=-1.d10
        ymin=1.d10
        ymax=-1.d10
        DO iwork = 1, tot_n
           IF (celldm_geo(1,iwork) > xmax) xmax=celldm_geo(1,iwork)  
           IF (celldm_geo(1,iwork) < xmin) xmin=celldm_geo(1,iwork)  
           IF (celldm_geo(2,iwork) > ymax) ymax=celldm_geo(2,iwork)  
           IF (celldm_geo(2,iwork) < ymin) ymin=celldm_geo(2,iwork)  
        END DO
        IF (ene_levels(1)==-1000._DP) THEN
           emin=1.d10
           emax=-1.d10
           DO iwork = 1, tot_n
              IF ( energy_geo(iwork) + pressure * omega_geo(iwork) > emax ) &
                 emax = energy_geo(iwork) + pressure * omega_geo(iwork)
              IF ( energy_geo(iwork) + pressure * omega_geo(iwork) < emin ) &
                 emin = energy_geo(iwork) + pressure * omega_geo(iwork)
           ENDDO
!
!    emax and emin are not used as contours 
!
           deltae = (emax-emin) / (ncontours+1)
           DO icont = 1, ncontours
              ene_levels_int(icont) = emin + icont * deltae
              color_levels(icont) = color(MOD((icont-1)/3,8)+1)
           END DO
        END IF

        WRITE(stdout,'(/,5x,"The plot will have ",i5," levels")') ncontours
        DO icont=1, ncontours
           WRITE(stdout,'(5x,"Level ",i5," Energy= ",f15.8," ", a)') icont, &
                          ene_levels_int(icont), TRIM(color_levels(icont))
        ENDDO
        DO ifile=1,ngeo(3)
           filenameps=TRIM(flpsenergy)//int_to_char(ifile)
           CALL add_pressure(filenameps)
           fileout='energy_files/'//TRIM(flenergy)//int_to_char(ifile)
           CALL add_pressure(fileout)
           filenameps=TRIM(filenameps)//TRIM(flext)
           CALL gnuplot_start_2dplot(ncontours, nx, ny)
           DO icont=1,ncontours
              CALL gnuplot_set_contour(fileout,ene_levels_int(icont), &
                                            color_levels(icont),tablefile)
           ENDDO
           CALL gnuplot_close_2dplot_prep()
           xlabel='a (a.u.)'
           ylabel='b/a '
           CALL gnuplot_do_2dplot(filenameps, xmin, xmax, ymin, ymax, xlabel, &
                                                  ylabel, tablefile, flext)
        ENDDO
     CASE DEFAULT
        RETURN
  END SELECT

  CALL gnuplot_end()

  IF (lgnuplot.AND.ionode) &
     ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!  IF (lgnuplot.AND.ionode) &
!     CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

  RETURN
  END SUBROUTINE plot_multi_energy_t
