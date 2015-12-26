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
  USE thermo_mod,           ONLY : ngeo, celldm_geo, energy_geo, omega_geo, &
                                   reduced_grid
  USE input_parameters,     ONLY : ibrav
  USE control_gnuplot,      ONLY : flgnuplot, gnuplot_command, lgnuplot
  USE data_files,           ONLY : flenergy, flevdat
  USE postscript_files,     ONLY : flpsenergy
  USE control_energy_plot,  ONLY : ncontours, ene_levels, color_levels
  USE control_quadratic_energy, ONLY : x_pos_min, hessian_v, degree, show_fit
  USE control_mur,          ONLY : nvol
  USE control_pressure,     ONLY : pressure, pressure_kb
  USE mp_images,            ONLY : my_image_id, root_image
  USE io_global,            ONLY : ionode
  USE gnuplot,              ONLY : gnuplot_start, gnuplot_end,             &
                                   gnuplot_start_2dplot,                   &
                                   gnuplot_set_contour, gnuplot_do_2dplot, &
                                   gnuplot_set_xticks, gnuplot_xlabel,     &
                                   gnuplot_write_header,                   &
                                   gnuplot_write_file_mul_data,            &
                                   gnuplot_write_file_mul_point,           &
                                   gnuplot_write_horizontal_line,          &
                                   gnuplot_ylabel, gnuplot_close_2dplot_prep, &
                                   gnuplot_line_v

  IMPLICIT NONE
  CHARACTER(LEN=256) :: gnu_filename, filename, filename1, label, filenameps
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=8) :: float_to_char
  CHARACTER(LEN=12) :: color(8), xlabel, ylabel
  REAL(DP) :: emax, emin, deltae, ene_levels_int(ncontours)
  REAL(DP) :: xmin, xmax, ymin, ymax
  INTEGER :: max_contours, nx, ny, icont, tot_n, iwork, ierr
  INTEGER :: system
  INTEGER :: compute_nwork

  IF ( my_image_id /= root_image ) RETURN

  gnu_filename=TRIM(flgnuplot)//'_energy'
  IF (reduced_grid.OR.show_fit) THEN
     filename=TRIM(flevdat)//'_quadratic'
  ELSE
     filename=TRIM(flenergy)//int_to_char(1)
  ENDIF
  filenameps=TRIM(flpsenergy)
  IF (pressure /= 0.0_DP) THEN
     gnu_filename=TRIM(gnu_filename)//'.'//TRIM(float_to_char(pressure_kb,1))
     filename=TRIM(filename)//'.'// TRIM(float_to_char(pressure_kb,1))
     filenameps=TRIM(filenameps)//'.'//TRIM(float_to_char(pressure_kb,1))
  END IF

  color(1)='color_red'
  color(2)='color_green'
  color(3)='color_blue'
  color(4)='color_yellow'
  color(5)='color_pink'
  color(6)='color_cyan'
  color(7)='color_orange'
  color(8)='color_black'
  CALL gnuplot_start(gnu_filename)

  SELECT CASE (ibrav) 
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
                                                                  1.0_DP)
        CALL gnuplot_xlabel('a (a.u.)',.FALSE.)

        IF (pressure /= 0.0_DP) THEN
           label='Enthalpy (Ry)    p= '&
                             &//TRIM(float_to_char(pressure_kb,1))//' kbar'
        ELSE
           label='Energy (Ry)'
        END IF
        CALL gnuplot_ylabel(TRIM(label),.FALSE.)
        filename1=TRIM(flevdat)//'_quadratic'
        IF (pressure /= 0.0_DP) filename1=TRIM(filename1)//'.'// &
                                   TRIM(float_to_char(pressure_kb,1))
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

     CASE(4,6,7)
        IF (ncontours==0) RETURN
        ene_levels_int(:)=ene_levels(:)
        IF (reduced_grid.OR.show_fit) THEN
           nx=nvol
           ny=nvol
        ELSE
           nx=ngeo(1)
           ny=ngeo(3)
        ENDIF
        tot_n=compute_nwork()
        xmin=1.d10
        xmax=-1.d10
        ymin=1.d10
        ymax=-1.d10
        DO iwork = 1, tot_n
           IF (celldm_geo(1,iwork) > xmax) xmax=celldm_geo(1,iwork)  
           IF (celldm_geo(1,iwork) < xmin) xmin=celldm_geo(1,iwork)  
           IF (celldm_geo(3,iwork) > ymax) ymax=celldm_geo(3,iwork)  
           IF (celldm_geo(3,iwork) < ymin) ymin=celldm_geo(3,iwork)  
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

        CALL gnuplot_start_2dplot(ncontours, nx, ny)
  
        DO icont=1,ncontours
           CALL gnuplot_set_contour(filename,ene_levels_int(icont), &
                                            color_levels(icont))
        ENDDO
        CALL gnuplot_close_2dplot_prep()
        xlabel='a (a.u.)'
        ylabel='c/a '

        CALL gnuplot_do_2dplot(filenameps, xmin, xmax, ymin, ymax, xlabel, &
                                                                   ylabel)

        IF (degree==2) THEN
           CALL gnuplot_line_v(hessian_v(1,1), hessian_v(2,1), x_pos_min(1),  &
                                       x_pos_min(2),.FALSE.,'color_blue')
           CALL gnuplot_line_v(hessian_v(1,2), hessian_v(2,2), x_pos_min(1),  &
                                       x_pos_min(2),.FALSE.,'color_blue')
        ENDIF
     CASE (8,9,91,10,11) 
        RETURN
     CASE DEFAULT
        RETURN
  END SELECT

  CALL gnuplot_end()

  IF (lgnuplot.AND.ionode) &
     ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

  RETURN
  END SUBROUTINE plot_multi_energy
