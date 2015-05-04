!
! Copyright (C) 2014 Andrea Dal Corso
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
  USE input_parameters,     ONLY : ibrav
  USE control_gnuplot,      ONLY : flgnuplot, gnuplot_command, lgnuplot
  USE data_files,           ONLY : flenergy 
  USE postscript_files,     ONLY : flpsenergy
  USE control_energy_plot,  ONLY : ncontours, ene_levels, color_levels
  USE control_mur,          ONLY : lmurn, celldm0
  USE control_quadratic_energy, ONLY : x_pos_min, hessian_v, degree
  USE mp_images,            ONLY : my_image_id, root_image
  USE io_global,            ONLY : ionode
  USE gnuplot,              ONLY : gnuplot_start, gnuplot_end,             &
                                   gnuplot_start_2dplot,                   &
                                   gnuplot_set_contour, gnuplot_do_2dplot, &
                                   gnuplot_set_xticks, gnuplot_xlabel,     &
                                   gnuplot_ylabel, gnuplot_close_2dplot_prep, &
                                   gnuplot_line_v

  IMPLICIT NONE
  CHARACTER(LEN=256) :: gnu_filename, fileout
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=12) :: color(8), xlabel, ylabel
  REAL(DP) :: emax, emin, deltae
  REAL(DP) :: xmin, xmax, ymin, ymax
  INTEGER :: max_contours, nx, ny, ifiles, icont, tot_n, iwork, ierr
  INTEGER :: system
  CHARACTER(LEN=8) :: float_to_char

  IF ( my_image_id /= root_image ) RETURN
  IF (ncontours==0) RETURN
  nx=ngeo(1)
  IF (lmurn) THEN
     ny=1
  ELSE
     ny=ngeo(3)
  ENDIF

  SELECT CASE (ibrav) 
     CASE(4,6,7)
        IF (lmurn) THEN
           tot_n=ngeo(1)
        ELSE
           tot_n=ngeo(1)*ngeo(3)
        ENDIF
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
        ifiles=1
     CASE DEFAULT
        RETURN
  END SELECT

  color(1)='color_red'
  color(2)='color_green'
  color(3)='color_blue'
  color(4)='color_yellow'
  color(5)='color_pink'
  color(6)='color_cyan'
  color(7)='color_orange'
  color(8)='color_black'

  IF (ene_levels(1)==-1000._DP) THEN
     emin=1.d10
     emax=-1.d10
     DO iwork = 1, tot_n
        IF ( energy_geo(iwork) > emax ) emax = energy_geo(iwork) 
        IF ( energy_geo(iwork) < emin ) emin = energy_geo(iwork)
     ENDDO
!
!    emax and emin are not used as contours 
!
     deltae = (emax-emin) / (ncontours+1)
     DO icont = 1, ncontours
        ene_levels(icont) = emin + icont * deltae
        color_levels(icont) = color(MOD((icont-1)/3,8)+1)
     END DO
  END IF

  gnu_filename=TRIM(flgnuplot)//'_energy'
  fileout=TRIM(flenergy)//int_to_char(ifiles)
  CALL gnuplot_start(gnu_filename)
  CALL gnuplot_start_2dplot(ncontours, nx, ny)
  
  DO icont=1,ncontours
     CALL gnuplot_set_contour(fileout,ene_levels(icont),color_levels(icont))
  ENDDO
  CALL gnuplot_close_2dplot_prep()

  SELECT CASE (ibrav) 
     CASE(4,6,7)
        xlabel='a (a.u.)'
        ylabel='c/a '
     CASE DEFAULT
        RETURN
  END SELECT

  fileout=flpsenergy
  CALL gnuplot_do_2dplot(fileout, xmin, xmax, ymin, ymax, xlabel, ylabel)

  IF (degree==2) THEN
     CALL gnuplot_line_v(hessian_v(1,1), hessian_v(2,1), x_pos_min(1),  &
                                       x_pos_min(2),.FALSE.,'color_blue')
     CALL gnuplot_line_v(hessian_v(1,2), hessian_v(2,2), x_pos_min(1),  &
                                       x_pos_min(2),.FALSE.,'color_blue')
  ENDIF
  CALL gnuplot_end()

  IF (lgnuplot.AND.ionode) &
     ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

  RETURN
  END SUBROUTINE plot_multi_energy
