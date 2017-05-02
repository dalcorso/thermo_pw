!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM test_colors
!
!  This is a simple program that tests the gnuplot_color module of the
!  thermo_pw library. It writes a gnuplot script in color.gnu.
!  Depending on the input, the script will plot:
!  1) The base colors
!  2) The dark and light colors
!  3) The red colors
!  4) The orange colors
!  5) The yellow colors
!  6) The green colors
!  7) The blue colors
!  8) The violet colors
!  9) The gray colors
! 10) All the colors available in gnuplot.
!
! To obtain just the gnuplot script and not the postscript file,
! set the variable lgnuplot to .FALSE.
!
USE kinds,           ONLY : DP
USE gnuplot,         ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                            gnuplot_rectangle, gnuplot_write_command, &
                            gnuplot_put_label
USE gnuplot_color,   ONLY : nbase_colors, base_colors, gnuplot_set_base_colors,&
                            gnuplot_set_dark_light_colors, ndark_light_colors, &
                            dark_light_colors, nred_colors, red_colors, &
                            gnuplot_set_reds, norange_colors, orange_colors, &
                            gnuplot_set_oranges, nyellow_colors, &
                            yellow_colors, gnuplot_set_yellows, &
                            ngreen_colors, green_colors, gnuplot_set_greens, &
                            nblue_colors, blue_colors, gnuplot_set_blues, &
                            nviolet_colors, violet_colors, &
                            gnuplot_set_violets, ngray_colors, gray_colors, &
                            gnuplot_set_grays, nall_colors, all_colors, &
                            gnuplot_set_all_colors, convert_color_name

USE mp_global,       ONLY : mp_startup, mp_global_end
USE environment,     ONLY : environment_start, environment_end
USE mp_images,       ONLY : root_image, my_image_id
USE io_global,       ONLY : ionode, stdout

IMPLICIT NONE
INTEGER :: ierr
CHARACTER(LEN=256) :: gnu_filename, filename, psfilename
CHARACTER(LEN=256) :: flgnuplot, gnuplot_command
LOGICAL :: lgnuplot
CHARACTER(LEN=20) :: color
REAL(DP) :: x(4), y(4), xmin, xmax, ymin, ymax, deltax, xm
INTEGER  :: icolor, system
CHARACTER(LEN=9) :: code='test_colors'

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

IF ( my_image_id /= root_image ) GOTO 100

lgnuplot=.TRUE.
gnuplot_command='gnuplot'
flgnuplot='color.gnu'

IF (ionode) THEN
   WRITE(stdout,'(/,5x,"1) Base colors ")')
   WRITE(stdout,'(5x,"2) Dark and light colors ")')
   WRITE(stdout,'(5x,"3) Reds ")')
   WRITE(stdout,'(5x,"4) Oranges ")')
   WRITE(stdout,'(5x,"5) Yellows ")')
   WRITE(stdout,'(5x,"6) Greens ")')
   WRITE(stdout,'(5x,"7) Blues ")')
   WRITE(stdout,'(5x,"8) Violets ")')
   WRITE(stdout,'(5x,"9) Grays ")')
   WRITE(stdout,'(5x,"10) All colors ")')
   WRITE(stdout,'(/,5x,"Your choice? ")')
   READ(5,*) icolor
ENDIF
gnu_filename=TRIM(flgnuplot)
CALL gnuplot_start(gnu_filename)

xmin=0.0_DP
xmax=1.0_DP
ymin=0.0_DP
ymax=1.0_DP
SELECT CASE (icolor)
   CASE(1)
     deltax=(xmax-xmin)/nbase_colors
     psfilename='color_base.ps'
   CASE(2)
     deltax=(xmax-xmin)/ndark_light_colors
     psfilename='color_dark_light.ps'
   CASE(3)
     deltax=(xmax-xmin)/nred_colors
     psfilename='color_reds.ps'
   CASE(4)
     deltax=(xmax-xmin)/norange_colors
     psfilename='color_oranges.ps'
   CASE(5)
     deltax=(xmax-xmin)/nyellow_colors
     psfilename='color_yellows.ps'
   CASE(6)
     deltax=(xmax-xmin)/ngreen_colors
     psfilename='color_greens.ps'
   CASE(7)
     deltax=(xmax-xmin)/nblue_colors
     psfilename='color_blues.ps'
   CASE(8)
     deltax=(xmax-xmin)/nviolet_colors
     psfilename='color_violets.ps'
   CASE(9)
     deltax=(xmax-xmin)/ngray_colors
     psfilename='color_grays.ps'
   CASE(10)
     deltax=(xmax-xmin)/nall_colors
     psfilename='color_all.ps'
   CASE DEFAULT
      WRITE(stdout,'(5x,"Option not available")')
      GOTO 100
END SELECT
filename=TRIM(psfilename)
CALL gnuplot_write_header(filename, xmin, xmax, ymin, ymax, 1.0_DP ) 
CALL gnuplot_write_command('unset border',.FALSE.)
CALL gnuplot_write_command('unset xtics',.FALSE.)
CALL gnuplot_write_command('unset ytics',.FALSE.)

x=0.0_DP
y=0.0_DP
y(1)=ymin
y(3)=0.3_DP
SELECT CASE (icolor)
   CASE(1)
      CALL gnuplot_set_base_colors()
      DO icolor=1, nbase_colors
         x(1) = xmin + (icolor-1) * deltax
         x(3) = x(1) + deltax
         CALL gnuplot_rectangle(x,y,'1.0','color_'//TRIM(base_colors(icolor)))
         CALL convert_color_name(base_colors(icolor), color)
         xm=(x(1)+x(3))*0.5_DP
         CALL gnuplot_put_label(xm,0.32_DP, icolor, TRIM(color), &
                             .FALSE.,'left rotate')
      ENDDO
   CASE(2)
      CALL gnuplot_set_dark_light_colors()
      DO icolor=1, ndark_light_colors
         x(1) = xmin + (icolor-1) * deltax
         x(3) = x(1) + deltax
         CALL gnuplot_rectangle(x,y,'1.0', &
                               'color_'//TRIM(dark_light_colors(icolor)))
         CALL convert_color_name(dark_light_colors(icolor), color)
         xm=(x(1)+x(3))*0.5_DP
         CALL gnuplot_put_label(xm,0.32_DP, icolor, TRIM(color), &
                             .FALSE.,'left rotate')
      ENDDO
   CASE(3)
      CALL gnuplot_set_reds()
      DO icolor=1, nred_colors
         x(1) = xmin + (icolor-1) * deltax
         x(3) = x(1) + deltax
         CALL gnuplot_rectangle(x,y,'1.0','color_'//TRIM(red_colors(icolor)))
         CALL convert_color_name(red_colors(icolor), color)
         xm=(x(1)+x(3))*0.5_DP
         CALL gnuplot_put_label(xm,0.32_DP, icolor, TRIM(color), &
                             .FALSE.,'left rotate')
      ENDDO
   CASE(4)
      CALL gnuplot_set_oranges()
      DO icolor=1, norange_colors
         x(1) = xmin + (icolor-1) * deltax
         x(3) = x(1) + deltax
         CALL gnuplot_rectangle(x,y,'1.0','color_'//TRIM(orange_colors(icolor)))
         CALL convert_color_name(orange_colors(icolor), color)
         xm=(x(1)+x(3))*0.5_DP
         CALL gnuplot_put_label(xm,0.32_DP, icolor, TRIM(color), &
                             .FALSE.,'left rotate')
      ENDDO
   CASE(5)
      CALL gnuplot_set_yellows()
      DO icolor=1, nyellow_colors
         x(1) = xmin + (icolor-1) * deltax
         x(3) = x(1) + deltax
         CALL gnuplot_rectangle(x,y,'1.0','color_'//TRIM(yellow_colors(icolor)))
         CALL convert_color_name(yellow_colors(icolor), color)
         xm=(x(1)+x(3))*0.5_DP
         CALL gnuplot_put_label(xm,0.32_DP, icolor, TRIM(color), &
                             .FALSE.,'left rotate')
      ENDDO
   CASE(6)
      CALL gnuplot_set_greens()
      DO icolor=1, ngreen_colors
         x(1) = xmin + (icolor-1) * deltax
         x(3) = x(1) + deltax
         CALL gnuplot_rectangle(x,y,'1.0','color_'//TRIM(green_colors(icolor)))
         CALL convert_color_name(green_colors(icolor), color)
         xm=(x(1)+x(3))*0.5_DP
         CALL gnuplot_put_label(xm,0.32_DP, icolor, TRIM(color), &
                             .FALSE.,'left rotate')
      ENDDO
   CASE(7)
      CALL gnuplot_set_blues()
      DO icolor=1, nblue_colors
         x(1) = xmin + (icolor-1) * deltax
         x(3) = x(1) + deltax
         CALL gnuplot_rectangle(x,y,'1.0','color_'//TRIM(blue_colors(icolor)))
         xm=(x(1)+x(3))*0.5_DP
         CALL convert_color_name(blue_colors(icolor), color)
         CALL gnuplot_put_label(xm,0.32_DP, icolor, TRIM(color), &
                             .FALSE.,'left rotate')
      ENDDO
   CASE(8)
      CALL gnuplot_set_violets()
      DO icolor=1, nviolet_colors
         x(1) = xmin + (icolor-1) * deltax
         x(3) = x(1) + deltax
         CALL gnuplot_rectangle(x,y,'1.0','color_'//TRIM(violet_colors(icolor)))
         xm=(x(1)+x(3))*0.5_DP
         CALL convert_color_name(violet_colors(icolor), color)
         CALL gnuplot_put_label(xm,0.32_DP, icolor, TRIM(color), &
                             .FALSE.,'left rotate')
      ENDDO
   CASE(9)
      CALL gnuplot_set_grays()
      DO icolor=1, ngray_colors
         x(1) = xmin + (icolor-1) * deltax
         x(3) = x(1) + deltax
         CALL gnuplot_rectangle(x,y,'1.0','color_'//TRIM(gray_colors(icolor)))
         CALL convert_color_name(gray_colors(icolor), color)
         xm=(x(1)+x(3))*0.5_DP
         CALL gnuplot_put_label(xm,0.32_DP,icolor,TRIM(color), &
                             .FALSE.,'left rotate')
      ENDDO
   CASE(10)
      CALL gnuplot_set_all_colors()
      DO icolor=1, nall_colors
         x(1) = xmin + (icolor-1) * deltax
         x(3) = x(1) + deltax
         CALL gnuplot_rectangle(x,y,'1.0','color_'//TRIM(all_colors(icolor)))
         CALL convert_color_name(all_colors(icolor), color)
         xm=(x(1)+x(3))*0.5_DP
         CALL gnuplot_put_label(xm,0.32_DP, icolor, TRIM(color), &
                             .FALSE.,'left rotate')
      ENDDO
END SELECT

CALL gnuplot_write_command('plot x+100',.FALSE.)
CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)
100 CONTINUE

CALL environment_end( code )
!
CALL mp_global_end ()

END PROGRAM test_colors
