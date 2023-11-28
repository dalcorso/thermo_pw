!
! Copyright (C) 2013-2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE gnuplot
!
    USE kinds, ONLY : DP
    USE io_global, ONLY : ionode
    IMPLICIT NONE
    PRIVATE
    SAVE
    INTEGER :: iun_gnuplot       ! unit where this module writes

    INTEGER :: contour_counter=0
    INTEGER :: contour_max
    INTEGER :: counterv=0
    CHARACTER(LEN=20), ALLOCATABLE :: contour_color(:)
    LOGICAL :: lbackspace=.TRUE.

    PUBLIC iun_gnuplot, gnuplot_write_header, &
           gnuplot_write_vertical_line, gnuplot_write_horizontal_line, &
           gnuplot_write_label, gnuplot_write_file_data, gnuplot_start, &
           gnuplot_set_eref, gnuplot_set_fact, gnuplot_set_gfact, &
           gnuplot_xlabel, gnuplot_ylabel, gnuplot_write_label_yl, &
           gnuplot_write_label_yl_bar, &
           gnuplot_unset_xticks, gnuplot_unset_yticks, &
           gnuplot_unset_border, gnuplot_write_horizontal_segment,&
           gnuplot_set_xticks, gnuplot_set_yticks,    &
           gnuplot_write_file_mul_data_minus, &
           gnuplot_write_file_mul_data_times, &
           gnuplot_write_file_mul_data, gnuplot_write_file_mul_point, &
           gnuplot_write_file_mul_line_point, &
           gnuplot_write_file_mul_data_log10, &
           gnuplot_write_file_mul_data_sum,   &
           gnuplot_write_file_mul_data_diff,   &
           gnuplot_write_file_mul_data_div,   &
           gnuplot_write_file_mul_data_linear,   &
           gnuplot_write_file_mul_point_sum, gnuplot_write_command, &
           gnuplot_end, gnuplot_do_2dplot, gnuplot_start_2dplot, &
           gnuplot_set_contour, gnuplot_rectangle, gnuplot_polygon, &
           gnuplot_put_label, gnuplot_put_label_yl, gnuplot_circle, &
           gnuplot_line_v, gnuplot_rectangle_yl,&
           gnuplot_line, gnuplot_print_objects, gnuplot_close_2dplot_prep, &
           determine_backspace

CONTAINS

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_header(filename, xmin, xmax, ymin, ymax, xscale, &
                                                                     flext)
!--------------------------------------------------------------------------
!
!  This routine sets the dimensions of the plot. It recives:
!  filename : the name of the postscript file where the gnuplot script writes
!  xmin, xmax : the minimum and maximum x values of the plot
!  ymin, ymax : the minimum and maximum y values of the plot
!  xscale : an optional multiplicative factor for all the x coordinate,
!           to change the units. 
!
IMPLICIT NONE
CHARACTER(LEN=*) :: filename
CHARACTER(LEN=*) :: flext
REAL(DP), INTENT(IN) :: xmin, xmax, ymin, ymax, xscale

IF (ionode) THEN
   WRITE(iun_gnuplot,'("set encoding iso_8859_15")')
   IF (flext=='.pdf') THEN
      WRITE(iun_gnuplot,'("set terminal pdf enhanced solid &
                        &color font ""AvantGarde-Book"" fontscale 0.5")')
   ELSE
      WRITE(iun_gnuplot,'("set terminal postscript enhanced solid &
                                             &color ""AvantGarde-Book"" 20")')
   ENDIF
   WRITE(iun_gnuplot,'("set output """, a, """")') TRIM(filename) 
   WRITE(iun_gnuplot,*)
   WRITE(iun_gnuplot,'("set key off")')
   WRITE(iun_gnuplot,'("xscale=",f15.6)') xscale
   WRITE(iun_gnuplot,'("xshift=0.0")') 

   IF (xmin /= xmax) THEN
      WRITE(iun_gnuplot,'("xmin=",f15.6)') xmin
      WRITE(iun_gnuplot,'("xmax=",f15.6)') xmax
      WRITE(iun_gnuplot,'("set xrange [xmin*xscale-xshift:xmax*xscale-xshift]")') 
   ENDIF
   IF (ymin /= ymax) THEN
      WRITE(iun_gnuplot,'("ymin=",f15.6)') ymin
      WRITE(iun_gnuplot,'("ymax=",f15.6)') ymax
      WRITE(iun_gnuplot,'("set yrange [ymin:ymax]")') 
   END IF

   WRITE(iun_gnuplot,'("set border lw 2")') 
   WRITE(iun_gnuplot,'("eref=0.0")') 
   WRITE(iun_gnuplot,'("fact=1.0")') 
   WRITE(iun_gnuplot,'("gfact=1.0")') 
   WRITE(iun_gnuplot,'("point_size=1.0")') 

   CALL gnuplot_write_command('color_red="red"',.FALSE.)
   CALL gnuplot_write_command('color_green="green"',.FALSE.)
   CALL gnuplot_write_command('color_blue="blue"',.FALSE.)
   CALL gnuplot_write_command('color_cyan="cyan"',.FALSE.)
   CALL gnuplot_write_command('color_magenta="magenta"',.FALSE.)
   CALL gnuplot_write_command('color_gold="gold"',.FALSE.)
   CALL gnuplot_write_command('color_pink="pink"',.FALSE.)
   CALL gnuplot_write_command('color_black="black"',.FALSE.)
   CALL gnuplot_write_command('color_olive="olive"',.FALSE.)
   CALL gnuplot_write_command('color_brown="brown"',.FALSE.)
   CALL gnuplot_write_command('color_gray="gray"',.FALSE.)
   CALL gnuplot_write_command('color_light_blue="light-blue"',.FALSE.)
   CALL gnuplot_write_command('color_orange="orange"',.FALSE.)
   CALL gnuplot_write_command('color_yellow="yellow"',.FALSE.)

ENDIF

RETURN
END SUBROUTINE gnuplot_write_header

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_vertical_line(xcoord, linewidth, pos, color, comment)
!--------------------------------------------------------------------------
!
!  This subroutine plots a vertical line
!
IMPLICIT NONE
REAL(DP) :: xcoord, ycoord
INTEGER  :: linewidth
CHARACTER(LEN=3) :: label
CHARACTER(LEN=*) :: color, pos
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt='set arrow from", f12.4, "*xscale-xshift,ymin to ", f12.4,&
                 & "*xscale-xshift,ymax nohead ",a," lw ",i3, " lc rgb ",a)'
IF (comment) frt = '# ' // TRIM(frt)
frt='("'//TRIM(frt)

IF (ionode) &
WRITE( iun_gnuplot, frt ) xcoord, xcoord, TRIM(pos), linewidth, TRIM(color)

RETURN
END SUBROUTINE gnuplot_write_vertical_line

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_horizontal_line(ycoord, linewidth, pos, color, &
                                                               comment)
!--------------------------------------------------------------------------
!  
!  This subroutine plots a horizontal line
!
IMPLICIT NONE
REAL(DP) :: ycoord
INTEGER  :: linewidth
CHARACTER(LEN=*) :: color, pos
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt='set arrow from xmin*xscale-xshift,", f12.4, " to xmax*xscale-xshift,", &
             & f12.4," nohead ",a," lw ",i3," lc rgb ",a)'
IF (comment) frt = '# ' // TRIM(frt)
frt='("'//TRIM(frt)
IF (ionode) &
WRITE(iun_gnuplot, frt) ycoord, ycoord, TRIM(pos), linewidth, TRIM(color)

RETURN
END SUBROUTINE gnuplot_write_horizontal_line

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_horizontal_segment(x, ycoord, linewidth, pos, &
                                                            color, comment)
!--------------------------------------------------------------------------
!
!   This routine is like gnuplot_write_horizontal_line but writes the
!   segment between x(1) and x(2)
!
IMPLICIT NONE
REAL(DP) :: x(2)
INTEGER  :: linewidth
CHARACTER(LEN=*) :: ycoord, color, pos
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt='set arrow from ",f12.4,"*xscale-xshift, '//TRIM(ycoord)//  &
      ' to ",f12.4,"*xscale-xshift,'//TRIM(ycoord) //' nohead ",a,&
      &" lw ",i3," lc rgb ",a)'
IF (comment) frt = '# ' // TRIM(frt)
frt='("'//TRIM(frt)
IF (ionode) &
WRITE(iun_gnuplot, frt) x(1), x(2), TRIM(pos), linewidth, TRIM(color)

RETURN
END SUBROUTINE gnuplot_write_horizontal_segment
 
!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_label(xcoord, ycoord, label, comment)
!--------------------------------------------------------------------------
!
!  This subroutine writes a label at the point (xcoord, ycoord)
!
IMPLICIT NONE
REAL(DP) :: xcoord, ycoord
CHARACTER(LEN=3) :: label
CHARACTER(LEN=20) :: ws
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

CALL set_ws_from_label(label, ws)

frt='set label """,a,""" at ", f12.4,"*xscale-xshift,",f12.4," center")'
IF (comment) frt = '# ' // TRIM(frt)
frt='("'//TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt)  TRIM(ws), xcoord, ycoord

RETURN
END SUBROUTINE gnuplot_write_label

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_put_label(xcoord, ycoord, tag, label, comment, where_lab)
!--------------------------------------------------------------------------
!
!  this routine writes only the label, without any transformation
!
IMPLICIT NONE
REAL(DP) :: xcoord, ycoord
INTEGER, INTENT(IN) :: tag
CHARACTER(LEN=*), OPTIONAL :: where_lab
CHARACTER(LEN=*) :: label
CHARACTER(LEN=256) :: frt, whel
INTEGER :: lens
LOGICAL :: comment

lens=LEN_TRIM(label)
IF (PRESENT(where_lab)) THEN
   whel=' '//TRIM(where_lab)
ELSE
   whel=' center'
ENDIF
frt='set label ",i7," """,a,""" at ", f12.4,"*xscale-xshift,",f12.4,a)'
IF (comment) frt = '# ' // TRIM(frt)
frt='("'//TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt)  tag, TRIM(label),  xcoord, ycoord, TRIM(whel)

RETURN
END SUBROUTINE gnuplot_put_label

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_put_label_yl(xcoord, wherey, tag, label, color, &
                                                        comment, where_lab)
!--------------------------------------------------------------------------
!
!  this routine writes only the label, without any transformation.
!  The y coordinate is calculated from the gnuplot script and here
!  the command to calculate this coordinate must be provided
!
IMPLICIT NONE
REAL(DP) :: xcoord
INTEGER, INTENT(IN) :: tag
CHARACTER(LEN=*), OPTIONAL :: where_lab
CHARACTER(LEN=*) :: label, wherey, color
CHARACTER(LEN=256) :: frt, whel
INTEGER :: lens
LOGICAL :: comment

lens=LEN_TRIM(label)
IF (PRESENT(where_lab)) THEN
   whel=' '//TRIM(where_lab)
ELSE
   whel=' center'
ENDIF
frt='set label ",i7," """,a,""" at ", f12.4,"*xscale-xshift,",a,a,&
                                                       &" tc rgb ",a)'
IF (comment) frt = '# ' // TRIM(frt)
frt='("'//TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt)  tag, TRIM(label),  xcoord, TRIM(wherey), &
                   TRIM(whel), TRIM(color)

RETURN
END SUBROUTINE gnuplot_put_label_yl

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_label_yl(xcoord, ylabel, label, comment)
!--------------------------------------------------------------------------
!
!   This subroutine puts a label in a y position calculated by the
!   gnuplot script
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: xcoord
CHARACTER(LEN=*), INTENT(IN) :: ylabel
CHARACTER(LEN=3) :: label
CHARACTER(LEN=20) :: ws
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

CALL set_ws_from_label(label, ws)

frt='set label """,a,""" at ", f12.4,"*xscale-xshift,",a," center")'
IF (comment) frt = '# ' // TRIM(frt)
frt='("'//TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt)  TRIM(ws), xcoord, TRIM(ylabel)

RETURN
END SUBROUTINE gnuplot_write_label_yl

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_label_yl_bar(xcoord, ylabel, label, comment)
!--------------------------------------------------------------------------
!
!   This subroutine puts a label in a y position calculated by the
!   gnuplot script
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: xcoord
CHARACTER(LEN=*), INTENT(IN) :: ylabel
CHARACTER(LEN=3), INTENT(IN) :: label
CHARACTER(LEN=20) :: ws
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

CALL set_ws_from_label(label, ws)

ws="~"//TRIM(ws)//"{0.6-}"

frt='set label """,a,""" at ", f12.4,"*xscale-xshift,",a," center")'
IF (comment) frt = '# ' // TRIM(frt)
frt='("'//TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt)  TRIM(ws), xcoord, TRIM(ylabel)

RETURN
END SUBROUTINE gnuplot_write_label_yl_bar

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_ylabel(label, comment)
!--------------------------------------------------------------------------
!
!   This subroutine sets the y label
!
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: label
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt = 'set ylabel """,a,"""")'
IF (comment) frt = '# ' // TRIM(frt)
frt='("'//TRIM(frt)

IF (ionode.AND. label /=' ') WRITE(iun_gnuplot, frt) TRIM(label)

RETURN
END SUBROUTINE gnuplot_ylabel

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_xlabel(label, comment)
!--------------------------------------------------------------------------
!
!   This subroutine sets the x label
!
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: label
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt='set xlabel """,a,"""")'
IF (comment) frt = '# ' // TRIM(frt)
frt='("'//TRIM(frt)

IF (ionode.AND. label /=' ') WRITE(iun_gnuplot, frt) TRIM(label)

RETURN
END SUBROUTINE gnuplot_xlabel

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_unset_xticks(comment)
!--------------------------------------------------------------------------
!
!   This subroutine hides the xticks on the x axis
!
IMPLICIT NONE
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt = 'unset xtics")'
IF (comment) frt = '# ' // TRIM(frt)
frt='("'//TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt) 

RETURN
END SUBROUTINE gnuplot_unset_xticks

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_unset_border(comment)
!--------------------------------------------------------------------------
!
!   This subroutine removes the border
!
IMPLICIT NONE
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt = 'unset border")'
IF (comment) frt = '# ' // TRIM(frt)
frt='("'//TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt) 

RETURN
END SUBROUTINE gnuplot_unset_border

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_set_xticks(xstart,delta,xend,comment)
!--------------------------------------------------------------------------
!
!   This routines sets the ticks on the x axis
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: xstart, delta, xend
LOGICAL, INTENT(IN) :: comment
CHARACTER(LEN=256) :: frt
WRITE(frt, '("set xtics ",f11.6,"*xscale-xshift,",f11.6,"*xscale,",f11.6,"*xscale-xshift")') xstart, delta, xend 
IF (comment) frt = '# ' // TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, '(a)') TRIM(frt) 

RETURN
END SUBROUTINE gnuplot_set_xticks

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_set_yticks(ystart,delta,yend,comment)
!--------------------------------------------------------------------------
!
!   This routines sets the ticks on the y axis
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: ystart, delta, yend
LOGICAL, INTENT(IN) :: comment
CHARACTER(LEN=256) :: frt

WRITE(frt, '("set ytics ",f10.5,",",f10.5,",",f10.5)') ystart, delta, yend 
IF (comment) frt = '# ' // TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, '(a)') TRIM(frt) 

RETURN
END SUBROUTINE gnuplot_set_yticks

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_unset_yticks(comment)
!--------------------------------------------------------------------------
!
!   This routines hides the ticks on the y axis
!
IMPLICIT NONE
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt = 'unset ytics")'
IF (comment) frt = '# ' // TRIM(frt)
frt='("'//TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt) 

RETURN
END SUBROUTINE gnuplot_unset_yticks

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_set_eref(eref, comment)
!--------------------------------------------------------------------------
!
!    This subroutine sets eref with a value which is subtracted to all
!    y values
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: eref
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt = 'eref=", e24.12)'
IF (comment) frt = '# ' // TRIM(frt)
frt='("'//TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt) eref

RETURN
END SUBROUTINE gnuplot_set_eref

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_set_gfact(gfact, comment)
!--------------------------------------------------------------------------
!
!    This subroutine sets gfact with a value which is multiplied to all
!    y values after subtraction of eref
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: gfact
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt = 'gfact=", e20.8)'
IF (comment) frt = '# ' // TRIM(frt)
frt='("'//TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt) gfact

RETURN
END SUBROUTINE gnuplot_set_gfact

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_set_fact(fact, comment)
!--------------------------------------------------------------------------
!
!    This subroutine sets fact with a value which is multiplied to all
!    y values before subtraction of eref
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: fact
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt = 'fact=", e20.8)'
IF (comment) frt = '#'//TRIM(frt)
frt='("'//TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt) fact

RETURN
END SUBROUTINE gnuplot_set_fact

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_file_data(data_file,lw,color,start,last, comment)
!--------------------------------------------------------------------------
!
!   This subroutine plots the data contained in a file
!   lw   line width
!   color color
!   start, last are logical variables set to .TRUE. if this is the first
!               line of a plot or the last line of a plot
!
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: data_file
CHARACTER(LEN=*), INTENT(IN) :: lw, color
LOGICAL, INTENT(IN) :: start, last

CHARACTER(LEN=256) :: string
LOGICAL :: comment

string=" """//TRIM(data_file)//""" u ($1*xscale-xshift):($2" &
            //"*fact-eref)*gfact w l lw "//TRIM(lw)//" lc rgb "//TRIM(color)

IF (start) string="plot "//TRIM(string)
IF (lbackspace) THEN
   IF (.NOT.last) string=TRIM(string)//", \ "
ELSE
   IF (.NOT.last) string=TRIM(string)//", \\"
ENDIF
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_data

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_file_mul_data_sum(data_file, col1, col2, col3,&
                     color, start, last, comment)
!--------------------------------------------------------------------------
!
!   This subroutine plots the data contained in a file. The plot is
!   col2+col3  versus col1
!
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: data_file
INTEGER, INTENT(IN) :: col1, col2, col3
CHARACTER(LEN=*), INTENT(IN) :: color
LOGICAL, INTENT(IN) :: start, last

CHARACTER(LEN=256) :: string
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: comment

string=" """//TRIM(data_file)//""" u ($"//TRIM(int_to_char(col1))//"*xscale-xshift):(($"// &
              TRIM(int_to_char(col2)) //"+ $"//TRIM(int_to_char(col3)) &
            //")*fact-eref)*gfact w l lw 3 lc rgb "//TRIM(color)

IF (start) string="plot "//TRIM(string)
IF (lbackspace) THEN
   IF (.NOT.last) string=TRIM(string)//", \ "
ELSE
   IF (.NOT.last) string=TRIM(string)//", \\"
ENDIF
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_data_sum

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_file_mul_data_diff(data_file, col1, col2, col3,&
                     color, start, last, comment)
!--------------------------------------------------------------------------
!
!   This subroutine plots the data contained in a file. The plot is
!   col2-col3 versus col1
!
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: data_file
INTEGER, INTENT(IN) :: col1, col2, col3
CHARACTER(LEN=*), INTENT(IN) :: color
LOGICAL, INTENT(IN) :: start, last

CHARACTER(LEN=256) :: string
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: comment

string=" """//TRIM(data_file)//""" u ($"//TRIM(int_to_char(col1))//"*xscale-xshift):(($"// &
              TRIM(int_to_char(col2)) //"- $"//TRIM(int_to_char(col3)) &
            //")*fact-eref)*gfact w l lw 3 lc rgb "//TRIM(color)

IF (start) string="plot "//TRIM(string)
IF (lbackspace) THEN
   IF (.NOT.last) string=TRIM(string)//", \ "
ELSE
   IF (.NOT.last) string=TRIM(string)//", \\"
ENDIF
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_data_diff
!
!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_file_mul_data_linear(data_file, col1, col2, col3,&
                     alpha, beta, color, start, last, comment)
!--------------------------------------------------------------------------
!
!   This subroutine plots the data contained in a file. The plot is
!   alpha*col2+beta*col3 versus col1
!
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: data_file
INTEGER, INTENT(IN) :: col1, col2, col3
REAL(DP), INTENT(IN) :: alpha, beta
CHARACTER(LEN=*), INTENT(IN) :: color
LOGICAL, INTENT(IN) :: start, last

CHARACTER(LEN=256) :: string, alphas, betas
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: comment

WRITE(alphas,'(f20.7)') alpha
WRITE(betas,'(f20.7)') beta
string=" """//TRIM(data_file)//""" u ($"//TRIM(int_to_char(col1))&
            //"*xscale-xshift):(($"//TRIM(int_to_char(col2))//"*"&
            //TRIM(alphas)//" +"//TRIM(betas)//"*$"&
            //TRIM(int_to_char(col3))//")*fact-eref)*gfact w l lw &
            &3 lc rgb "//TRIM(color)

IF (start) string="plot "//TRIM(string)
IF (lbackspace) THEN
   IF (.NOT.last) string=TRIM(string)//", \ "
ELSE
   IF (.NOT.last) string=TRIM(string)//", \\"
ENDIF
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_data_linear
!
!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_file_mul_data_div(data_file, col1, col2, col3,&
                     color, start, last, comment)
!--------------------------------------------------------------------------
!
!   This subroutine plots the data contained in a file. The plot is
!   col2/col3 versus col1
!
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: data_file
INTEGER, INTENT(IN) :: col1, col2, col3
CHARACTER(LEN=*), INTENT(IN) :: color
LOGICAL, INTENT(IN) :: start, last

CHARACTER(LEN=256) :: string
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: comment

string=" """//TRIM(data_file)//""" u ($"//TRIM(int_to_char(col1))//"*xscale-xshift):(($"// &
              TRIM(int_to_char(col2)) //"/ $"//TRIM(int_to_char(col3)) &
            //")*fact-eref)*gfact w l lw 3 lc rgb "//TRIM(color)

IF (start) string="plot "//TRIM(string)
IF (lbackspace) THEN
   IF (.NOT.last) string=TRIM(string)//", \ "
ELSE
   IF (.NOT.last) string=TRIM(string)//", \\"
ENDIF
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_data_div

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_file_mul_point_sum(data_file, col1, col2, col3,&
                     color, start, last, comment)
!--------------------------------------------------------------------------
!
!   This subroutine plots the data contained in a file. The plot is
!   col2+col3  versus col1 and points are used
!
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: data_file
INTEGER, INTENT(IN) :: col1, col2, col3
CHARACTER(LEN=*), INTENT(IN) :: color
LOGICAL, INTENT(IN) :: start, last

CHARACTER(LEN=256) :: string
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: comment

string=" """//TRIM(data_file)//""" u ($"//TRIM(int_to_char(col1))//"*xscale-xshift):(($"// &
              TRIM(int_to_char(col2)) //"+ $"//TRIM(int_to_char(col3)) &
      //")*fact-eref)*gfact w p pt 82 ps point_size lc rgb "//TRIM(color)

IF (start) string="plot "//TRIM(string)
IF (lbackspace) THEN
   IF (.NOT.last) string=TRIM(string)//", \ "
ELSE
   IF (.NOT.last) string=TRIM(string)//", \\"
ENDIF
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_point_sum

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_file_mul_data_times(data_file, col1, col2, col3,&
                     color, start, last, comment)
!--------------------------------------------------------------------------
!
!   This subroutine plots the data contained in a file. The plot is
!   col2*col3  versus col1
!
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: data_file
INTEGER, INTENT(IN) :: col1, col2, col3
CHARACTER(LEN=*), INTENT(IN) :: color
LOGICAL, INTENT(IN) :: start, last

CHARACTER(LEN=256) :: string
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: comment

string=" """//TRIM(data_file)//""" u ($"//TRIM(int_to_char(col1))//"*xscale-xshift):(($"// &
              TRIM(int_to_char(col2)) //"* $"//TRIM(int_to_char(col3)) &
            //")*fact-eref)*gfact w l lw 3 lc rgb "//TRIM(color)

IF (start) string="plot "//TRIM(string)
IF (lbackspace) THEN
   IF (.NOT.last) string=TRIM(string)//", \ "
ELSE
   IF (.NOT.last) string=TRIM(string)//", \\"
ENDIF
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_data_times

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_file_mul_data(data_file, col1, col2, color, start, &
                                       last, comment)
!--------------------------------------------------------------------------
!
!   This subroutine plots the data contained in a file. The plot is
!   col2 versus col1
!
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: data_file
INTEGER, INTENT(IN) :: col1, col2
CHARACTER(LEN=*), INTENT(IN) :: color
LOGICAL, INTENT(IN) :: start, last

CHARACTER(LEN=256) :: string
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: comment

string=" """//TRIM(data_file)//""" u ($"//TRIM(int_to_char(col1))//&
                     "*xscale-xshift):($"//TRIM(int_to_char(col2)) &
            //"*fact-eref)*gfact w l lw 3 lc rgb "//TRIM(color)

IF (start) string="plot "//TRIM(string)
IF (lbackspace) THEN
   IF (.NOT.last) string=TRIM(string)//", \ "
ELSE
   IF (.NOT.last) string=TRIM(string)//", \\"
ENDIF
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_data

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_file_mul_data_log10(data_file, col1, col2, color, &
                                                       start, last, comment)
!--------------------------------------------------------------------------
!
!   This subroutine plots the logarithm with basis 10 of the data contained 
!   in a file. The plot is log10 (col2) versus col1
!
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: data_file
INTEGER, INTENT(IN) :: col1, col2
CHARACTER(LEN=*), INTENT(IN) :: color
LOGICAL, INTENT(IN) :: start, last

CHARACTER(LEN=256) :: string
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: comment

string=" """//TRIM(data_file)//""" u ($"//TRIM(int_to_char(col1))//&
                     "*xscale-xshift):(log10($"//TRIM(int_to_char(col2)) &
            //"*fact-eref)*gfact) w l lw 3 lc rgb "//TRIM(color)

IF (start) string="plot "//TRIM(string)
IF (lbackspace) THEN
   IF (.NOT.last) string=TRIM(string)//", \ "
ELSE
   IF (.NOT.last) string=TRIM(string)//", \\"
ENDIF
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_data_log10

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_file_mul_data_minus(data_file, col1, col2, &
                                      color, start, last, comment)
!--------------------------------------------------------------------------
!
!   This subroutine plots the data contained in a file. The plot is
!   -col2 versus col1
!
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: data_file
INTEGER, INTENT(IN) :: col1, col2
CHARACTER(LEN=*), INTENT(IN) :: color
LOGICAL, INTENT(IN) :: start, last

CHARACTER(LEN=256) :: string
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: comment

string=" """//TRIM(data_file)//""" u ($"//TRIM(int_to_char(col1))//&
                     "*xscale-xshift):(-$"//TRIM(int_to_char(col2)) &
            //"*fact-eref)*gfact w l lw 3 lc rgb "//TRIM(color)

IF (start) string="plot "//TRIM(string)
IF (lbackspace) THEN
   IF (.NOT.last) string=TRIM(string)//", \ "
ELSE
   IF (.NOT.last) string=TRIM(string)//", \\"
ENDIF
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_data_minus

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_file_mul_point(data_file, col1, col2, color, start, &
                                        last, comment)
!--------------------------------------------------------------------------
!
!   This subroutine plots the data contained in a file with points. The plot is
!   col2 versus col1
!
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: data_file
INTEGER, INTENT(IN) :: col1, col2
CHARACTER(LEN=*), INTENT(IN) :: color
LOGICAL, INTENT(IN) :: start, last

CHARACTER(LEN=256) :: string
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: comment

string=" """//TRIM(data_file)//""" u ($"//TRIM(int_to_char(col1))//"*xscale-xshift):($"// &
              TRIM(int_to_char(col2)) &
            //"*fact-eref)*gfact w p pt 82 ps point_size lc rgb "//TRIM(color)

IF (start) string="plot "//TRIM(string)
IF (lbackspace) THEN
   IF (.NOT.last) string=TRIM(string)//", \ "
ELSE
   IF (.NOT.last) string=TRIM(string)//", \\"
ENDIF
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_point

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_file_mul_line_point(data_file, col1, col2, color, &
                                  start, continue, last, comment)
!--------------------------------------------------------------------------
!
!   This subroutine plots the data contained in a file with line and
!   points. The plot is col2 versus col1
!
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: data_file
INTEGER, INTENT(IN) :: col1, col2
CHARACTER(LEN=*), INTENT(IN) :: color
LOGICAL, INTENT(IN) :: start, continue, last

CHARACTER(LEN=256) :: string
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: comment

string=" """//TRIM(data_file)//""" u ($"//TRIM(int_to_char(col1))//"*xscale-xshift):($"// &
              TRIM(int_to_char(col2)) &
            //"*fact-eref)*gfact w lp lw 3 pt 82 ps point_size lc rgb "&
            //TRIM(color)

IF (start) string="plot "//TRIM(string)
IF (continue) string="replot "//TRIM(string)
IF (lbackspace) THEN
   IF (.NOT.last) string=TRIM(string)//", \ "
ELSE
   IF (.NOT.last) string=TRIM(string)//", \\"
ENDIF
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_line_point

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_write_command(command, comment)
!--------------------------------------------------------------------------
!
!   this subroutine writes a command in a gnuplot script
!
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: command
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt = TRIM(command)
IF (comment) frt = '# ' // TRIM(frt)

IF (ionode.AND. command /=' ') WRITE(iun_gnuplot, '(a)') TRIM(frt) 

RETURN
END SUBROUTINE gnuplot_write_command

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_initialize_contour_counter(max_contours)
!--------------------------------------------------------------------------
!
!   This routine initializes the contour counter
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: max_contours

contour_counter=0
contour_max=max_contours
ALLOCATE(contour_color(contour_max))

RETURN
END SUBROUTINE gnuplot_initialize_contour_counter

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_set_contour(file2d_dat,levr,color,tablefile)
!--------------------------------------------------------------------------
!
!  This routine gives the commands to search a contour level in a file 
!  with a 2d function and write it in a table file. It assumes that 
!  gnuplot_start_2dplot has been already called
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: levr
CHARACTER(LEN=*),INTENT(IN) :: file2d_dat, color
CHARACTER(LEN=256) :: filename, tablefile
CHARACTER(LEN=6) :: int_to_char

contour_counter=contour_counter+1
IF (contour_counter > contour_max) CALL errore('gnuplot_set_contour',&
                                                 'too many contours',1)
contour_color(contour_counter) = color
IF (ionode) THEN
   WRITE(iun_gnuplot,'("set cntrparam levels discrete",f12.6)') levr
   filename=TRIM(tablefile)//'_'//TRIM(int_to_char(contour_counter))//'.dat'
   WRITE(iun_gnuplot,'("set output """,a,"""")') TRIM(filename)
   WRITE(iun_gnuplot,'("splot ''",a,"'' using 1:2:3 w l")') TRIM(file2d_dat)
ENDIF

END SUBROUTINE gnuplot_set_contour

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_start_2dplot(max_contours, nx, ny)
!--------------------------------------------------------------------------
!
!  this soubroutine initialize a 2d plot on a file with nx, ny grid of points
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: max_contours
INTEGER, INTENT(IN) :: nx, ny

IF (ionode) THEN
   WRITE(iun_gnuplot, '("nyp=",i5)') ny
   WRITE(iun_gnuplot, '("nxp=",i5)') nx
   WRITE(iun_gnuplot, '("set view map")')
   WRITE(iun_gnuplot, '("set size square")')
   WRITE(iun_gnuplot, '("unset surface")')
   WRITE(iun_gnuplot, '("unset clabel")')
   WRITE(iun_gnuplot, '("set contour")')
   WRITE(iun_gnuplot, '("set dgrid3d nyp,nxp")')
   WRITE(iun_gnuplot, '("set cntrparam cubicspline")')
   WRITE(iun_gnuplot, '("set table")')
ENDIF

CALL gnuplot_initialize_contour_counter(max_contours)

RETURN
END SUBROUTINE gnuplot_start_2dplot

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_line_v( v_x, v_y, x_0, y_0, start, color)
!--------------------------------------------------------------------------
!
!
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP), INTENT(IN) :: v_x, v_y, x_0, y_0
LOGICAL, INTENT(IN) :: start
CHARACTER(LEN=*),INTENT(IN) :: color
CHARACTER(LEN=10) :: avx, avy, ax0, ay0
CHARACTER(LEN=6) :: int_to_char

counterv=counterv+1

avx='v_x_'//TRIM(int_to_char(counterv))
avy='v_y_'//TRIM(int_to_char(counterv))
ax0='x_0_'//TRIM(int_to_char(counterv))
ay0='y_0_'//TRIM(int_to_char(counterv))
IF (ionode) THEN
   WRITE(iun_gnuplot,'(a,"=",f15.8)') avx, v_x
   WRITE(iun_gnuplot,'(a,"=",f15.8)') avy, v_y
   WRITE(iun_gnuplot,'(a,"=",f15.8)') ax0, x_0
   WRITE(iun_gnuplot,'(a,"=",f15.8)') ay0, y_0

   IF (start) THEN
      WRITE(iun_gnuplot,'("plot ",a,"/",a,"*(x-",a,")+",a," lw 2 lc rgb ",a)')&
             avy, avx, ax0, ay0, TRIM(color)
   ELSE
      WRITE(iun_gnuplot,'("replot ",a,"/",a,"*(x-",a,")+",a," lw 2 lc rgb ",a)')&
            avy, avx, ax0, ay0, TRIM(color)
   ENDIF
ENDIF
RETURN
END SUBROUTINE gnuplot_line_v

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_close_2dplot_prep
!--------------------------------------------------------------------------
!
IMPLICIT NONE
IF (ionode) WRITE(iun_gnuplot,'("unset table")')

RETURN
END SUBROUTINE gnuplot_close_2dplot_prep

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_do_2dplot(filename, xmin, xmax, ymin, ymax, xlabel, &
                                       ylabel, tablefile, flext)
!--------------------------------------------------------------------------
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: xmin, xmax, ymin, ymax
CHARACTER(LEN=*), INTENT(IN) :: filename, flext, xlabel, ylabel
CHARACTER(LEN=256) :: filename1, tablefile
CHARACTER(LEN=6) :: int_to_char
INTEGER :: iplot

CALL gnuplot_write_header(filename, xmin, xmax, ymin, ymax, 1.0_DP, flext)
CALL gnuplot_xlabel(xlabel,.FALSE.)
CALL gnuplot_ylabel(ylabel,.FALSE.)

IF (ionode) THEN
   IF ((xmax-xmin) < 2.0_DP*(ymax-ymin)) &
     CALL gnuplot_set_xticks(xmin,(xmax-xmin)/4.0_DP,xmax,.FALSE.)
   IF ((ymax-ymin) < 2.0_DP*(xmax-xmin)) &
     CALL gnuplot_set_yticks(ymin,(ymax-ymin)/4.0_DP,ymax,.FALSE.)
   WRITE(iun_gnuplot,'("set size ratio",f15.7)') (ymax-ymin)/(xmax-xmin)
   DO iplot=1, contour_counter
      filename1=TRIM(tablefile)//'_'//TRIM(int_to_char(iplot))//'.dat'
      IF (iplot==1.AND.contour_counter>1) THEN
         WRITE(iun_gnuplot,'("plot """,a,""" u 1:2 w l lw 4 lc rgb ",a,",\")') & 
                  TRIM(filename1), TRIM(contour_color(iplot))
      ELSEIF (iplot==contour_counter.AND.contour_counter>1) THEN
         WRITE(iun_gnuplot,'("""",a,""" u 1:2 w l lw 4 lc rgb ",a)') & 
                  TRIM(filename1), TRIM(contour_color(iplot))
      ELSEIF (contour_counter==1) THEN
         WRITE(iun_gnuplot,'("plot """,a,""" u 1:2 w l lw 4 lc rgb ",a)') & 
                  TRIM(filename1), TRIM(contour_color(iplot))
      ELSE
         WRITE(iun_gnuplot,'("""",a,""" u 1:2 w l lw 4 lc rgb ",a,",\")') & 
                  TRIM(filename1), TRIM(contour_color(iplot))
      ENDIF
   ENDDO
ENDIF

CALL gnuplot_stop_2dplot()

RETURN
END SUBROUTINE gnuplot_do_2dplot

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_line(x, y, lw, front, color)
!--------------------------------------------------------------------------
!
!  write a line from (x(1), y(1)) to (x(2),y(2))
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: x(2), y(2)
CHARACTER(LEN=*), INTENT(IN) :: color, lw, front

IF (ionode) THEN
   IF (TRIM(front)=='front') THEN
      WRITE(iun_gnuplot, &
        '("set arrow from ",f12.6,"*xscale-xshift,",f12.6," to ",f12.6,"*xscale-xshift,",f12.6,&
            &" nohead front lw ",a," lc rgb ",a)') x(1), y(1), x(2), y(2), &
                              TRIM(lw), TRIM(color)
   ELSE
      WRITE(iun_gnuplot, &
     '("set arrow from ",f12.6,"*xscale-xshift,",f12.6," to ",f12.6,"*xscale-xshift,",f12.6,&
            &" nohead back lw ",a," lc rgb ",a)') x(1), y(1), x(2), y(2), &
                               TRIM(lw), TRIM(color)
   ENDIF
ENDIF
RETURN
END SUBROUTINE gnuplot_line

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_rectangle(x, y, opacity, color)
!--------------------------------------------------------------------------
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: x(4), y(4)
CHARACTER(LEN=*), INTENT(IN) :: opacity, color

IF (ionode) &
   WRITE(iun_gnuplot, &
     '("set obj rect from ",f12.6,"*xscale-xshift,",f12.6," to ",f12.6,"*xscale-xshift,",f12.6,&
                     &" behind fs solid ",a," noborder fc rgb ",a)') x(1), y(1), &
                                     x(3), y(3), TRIM(opacity), TRIM(color)
RETURN
END SUBROUTINE gnuplot_rectangle

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_rectangle_yl(x, opacity, color)
!--------------------------------------------------------------------------
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: x(2)
CHARACTER(LEN=*), INTENT(IN) :: opacity, color

IF (ionode) &
   WRITE(iun_gnuplot, &
     '("set obj rect from ",f12.6,"*xscale-xshift, ymin  to ",f12.6,&
       &"*xscale-xshift, ymax behind fs solid ",a," noborder fc rgb ",a)')&
                x(1), x(2), TRIM(opacity), TRIM(color)
RETURN
END SUBROUTINE gnuplot_rectangle_yl

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_circle(x, y, radius, opacity, color)
!--------------------------------------------------------------------------
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: x, y
CHARACTER(LEN=*), INTENT(IN) :: radius, opacity, color

IF (ionode) &
   WRITE(iun_gnuplot, &
     '("set obj circle at ",f12.6,"*xscale-xshift,",f12.6, " size ", a, &
                     &" front fs solid ",a," noborder fc rgb ",a)') x, y, &
                                        TRIM(radius), TRIM(opacity), TRIM(color)
RETURN
END SUBROUTINE gnuplot_circle

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_polygon(n, x, y, opacity, color)
!--------------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
REAL(DP), INTENT(IN) :: x(n), y(n)
CHARACTER(LEN=*), INTENT(IN) :: opacity, color
CHARACTER(LEN=512) :: fmt_str, aux
INTEGER :: i
fmt_str="set obj polygon from "
DO i=1,n
   WRITE(aux,'(f12.6,"*xscale-xshift,",f12.6," to")') x(i), y(i)
   fmt_str=TRIM(fmt_str) // TRIM(aux)
ENDDO
WRITE(aux,'(f12.6,"*xscale-xshift,",f12.6, " behind fs solid ",a," noborder fc rgb ", a)') &
                                             x(1), y(1), TRIM(opacity), TRIM(color)
fmt_str=TRIM(fmt_str) // TRIM(aux)

IF (ionode)  WRITE(iun_gnuplot, '(a)') TRIM(fmt_str)

RETURN
END SUBROUTINE gnuplot_polygon

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_print_objects()
!--------------------------------------------------------------------------
!
IMPLICIT NONE

IF (ionode)  WRITE(iun_gnuplot, '("plot x+1e6" )') 

RETURN
END SUBROUTINE gnuplot_print_objects

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_stop_2dplot()
!--------------------------------------------------------------------------

IMPLICIT NONE

DEALLOCATE(contour_color)

RETURN
END SUBROUTINE gnuplot_stop_2dplot

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_start(filename_gnu)
!--------------------------------------------------------------------------
!
!  This routine opens the gnuplot script file
!
IMPLICIT NONE
CHARACTER(LEN=*) :: filename_gnu
INTEGER :: find_free_unit

IF (ionode) THEN
   iun_gnuplot=find_free_unit()
   OPEN(UNIT=iun_gnuplot, FILE=TRIM(filename_gnu), STATUS='unknown', &
                                                           FORM='formatted')
ENDIF

RETURN
END SUBROUTINE gnuplot_start

!--------------------------------------------------------------------------
SUBROUTINE gnuplot_end
!--------------------------------------------------------------------------
!
!  This routine closes the gnuplot script file
!
IMPLICIT NONE

IF (ionode) CLOSE(UNIT=iun_gnuplot, STATUS='KEEP') 

RETURN
END SUBROUTINE gnuplot_end

!--------------------------------------------------------------------------
SUBROUTINE set_ws_from_label(label, ws)
!--------------------------------------------------------------------------
!
IMPLICIT NONE
CHARACTER(LEN=3), INTENT(IN) :: label
CHARACTER(LEN=20), INTENT(OUT) :: ws
INTEGER :: lens

lens=LEN_TRIM(label)
IF (label=='gG ') THEN
   ws="{/Symbol G}"
ELSEIF (label=='gS ') THEN
   ws="{/Symbol S}"
ELSEIF (label=='gS0') THEN
   ws="{/Symbol S_0}"
ELSEIF (label=='gS1') THEN
   ws="{/Symbol S_1}"
ELSEIF (label=='gD0') THEN
   ws="{/Symbol D_0}"
ELSEIF (label=='gL0') THEN
   ws="{/Symbol L_0}"
ELSEIF (label=='A1') THEN
   ws="A_1"
ELSEIF (label=='B1') THEN
   ws="B_1"
ELSEIF (label=='L1') THEN
   ws="L_1"
ELSEIF (label=='P1') THEN
   ws="P_1"
ELSEIF (label=='P2') THEN
   ws="P_2"
ELSEIF (label=='X1') THEN
   ws="X_1"
ELSEIF (label=='Y1') THEN
   ws="Y_1"
ELSEIF (label=='Z1') THEN
   ws="Z_1"
ELSEIF (label(1:1)=='g') THEN
    ws="{/Symbol "//label(2:3)//"}"
ELSEIF (lens>1.AND.TRIM(label(lens:lens))/=' ') THEN
   ws=label(lens-1:lens-1)//"_"//label(lens:lens)
ELSE
    ws=label
ENDIF
RETURN
END SUBROUTINE set_ws_from_label

!----------------------------------------------------------------------------
  SUBROUTINE determine_backspace()
!----------------------------------------------------------------------------
  !
  ! This routine has to be called once. It determines if the 
  ! Fortran treats the backspace character as an escape (as in C)
  ! or as a character as is usually the case.
  ! lbackspace is set to .TRUE. in the second case.
  ! This routine must be called by all processors.
  !
  ! If the routine is not called this module sets lbackspace to .TRUE..
  !

  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  USE io_global, ONLY : meta_ionode, meta_ionode_id
  IMPLICIT NONE
  INTEGER :: outunit, filesize
  INTEGER :: find_free_unit

  IF (meta_ionode) THEN 
     outunit=find_free_unit()
     OPEN(UNIT=outunit, FILE='test.txt', STATUS='unknown')
     WRITE(outunit,'(a)') '\\'
     CLOSE(UNIT=outunit,STATUS='KEEP')
     OPEN(UNIT=outunit, FILE='test.txt', STATUS='old')
     INQUIRE(UNIT=outunit, SIZE=filesize)
     CLOSE(UNIT=outunit,STATUS='DELETE')
  ENDIF
  CALL mp_bcast(filesize, meta_ionode_id, world_comm)
  IF (filesize==3) THEN
     lbackspace=.TRUE.
  ELSE
     lbackspace=.FALSE.
  ENDIF
   
  END SUBROUTINE determine_backspace

END MODULE gnuplot
