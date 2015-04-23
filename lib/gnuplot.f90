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

    PUBLIC iun_gnuplot, gnuplot_write_header, &
           gnuplot_write_vertical_line, gnuplot_write_horizontal_line, &
           gnuplot_write_label, gnuplot_write_file_data, gnuplot_start, &
           gnuplot_set_eref, gnuplot_set_fact, gnuplot_set_gfact, &
           gnuplot_xlabel, gnuplot_ylabel, gnuplot_write_label_yl, &
           gnuplot_unset_xticks, gnuplot_unset_yticks, &
           gnuplot_set_xticks, gnuplot_set_yticks,    &
           gnuplot_write_file_mul_data, gnuplot_write_file_mul_point, &
           gnuplot_write_file_mul_data_sum, gnuplot_write_command, &
           gnuplot_end, gnuplot_do_2dplot, gnuplot_start_2dplot, &
           gnuplot_set_contour, gnuplot_rectangle, gnuplot_polygon, &
           gnuplot_put_label, gnuplot_circle, gnuplot_line_v, &
           gnuplot_line, gnuplot_print_objects, gnuplot_close_2dplot_prep

CONTAINS

SUBROUTINE gnuplot_write_header(filename, xmin, xmax, ymin, ymax, xscale)
!
!  This routine sets the dimensions of the plot. It recives:
!  filename : the name of the postscript file where the gnuplot script writes
!  xmin, xmax : the minimum and maximum x values of the plot
!  ymin, ymax : the minimum and maximum y values of the plot
!  xscale : an optional multiplicative factor for all the x coordinate,
!           to change the units. It is a comment, must be activated in the
!           script 
!
IMPLICIT NONE
CHARACTER(LEN=*) :: filename
REAL(DP), INTENT(IN) :: xmin, xmax, ymin, ymax, xscale

IF (ionode) THEN
   WRITE(iun_gnuplot,'("set encoding iso_8859_15")')
   WRITE(iun_gnuplot,'("set terminal postscript enhanced solid color ""AvantGarde-Book"" 20")')
   WRITE(iun_gnuplot,'("set output """, a, """")') TRIM(filename) 
   WRITE(iun_gnuplot,*)
   WRITE(iun_gnuplot,'("set key off")')
   WRITE(iun_gnuplot,'("#xscale=",f15.6)') xscale
   WRITE(iun_gnuplot,'("xscale=1.0")') 
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

ENDIF

RETURN
END SUBROUTINE gnuplot_write_header

SUBROUTINE gnuplot_write_vertical_line(xcoord, linewidth, pos, color, comment)
IMPLICIT NONE
REAL(DP) :: xcoord, ycoord
INTEGER  :: linewidth
CHARACTER(LEN=3) :: label
CHARACTER(LEN=*) :: color, pos
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt='("set arrow from", f12.4, "*xscale-xshift,ymin to ", f12.4,&
                 & "*xscale-xshift,ymax nohead ",a," lw ",i3, " lc rgb ",a)'
IF (comment) frt = '# ' // TRIM(frt)
IF (ionode) &
WRITE( iun_gnuplot, frt ) xcoord, xcoord, TRIM(pos), linewidth, TRIM(color)

RETURN
END SUBROUTINE gnuplot_write_vertical_line

SUBROUTINE gnuplot_write_horizontal_line(ycoord, linewidth, pos, color, comment)
IMPLICIT NONE
REAL(DP) :: ycoord
INTEGER  :: linewidth
CHARACTER(LEN=*) :: color, pos
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt='("set arrow from xmin*xscale-xshift,", f12.4, " to xmax*xscale-xshift,", &
             & f12.4," nohead ",a," lw ",i3," lc rgb ",a)'
IF (comment) frt = '# ' // TRIM(frt)
IF (ionode) &
WRITE(iun_gnuplot, frt) ycoord, ycoord, TRIM(pos), linewidth, TRIM(color)

RETURN
END SUBROUTINE gnuplot_write_horizontal_line
 
SUBROUTINE gnuplot_write_label(xcoord, ycoord, label, comment)
IMPLICIT NONE
REAL(DP) :: xcoord, ycoord
CHARACTER(LEN=3) :: label
CHARACTER(LEN=20) :: ws
CHARACTER(LEN=256) :: frt
INTEGER :: lens
LOGICAL :: comment

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
ELSEIF (label=='Y1') THEN
   ws="Y_1"
ELSEIF (label=='Z1') THEN
   ws="Z_1"
ELSEIF (lens>1.AND.TRIM(label(lens:lens))/=' ') THEN
   ws=label(lens-1:lens-1)//"_"//label(lens:lens)
ELSE
   ws=label
ENDIF

frt='("set label """,a,""" at ", f12.4,"*xscale-xshift,",f12.4," center")'
IF (comment) frt = '# ' // TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt)  TRIM(ws), xcoord, ycoord

RETURN
END SUBROUTINE gnuplot_write_label

SUBROUTINE gnuplot_put_label(xcoord, ycoord, tag, label, comment, where_lab)
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
frt='("set label ",i7," """,a,""" at ", f12.4,"*xscale-xshift,",f12.4,a)'
IF (comment) frt = '# ' // TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt)  tag, TRIM(label),  xcoord, ycoord, TRIM(whel)

RETURN
END SUBROUTINE gnuplot_put_label

SUBROUTINE gnuplot_write_label_yl(xcoord, ylabel, label, comment)
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
ELSE
    ws=label
ENDIF

frt='("set label """,a,""" at ", f12.4,",",a," center")'
IF (comment) frt = '# ' // TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt)  TRIM(ws), xcoord, TRIM(ylabel)

RETURN
END SUBROUTINE gnuplot_write_label_yl

SUBROUTINE gnuplot_ylabel(label, comment)
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: label
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt = '("set ylabel """,a,"""")'
IF (comment) frt = '# ' // TRIM(frt)

IF (ionode.AND. label /=' ') WRITE(iun_gnuplot, frt) TRIM(label)

RETURN
END SUBROUTINE gnuplot_ylabel


SUBROUTINE gnuplot_xlabel(label, comment)
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: label
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt='set xlabel """,a,"""")'
IF (comment) THEN
   frt = '("# ' // TRIM(frt)
ELSE
   frt = '(" ' // TRIM(frt)
ENDIF

IF (ionode.AND. label /=' ') WRITE(iun_gnuplot, frt) TRIM(label)

RETURN
END SUBROUTINE gnuplot_xlabel

SUBROUTINE gnuplot_unset_xticks(comment)
IMPLICIT NONE
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt = '("unset xtics")'
IF (comment) frt = '# ' // TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt) 

RETURN
END SUBROUTINE gnuplot_unset_xticks

SUBROUTINE gnuplot_set_xticks(xstart,delta,xend,comment)
IMPLICIT NONE
REAL(DP), INTENT(IN) :: xstart, delta, xend
LOGICAL, INTENT(IN) :: comment
CHARACTER(LEN=256) :: frt
WRITE(frt, '("set xtics ",f11.6,"*xscale-xshift,",f11.6,"*xscale,",f11.6,"*xscale-xshift")') xstart, delta, xend 
IF (comment) frt = '# ' // TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, '(a)') TRIM(frt) 

RETURN
END SUBROUTINE gnuplot_set_xticks

SUBROUTINE gnuplot_set_yticks(ystart,delta,yend,comment)
IMPLICIT NONE
REAL(DP), INTENT(IN) :: ystart, delta, yend
LOGICAL, INTENT(IN) :: comment
CHARACTER(LEN=256) :: frt

WRITE(frt, '("set ytics ",f10.5,",",f10.5,",",f10.5)') ystart, delta, yend 
IF (comment) frt = '# ' // TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, '(a)') TRIM(frt) 

RETURN
END SUBROUTINE gnuplot_set_yticks

SUBROUTINE gnuplot_unset_yticks(comment)
IMPLICIT NONE
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt = '("unset ytics")'
IF (comment) frt = '# ' // TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt) 

RETURN
END SUBROUTINE gnuplot_unset_yticks

SUBROUTINE gnuplot_set_eref(eref, comment)
IMPLICIT NONE
REAL(DP), INTENT(IN) :: eref
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt = '("eref=", e20.8)'
IF (comment) frt = '# ' // TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt) eref

RETURN
END SUBROUTINE gnuplot_set_eref

SUBROUTINE gnuplot_set_gfact(gfact, comment)
IMPLICIT NONE
REAL(DP), INTENT(IN) :: gfact
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt = '("gfact=", e20.8)'
IF (comment) frt = '# ' // TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt) gfact

RETURN
END SUBROUTINE gnuplot_set_gfact

SUBROUTINE gnuplot_set_fact(fact, comment)
IMPLICIT NONE
REAL(DP), INTENT(IN) :: fact
CHARACTER(LEN=256) :: frt
LOGICAL :: comment

frt = '("fact=", e20.8)'
IF (comment) frt = '# ' // TRIM(frt)

IF (ionode) WRITE(iun_gnuplot, frt) fact

RETURN
END SUBROUTINE gnuplot_set_fact

SUBROUTINE gnuplot_write_file_data(data_file,lw,color,start,last, comment)
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: data_file
CHARACTER(LEN=*), INTENT(IN) :: lw, color
LOGICAL, INTENT(IN) :: start, last

CHARACTER(LEN=256) :: string
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: comment

string=" """//TRIM(data_file)//""" u ($1*xscale-xshift):($2" &
            //"*fact-eref)*gfact w l lw "//TRIM(lw)//" lc rgb "//TRIM(color)

IF (start) string="plot "//TRIM(string)
IF (.NOT.last) string=TRIM(string)//", \"
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_data

SUBROUTINE gnuplot_write_file_mul_data_sum(data_file, col1, col2, col3,&
                     color, start, last, comment)
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
IF (.NOT.last) string=TRIM(string)//", \"
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_data_sum

SUBROUTINE gnuplot_write_file_mul_data(data_file, col1, col2, color, start, &
                                       last, comment)
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
IF (.NOT.last) string=TRIM(string)//", \"
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_data

SUBROUTINE gnuplot_write_file_mul_point(data_file, col1, col2, color, start, &
                                        last, comment)
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
            //"*fact-eref)*gfact w p pt 82 lc rgb "//TRIM(color)

IF (start) string="plot "//TRIM(string)
IF (.NOT.last) string=TRIM(string)//", \"
IF (comment) string = '# ' // TRIM(string)

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_point

SUBROUTINE gnuplot_write_command(command, comment)
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

SUBROUTINE gnuplot_initialize_contour_counter(max_contours)

IMPLICIT NONE
INTEGER, INTENT(IN) :: max_contours

contour_counter=0
contour_max=max_contours
ALLOCATE(contour_color(contour_max))

RETURN
END SUBROUTINE gnuplot_initialize_contour_counter

SUBROUTINE gnuplot_set_contour(file2d_dat,levr,color)
!
!  This routine gives the commands to search a contour levr in a file with a 2d
!  function and write it in a table file. It assumes that gnuplot_start_2dplot
!  has been already called
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: levr
CHARACTER(LEN=*),INTENT(IN) :: file2d_dat, color
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=6) :: int_to_char

contour_counter=contour_counter+1
IF (contour_counter > contour_max) CALL errore('gnuplot_set_contour',&
                                                 'too many contours',1)
contour_color(contour_counter) = color
IF (ionode) THEN
   WRITE(iun_gnuplot,'("set cntrparam levels discrete",f12.6)') levr
   filename='table_'//TRIM(int_to_char(contour_counter))//'.dat'
   WRITE(iun_gnuplot,'("set output """,a,"""")') TRIM(filename)
   WRITE(iun_gnuplot,'("splot ''",a,"'' using 1:2:3 w l")') TRIM(file2d_dat)
ENDIF

END SUBROUTINE gnuplot_set_contour

SUBROUTINE gnuplot_start_2dplot(max_contours, nx, ny)
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

SUBROUTINE gnuplot_line_v( v_x, v_y, x_0, y_0, start, color)
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

SUBROUTINE gnuplot_close_2dplot_prep

IMPLICIT NONE
IF (ionode) WRITE(iun_gnuplot,'("unset table")')

RETURN
END SUBROUTINE gnuplot_close_2dplot_prep

SUBROUTINE gnuplot_do_2dplot(filename, xmin, xmax, ymin, ymax, xlabel, ylabel)

IMPLICIT NONE
REAL(DP), INTENT(IN) :: xmin, xmax, ymin, ymax
CHARACTER(LEN=*), INTENT(IN) :: filename, xlabel, ylabel
CHARACTER(LEN=256) :: filename1
CHARACTER(LEN=6) :: int_to_char
INTEGER :: iplot

CALL gnuplot_write_header(filename, xmin, xmax, ymin, ymax, 1.0_DP)
CALL gnuplot_xlabel(xlabel,.FALSE.)
CALL gnuplot_ylabel(ylabel,.FALSE.)

IF (ionode) THEN
   IF ((xmax-xmin) < 2.0_DP*(ymax-ymin)) &
     CALL gnuplot_set_xticks(xmin,(xmax-xmin)/4.0_DP,xmax,.FALSE.)
   IF ((ymax-ymin) < 2.0_DP*(xmax-xmin)) &
     CALL gnuplot_set_yticks(ymin,(ymax-ymin)/4.0_DP,ymax,.FALSE.)
   WRITE(iun_gnuplot,'("set size ratio",f15.7)') (ymax-ymin)/(xmax-xmin)
   DO iplot=1, contour_counter
      filename1='table_'//TRIM(int_to_char(iplot))//'.dat'
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

SUBROUTINE gnuplot_line(x, y, lw, front, color)
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

SUBROUTINE gnuplot_rectangle(x, y, opacity, color)
IMPLICIT NONE
REAL(DP), INTENT(IN) :: x(4), y(4)
CHARACTER(LEN=*), INTENT(IN) :: opacity, color

IF (ionode) &
   WRITE(iun_gnuplot, &
     '("set obj rect from ",f12.6,"*xscale-xshift,",f12.6," to ",f12.6,"*xscale-shift,",f12.6,&
                     &" behind fs solid ",a," noborder fc rgb ",a)') x(1), y(1), &
                                     x(3), y(3), TRIM(opacity), TRIM(color)
RETURN
END SUBROUTINE gnuplot_rectangle

SUBROUTINE gnuplot_circle(x, y, radius, opacity, color)
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


SUBROUTINE gnuplot_polygon(n, x, y, opacity, color)
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

SUBROUTINE gnuplot_print_objects()

IMPLICIT NONE

IF (ionode)  WRITE(iun_gnuplot, '("plot x+1e6" )') 

RETURN
END SUBROUTINE gnuplot_print_objects

SUBROUTINE gnuplot_stop_2dplot()

IMPLICIT NONE

DEALLOCATE(contour_color)

RETURN
END SUBROUTINE gnuplot_stop_2dplot


SUBROUTINE gnuplot_start(filename_gnu)
!
!  This routine opens the gnuplot script file
!
IMPLICIT NONE
CHARACTER(LEN=*) :: filename_gnu

iun_gnuplot=55

IF (ionode) OPEN(UNIT=iun_gnuplot, FILE=TRIM(filename_gnu), &
                 STATUS='unknown', FORM='formatted')

RETURN
END SUBROUTINE gnuplot_start

SUBROUTINE gnuplot_end
!
!  This routine closes the gnuplot script file
!
IMPLICIT NONE

IF (ionode) CLOSE(UNIT=iun_gnuplot, STATUS='KEEP') 

RETURN
END SUBROUTINE gnuplot_end

END MODULE gnuplot
