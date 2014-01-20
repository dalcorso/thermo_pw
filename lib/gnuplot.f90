!
! Copyright (C) 2013-2014 A. Dal Corso 
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

    PUBLIC iun_gnuplot, gnuplot_write_header, &
           gnuplot_write_vertical_line, gnuplot_write_horizontal_line, &
           gnuplot_write_label, gnuplot_write_file_data, gnuplot_start, &
           gnuplot_set_eref, gnuplot_set_fact, gnuplot_xlabel, gnuplot_ylabel, &
           gnuplot_unset_xticks, gnuplot_unset_yticks, &
           gnuplot_write_file_mul_data, gnuplot_write_file_mul_point, &
           gnuplot_end

CONTAINS

SUBROUTINE gnuplot_write_header(filename, xmin, xmax, ymin, ymax)
!
!  This routine sets the dimensions of the plot. It recives:
!  filename : the name of the postscript file where the gnuplot script writes
!  xmin, xmax : the minimum and maximum x values of the plot
!  ymin, ymax : the minimum and maximum y values of the plot
!  
!
IMPLICIT NONE
CHARACTER(LEN=*) :: filename
REAL(DP), INTENT(IN) :: xmin, xmax, ymin, ymax

IF (ionode) THEN
   WRITE(iun_gnuplot,'("set encoding iso_8859_15")')
   WRITE(iun_gnuplot,'("set terminal postscript enhanced solid color ""Helvetica"" 20")')
   WRITE(iun_gnuplot,'("set output """, a, """")') TRIM(filename) 
   WRITE(iun_gnuplot,*)
   WRITE(iun_gnuplot,'("set key off")')

   IF (xmin /= xmax) THEN
      WRITE(iun_gnuplot,'("xmin=",f15.6)') xmin
      WRITE(iun_gnuplot,'("xmax=",f15.6)') xmax
      WRITE(iun_gnuplot,'("set xrange [xmin:xmax]")') 
   ENDIF
   IF (ymin /= ymax) THEN
      WRITE(iun_gnuplot,'("ymin=",f15.6)') ymin
      WRITE(iun_gnuplot,'("ymax=",f15.6)') ymax
      WRITE(iun_gnuplot,'("set yrange [ymin:ymax]")') 
   END IF

   WRITE(iun_gnuplot,'("set border lw 2")') 
   WRITE(iun_gnuplot,'("eref=0.0")') 
   WRITE(iun_gnuplot,'("fact=1.0")') 
ENDIF

RETURN
END SUBROUTINE gnuplot_write_header

SUBROUTINE gnuplot_write_vertical_line(xcoord, linewidth, pos, color)
IMPLICIT NONE
REAL(DP) :: xcoord, ycoord
INTEGER  :: linewidth
CHARACTER(LEN=3) :: label
CHARACTER(LEN=*) :: color, pos

IF (ionode) &
WRITE(iun_gnuplot, '("set arrow from", f12.4, ",ymin to ", f12.4,&
                 & ",ymax nohead ",a," lw ",i3, " lc rgb """,a,"""")') &
                          xcoord, xcoord, TRIM(pos), linewidth, TRIM(color)

RETURN
END SUBROUTINE gnuplot_write_vertical_line

SUBROUTINE gnuplot_write_horizontal_line(ycoord, linewidth, pos, color)
IMPLICIT NONE
REAL(DP) :: ycoord
INTEGER  :: linewidth
CHARACTER(LEN=*) :: color, pos

IF (ionode) &
WRITE(iun_gnuplot,'("set arrow from xmin,", f12.4, " to xmax,", &
                    & f12.4," nohead ",a," lw ",i3," lc rgb """,a,"""")') &
                            ycoord, ycoord, TRIM(pos), linewidth, TRIM(color)

RETURN
END SUBROUTINE gnuplot_write_horizontal_line
 
SUBROUTINE gnuplot_write_label(xcoord, ycoord, label)
IMPLICIT NONE
REAL(DP) :: xcoord, ycoord
CHARACTER(LEN=3) :: label
CHARACTER(LEN=20) :: ws

IF (label(1:1)=='g') THEN
    ws="{/Symbol "//label(2:3)//"}"
ELSE
    ws=label
ENDIF
IF (ionode) &
   WRITE(iun_gnuplot,'("set label """,a,""" at ", f12.4,",",f12.4," center")') &
                     TRIM(ws), xcoord, ycoord

RETURN
END SUBROUTINE gnuplot_write_label

SUBROUTINE gnuplot_ylabel(label)
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: label

IF (ionode.AND. label /=' ') &
   WRITE(iun_gnuplot,'("set ylabel """,a,"""")') TRIM(label)

RETURN
END SUBROUTINE gnuplot_ylabel

SUBROUTINE gnuplot_xlabel(label)
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: label

IF (ionode.AND. label /=' ') &
   WRITE(iun_gnuplot,'("set xlabel """,a,"""")') TRIM(label)

RETURN
END SUBROUTINE gnuplot_xlabel

SUBROUTINE gnuplot_unset_xticks()
IMPLICIT NONE
   IF (ionode) WRITE(iun_gnuplot,'("unset xtics")') 
RETURN
END SUBROUTINE gnuplot_unset_xticks

SUBROUTINE gnuplot_unset_yticks()
IMPLICIT NONE
   WRITE(iun_gnuplot,'("unset ytics")') 
RETURN
END SUBROUTINE gnuplot_unset_yticks


SUBROUTINE gnuplot_set_eref(eref)
IMPLICIT NONE
REAL(DP), INTENT(IN) :: eref

IF (ionode) &
   WRITE(iun_gnuplot,'("eref=", e20.8)') eref

RETURN
END SUBROUTINE gnuplot_set_eref

SUBROUTINE gnuplot_set_fact(fact)
IMPLICIT NONE
REAL(DP), INTENT(IN) :: fact

IF (ionode) &
   WRITE(iun_gnuplot,'("fact=", e20.8)') fact

RETURN
END SUBROUTINE gnuplot_set_fact


SUBROUTINE gnuplot_write_file_data(data_file,color,start,last)
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: color
LOGICAL, INTENT(IN) :: start, last

CHARACTER(LEN=*) :: data_file
CHARACTER(LEN=256) :: string
CHARACTER(LEN=6) :: int_to_char

string=" """//TRIM(data_file)//""" u 1:($2" &
            //"*fact-eref) w l lw 3 lc rgb """//TRIM(color)//""""

IF (start) string="plot "//TRIM(string)
IF (.NOT.last) string=TRIM(string)//", \"

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_data

SUBROUTINE gnuplot_write_file_mul_data(data_file,col1,col2,color,start,last)
IMPLICIT NONE

INTEGER, INTENT(IN) :: col1, col2
CHARACTER(LEN=*), INTENT(IN) :: color
LOGICAL, INTENT(IN) :: start, last

CHARACTER(LEN=*) :: data_file
CHARACTER(LEN=256) :: string
CHARACTER(LEN=6) :: int_to_char

string=" """//TRIM(data_file)//""" u ($"//TRIM(int_to_char(col1))//"):($"// &
              TRIM(int_to_char(col2)) &
            //"*fact-eref) w l lw 3 lc rgb """//TRIM(color)//""""

IF (start) string="plot "//TRIM(string)
IF (.NOT.last) string=TRIM(string)//", \"

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_data

SUBROUTINE gnuplot_write_file_mul_point(data_file,col1,col2,color,start,last)
IMPLICIT NONE

INTEGER, INTENT(IN) :: col1, col2
CHARACTER(LEN=*), INTENT(IN) :: color
LOGICAL, INTENT(IN) :: start, last

CHARACTER(LEN=*) :: data_file
CHARACTER(LEN=256) :: string
CHARACTER(LEN=6) :: int_to_char

string=" """//TRIM(data_file)//""" u ($"//TRIM(int_to_char(col1))//"):($"// &
              TRIM(int_to_char(col2)) &
            //"*fact-eref) w p pt 82 lc rgb """//TRIM(color)//""""

IF (start) string="plot "//TRIM(string)
IF (.NOT.last) string=TRIM(string)//", \"

IF (ionode) &
   WRITE(iun_gnuplot,'(a)') TRIM(string)

RETURN
END SUBROUTINE gnuplot_write_file_mul_point

SUBROUTINE gnuplot_start(filename_gnu)
IMPLICIT NONE
CHARACTER(LEN=*) :: filename_gnu

iun_gnuplot=55

IF (ionode) OPEN(UNIT=iun_gnuplot, FILE=TRIM(filename_gnu), &
                 STATUS='unknown', FORM='formatted')

RETURN
END SUBROUTINE gnuplot_start

SUBROUTINE gnuplot_end
IMPLICIT NONE

IF (ionode) CLOSE(UNIT=iun_gnuplot, STATUS='KEEP') 

RETURN
END SUBROUTINE gnuplot_end

END MODULE gnuplot
