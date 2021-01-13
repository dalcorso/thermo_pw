!
! Copyright (C) 2013-2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE asy
!
!  This module provides subroutines to write a script for the asymptote
!  code to plot the Brillouin zone.
!  It offers the following subroutines:
!
!  asy_writepoint   writes the 3d coordinates of a point and gives a
!                   label to the point.
!
!  asy_write2dpoint writes the 2d coordinates of a point and gives a
!                   label to the point.
!
!  asy_putlabel     writes a label in a 3d position and gives it a small
!                   shift with respect to the 3d point. This routine
!                   can also receive a label gG to write the
!                   greek letter gamma. Only the labels of the different
!                   Brillouin zone are supported.
!
!  asy_put2dlabel   writes a label in a 2d position and gives it a small
!                   shift. Also this routine supports the gG label.
!  
!  asy_writesurface writes the commands to plot a 3d surface defined by 
!                   points given before in the script. The input is the
!                   number of points and their numbers in the list.
!
!  asy_write2dsurface writes the commands to plot a 2d surface defined by 
!                   points given before in the script. The input is the
!                   number of points and their numbers in the list.
!
!  asy_openplot     writes the commands to start an asymptote script and
!                   make a 3d plot.
!                   It also opens the file for the script.
!   
!  asy_open2dplot   writes the commands to start an asymptote script and
!                   make a 2d plot.
!                   It also opens the file for the script.
! 
!  asy_closeplot    closes the file with the script.
!
!  asy_plotaxis     plots the k_x, k_y, and k_z axis. It receives the length
!                   of each axis inside the Brillouin zone and plots an
!                   asis that goes outside the BZ of 0.55, 0.45, and 0.55 times
!                   the length inside. The part inside the BZ is dashed.
!                   It puts also the labels k_x, k_y, and k_z.
!
!  asy_2d_plotaxis  plots the k_x and k_y axis. It receives the length
!                   of each axis inside the Brillouin zone and plots an
!                   asis that goes outside the BZ of 0.55, 0.45 times
!                   the length inside. The part inside the BZ is dashed.
!                   It puts also the labels k_x, and k_y.
!
!  asy_join         writes the commands to draw a line between two 3d 
!                   points. The style and color of the line can be specified
!                   in input.      
! 
!  asy_2d_join      writes the commands to draw a line between two 2d 
!                   points. The style and color of the line can be specified
!                   in input.
!
!
USE kinds, ONLY : DP
USE io_global, ONLY : ionode, ionode_id
USE mp_images, ONLY : intra_image_comm
IMPLICIT NONE
PRIVATE
SAVE

INTEGER :: asyu = 56
REAL(DP) :: asy_proj(3) = (/ 5.0_DP, 2.0_DP, 1.0_DP /)
PUBLIC asy_writepoint, asy_putlabel, asy_writesurface, asy_openplot, &
       asy_closeplot, asy_join, asy_plotaxis, asy_proj, asy_open2dplot, &
       asy_write_2d_point, asy_put2dlabel, asy_2d_join, asy_write_2d_surface, &
       asy_2d_plotaxis

CONTAINS
!
!------------------------------------------------------------------
SUBROUTINE asy_writepoint(point,label)
!------------------------------------------------------------------
!
IMPLICIT NONE

CHARACTER(LEN=5) :: label
REAL(DP) :: point(3)

IF (ionode) &
   WRITE(asyu,'("triple ",a5,"=(",f10.6,",",f10.6,",",f10.6,");")') &
          TRIM(label), point

RETURN
END SUBROUTINE asy_writepoint

!------------------------------------------------------------------
SUBROUTINE asy_write_2d_point(point,label)
!------------------------------------------------------------------
!
IMPLICIT NONE

CHARACTER(LEN=5) :: label
REAL(DP) :: point(2)

IF (ionode) &
   WRITE(asyu,'("pair ",a5,"=(",f10.6,",",f10.6,");")') &
          TRIM(label), point

RETURN
END SUBROUTINE asy_write_2d_point

!------------------------------------------------------------------
SUBROUTINE asy_putlabel(label, a, pos)
!------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER, SAVE :: counter=0
CHARACTER(LEN=3), INTENT(IN) :: label, pos
REAL(DP) :: a(3)
INTEGER :: lens
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=5) :: pointl
CHARACTER(LEN=3) :: label_1

counter=counter+1
pointl="PP"//TRIM(int_to_char(counter))
CALL asy_writepoint(a,pointl)

IF (ionode) THEN
   lens=LEN_TRIM(label)
   IF (label=='gG ') THEN
      WRITE(asyu,'("label(scale(1.6)*""$\Gamma$"",",a,",",a,",red);")') &
                 TRIM(pointl), TRIM(pos)
   ELSEIF (label=='gS ') THEN
      WRITE(asyu,'("label(scale(1.6)*""$\Sigma$"",",a,",",a,",red);")') &
                 TRIM(pointl), TRIM(pos)
   ELSEIF (label=='gS0') THEN
      WRITE(asyu,'("label(scale(1.6)*""$\Sigma_0$"",",a,",",a,",red);")') &
                 TRIM(pointl), TRIM(pos)
   ELSEIF (label=='gS1') THEN
      WRITE(asyu,'("label(scale(1.6)*""$\Sigma_1$"",",a,",",a,",red);")') &
                 TRIM(pointl), TRIM(pos)
   ELSEIF (label=='gD0') THEN
      WRITE(asyu,'("label(scale(1.6)*""$\Delta_0$"",",a,",",a,",red);")') &
                 TRIM(pointl), TRIM(pos)
   ELSEIF (label=='gL0') THEN
      WRITE(asyu,'("label(scale(1.6)*""$\Lambda_0$"",",a,",",a,",red);")')  &
                 TRIM(pointl), TRIM(pos)
   ELSEIF (lens>1.AND.TRIM(label(lens:lens))/=' ') THEN
      label_1=label(lens-1:lens-1)//"_"//label(lens:lens)
      WRITE(asyu,'("label(scale(1.6)*""$",a,"$"",",a,",",a,",red);")') &
                         TRIM(label_1), TRIM(pointl), TRIM(pos)
   ELSE
      WRITE(asyu,'("label(scale(1.6)*""",a,""",",a,",",a,",red);")')  &
                          TRIM(label), TRIM(pointl), TRIM(pos)
   ENDIF
ENDIF

RETURN
END SUBROUTINE asy_putlabel

!------------------------------------------------------------------
SUBROUTINE asy_put2dlabel(label, a, pos)
!------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER, SAVE :: counter=0
CHARACTER(LEN=3), INTENT(IN) :: label, pos
REAL(DP) :: a(2)
INTEGER :: lens
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=5) :: pointl
CHARACTER(LEN=3) :: label_1

counter=counter+1
pointl="PF"//TRIM(int_to_char(counter))
CALL asy_write_2d_point(a,pointl)

IF (ionode) THEN
   lens=LEN_TRIM(label)
   IF (label=='gG ') THEN
      WRITE(asyu,'("label(scale(2.2)*""$\Gamma$"",",a,",",a,",blue);")') &
                 TRIM(pointl), TRIM(pos)
   ELSEIF (lens>1.AND.TRIM(label(lens:lens))/=' ') THEN
      label_1=label(lens-1:lens-1)//"_"//label(lens:lens)
      WRITE(asyu,'("label(scale(2.2)*""$",a,"$"",",a,",",a,",blue);")') &
                         TRIM(label_1), TRIM(pointl), TRIM(pos)
   ELSE
      WRITE(asyu,'("label(scale(2.2)*""",a,""",",a,",",a,",blue);")')  &
                          TRIM(label), TRIM(pointl), TRIM(pos)
   ENDIF
ENDIF

RETURN
END SUBROUTINE asy_put2dlabel

!------------------------------------------------------------------
SUBROUTINE asy_writesurface(indeces)
!------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER :: indeces(9)
CHARACTER(LEN=100) :: linesur
CHARACTER(LEN=6) :: int_to_char
INTEGER :: i

linesur="path3 yp="
DO i=1,indeces(1)
   linesur=TRIM(linesur)//"F"//TRIM(int_to_char(indeces(i+1)))//"--"
ENDDO
linesur=TRIM(linesur)//"cycle;"

IF (ionode) THEN
   WRITE(asyu,*)

   WRITE(asyu,"(a)") TRIM(linesur)
   WRITE(asyu,'("surface y=surface(yp);")')
   WRITE(asyu,'("draw(y,blue+opacity(0.7));")')
   WRITE(asyu,'("draw(yp);")')
ENDIF

RETURN
END SUBROUTINE asy_writesurface

!------------------------------------------------------------------
SUBROUTINE asy_write_2d_surface(nvertices)
!------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER :: nvertices
CHARACTER(LEN=100) :: linesur
CHARACTER(LEN=6) :: int_to_char
INTEGER :: i

linesur="path yp="
DO i=1,nvertices
   linesur=TRIM(linesur)//"G"//TRIM(int_to_char(i))//"--"
ENDDO
linesur=TRIM(linesur)//"cycle;"

IF (ionode) THEN
   WRITE(asyu,*)

   WRITE(asyu,"(a)") TRIM(linesur)
   WRITE(asyu,'("filldraw(yp,yellow);")')
ENDIF

RETURN
END SUBROUTINE asy_write_2d_surface

!------------------------------------------------------------------
SUBROUTINE asy_openplot(filename_asy)
!------------------------------------------------------------------
!
USE mp, ONLY : mp_bcast
IMPLICIT NONE
CHARACTER(LEN=*) :: filename_asy
INTEGER :: find_free_unit
INTEGER :: ios

IF (ionode) THEN
   asyu=find_free_unit()
   OPEN(unit=asyu, file=TRIM(filename_asy), status='unknown', &
             form='formatted', err=30, iostat=ios)
ENDIF
CALL mp_bcast(ios, ionode_id, intra_image_comm)
30 CALL errore('asy_openplot','opening asy plotting file',ABS(ios))

IF (ionode) THEN
   WRITE(asyu,'("import settings;")')
   WRITE(asyu,'("defaultpen(1.6);")')

   WRITE(asyu,'("import three;")')
   WRITE(asyu,'("currentprojection=orthographic(",2(f8.3,","),f8.3,");")') &
              asy_proj(1), asy_proj(2), asy_proj(3)
   WRITE(asyu,'("size(15cm);")')
   WRITE(asyu,'("size3(40cm,40cm,40cm);")')
   WRITE(asyu,*)
ENDIF

RETURN
END SUBROUTINE asy_openplot

!------------------------------------------------------------------
SUBROUTINE asy_open2dplot(filename_asy)
!------------------------------------------------------------------
!
USE mp, ONLY : mp_bcast
IMPLICIT NONE
CHARACTER(LEN=*) :: filename_asy
INTEGER :: find_free_unit
INTEGER :: ios

IF (ionode) THEN
   asyu=find_free_unit()
   OPEN(unit=asyu, file=TRIM(filename_asy), status='unknown', &
             form='formatted', err=30, iostat=ios)
ENDIF
CALL mp_bcast(ios, ionode_id, intra_image_comm)
30 CALL errore('asy_openplot','opening asy plotting file',ABS(ios))

IF (ionode) THEN
   WRITE(asyu,'("import settings;")')
   WRITE(asyu,'("defaultpen(1.6);")')

   WRITE(asyu,'("size(15cm);")')
   WRITE(asyu,*)
ENDIF

RETURN
END SUBROUTINE asy_open2dplot

!------------------------------------------------------------------
SUBROUTINE asy_closeplot()
!------------------------------------------------------------------
!
IMPLICIT NONE
IF (ionode) &
   CLOSE(UNIT=asyu, STATUS='KEEP')

RETURN
END SUBROUTINE asy_closeplot

!------------------------------------------------------------------
SUBROUTINE asy_plotaxis(xk)
!------------------------------------------------------------------
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: xk(3)
REAL(DP) :: maxdx

IF (ionode) THEN
   WRITE(asyu,'("triple G=(0.0,0.0,0.0);")')
   WRITE(asyu,'("triple M1=(",f10.6,",0.0,0.0);")') xk(1)
   WRITE(asyu,'("triple M2=(0.0,",f10.6,",0.0);")') xk(2)
   WRITE(asyu,'("triple M3=(0.0,0.0,",f10.6,");")') xk(3)

   maxdx=max(xk(1),xk(2))
   maxdx=max(maxdx,xk(3))
   WRITE(asyu,'("draw(G--M1,dotted);")')
   WRITE(asyu,'("draw(M1--M1+(",f10.6,",0.0,0.0),Arrow3);")') 0.55_DP*maxdx

   WRITE(asyu,'("draw(G--M2,dotted);")')
   WRITE(asyu,'("draw(M2--M2+(0.0,",f10.6,",0.0),Arrow3);")') 0.45_DP*maxdx

   WRITE(asyu,'("draw(G--M3,dotted);")')
   WRITE(asyu,'("draw(M3--M3+(0.0,0.0,",f10.6,"),Arrow3);")') 0.55_DP*maxdx
   WRITE(asyu,*)
   WRITE(asyu,'("label(scale(1.9)*""$k_x$"",M1+(",f10.6,",0.0,0.0),NW);")') &
                                                         0.55_DP*maxdx
   WRITE(asyu,'("label(scale(1.9)*""$k_y$"",M2+(0.0,",f10.6,",0.0),N);")') &
                                                         0.45_DP*maxdx
   WRITE(asyu,'("label(scale(1.9)*""$k_z$"",M3+(0.0,0.02,",f10.6,"),SE);")') &
                                                         0.55_DP*maxdx
   WRITE(asyu,*)
ENDIF

RETURN
END SUBROUTINE asy_plotaxis

!------------------------------------------------------------------
SUBROUTINE asy_2d_plotaxis(xk)
!------------------------------------------------------------------
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: xk(3)
REAL(DP) :: maxdx

IF (ionode) THEN
   WRITE(asyu,'("triple G=(0.0,0.0,0.0);")')
   WRITE(asyu,'("triple M1=(",f10.6,",0.0,0.0);")') xk(1)
   WRITE(asyu,'("triple M2=(0.0,",f10.6,",0.0);")') xk(2)

   maxdx=max(xk(1),xk(2))
   WRITE(asyu,'("draw(G--M1,dotted);")')
   WRITE(asyu,'("draw(M1--M1+(",f10.6,",0.0,0.0),Arrow3);")') 0.55_DP*maxdx

   WRITE(asyu,'("draw(G--M2,dotted);")')
   WRITE(asyu,'("draw(M2--M2+(0.0,",f10.6,",0.0),Arrow3);")') 0.45_DP*maxdx

   WRITE(asyu,*)
   WRITE(asyu,'("label(scale(1.9)*""$k_x$"",M1+(",f10.6,",0.0,0.0),N);")') &
                                                         0.55_DP*maxdx
   WRITE(asyu,'("label(scale(1.9)*""$k_y$"",M2+(0.0,",f10.6,",0.0),NW);")') &
                                                         0.45_DP*maxdx
   WRITE(asyu,*)
ENDIF

RETURN
END SUBROUTINE asy_2d_plotaxis

!------------------------------------------------------------------
SUBROUTINE asy_join(a,b,pen)
!------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER, SAVE :: counter=0
CHARACTER(LEN=*) :: pen
REAL(DP) :: a(3), b(3)
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=8) :: label, label1
CHARACTER(LEN=256) :: string

counter=counter+1
label='MM'//TRIM(int_to_char(counter))
IF (ionode) &
WRITE(asyu,'("triple ",a8,"=(",f10.6,",",f10.6,",",f10.6,");")') TRIM(label), a

counter=counter+1
label1='MM'//TRIM(int_to_char(counter))
IF (ionode) &
WRITE(asyu,'("triple ",a8,"=(",f10.6,",",f10.6,",",f10.6,");")') TRIM(label1), b

string="draw("//TRIM(label)//"--"//TRIM(label1)//","//TRIM(pen)//");"
IF (ionode) &
WRITE(asyu,'(a)') TRIM(string)

RETURN
END SUBROUTINE asy_join

!------------------------------------------------------------------
SUBROUTINE asy_2d_join(a,b,pen)
!------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER, SAVE :: counter=0
CHARACTER(LEN=*) :: pen
REAL(DP) :: a(2), b(2)
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=8) :: label, label1
CHARACTER(LEN=256) :: string

counter=counter+1
label='MF'//TRIM(int_to_char(counter))
IF (ionode) &
   WRITE(asyu,'("pair ",a8,"=(",f10.6,",",f10.6,");")') TRIM(label), a

counter=counter+1
label1='MF'//TRIM(int_to_char(counter))
IF (ionode) &
   WRITE(asyu,'("pair ",a8,"=(",f10.6,",",f10.6,");")') TRIM(label1), b

string="draw("//TRIM(label)//"--"//TRIM(label1)//","//TRIM(pen)//");"
IF (ionode) &
   WRITE(asyu,'(a)') TRIM(string)

RETURN
END SUBROUTINE asy_2d_join

END MODULE asy
