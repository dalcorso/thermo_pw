!
! Copyright (C) 2020 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE freecad
!
!   This module contains the routines to write a freecad script to plot
!   the Brillouin zone (3d and 2d). The script is optimized for freecad 0.18.
!   The freecad script is read opening freecad and choosing 
!   macro/macros/execute. 
!   In the case of the 3d Brillouin zone, a standard view can be obtained 
!   with the commands
!   view/freeze display/load view (reading the file brillouin_zone.cam)
!   view/freeze display/restore view 1.
!
!   The module contains the following routines:
!
!   freecad_openplot    opens the file with the scripts and writes on it
!                       a few commands to start the script
!   freecad_closeplot   closes the file with the script
!   freecad_centerview  gives a command to put brillouin zone in the center
!                       of the plot
!   freecad_isoview     gives the commands to choose the standard isometric 
!                       view of freecad
!   freecad_plotaxis    plot the axis k_x, k_y, k_z and puts a label close
!                       to each axis
!   freecad_2d_plotaxis plot the axis k_x and k_y and puts a label close to
!                       each axis.
!   freecad_writepoint  writes the coordinates of a point on the script.
!                       The point is not visualized.
!   freecad_join        put an edge between two 3d points. The edge become
!                       visible. It controls also the color of the edge.
!   freecad_setcolor    receives an object and a color in the form of a 
!                       rgb array and sets the color of the object
!   freecad_putlabel    put a label in a given position in 3d or 2d
!                       The routine can transform subscript and superscript
!                       and deal with greek letters as specified in pw input
!   freecad_writesurface gives the commands to create a surface from a
!                       set of points that must be given previously.
!                       no surface is plotted on the screen
!   freecad_createsolid receives the number of surfaces and creates 
!                       a solid inside all the surfaces. The solid is
!                       plotted on the screen
!   freecad_createshell this is used for the 2d brillouin zone. The routine
!                       makes the shell visible on the screen. Only one
!                       shell must have been defined on the script.
!   freecad_setfcfact   set the factor that converts from units 2\pi/a
!                       of the Brillouin zone to mm as in freecad 
!                       conventions
!   freecad_setfontsize set the size of the font of freecad
!                      
!   freecad_createpdf   tries to create a pdf with the TechDraw workbench.
!                       It is usually easier to give the few necessary 
!                       commands inside freecad
!
! The module has also a few routines internal to the module
! 
!   add_group           gives the commands to put an object in a group
!   create_index        receives an integer and creates a string with three
!                       characters with leading zeros
!   convert_label       converts a label of the BZ in the form gG or A1
!                       in a unicode command that prints the greek letter
!                       or the subscript.
!   convert_pos         takes the string N,S,E,W,NE,NW,SE,SW and converts
!                       it in a small 3d shift.
!   convert_2d_pos      takes the string N,S,E,W,NE,NW,SE,SW and converts
!                       it in a small 2d shift in the xy plane.
!   hide_object         set to false the Visibility of an object. Press
!                       a space on the name of the object to show it.
!   put_text            Put a 2d text on a given 3d point. No conversion 
!                       is made of the text.
!   initialize_text     writes the command to put a 2d text, set its size,
!                       and color. It is called after put_text.
!   create_cone         gives the commands to create a cone in a given 
!                       position and with a given orientation.

USE kinds, ONLY : DP
USE io_global, ONLY : ionode, ionode_id, stdout
USE mp_images, ONLY : intra_image_comm
IMPLICIT NONE
PRIVATE
SAVE

INTEGER :: fcu = 56        ! unit of the script file

REAL(DP) :: fcfact=100     ! This factor converts the units of the BZ into
                           ! freecad units
!  freecad uses mm units while the point coordinates are in units 2pi/a. We
!  multiply by fcfact to have a simply movable BZ and a reasonable size on
!  on TechDraw without rescaling. So 1 in 2\pi/a units corresponds to 
!  10 cm in freecad units. This default can be changed from input.
REAL(DP) :: ftsize=35.0_DP,    & ! size of the fonts.
            ftsizesub=30.0_DP    ! size of the fonts subscript

INTEGER :: counter_text=-1       ! count the number of text written to give
                                 ! them unique names
!

PUBLIC freecad_writepoint, freecad_openplot, freecad_closeplot,  &
       freecad_writesurface, freecad_createsolid,  &
       freecad_centerview, freecad_setcolor, freecad_plotaxis,   &
       freecad_join, freecad_putlabel, freecad_setfcfact,        &
       freecad_createshell, freecad_createpdf, freecad_2d_plotaxis, &
       freecad_isoview, freecad_setfontsize

CONTAINS
!
!--------------------------------------------------------------------
SUBROUTINE freecad_openplot(filename_freecad)
!--------------------------------------------------------------------
!
!   This routine opens the freecad script. It must be the first to be
!   called to write a script.
!
USE mp,        ONLY : mp_bcast
IMPLICIT NONE
CHARACTER(LEN=*) :: filename_freecad
INTEGER :: find_free_unit

INTEGER :: ios

IF (ionode) THEN
   fcu=find_free_unit()
   OPEN(UNIT=fcu, FILE=TRIM(filename_freecad)//'.FCMacro', STATUS='unknown', &
             FORM='formatted', ERR=30, IOSTAT=ios)
ENDIF
CALL mp_bcast(ios, ionode_id, intra_image_comm)
30 CALL errore('freecad_openplot','opening freecad plotting file',ABS(ios))

IF (ionode) THEN
   WRITE(fcu,'("import FreeCAD")') 
   WRITE(fcu,'("import Draft")')
   WRITE(fcu,'("import Part")')
   WRITE(fcu,'("import TechDraw")')
   WRITE(fcu,'("doc=App.newDocument(""Brillouin"")")')
!
!   set a global group called Brillouin zone. We put everything in this 
!   group, except the k_x, k_y and k_z axis. It can be plot  
!   with or without the axis using TechDraw.
!
   WRITE(fcu,'("App.activeDocument().Tip = App.activeDocument().&
                &addObject(''App::DocumentObjectGroup'',''Group'')")')
   WRITE(fcu,'("App.activeDocument().Group.Label = ''Brillouin zone''")')
ENDIF

RETURN
END SUBROUTINE freecad_openplot

!----------------------------------------------------------------------------
SUBROUTINE freecad_centerview()
!----------------------------------------------------------------------------
!
!  This routine sends the command to center the Brillouin zone on the screen
!
IMPLICIT NONE

IF (ionode) &
   WRITE(fcu,'("Gui.SendMsgToActiveView(""ViewFit"")")')

RETURN
END SUBROUTINE freecad_centerview

!----------------------------------------------------------------------------
SUBROUTINE freecad_isoview()
!----------------------------------------------------------------------------
!
!  This routine sends the commands to set an isometric view 
!
IMPLICIT NONE

IF (ionode) &
   WRITE(fcu,'("Gui.activeDocument().activeView().viewIsometric()")')

RETURN
END SUBROUTINE freecad_isoview

!-----------------------------------------------------------------------
SUBROUTINE freecad_setcolor(object,rgb,transparency)
!----------------------------------------------------------------------
!
!  The color is given in r,g,b form (from 0.0 to 1.0). The routine sets
!  also the transparency of the solid and the width of the lines.
!  The routine must be called after the object has been created.
!
IMPLICIT NONE
CHARACTER(LEN=*) :: object
REAL(DP) :: rgb(3)
INTEGER  :: transparency

IF (ionode) THEN
   WRITE(fcu,'("FreeCADGui.getDocument(""Brillouin"").&
    &getObject(""",a,""").ShapeColor = (",f4.2,",",f4.2,",",f4.2,")")') &
                            TRIM(object), rgb 
   WRITE(fcu,'("FreeCADGui.ActiveDocument.getObject(""",a,""").&
                        &Transparency=",i3)') TRIM(object), transparency
   WRITE(fcu,'("FreeCADGui.ActiveDocument.getObject(""",a,""").&
                                       &LineWidth = 3.00")') TRIM(object)
ENDIF

RETURN
END SUBROUTINE freecad_setcolor

!--------------------------------------------------------------------------
SUBROUTINE freecad_plotaxis(xk)
!--------------------------------------------------------------------------
!
!  The routine receives the intecepts of the axis with the BZ along 
!  k_x, k_y, and k_z and generates a set of axis from this point 
!  to outside the BZ. These axis are an edge and a small cone.
!  The routine adds also the labels k_x, k_y, and k_z close to the axis.
!
IMPLICIT NONE

REAL(DP) :: xk(3)
REAL(DP) :: maxdx, x1(3), x2(3), r(3), ang, axis(3), rgb(3)
CHARACTER(LEN=256) :: strin
CHARACTER(LEN=3) :: al

maxdx=MAX(xk(1),xk(2))
maxdx=MAX(maxdx,xk(3))
IF (ionode) THEN
!
!   Set a group to collect all the axis, so they can he hidden all together.
!   Note that by default the axis are hidden. To show them click on 
!   axis on the application tree and then press space twice.
!
   WRITE(fcu,'("App.activeDocument().Tip = App.activeDocument().&
                &addObject(''App::DocumentObjectGroup'',''Group001'')")')
   WRITE(fcu,'("App.activeDocument().Group001.Label = ''axis''")')
!
!  Axis z, an edge from the intersection with the bz to 0.55 the
!  maximum size of the intercepts with the BZ
!
   x1=0.0_DP
   rgb=0.0_DP
   x1(3)=xk(3)
   x2=0.0_DP
   x2(3)=xk(3)+0.55_DP*maxdx
   CALL freecad_join(x1,x2,rgb,al)
   CALL add_group('Group001','Edge'//TRIM(al))
   CALL hide_object('Edge'//TRIM(al))
!
!  Axis x
!
   x1=0.0_DP
   rgb=0.0_DP
   x1(1)=xk(1)
   x2=0.0_DP
   x2(1)=xk(1)+0.55_DP*maxdx
   CALL freecad_join(x1,x2,rgb,al)
   CALL add_group('Group001','Edge'//TRIM(al))
   CALL hide_object('Edge'//TRIM(al))
!
!  Axis y
!
   x1=0.0_DP
   rgb=0.0_DP
   x1(2)=xk(2)
   x2=0.0_DP
   x2(2)=xk(2)+0.45_DP*maxdx
   CALL freecad_join(x1,x2,rgb,al)
   CALL add_group('Group001','Edge'//TRIM(al))
   CALL hide_object('Edge'//TRIM(al))
!
!  Tip of the z axis.
!
   r=0.0_DP
   r(3)=(xk(3)+maxdx*0.45_DP)*fcfact
   axis=0.0_DP
   axis(3)=1.0_DP
   ang=0.0_DP
   CALL create_cone('Cone', fcfact/25.0_DP, r, axis, ang)
   CALL add_group('Group001','Cone')
   CALL hide_object('Cone')

!
!  Tip of the x axis.
!
   r=0.0_DP
   r(1)=(xk(1)+maxdx*0.55_DP)*fcfact
   axis=0.0_DP
   axis(2)=1.0_DP
   ang=90.0_DP
   CALL create_cone('Cone001', fcfact/25.0_DP, r, axis, ang)
   CALL add_group('Group001','Cone001')
   CALL hide_object('Cone001')
!
!  Tip of the y axis.
!
   r=0.0_DP
   r(2)=(xk(2)+maxdx*0.45_DP)*fcfact
   axis=0.0_DP
   axis(1)=1.0_DP
   ang=-90.0_DP
   CALL create_cone('Cone002', fcfact/25.0_DP, r, axis, ang)
   CALL add_group('Group001','Cone002')
   CALL hide_object('Cone002')
!
!  The label k of k_x. The text here is red
!
   rgb=0.0_DP
   rgb(1)=1.0_DP

   counter_text=counter_text+1
   CALL create_index(counter_text, al)
   r(1)=(xk(1)+maxdx*0.55_DP)*fcfact
   r(2)=0.10_DP*fcfact*maxdx
   r(3)=-0.12_DP*fcfact*maxdx
   CALL put_text('k',r)
   CALL initialize_text('Text'//TRIM(al), ftsize, rgb)
   CALL add_group('Group001','Text'//TRIM(al))
   CALL hide_object('Text'//TRIM(al))
!
!  The label x of k_x.
!
   counter_text=counter_text+1
   CALL create_index(counter_text, al)
   r(1)=(xk(1)+maxdx*0.55_DP)*fcfact
   r(2)=(0.10_DP+0.13_DP)*fcfact*maxdx
   r(3)=(-0.12_DP-0.06_DP)*fcfact*maxdx
   CALL put_text('x',r)
   CALL initialize_text('Text'//TRIM(al), ftsizesub, rgb)
   CALL add_group('Group001','Text'//TRIM(al))
   CALL hide_object('Text'//TRIM(al))

!
!  The label k of k_y.
!
   counter_text=counter_text+1
   CALL create_index(counter_text, al)
   r(1)=0.20_DP*fcfact*maxdx
   r(2)=(xk(2)+maxdx*0.45_DP)*fcfact
   r(3)=0.20_DP*fcfact*maxdx
   CALL put_text('k',r)
   CALL initialize_text('Text'//TRIM(al), ftsize, rgb)
   CALL add_group('Group001','Text'//TRIM(al))
   CALL hide_object('Text'//TRIM(al))
!
!  The label y of k_y.
!
   counter_text=counter_text+1
   CALL create_index(counter_text, al)
   r(1)=(0.20_DP)*fcfact*maxdx
   r(2)=(xk(2)+maxdx*0.45_DP+0.12_DP)*fcfact
   r(3)=(0.20_DP-0.06)*fcfact*maxdx
   CALL put_text('y',r)
   CALL initialize_text('Text'//TRIM(al), ftsizesub, rgb)
   CALL add_group('Group001','Text'//TRIM(al))
   CALL hide_object('Text'//TRIM(al))
!
!  The label k of k_z. 
!
   counter_text=counter_text+1
   CALL create_index(counter_text, al)
   r(1)=0.20_DP*fcfact*maxdx
   r(2)=0.20_DP*fcfact*maxdx
   r(3)=(xk(3)+maxdx*0.55_DP)*fcfact
   CALL put_text('k',r)
   CALL initialize_text('Text'//TRIM(al), ftsize, rgb)
   CALL add_group('Group001','Text'//TRIM(al))
   CALL hide_object('Text'//TRIM(al))
!
!  The label z of k_z
!
   counter_text=counter_text+1
   CALL create_index(counter_text, al)
   r(1)=0.20_DP*fcfact*maxdx
   r(2)=(0.20_DP+0.12_DP)*fcfact*maxdx
   r(3)=(xk(3)+maxdx*0.45_DP)*fcfact
   CALL put_text('z',r)
   CALL initialize_text('Text'//TRIM(al), ftsizesub, rgb)
   CALL add_group('Group001','Text'//TRIM(al))
   CALL hide_object('Text'//TRIM(al))

ENDIF

RETURN
END SUBROUTINE freecad_plotaxis

!--------------------------------------------------------------------------
SUBROUTINE freecad_2d_plotaxis(xk)
!--------------------------------------------------------------------------
!
!  The routine receives the intecept of the axis with the BZ along 
!  k_x and k_y and generates a set of axis from this point to outside the BZ.
!  These axis are an edge and a small cone.
!  It adds also the labels k_x and k_y.
!
IMPLICIT NONE

REAL(DP) :: xk(3)
REAL(DP) :: x1(3), x2(3), r(3), axis(3), ang, rgb(3), maxdx
CHARACTER(LEN=256) :: strin
CHARACTER(LEN=3) :: al

maxdx=MAX(xk(1),xk(2))
!
!  Axis x, an edge from the intersection with the bz to 0.55 the
!  maximum intersection with the BZ
!
IF (ionode) THEN
!
!   Open a group called axis where all the object created here are
!   collected. In this way the axis can be hidden with a single command
!
   WRITE(fcu,'("App.activeDocument().Tip = App.activeDocument().&
                &addObject(''App::DocumentObjectGroup'',''Group001'')")')
   WRITE(fcu,'("App.activeDocument().Group001.Label = ''axis''")')
!
!  Axis x
!
   x1=0.0_DP
   rgb=0.0_DP
   x1(1)=xk(1)
   x2=0.0_DP
   x2(1)=xk(1)+0.55_DP*maxdx
   CALL freecad_join(x1,x2,rgb,al)
   CALL add_group('Group001','Edge'//TRIM(al))
   CALL hide_object('Edge'//TRIM(al))
!
!  Axis y
!
   x1=0.0_DP
   rgb=0.0_DP
   x1(2)=xk(2)
   x2=0.0_DP
   x2(2)=xk(2)+0.55_DP*maxdx
   CALL freecad_join(x1,x2,rgb,al)
   CALL add_group('Group001','Edge'//TRIM(al))
   CALL hide_object('Edge'//TRIM(al))
!
!  Tip of the x axis
!
   r=0.0_DP
   r(1)=(xk(1)+maxdx*0.55_DP)*fcfact
   axis=0.0_DP
   axis(2)=1.0_DP
   ang=90.0_DP
   CALL create_cone('Cone', fcfact/25.0_DP, r, axis, ang)
   CALL add_group('Group001','Cone')
   CALL hide_object('Cone')
!
!  Tip of the y axis
!
   r=0.0_DP
   r(2)=(xk(2)+maxdx*0.55_DP)*fcfact
   axis=0.0_DP
   axis(1)=1.0_DP
   ang=-90.0_DP
   CALL create_cone('Cone001', fcfact/25.0_DP, r, axis, ang)
   CALL add_group('Group001','Cone001')
   CALL hide_object('Cone001')
!
!  The label k of k_x. The text here is red
!
   rgb=0.0_DP
   rgb(1)=1.0_DP

   counter_text=counter_text+1
   CALL create_index(counter_text, al)
   r(1)=(xk(1)+maxdx*0.55_DP)*fcfact
   r(2)=-0.20_DP*fcfact*maxdx
   r(3)=0.0_DP
   CALL put_text('k',r)
   CALL initialize_text('Text'//TRIM(al), ftsize, rgb)
   CALL add_group('Group001','Text'//TRIM(al))
   CALL hide_object('Text'//TRIM(al))
!
!  The label x of k_x
!
   counter_text=counter_text+1
   CALL create_index(counter_text, al)
   r(1)=(xk(1)+maxdx*0.65_DP)*fcfact
   r(2)=(-0.20_DP-0.06_DP)*fcfact*maxdx
   r(3)=0.0_DP
   CALL put_text('x',r)
   CALL initialize_text('Text'//TRIM(al), ftsizesub, rgb)
   CALL add_group('Group001','Text'//TRIM(al))
   CALL hide_object('Text'//TRIM(al))
!
!  The label k of k_y
!
   counter_text=counter_text+1
   CALL create_index(counter_text, al)
   r(1)=0.1_DP*maxdx*fcfact
   r(2)=(xk(2)+maxdx*0.55_DP)*fcfact
   r(3)=0.0_DP
   CALL put_text('k',r)
   CALL initialize_text('Text'//TRIM(al), ftsize, rgb)
   CALL add_group('Group001','Text'//TRIM(al))
   CALL hide_object('Text'//TRIM(al))
!
!  The label y of k_y
!
   counter_text=counter_text+1
   CALL create_index(counter_text, al)
   r(1)=.20_DP*maxdx*fcfact
   r(2)=(xk(2)+maxdx*0.50_DP)*fcfact
   r(3)=0.0_DP
   CALL put_text('y',r)
   CALL initialize_text('Text'//TRIM(al), ftsizesub, rgb)
   CALL add_group('Group001','Text'//TRIM(al))
   CALL hide_object('Text'//TRIM(al))
ENDIF

RETURN
END SUBROUTINE freecad_2d_plotaxis

!----------------------------------------------------------------------
SUBROUTINE freecad_writepoint(xk,npoint)
!----------------------------------------------------------------------
!
!  This routine receives the coordinates of a point and write them
!  as a FreeCAD vector on the script. The name of the vector is point#npoint
!
IMPLICIT NONE
REAL(DP) :: xk(3)
INTEGER  :: npoint

CHARACTER(LEN=256) :: strin
CHARACTER(LEN=3)   :: al
INTEGER :: ipol
!
!
IF (ionode) THEN
   CALL create_index(npoint, al)
   WRITE(strin,'("point",a,"=FreeCAD.Vector(",f13.8,",",f13.8,",&
                           &",f13.8,")")') TRIM(al), (xk(ipol)*fcfact,ipol=1,3)
   WRITE(fcu,'(a)') TRIM(strin)
ENDIF

RETURN
END SUBROUTINE freecad_writepoint

!----------------------------------------------------------------------
SUBROUTINE freecad_join(x1, x2, rgb, al)
!----------------------------------------------------------------------
!
!  This routine plots a line between the point x1 and x2
!
IMPLICIT NONE
CHARACTER(LEN=3), INTENT(OUT) :: al
INTEGER, SAVE  :: counter=-1
REAL(DP), INTENT(IN) :: x1(3), x2(3), rgb(3)
REAL(DP) :: dx(3), distance

counter=counter+1
CALL freecad_writepoint(x1,999)
CALL freecad_writepoint(x2,998)
CALL create_index(counter, al)
IF (ionode) THEN
   WRITE(fcu,'("_=Part.Edge(Part.LineSegment(point999, point998))")')
   WRITE(fcu,'("App.ActiveDocument.addObject(''Part::Feature'',''Edge'').&
                                                                 &Shape=_")')
   WRITE(fcu,'("del _")')
   WRITE(fcu,'("FreeCADGui.getDocument(""Brillouin"").getObject(""Edge",a,""").&
                 &LineColor = (",f7.2,",",f7.2,",",f7.2,")")') TRIM(al), rgb
   WRITE(fcu,'("FreeCADGui.getDocument(""Brillouin"").getObject(""Edge",a,""").&
                       &LineWidth = 3")') TRIM(al)
   CALL add_group('Group','Edge'//TRIM(al))
   WRITE(fcu,'("App.activeDocument().recompute()")')
ENDIF

RETURN
END SUBROUTINE freecad_join

!------------------------------------------------------------------
SUBROUTINE freecad_putlabel(label, a, pos, l2d)
!------------------------------------------------------------------
!
IMPLICIT NONE
CHARACTER(LEN=3), INTENT(IN) :: label, pos
LOGICAL, INTENT(IN) :: l2d      ! if true the pos are interpreted in 2d
REAL(DP) :: a(3)
INTEGER  :: lens
CHARACTER(LEN=3) :: al
CHARACTER(LEN=12) :: labelo
REAL(DP) :: r(3), da(3), rgb(3)

IF (l2d) THEN
   CALL convert_2d_pos(da, pos)
ELSE
   CALL convert_pos(da,pos)
ENDIF
CALL convert_label(label,labelo)
r=a+da
counter_text=counter_text+1
CALL create_index(counter_text, al)

IF (ionode) THEN
   WRITE(fcu,'("text = Draft.makeText([""",a,"""],point=FreeCAD.Vector(",f7.2,&
                     &",",f7.2,",",f7.2,"))")') TRIM(labelo), r*fcfact
   WRITE(fcu,'("Draft.autogroup(text)")')
   rgb=0.0_DP
   rgb(1)=1.0_DP
   CALL initialize_text('Text'//TRIM(al), ftsize, rgb)

   CALL add_group('Group','Text'//TRIM(al))
ENDIF

RETURN
END SUBROUTINE freecad_putlabel

!------------------------------------------------------------------
SUBROUTINE convert_label(label, labelo)
!------------------------------------------------------------------
!
!  This subroutine converts the labels in the form gLetter in the
!  corresponding unicode code to print the greek letter in freecad
!
IMPLICIT NONE
CHARACTER(LEN=3), INTENT(IN) :: label
CHARACTER(LEN=12), INTENT(OUT) :: labelo

IF (ionode) THEN
   IF (label=='gG ') THEN
      labelo='\u0393'
   ELSEIF (label=='gS ') THEN
      labelo='\u03a3'
   ELSEIF (label=='gS0') THEN
      labelo='\u03a3\u2080'
   ELSEIF (label=='gS1') THEN
      labelo='\u03a3\u2081'
   ELSEIF (label=='gD0') THEN
      labelo='\u0394\u2080'
   ELSEIF (label=='gL0') THEN
      labelo='\u039b\u2080'
   ELSEIF (TRIM(label)=='A1') THEN
      labelo='A\u2081'
   ELSEIF (TRIM(label)=='B1') THEN
      labelo='B\u2081'
   ELSEIF (TRIM(label)=='C1') THEN
      labelo='C\u2081'
   ELSEIF (TRIM(label)=='D1') THEN
      labelo='D\u2081'
   ELSEIF (TRIM(label)=='F0') THEN
      labelo='F\u2080'
   ELSEIF (TRIM(label)=='G0') THEN
      labelo='G\u2080'
   ELSEIF (TRIM(label)=='H1') THEN
      labelo='H\u2081'
   ELSEIF (TRIM(label)=='K1') THEN
      labelo='K\u2081'
   ELSEIF (TRIM(label)=='L1') THEN
      labelo='L\u2081'
   ELSEIF (TRIM(label)=='L2') THEN
      labelo='L\u2082'
   ELSEIF (TRIM(label)=='M2') THEN
      labelo='M\u2082'
   ELSEIF (TRIM(label)=='N0') THEN
      labelo='N\u2080'
   ELSEIF (TRIM(label)=='P1') THEN
      labelo='P\u2081'
   ELSEIF (TRIM(label)=='P2') THEN
      labelo='P\u2082'
   ELSEIF (TRIM(label)=='Q1') THEN
      labelo='Q\u2081'
   ELSEIF (TRIM(label)=='R0') THEN
      labelo='R\u2080'
   ELSEIF (TRIM(label)=='S0') THEN
      labelo='S\u2080'
   ELSEIF (TRIM(label)=='S2') THEN
      labelo='S\u2082'
   ELSEIF (TRIM(label)=='T4') THEN
      labelo='T\u2084'
   ELSEIF (TRIM(label)=='U1') THEN
      labelo='U\u2081'
   ELSEIF (TRIM(label)=='W1') THEN
      labelo='W\u2081'
   ELSEIF (TRIM(label)=='W3') THEN
      labelo='W\u2083'
   ELSEIF (TRIM(label)=='X1') THEN
      labelo='X\u2081'
   ELSEIF (TRIM(label)=='X3') THEN
      labelo='X\u2083'
   ELSEIF (TRIM(label)=='Y1') THEN
      labelo='Y\u2081'
   ELSEIF (TRIM(label)=='Z1') THEN
      labelo='Z\u2081'
   ELSE
      labelo=label
   ENDIF
ENDIF

RETURN
END SUBROUTINE convert_label

!-------------------------------------------------------------------
SUBROUTINE convert_pos(dx,pos)
!-------------------------------------------------------------------
!
!   The BZ letters have some standard position (N,S,E,W,NE,NW,SE,SW)
!   This routine converts this position in a shift of the position of 
!   the letter.
!
IMPLICIT NONE
CHARACTER(LEN=3) :: pos
REAL(DP) :: dx(3)

dx=0.0_DP

SELECT CASE (TRIM(pos))
   CASE('E')          ! East
       dx(1)=0.06_DP
       dx(2)=0.06_DP
   CASE('W')          ! West
       dx(1)=-0.06_DP
       dx(2)=-0.06_DP
   CASE('N')          ! North
       dx(3)=0.08_DP
   CASE('S')          ! South
       dx(3)=-0.12_DP
   CASE('NE')         ! North-East
       dx(1)=0.06_DP
       dx(2)=0.06_DP
       dx(3)=0.08_DP
   CASE('NW')         ! North-West
       dx(1)=-0.06_DP
       dx(2)=-0.06_DP
       dx(3)=0.08_DP
   CASE('SE')         ! South-East
       dx(2)=0.06_DP
       dx(2)=0.06_DP
       dx(3)=-0.12_DP
   CASE('SW')         ! South-West
       dx(1)=-0.06_DP
       dx(2)=-0.06_DP
       dx(3)=-0.12_DP
   CASE('')
   CASE DEFAULT
      CALL errore('convert_pos','unknown pos',1)
END SELECT  

RETURN
END SUBROUTINE convert_pos

!-------------------------------------------------------------------
SUBROUTINE convert_2d_pos(dx,pos)
!-------------------------------------------------------------------
!
!   The BZ letters have some standard position (N,S,E,W,NE,NW,SE,SW)
!   This routine converts this position in a shift of the position of 
!   the letter. This routine makes the shifts in the xy plane.
!   It puts also a small shift on the z direction to put the letters
!   above the figure.
!   
!
IMPLICIT NONE
CHARACTER(LEN=3) :: pos
REAL(DP) :: dx(3)

dx=0.0_DP
dx(3)=0.05_DP

SELECT CASE (TRIM(pos))
   CASE('E')          ! East
       dx(1)=0.06_DP
   CASE('W')          ! West
       dx(1)=-0.06_DP
   CASE('N')          ! North
       dx(2)=0.1_DP
   CASE('S')          ! South
       dx(2)=-0.1_DP
   CASE('NE')         ! North-East
       dx(1)=0.06_DP
       dx(2)=0.1_DP
   CASE('NW')         ! North-West
       dx(1)=-0.06_DP
       dx(2)=0.1_DP
   CASE('SE')         ! South-East
       dx(1)=0.06_DP
       dx(2)=-0.1_DP
   CASE('SW')         ! South-West
       dx(1)=-0.06_DP
       dx(2)=-0.1_DP
   CASE('')
   CASE DEFAULT
      CALL errore('convert_2d_pos','unknown pos',1)
END SELECT  

RETURN
END SUBROUTINE convert_2d_pos
!
!----------------------------------------------------------------------
SUBROUTINE freecad_writesurface(indeces, nface)
!----------------------------------------------------------------------
!
!  This routine assumes that the coordinates of all the vertices of the
!  BZ have been given. It receives the indices of the vertices of a
!  face of the BZ and write on the script the commands to define the surface.
!
IMPLICIT NONE
REAL(DP) :: xk(3)

INTEGER :: indeces(9), nface
CHARACTER(LEN=1000) :: linesur
CHARACTER(LEN=3)    :: al, al1, al2
INTEGER :: i, counter

IF (ionode) THEN
   counter=0
   DO i=1,indeces(1)-1
      counter=counter+1
      CALL create_index(indeces(i+1), al)
      CALL create_index(indeces(i+2), al1)
      CALL create_index(counter, al2)
      WRITE(fcu,'("L",a,"=Part.LineSegment(point",a,",point",a,")")') al2,al,al1
   ENDDO
   counter=counter+1
   CALL create_index(indeces(indeces(1)+1), al)
   CALL create_index(indeces(2), al1)
   CALL create_index(counter, al2)
   WRITE(fcu,'("L",a,"=Part.LineSegment(point",a,",point",a,")")') al2,al,al1

   linesur="W=Part.Wire(["
   DO i=1,counter
      CALL create_index(i, al2)
      WRITE(fcu,'("E",a,"=Part.Edge(L",a,")")') al2, al2
      linesur=TRIM(linesur)//"E"//TRIM(al2)//","
   ENDDO
   linesur=TRIM(linesur)//"])"
   WRITE(fcu,'(a)') TRIM(linesur) 

   CALL create_index(nface, al)
   WRITE(fcu,'("F",a,"=Part.Face(W)")') al
ENDIF

RETURN
END SUBROUTINE freecad_writesurface

!----------------------------------------------------------------------
SUBROUTINE freecad_createsolid(nfaces)
!----------------------------------------------------------------------
!
!  This routine assumes that all the faces of the solid have been created
!  and uses them to build a shell and then to transform the shell in 
!  a solid. Only the solid is added to the FreeCAD Application.
!
IMPLICIT NONE
INTEGER :: nfaces

CHARACTER(LEN=1000) :: linesur
CHARACTER(LEN=3)    :: al
INTEGER :: i

IF (ionode) THEN
   linesur="S=Part.Shell(["
   DO i=1,nfaces
      CALL create_index(i, al)
      linesur=TRIM(linesur)//"F"//TRIM(al)//","
   ENDDO
   linesur=TRIM(linesur)//" ])"
   WRITE(fcu,'(a)') TRIM(linesur)

   WRITE(fcu,'("P=Part.Solid(S)")')
   WRITE(fcu,'("App.ActiveDocument.addObject(''Part::Feature'',&
                                                      &''Solid'').Shape=P")')
   CALL add_group('Group','Solid')
ENDIF

RETURN
END SUBROUTINE freecad_createsolid

!----------------------------------------------------------------------
SUBROUTINE freecad_createshell()
!----------------------------------------------------------------------
!
!  This routine assumes that a 2d BZ has been created and called F001.
!  It associates it to the Active Document to make it visible.
!  
!
IMPLICIT NONE

IF (ionode) THEN
   WRITE(fcu,'("App.ActiveDocument.addObject(''Part::Feature'',&
                                                 &''Shell'').Shape=F001")')
   CALL add_group('Group','Shell')
ENDIF

RETURN
END SUBROUTINE freecad_createshell

!----------------------------------------------------------------
SUBROUTINE create_index(intinp, al)
!----------------------------------------------------------------
!
!   This routine receives a positive integer and creates a string 
!   with three characters with a sufficient number of leading 0s.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: intinp
CHARACTER(LEN=3),INTENT(OUT) :: al
CHARACTER(LEN=6) :: int_to_char

IF (intinp==0) THEN 
   al='    '
ELSEIF (intinp>0 .AND. intinp<10) THEN
   al='00'//TRIM(int_to_char(intinp))
ELSEIF (intinp<100) THEN
   al='0'//TRIM(int_to_char(intinp))
ELSEIF (intinp<1000) THEN
   al=TRIM(int_to_char(intinp))
ELSE
   CALL errore('create_index','index is too large',1)
ENDIF

RETURN
END SUBROUTINE create_index

!--------------------------------------------------------------------
SUBROUTINE freecad_setfcfact(fact)
!--------------------------------------------------------------------
!
!  This routine changes the default for the factor that converts 
!  BZ units in mm used in freecad. The default is 100., meaning that
!  2 \pi / a corresponds to 10 cm in final plot. 
!  Using here for instance
!  100 2 \pi /a and a in Angstrom, one has that 10 cm in the freecad
!  plot corresponds to an A^{-1}. 
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: fact

fcfact=fact

RETURN
END SUBROUTINE freecad_setfcfact
!
!------------------------------------------------------------------------
SUBROUTINE freecad_setfontsize(ftsize_in)
!------------------------------------------------------------------------
!
IMPLICIT NONE
REAL(DP) :: ftsize_in

ftsize=ftsize_in * fcfact
ftsizesub=ftsize * 5.0_DP / 6.0_DP

RETURN
END SUBROUTINE freecad_setfontsize
!
!--------------------------------------------------------------------
SUBROUTINE freecad_closeplot()
!--------------------------------------------------------------------
!
!   This routine closes the file with the freecad script
!
IMPLICIT NONE

IF (ionode) &
   CLOSE(UNIT=fcu,STATUS='KEEP')

RETURN
END SUBROUTINE freecad_closeplot
!
!--------------------------------------------------------------------------
SUBROUTINE freecad_createpdf(filename)
!--------------------------------------------------------------------------
!
!  This subroutine creates a view of the Brillouin zone, reads a template
!  for the paper (the default is a A4 paper in landscape mode) and
!  generates the pdf file with the Brillouin zone figure
!
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename

IF (ionode) THEN
   WRITE(fcu,'("page=App.activeDocument().addObject(''TechDraw:&
                   &:DrawPage'',''Page'')")')
   WRITE(fcu,'("App.activeDocument().addObject(''TechDraw:&
                                  &:DrawSVGTemplate'',''Template'')")')
   WRITE(fcu,'("App.activeDocument().Template.Template = ''./tpw_bz.svg''")')
   WRITE(fcu,'("App.activeDocument().Page.Template = App.activeDocument().&
                                                           &Template")')
   WRITE(fcu,'("page.ViewObject.show()")')
   WRITE(fcu,'("App.activeDocument().addObject(''TechDraw::DrawViewPart'',&
                                                &''View'')")')
   WRITE(fcu,'("App.activeDocument().Page.addView(App.activeDocument().View)")')
   !
   !   Group001 contains the BZ and the axis. Use Group here to remove 
   !   the axis
   !
   WRITE(fcu,'("App.ActiveDocument.View.Source = [App.ActiveDocument.&
                                &getObject(""Group"")]")')
   ! 
   !    This is the view direction that produces figures similar to those
   !    given by the asy module. 
   !
   WRITE(fcu,'("App.ActiveDocument.View.Direction = (0.88,0.47,0.12)")')
   WRITE(fcu,'("App.ActiveDocument.View.Scale = 0.65")')
   !
   !    These two commands plot the hidden lines. Set to False to remove them
   !
   WRITE(fcu,'("App.ActiveDocument.View.HardHidden = True")')
   WRITE(fcu,'("FreeCADGui.ActiveDocument.View.HiddenWidth = 1.40")')
   !
   !  Positioning of the BZ in the page. Might need some tuning depending
   !  on the size of the BZ. Can be done from inside freecad
   !
   WRITE(fcu,'("App.ActiveDocument.View.X = 95")')
   WRITE(fcu,'("App.ActiveDocument.View.Y = 115")')
   !
   !   Here you can change the name of the file with the BZ pdf figure
   !
   WRITE(fcu,'("App.activeDocument().recompute()")')
   WRITE(fcu,'("obj = FreeCAD.ActiveDocument.getObject(""Page"")")')
   WRITE(fcu,'("App.activeDocument().View.recompute()")')
   WRITE(fcu,'("FreeCADGui.export([obj], u""./",a,""")")') TRIM(filename)
   WRITE(fcu,'("App.activeDocument().recompute()")')
ENDIF

RETURN
END SUBROUTINE freecad_createpdf

!-----------------------------------------------------------------------
SUBROUTINE add_group(group,object)
!-----------------------------------------------------------------------

IMPLICIT NONE
CHARACTER(LEN=*) :: group
CHARACTER(LEN=*) :: object

IF (ionode) &
   WRITE(fcu,'("App.getDocument(""Brillouin"").getObject(""",a,""").&
     &addObject(App.getDocument(""Brillouin"").getObject(""",a,"""))")') &
             TRIM(group), TRIM(object)

RETURN
END SUBROUTINE add_group

!--------------------------------------------------------------------------
SUBROUTINE hide_object(object)
!--------------------------------------------------------------------------

IMPLICIT NONE
CHARACTER(LEN=*) :: object

IF (ionode) &
   WRITE(fcu,'("Gui.getDocument(""Brillouin"").getObject(""",a,""").&
                                &Visibility=False")') TRIM(object)

RETURN
END SUBROUTINE hide_object

!--------------------------------------------------------------------------
SUBROUTINE put_text(text,r)
!--------------------------------------------------------------------------
!
IMPLICIT NONE
CHARACTER(LEN=*) :: text
REAL(DP) :: r(3)

IF (ionode) THEN
   WRITE(fcu,'("text = Draft.makeText([""",a,"""],point=FreeCAD.Vector(",f7.2,&
  &",",f7.2,",",f7.2,"))")') TRIM(text), r
   WRITE(fcu,'("Draft.autogroup(text)")')
ENDIF

RETURN
END SUBROUTINE put_text

!--------------------------------------------------------------------------
SUBROUTINE initialize_text(object, fontsize, rgb)
!--------------------------------------------------------------------------

IMPLICIT NONE
CHARACTER(LEN=*) :: object
REAL(DP) :: rgb(3), fontsize

IF (ionode) THEN
   WRITE(fcu,'("FreeCADGui.getDocument(""Brillouin"").getObject(""",a,""").&
                          &DisplayMode = u""2D text""")') TRIM(object)
   WRITE(fcu,'("FreeCADGui.getDocument(""Brillouin"").getObject(""",a,""").&
                        &FontSize = ''",f4.0," mm''")') TRIM(object), fontsize
   WRITE(fcu,'("FreeCADGui.getDocument(""Brillouin"").getObject(""",a,""").&
                          &TextColor = (",f4.0,",",f4.0,",",f4.0,")")') &
                                                          TRIM(object), rgb
ENDIF
RETURN
END SUBROUTINE initialize_text

!------------------------------------------------------------------------
SUBROUTINE create_cone(label, radius, pos, rot, angle)
!------------------------------------------------------------------------

IMPLICIT NONE
REAL(DP) :: radius, pos(3), rot(3), angle
CHARACTER(LEN=*) :: label

IF (ionode) THEN
   WRITE(fcu,'("App.ActiveDocument.addObject(""Part::Cone"",""",a,""")")') &
                                                                 TRIM(label)
   WRITE(fcu,'("App.ActiveDocument.ActiveObject.Label = """,a,"""")') &
                                                                 TRIM(label)
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""",a,""").&
                     &Radius1 = ''",f7.2," mm''")') TRIM(label), radius
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""",a,""").&
                                 &Radius2= ''0 mm''")') TRIM(label)
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""",a,""").      &
           &Placement = App.Placement(App.Vector(",f7.2,",",f7.2,",",f7.2,"), &
           &App.Rotation(App.Vector(",f7.2,",",f7.2,",",f7.2,"),",f7.2,"))")')&
                                             TRIM(label), pos, rot, angle
ENDIF

RETURN
END SUBROUTINE create_cone

END MODULE freecad
