!
! Copyright (C) 2020 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE freecad
!
USE kinds, ONLY : DP
USE io_global, ONLY : ionode, ionode_id, stdout
USE mp_images, ONLY : intra_image_comm
IMPLICIT NONE
PRIVATE
SAVE

INTEGER :: fcu = 56        ! unit of the script file

REAL(DP) :: fcfact=100     ! This factor converts the units of the BZ into
                           ! freecad units
!
!  freecad uses mm units while the point coordinates are in units 2pi/a. We
!  multiply by fcfact to have a simply movable BZ and a reasonable size on
!  on techdraw without rescaling. So 1 in 2\pi/a units corresponds to 
!  10 cm in freecad units. This default can be changed from input.

PUBLIC freecad_writepoint, freecad_openplot, freecad_closeplot,  &
       freecad_writesurface, freecad_createsolid,  &
       freecad_centerview, freecad_setcolor, freecad_plotaxis,   &
       freecad_join, freecad_putlabel, freecad_setfcfact, freecad_createpdf

CONTAINS
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
ENDIF

RETURN
END SUBROUTINE freecad_openplot

!----------------------------------------------------------------------------
SUBROUTINE freecad_centerview()
!----------------------------------------------------------------------------
!
!  This routine sends the commands to set an isometric view and center
!  the Brillouin zone on the screen
!
IMPLICIT NONE

IF (ionode) THEN
   WRITE(fcu,'("Gui.SendMsgToActiveView(""ViewFit"")")')
   WRITE(fcu,'("Gui.activeDocument().activeView().viewIsometric()")')
ENDIF

RETURN
END SUBROUTINE freecad_centerview

!-----------------------------------------------------------------------
SUBROUTINE freecad_setcolor(r,g,b)
!----------------------------------------------------------------------
!
!  The color is given in r,g,b form (from 0 to 1). The routine must
!  be called after the BZ solid has been created.
!
IMPLICIT NONE
REAL(DP) :: r, g, b

IF (ionode) &
   WRITE(fcu,'("FreeCADGui.getDocument(""Brillouin"").&
    &getObject(""Solid"").ShapeColor = (",f4.2,",",f4.2,",",f4.2,")")') r, g, b

RETURN
END SUBROUTINE freecad_setcolor

!--------------------------------------------------------------------------
SUBROUTINE freecad_plotaxis(xk)
!--------------------------------------------------------------------------
!
!  The routine receives the intecept of the axis with the BZ along x, y, and
!  z and generates a set of axis from this point to outside the BZ.
!  These axis are a small cylinder and a small cone.
!
IMPLICIT NONE

REAL(DP) :: xk(3)
REAL(DP) :: maxdx
CHARACTER(LEN=256) :: strin

maxdx=MAX(xk(1),xk(2))
maxdx=MAX(maxdx,xk(3))

!
!  Axis z, a cilinder from the intersection with the bz to 0.55 the
!  maximum dimension of the BZ
!

IF (ionode) THEN
!
!  Axis z
!
   WRITE(fcu,'("App.ActiveDocument.addObject(""Part::Cylinder"",&
                                                &""Cylinder"")")')
   WRITE(fcu,'("App.ActiveDocument.ActiveObject.Label = ""Cylinder""")')
   WRITE(strin,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cylinder"")&
                    &.Height = ''",f7.2," mm''")') maxdx*0.45_DP*fcfact
   WRITE(fcu,'(a)') TRIM(strin)
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cylinder"").&
                               &Radius = ''",f7.2," mm''")') fcfact/200._DP
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cylinder"").&
           &Placement = App.Placement(App.Vector(0,0,",f7.2,"),App.&
           &Rotation(App.Vector(0,0,1),0))")') xk(3)*fcfact
!
!  Axis x
!
   WRITE(fcu,'("App.ActiveDocument.addObject(""Part::Cylinder"",&
                                                &""Cylinder001"")")')
   WRITE(fcu,'("App.ActiveDocument.ActiveObject.Label = ""Cylinder001""")')
   WRITE(strin,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cylinder001"")&
                    &.Height = ''",f7.2," mm''")') maxdx*0.55_DP*fcfact
   WRITE(fcu,'(a)') TRIM(strin)
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cylinder001"").&
                              &Radius = ''",f7.2," mm''")') fcfact/200.0_DP
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cylinder001"").&
           &Placement = App.Placement(App.Vector("f7.2,",0,0),App.&
              &Rotation(App.Vector(0,1,0),90))")') xk(1)*fcfact
!
!  Axis y
!

   WRITE(fcu,'("App.ActiveDocument.addObject(""Part::Cylinder"",&
                                                &""Cylinder002"")")')
   WRITE(fcu,'("App.ActiveDocument.ActiveObject.Label = ""Cylinder002""")')
   WRITE(strin,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cylinder002"")&
                    &.Height = ''",f7.2," mm''")') maxdx*0.45_DP*fcfact
   WRITE(fcu,'(a)') TRIM(strin)
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cylinder002"").&
                                &Radius = ''",f7.2," mm''")') fcfact/200.0_DP
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cylinder002"").&
           &Placement = App.Placement(App.Vector(0,",f7.2,",0),App.&
              &Rotation(App.Vector(1,0,0),-90))")') xk(2)*fcfact
!
!  Tip of the z axis
!

   WRITE(fcu,'("App.ActiveDocument.addObject(""Part::Cone"",""Cone"")")')
   WRITE(fcu,'("App.ActiveDocument.ActiveObject.Label = ""Cone""")')
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cone"").&
                                 &Radius1 = ''",f7.2," mm''")') fcfact/25.0_DP
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cone"").Radius2 &
                                                &= ''0 mm''")')
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cone"").&
           &Placement = App.Placement(App.Vector(0,0,",f7.2,"),App.&
           &Rotation(App.Vector(0,0,1),0))")') (xk(3)+maxdx*0.45_DP)*fcfact

!
!  Tip of the x axis
!

   WRITE(fcu,'("App.ActiveDocument.addObject(""Part::Cone"",""Cone001"")")')
   WRITE(fcu,'("App.ActiveDocument.ActiveObject.Label = ""Cone001""")')
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cone001"").&
                                &Radius1 = ''",f7.2," mm''")') fcfact/25.0_DP
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cone001"").&
                                                &Radius2 = ''0 mm''")')
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cone001"").&
           &Placement = App.Placement(App.Vector(",f7.2,",0,0),App.&
           &Rotation(App.Vector(0,1,0),90))")') (xk(1)+maxdx*0.55_DP)*fcfact
!
!  Tip of the y axis
!

   WRITE(fcu,'("App.ActiveDocument.addObject(""Part::Cone"",""Cone002"")")')
   WRITE(fcu,'("App.ActiveDocument.ActiveObject.Label = ""Cone002""")')
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cone002"").&
                              &Radius1 = ''",f7.2," mm''")') fcfact/25.0_DP
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cone002"").&
                              &Radius2 = ''0 mm''")')
   WRITE(fcu,'("FreeCAD.getDocument(""Brillouin"").getObject(""Cone002"").&
           &Placement = App.Placement(App.Vector(0,",f7.2,",0),App.&
           &Rotation(App.Vector(1,0,0),-90))")') (xk(2)+maxdx*0.45_DP)*fcfact
!
!  Put all the axis into a group, so they can be removed all together
!
   WRITE(fcu,'("App.activeDocument().Tip = App.activeDocument().&
                &addObject(''App::DocumentObjectGroup'',''Group'')")')
   WRITE(fcu,'("App.activeDocument().Group.Label = ''axis''")')
   WRITE(fcu,'("App.getDocument(""Brillouin"").getObject(""Group"").&
          &addObject(App.getDocument(""Brillouin"").getObject(""Cone""))")')
   WRITE(fcu,'("App.getDocument(""Brillouin"").getObject(""Group"").&
          &addObject(App.getDocument(""Brillouin"").getObject(""Cone001""))")')
   WRITE(fcu,'("App.getDocument(""Brillouin"").getObject(""Group"").&
          &addObject(App.getDocument(""Brillouin"").getObject(""Cone002""))")')
   WRITE(fcu,'("App.getDocument(""Brillouin"").getObject(""Group"").&
       &addObject(App.getDocument(""Brillouin"").getObject(""Cylinder""))")')
   WRITE(fcu,'("App.getDocument(""Brillouin"").getObject(""Group"").&
       &addObject(App.getDocument(""Brillouin"").getObject(""Cylinder001""))")')
   WRITE(fcu,'("App.getDocument(""Brillouin"").getObject(""Group"").&
       &addObject(App.getDocument(""Brillouin"").getObject(""Cylinder002""))")')

   WRITE(fcu,'("Gui.getDocument(""Brillouin"").getObject(""Cone"").&
                                                          &Visibility=False")')
   WRITE(fcu,'("Gui.getDocument(""Brillouin"").getObject(""Cone001"").&
                                                          &Visibility=False")')
   WRITE(fcu,'("Gui.getDocument(""Brillouin"").getObject(""Cone002"").&
                                                          &Visibility=False")')
   WRITE(fcu,'("Gui.getDocument(""Brillouin"").getObject(""Cylinder"").&
                                                          &Visibility=False")')
   WRITE(fcu,'("Gui.getDocument(""Brillouin"").getObject(""Cylinder001"").&
                                                          &Visibility=False")')
   WRITE(fcu,'("Gui.getDocument(""Brillouin"").getObject(""Cylinder002"").&
                                                          &Visibility=False")')

ENDIF

RETURN
END SUBROUTINE freecad_plotaxis


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
                           ",f13.8,")")') TRIM(al), (xk(ipol)*fcfact,ipol=1,3)
   WRITE(fcu,'(a)') TRIM(strin)
ENDIF

RETURN
END SUBROUTINE freecad_writepoint

!----------------------------------------------------------------------
SUBROUTINE freecad_join(x1, x2)
!----------------------------------------------------------------------
!
!  This routine plots a cylinder between the point x1 and x2
!
IMPLICIT NONE
INTEGER, SAVE  :: counter=-1
REAL(DP), INTENT(IN) :: x1(3), x2(3)
REAL(DP) :: dx(3), distance
CHARACTER(LEN=3) :: al

counter=counter+1
CALL freecad_writepoint(x1,999)
CALL freecad_writepoint(x2,998)
CALL create_index(counter, al)
IF (ionode) THEN
   WRITE(fcu,'("_=Part.Edge(Part.LineSegment(point999, point998))")')
   WRITE(fcu,'("App.ActiveDocument.addObject(''Part::Feature'',''Edge'').&
                                                                 &Shape=_")')
   WRITE(fcu,'("del _")')
   WRITE(fcu,'("App.getDocument(""Brillouin"").getObject(""Group"").&
     &addObject(App.getDocument(""Brillouin"").getObject(""Edge",a,"""))")') &
             TRIM(al)
   WRITE(fcu,'("FreeCADGui.getDocument(""Brillouin"").getObject(""Edge",a,""").&
                       &LineColor = (1.00,0.00,0.00)")') TRIM(al)
   WRITE(fcu,'("FreeCADGui.getDocument(""Brillouin"").getObject(""Edge",a,""").&
                       &LineWidth = 3")') TRIM(al)
   WRITE(fcu,'("App.activeDocument().recompute()")')
ENDIF

RETURN
END SUBROUTINE freecad_join

!------------------------------------------------------------------
SUBROUTINE freecad_putlabel(label, a, pos)
!------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER, SAVE :: counter=-1
CHARACTER(LEN=3), INTENT(IN) :: label, pos
REAL(DP) :: a(3)
INTEGER :: lens
CHARACTER(LEN=3) :: al
CHARACTER(LEN=7) :: labelo
REAL(DP) :: r(3), da(3)

CALL convert_pos(da,pos)
CALL convert_label(label,labelo)
r=a+da
counter=counter+1
CALL create_index(counter, al)

IF (ionode) THEN
   WRITE(fcu,'("text = Draft.makeText([""",a,"""],point=FreeCAD.Vector(",f7.2,&
                     &",",f7.2,",",f7.2,"))")') TRIM(labelo), r*fcfact
   WRITE(fcu,'("Draft.autogroup(text)")')
   WRITE(fcu,'("FreeCADGui.getDocument(""Brillouin"").getObject(""Text",a,""").&
                          &DisplayMode = u""2D text""")') TRIM(al)
   WRITE(fcu,'("FreeCADGui.getDocument(""Brillouin"").getObject(""Text",a,""").&
                          &FontSize = ''35 mm''")') TRIM(al)
   WRITE(fcu,'("FreeCADGui.getDocument(""Brillouin"").getObject(""Text",a,""").&
                          &TextColor = (1.00,0.00,0.00)")') TRIM(al)
ENDIF

!WRITE(fcu,'("ss=Draft.makeShapeString(String=""",a,""",&
!            &FontFile=""/usr/share/fonts/opentype/cantarell/&
!            &Cantarell-Regular.otf"",Size=5.0,Tracking=0.0)")') TRIM(label)
!WRITE(fcu,'("plm=FreeCAD.Placement()")')
!WRITE(fcu,'("plm.Base=FreeCAD.Vector(",f7.2,",",f7.2,",",f7.2,")")') r*fcfact
!WRITE(fcu,'("plm.Rotation.Q=(0.,0.,0.,0.)")')
!WRITE(fcu,'("ss.Placement=plm")')
!WRITE(fcu,'("ss.Support=None")')
!WRITE(fcu,'("Draft.autogroup(ss)")')

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
CHARACTER(LEN=7), INTENT(OUT) :: labelo

IF (ionode) THEN
   IF (label=='gG ') THEN
      labelo='\u0393'
   ELSEIF (label=='gS ') THEN
      labelo='\u03a3'
   ELSEIF (label=='gS0') THEN
      labelo='\u03a30'
   ELSEIF (label=='gS1') THEN
      labelo='\u03a31'
   ELSEIF (label=='gD0') THEN
      labelo='\u03940'
   ELSEIF (label=='gL0') THEN
      labelo='\u039b0'
   ELSE
      labelo=label
   ENDIF
ENDIF

RETURN
END SUBROUTINE convert_label

SUBROUTINE convert_pos(dx,pos)

IMPLICIT NONE
CHARACTER(LEN=3) :: pos
REAL(DP) :: dx(3)

dx=0.0_DP

SELECT CASE (TRIM(pos))
   CASE('E')
       dx(1)=0.06_DP
       dx(2)=0.06_DP
   CASE('W')
       dx(1)=-0.06_DP
       dx(2)=-0.06_DP
   CASE('N')
       dx(3)=0.08_DP
   CASE('S')
       dx(3)=-0.12_DP
   CASE('NE')
       dx(1)=0.06_DP
       dx(2)=0.06_DP
       dx(3)=0.08_DP
   CASE('NW')
       dx(1)=-0.06_DP
       dx(2)=-0.06_DP
       dx(3)=0.08_DP
   CASE('SE')
       dx(2)=0.06_DP
       dx(2)=0.06_DP
       dx(3)=-0.12_DP
   CASE('SW')
       dx(1)=-0.06_DP
       dx(2)=-0.06_DP
       dx(3)=-0.12_DP
   CASE DEFAULT
      CALL errore('convert_pos','unknown pos',1)
END SELECT  

RETURN
END SUBROUTINE convert_pos
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
   WRITE(fcu,'("FreeCADGui.ActiveDocument.getObject(""Solid"").&
                                                      &Transparency=15")')
   WRITE(fcu,'("FreeCADGui.ActiveDocument.getObject(""Solid"").&
                                                      &LineWidth = 3.00")')
ENDIF

RETURN
END SUBROUTINE freecad_createsolid

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
!---------------------------------------------------------------------------
SUBROUTINE freecad_create_bz()
!---------------------------------------------------------------------------
!
!  This routine puts in a single group, the BZ itself and the axis
!  It is then possible to choose what to plot.
!
IMPLICIT NONE

IF (ionode) THEN
   WRITE(fcu, '("App.activeDocument().Tip = App.activeDocument().&
                 &addObject(''App::DocumentObjectGroup'',''Group001'')")')
   WRITE(fcu, '("App.activeDocument().Group001.Label = ''Brillouin zone''")')
   WRITE(fcu, '("App.getDocument(""Brillouin"").getObject(""Group001"").&
           &addObject(App.getDocument(""Brillouin"").getObject(""Group""))")')
   WRITE(fcu, '("App.getDocument(""Brillouin"").getObject(""Group001"").&
           &addObject(App.getDocument(""Brillouin"").getObject(""Solid""))")')
ENDIF

RETURN
END SUBROUTINE freecad_create_bz

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

CALL freecad_create_bz()
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
                                &getObject(""Group001"")]")')
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

END MODULE freecad
