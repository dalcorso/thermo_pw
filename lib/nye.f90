!
! Copyright (C) 2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE nye
!
!   this module contains the support routines for the shape of tensors
!   conventions, definitions, etc.
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  CHARACTER(LEN=16) :: thermal_gnuplot_name(6)

  DATA  thermal_gnuplot_name / '{/Symbol a}_{xx}', '{/Symbol a}_{yy}',  &
                               '{/Symbol a}_{zz}', '{/Symbol a}_{yz}',  &
                               '{/Symbol a}_{xz}', '{/Symbol a}_{xy}'   /

  PUBLIC print_vectors_shape, print_tensor2_shape, print_piezo_shape, &
         print_el_cons_shape, print_b_fact_shape, thermal_gnuplot_name, &
         needed_tensor2

CONTAINS

!--------------------------------------------------------------------
SUBROUTINE print_vectors_shape(code_group, ibrav)
!--------------------------------------------------------------------

USE io_global, ONLY : stdout

IMPLICIT NONE
INTEGER, INTENT(IN) :: code_group, ibrav

SELECT CASE (ibrav)
   CASE(1,2,3)  
!
!   cubic
!
      IF (code_group==29.OR. code_group==32) THEN
         WRITE(stdout,'(/,5x, "This solid has inversion symmetry.")')
         WRITE(stdout,'(5x, "First-rank tensors, such as the &
                           &spontaneous polarization, vanish.")')
      ELSEIF (code_group==28.OR.code_group==30.OR.code_group==31) THEN
         WRITE(stdout,'(/,5x, "This solid has not inversion but, &
                                 &first-rank tensors, such as ")')
         WRITE(stdout,'(5x, "the spontaneous polarization, vanish.")')
      ENDIF
   CASE(4,5,6,7)  
!
!  hexagonal, trigonal, tetragonal
!
      IF (code_group==18.OR.code_group==22.OR.code_group==23&
                                                 .OR.code_group==25) THEN
         WRITE(stdout,'(/,5x, "This solid has inversion symmetry.")')
         WRITE(stdout,'(5x, "First-rank tensors, such as the &
                                 &spontaneous polarization, vanish.")')
      ELSEIF (code_group==9.OR.code_group==10.OR.code_group==17.OR. &
                      code_group==21.OR.code_group==24.OR.code_group==26 ) THEN
         WRITE(stdout,'(/,5x, "This solid has not inversion but, &
                                      &first-rank tensors, such as ")')
         WRITE(stdout,'(5x, "the spontaneous polarization, vanish.")')
      ELSEIF (code_group==5.OR.code_group==6.OR.code_group==7.OR. &
                  code_group==13.OR.code_group==14.OR.code_group==15) THEN
         WRITE(stdout,'(/,5x, "First-rank tensors, such as the &
                        &spontaneous polarization, have the form:")')
         WRITE(stdout,'(/,5x, "(  .   .   p3 )")')
      ENDIF
   CASE(8,9,10,11)  
!
!  orthorhombic
!
      IF (code_group==20) THEN
         WRITE(stdout,'(/,5x, "This solid has inversion symmetry.")')
         WRITE(stdout,'(5x, "First-rank tensors, such as the &
                            &spontaneous polarization, vanish.")')
      ELSEIF (code_group==8) THEN
         WRITE(stdout,'(/,5x, "This solid has not inversion but, &
                                   &first-rank tensors, such as ")')
         WRITE(stdout,'(5x, "the spontaneous polarization, vanish.")')
      ELSEIF (code_group==12) THEN
         WRITE(stdout,'(/,5x, "First-rank tensors, such as the &
                      &spontaneous polarization, have the form:")')
         WRITE(stdout,'(/,5x, "(  .   .   p3 )")')
      ENDIF
   CASE(12,13)  
!
!   monoclinic c orientation
!
      IF (code_group==16) THEN
         WRITE(stdout,'(/,5x, "This solid has inversion symmetry.")')
         WRITE(stdout,'(5x, "First-rank tensors, such as the &
                                 &spontaneous polarization, vanish.")')
      ELSEIF (code_group==4) THEN
         WRITE(stdout,'(/,5x, "First-rank tensors, such as the &
                               &spontaneous polarization, have the form:")')
         WRITE(stdout,'(/,5x, "(  .   .   p3 )")')
      ELSEIF (code_group==3) THEN
         WRITE(stdout,'(/,5x, "First-rank tensors, such as the &
                      &spontaneous polarization, have the form:")')
         WRITE(stdout,'(/,5x, "( p1   p2   . )")')
      ENDIF
   CASE(-12,-13)  
!
!   monoclinic  unique b
!
      IF (code_group==16) THEN
         WRITE(stdout,'(/,5x, "This solid has inversion symmetry.")')
         WRITE(stdout,'(5x, "First-rank tensors, such as the &
                                 &spontaneous polarization, vanish.")')
      ELSEIF (code_group==4) THEN
         WRITE(stdout,'(/,5x, "First-rank tensors, such as the &
                               &spontaneous polarization, have the form:")')
         WRITE(stdout,'(/,5x, "(  .   p2   .  )")')
      ELSEIF (code_group==3) THEN
         WRITE(stdout,'(/,5x, "First-rank tensors, such as the &
                               &spontaneous polarization, have the form:")')
         WRITE(stdout,'(/,5x, "(  p1   .   p3 )")')
      ENDIF
   CASE(14)  
!
!  triclinc 
!
      IF (code_group==2) THEN
         WRITE(stdout,'(/,5x, "This solid has inversion symmetry.")')
         WRITE(stdout,'(5x, "First-rank tensors, such as the &
                                 &spontaneous polarization, vanish.")')
      ELSEIF (code_group==1) THEN
         WRITE(stdout,'(/,5x, "First-rank tensors, such as the &
                               &spontaneous polarization, have the form:")')
         WRITE(stdout,'(/,5x, "( p1   p2   p3 )")')
      ENDIF
   CASE DEFAULT 
END SELECT

RETURN
END SUBROUTINE print_vectors_shape

!--------------------------------------------------------------------
SUBROUTINE print_tensor2_shape(ibrav)
!--------------------------------------------------------------------

USE io_global, ONLY : stdout

IMPLICIT NONE

INTEGER, INTENT(IN) :: ibrav

SELECT CASE (ibrav)
   CASE(1,2,3)  
!
!   cubic
!
      WRITE(stdout,'(/,5x, "( e11   .    .  )")')
      WRITE(stdout,'(  5x, "(  .   e11   .  )")')
      WRITE(stdout,'(  5x, "(  .    .   e11 )")')
   CASE(4,5,6,7)  
!
!  hexagonal, trigonal, tetragonal
!
      WRITE(stdout,'(/,5x, "( e11   .    .  )")')
      WRITE(stdout,'(  5x, "(  .   e11   .  )")')
      WRITE(stdout,'(  5x, "(  .    .   e33 )")')
   CASE(8,9,10,11)  
!
!  orthorhombic
!
      WRITE(stdout,'(/,5x, "( e11   .    .  )")')
      WRITE(stdout,'(  5x, "(  .   e22   .  )")')
      WRITE(stdout,'(  5x, "(  .    .   e33 )")')
   CASE(12,13)  
!
!   monoclinic unique c
!
      WRITE(stdout,'(/,5x, "( e11  e12   .  )")')
      WRITE(stdout,'(  5x, "( e12  e22   .  )")')
      WRITE(stdout,'(  5x, "(  .    .   e33 )")')
   CASE(-12,-13)  
!
!   monoclinic unique b
!
      WRITE(stdout,'(/,5x, "( e11   .   e13 )")')
      WRITE(stdout,'(  5x, "(  .   e22   .  )")')
      WRITE(stdout,'(  5x, "( e13   .   e33 )")')
   CASE(14,0)  
!
!  triclinc 
!
      WRITE(stdout,'(/,5x, "( e11  e12  e13 )")')
      WRITE(stdout,'(  5x, "( e12  e22  e23 )")')
      WRITE(stdout,'(  5x, "( e13  e23  e33 )")')
   CASE DEFAULT 
END SELECT

RETURN
END SUBROUTINE print_tensor2_shape

!--------------------------------------------------------------------
SUBROUTINE needed_tensor2(ibrav, tensor_in_use)
!--------------------------------------------------------------------

USE io_global, ONLY : stdout

IMPLICIT NONE

INTEGER, INTENT(IN) :: ibrav
LOGICAL, INTENT(INOUT) :: tensor_in_use(6)   ! use Voigt notation

tensor_in_use=.FALSE.
tensor_in_use(1)=.TRUE.
SELECT CASE (ibrav)
   CASE(1,2,3)
!
!  cubic case, only alpha_xx is printed
!
   CASE(4,5,6,7)  
!
!  hexagonal, trigonal, tetragonal
!
      tensor_in_use(3)=.TRUE.
   CASE(8,9,10,11)  
!
!  orthorhombic
!
      tensor_in_use(2)=.TRUE.
      tensor_in_use(3)=.TRUE.
   CASE(12,13)  
!
!   monoclinic unique c
!
      tensor_in_use(2)=.TRUE.
      tensor_in_use(3)=.TRUE.
      tensor_in_use(6)=.TRUE.
   CASE(-12,-13)  
!
!   monoclinic unique b
!
      tensor_in_use(2)=.TRUE.
      tensor_in_use(3)=.TRUE.
      tensor_in_use(5)=.TRUE.
   CASE(14,0)  
!
!  triclinc, all printed 
!
      tensor_in_use=.TRUE.
   CASE DEFAULT 
END SELECT

RETURN
END SUBROUTINE needed_tensor2


!--------------------------------------------------------------------
SUBROUTINE print_piezo_shape(code_group, ibrav)
!--------------------------------------------------------------------

USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: code_group, ibrav

INTEGER :: code_group_

code_group_=code_group
IF (ibrav==0) code_group_=0

SELECT CASE (code_group_) 

   CASE(2,16,18,19,20,22,23,25,27,29,32) 
       WRITE(stdout,'(/,5x,"Solid with inversion symmetry.")') 
       WRITE(stdout,'(5x,"Third-rank tensors, such as the piezoelectic&
                     & tensor, vanish.")')
   CASE(3)
!
!  C_s   Monoclinic
!
      IF (ibrav==-12.OR.ibrav==-13) THEN
         WRITE(stdout,'(/,5x,"( g11  g12  g13   .   g15   .  )")') 
         WRITE(stdout,'(  5x,"(  .    .    .   g24   .   g26 )")') 
         WRITE(stdout,'(  5x,"( g31  g32  g33   .   g35   .  )")') 
      ELSE
         WRITE(stdout,'(/,5x,"( g11  g12  g13   .    .   g16 )")') 
         WRITE(stdout,'(  5x,"( g21  g22  g23   .    .   g26 )")') 
         WRITE(stdout,'(  5x,"(  .    .    .   g16  g26   .  )")') 
      ENDIF

   CASE(4)
!
!  C_2   Monoclinic
!
      IF (ibrav==-12.OR.ibrav==-13) THEN
         WRITE(stdout,'(/,5x,"(  .    .    .   g14   .   g16 )")') 
         WRITE(stdout,'(  5x,"( g21  g22  g23   .   g25   .  )")') 
         WRITE(stdout,'(  5x,"(  .    .    .   g34   .   g36 )")') 
      ELSE
         WRITE(stdout,'(/,5x,"(  .    .    .   g14  g15   .  )")') 
         WRITE(stdout,'(  5x,"(  .    .    .   g24  g25   .  )")') 
         WRITE(stdout,'(  5x,"( g31  g32  g33   .    .   g66 )")') 
      ENDIF

   CASE(6,7)
!
!  C_4, tetragonal, C_6 hexagonal
!
      WRITE(stdout,'(/,5x,"(  .    .    .   g14  g15   .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .   g24 -g14   .  )")') 
      WRITE(stdout,'(  5x,"( g31  g31  g33   .    .    .  )")') 

   CASE(8)
!
!  D_2 (222) Orthorombic
!
      WRITE(stdout,'(/,5x,"(  .    .    .   g14   .    .  )")') 
      WRITE(stdout,'(5x,"(  .    .    .    .   g25   .  )")') 
      WRITE(stdout,'(5x,"(  .    .    .    .    .   g36 )")') 

   CASE(9)
!
! D_3  Trigonal 
!
      WRITE(stdout,'(/,5x,"( g11 -g11   .   g14   .    .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .  -g14 2g11 )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .    .    .  )")') 

   CASE(10,11)
!
! D_4  tetragonal, D_6 hexagonal
!
      WRITE(stdout,'(/,5x,"(  .    .    .   g14   .    .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .  -g14   .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .    .    .  )")') 

   CASE(12)
!
! C_2v  Orthorombic
!
      WRITE(stdout,'(/,5x,"(  .    .    .    .   g15   .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .   g24   .    .  )")') 
      WRITE(stdout,'(  5x,"( g31  g32  g33   .    .    .  )")') 

   CASE(13)
!
! C_3v  Trigonal. Assuming m perpendicular to x1
!
      WRITE(stdout,'(/,5x,"(  .    .    .    .   g15 -g21 )")') 
      WRITE(stdout,'(  5x,"( g21 -g21   .   g15   .    .  )")') 
      WRITE(stdout,'(  5x,"( g31  g31  g33   .    .    .  )")') 

   CASE(14,15)
!
! C_4v tetragonal, C_6v hexagonal
!
      WRITE(stdout,'(/,5x,"(  .    .    .    .   g15   .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .   g15   .    .  )")') 
      WRITE(stdout,'(  5x,"( g31  g31  g33   .    .    .  )")') 

   CASE(17)
!
! C_3h hexagonal
!
      WRITE(stdout,'(/,5x,"( g11 -g11   .    .    .  -g12 )")') 
      WRITE(stdout,'(  5x,"( g12 -g12   .    .    .   g11 )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .    .    .  )")') 

   CASE(21)
!
! D_3h hexagonal
!
      WRITE(stdout,'(/,5x,"(  .    .    .    .    .  -g12 )")') 
      WRITE(stdout,'(  5x,"( g12 -g12   .    .    .    .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .    .    .  )")') 

   CASE(24)
!
! D_2d tetragonal: axis 2 || x1
!
      WRITE(stdout,'(/,5x,"(  .    .    .   g14   .    .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .   g14   .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .    .   g34 )")') 

   CASE(26)
!
! S_4 tetragonal
!
      WRITE(stdout,'(/,5x,"(  .    .    .   g14  g15   .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .  -g15  g14   .  )")') 
      WRITE(stdout,'(  5x,"( g31 -g31   .    .    .   g36 )")') 

   CASE(28,30)
!
! T, T_d cubic
!
      WRITE(stdout,'(/,5x,"(  .    .    .   g14   .    .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .   g14   .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .    .   g14 )")') 

   CASE(31)
      WRITE(stdout,'(/,5x,"Solid with O symmetry. The &
                             &piezoelectic tensor vanishes")')
   CASE DEFAULT
!
!  C_1 
!
      WRITE(stdout,'(/,5x,"No symmetry or ibrav=0")')
      WRITE(stdout,'(/,5x,"( g11  g12  g13  g14  g15  g16 )")') 
      WRITE(stdout,'(  5x,"( g21  g22  g23  g24  g25  g26 )")') 
      WRITE(stdout,'(  5x,"( g31  g32  g33  g34  g35  g36 )")') 
END SELECT

RETURN
END SUBROUTINE print_piezo_shape

!--------------------------------------------------------------------
SUBROUTINE print_el_cons_shape(laue, ibrav)
!--------------------------------------------------------------------
!
!  This routine prints on output the chape of the elastic constants
!  tensor
!
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: laue, ibrav

INTEGER :: laue_

laue_=laue
IF (ibrav==0) laue_=0


SELECT CASE (laue_) 
   CASE (16)
!
!    monoclinic case, class C_2h (2/m) 
!
      IF (ibrav==-12.OR.ibrav==-13) THEN
!
!    unique axis b
!
         WRITE(stdout,'(/,5x,"( c11  c12  c13   .   c15   .  )")') 
         WRITE(stdout,'(  5x,"( c12  c22  c23   .   c25   .  )")') 
         WRITE(stdout,'(  5x,"( c13  c23  c33   .   c35   .  )")') 
         WRITE(stdout,'(  5x,"(  .    .    .   c44   .   c46 )")') 
         WRITE(stdout,'(  5x,"( c15  c25  c35   .   c55   .  )")') 
         WRITE(stdout,'(  5x,"(  .    .    .   c46   .   c66 )")') 
      ELSE
!
!   unique axis c
!
         WRITE(stdout,'(/,5x,"( c11  c12  c13   .    .   c16 )")') 
         WRITE(stdout,'(  5x,"( c12  c22  c23   .    .   c26 )")') 
         WRITE(stdout,'(  5x,"( c13  c23  c33   .    .   c36 )")') 
         WRITE(stdout,'(  5x,"(  .    .    .   c44  c45   .  )")') 
         WRITE(stdout,'(  5x,"(  .    .    .   c45  c55   .  )")') 
         WRITE(stdout,'(  5x,"( c16  c26  c36   .    .   c66 )")') 
      ENDIF 
   CASE (20)
!
!  orthorhombic D_2h (mmm)
!
      WRITE(stdout,'(/,5x,"( c11  c12  c13   .    .    .  )")') 
      WRITE(stdout,'(  5x,"( c12  c22  c23   .    .    .  )")') 
      WRITE(stdout,'(  5x,"( c13  c23  c33   .    .    .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .   c44   .    .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .   c55   .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .    .   c66 )")') 
   CASE (18)
!
!  tetragonal C_4h (4/m)
!
      WRITE(stdout,'(/,5x,"( c11  c12  c13   .    .   c16 )")') 
      WRITE(stdout,'(  5x,"( c12  c11  c13   .    .  -c16 )")') 
      WRITE(stdout,'(  5x,"( c13  c13  c33   .    .    .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .   c44   .    .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .   c44   .  )")') 
      WRITE(stdout,'(  5x,"( c16 -c16   .    .    .   c66 )")') 
   CASE (22)
!
!  tetragonal D_4h (4/mmm)
!
      WRITE(stdout,'(/,5x,"( c11  c12  c13   .    .    .  )")') 
      WRITE(stdout,'(  5x,"( c12  c11  c13   .    .    .  )")') 
      WRITE(stdout,'(  5x,"( c13  c13  c33   .    .    .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .   c44   .    .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .   c44   .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .    .   c66 )")') 
   CASE (27)
!
!  trigonal S_6 (-3)
!
      WRITE(stdout,'(/,5x,"( c11  c12  c13  c14  c15   .  )")') 
      WRITE(stdout,'(  5x,"( c12  c11  c13 -c14 -c15   .  )")') 
      WRITE(stdout,'(  5x,"( c13  c13  c33   .    .    .  )")') 
      WRITE(stdout,'(  5x,"( c14 -c14   .   c44   .  -c15 )")') 
      WRITE(stdout,'(  5x,"( c15 -c15   .    .   c44  c14 )")') 
      WRITE(stdout,'(  5x,"(  .    .    .  -c15  c14   X  )")') 
      WRITE(stdout,'(  5x,"X=(c11-c12)/2")') 
   CASE (25)
!
!  trigonal D_3d (-3m)
!
      WRITE(stdout,'(/,5x,"( c11  c12  c13  c14   .    .  )")') 
      WRITE(stdout,'(  5x,"( c12  c11  c13 -c14   .    .  )")') 
      WRITE(stdout,'(  5x,"( c13  c13  c33   .    .    .  )")') 
      WRITE(stdout,'(  5x,"( c14 -c14   .   c44   .    .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .   c44  c14 )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .   c14   X  )")') 
      WRITE(stdout,'(  5x,"X=(c11-c12)/2")') 
   CASE (19,23)
!
!  hexagonal C_6h (6/m), D_6h (6/mmm)
!
      WRITE(stdout,'(/,5x,"( c11  c12  c13   .    .    . )")') 
      WRITE(stdout,'(  5x,"( c12  c11  c13   .    .    . )")') 
      WRITE(stdout,'(  5x,"( c13  c13  c33   .    .    . )")') 
      WRITE(stdout,'(  5x,"(  .    .    .   c44   .    . )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .   c44   . )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .    .    X )")') 
      WRITE(stdout,'(  5x,"X=(c11-c12)/2")') 
   CASE (29,32)
!
!  cubic T_h (m-3), O_h (m-3m)
!
      WRITE(stdout,'(/,5x,"( c11  c12  c12   .    .    .  )")') 
      WRITE(stdout,'(  5x,"( c12  c11  c12   .    .    .  )")') 
      WRITE(stdout,'(  5x,"( c12  c12  c11   .    .    .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .   c44   .    .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .   c44   .  )")') 
      WRITE(stdout,'(  5x,"(  .    .    .    .    .   c44 )")') 
   CASE DEFAULT
!
!  Triclinic or generic
!
      IF (laue /=2) &
      WRITE(stdout,'(5x,"Laue class not programmed or ibrav=0 using C_i")') 
      WRITE(stdout,'(/,5x,"( c11  c12  c13  c14  c15  c16 )")') 
      WRITE(stdout,'(  5x,"( c12  c22  c23  c24  c25  c26 )")') 
      WRITE(stdout,'(  5x,"( c13  c23  c33  c34  c35  c36 )")') 
      WRITE(stdout,'(  5x,"( c14  c24  c34  c44  c45  c46 )")') 
      WRITE(stdout,'(  5x,"( c15  c25  c35  c45  c55  c56 )")') 
      WRITE(stdout,'(  5x,"( c16  c26  c36  c46  c56  c66 )")') 
END SELECT

RETURN
END SUBROUTINE print_el_cons_shape

!--------------------------------------------------------------------
SUBROUTINE print_b_fact_shape(b_rest)
!--------------------------------------------------------------------
!
!  This routine prints on output the shape of the B factor tensor. 
!
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: b_rest

INTEGER :: b_rest_

b_rest_=b_rest

!WRITE(stdout,'(/,5x,"case n.",i3)') b_rest_

SELECT CASE (b_rest_)

   CASE (1)
      WRITE(stdout,'(/,5x,"( b11   .   b13 )")')
      WRITE(stdout,'(  5x,"(  .   b22   .  )")')
      WRITE(stdout,'(  5x,"( b13   .   b11 )")')

   CASE (2)
      WRITE(stdout,'(/,5x,"( b11  b12   .  )")')
      WRITE(stdout,'(  5x,"( b12  b22   .  )")')
      WRITE(stdout,'(  5x,"(  .    .   b33 )")')

   CASE (3)
      WRITE(stdout,'(/,5x,"( b11   .    .  )")')
      WRITE(stdout,'(  5x,"(  .   b22  b23 )")')
      WRITE(stdout,'(  5x,"(  .   b23  b33 )")')

   CASE (4)
      WRITE(stdout,'(/,5x,"( b11   .    .  )")')
      WRITE(stdout,'(  5x,"(  .   b22   .  )")')
      WRITE(stdout,'(  5x,"(  .    .   b33 )")')

   CASE (5)
      WRITE(stdout,'(/,5x,"( b11  b12   .  )")')
      WRITE(stdout,'(  5x,"( b12  b11   .  )")')
      WRITE(stdout,'(  5x,"(  .    .   b33 )")')

   CASE (6)
      WRITE(stdout,'(/,5x,"( b11  b12  b13 )")')
      WRITE(stdout,'(  5x,"( b12  b11  b13 )")')
      WRITE(stdout,'(  5x,"( b13  b13  b33 )")')

   CASE (7)
      WRITE(stdout,'(/,5x,"( b11   b12   b13 )")')
      WRITE(stdout,'(  5x,"( b12   b11  -b13 )")')
      WRITE(stdout,'(  5x,"( b13  -b13   b33 )")')

   CASE (8)
      WRITE(stdout,'(/,5x,"( b11   .    .  )")')
      WRITE(stdout,'(  5x,"(  .   b11   .  )")')
      WRITE(stdout,'(  5x,"(  .    .   b33 )")')

   CASE (9)
      WRITE(stdout,'(/,5x,"( b11   .    .  )")')
      WRITE(stdout,'(  5x,"(  .   b22  b23 )")')
      WRITE(stdout,'(  5x,"(  .   b23  b22 )")')

   CASE (10)
      WRITE(stdout,'(/,5x,"( b11  b12  b12 )")')
      WRITE(stdout,'(  5x,"( b12  b22  b23 )")')
      WRITE(stdout,'(  5x,"( b12  b23  b22 )")')

   CASE (11)
      WRITE(stdout,'(/,5x,"(  b11  b12  -b12 )")')
      WRITE(stdout,'(  5x,"(  b12  b22   b23 )")')
      WRITE(stdout,'(  5x,"( -b12  b23   b22 )")')

   CASE (12)
      WRITE(stdout,'(/,5x,"( b11   .    .  )")')
      WRITE(stdout,'(  5x,"(  .   b22   .  )")')
      WRITE(stdout,'(  5x,"(  .    .   b22 )")')

   CASE (13)
      WRITE(stdout,'(/,5x,"( b11  b22  b13 )")')
      WRITE(stdout,'(  5x,"( b22  b22   .  )")')
      WRITE(stdout,'(  5x,"( b13   .   b33 )")')

   CASE (14)
      WRITE(stdout,'(/,5x,"( b11  b22   .  )")')
      WRITE(stdout,'(  5x,"( b22  b22   .  )")')
      WRITE(stdout,'(  5x,"(  .    .   b33 )")')

   CASE (15)
      WRITE(stdout,'(/,5x,"( b11   b22   b13 )")')
      WRITE(stdout,'(  5x,"( b22   b22  2b13 )")')
      WRITE(stdout,'(  5x,"( b13  2b13   b33 )")')

   CASE (16)
      WRITE(stdout,'(/,5x,"( b11  b11   .  )")')
      WRITE(stdout,'(  5x,"( b11  b11   .  )")')
      WRITE(stdout,'(  5x,"(  .    .   b33 )")')

   CASE (17)
      WRITE(stdout,'(/,5x,"( b11   .    .  )")')
      WRITE(stdout,'(  5x,"(  .   b11   .  )")')
      WRITE(stdout,'(  5x,"(  .    .   b11 )")')

   CASE (18)
      WRITE(stdout,'(/,5x,"( b11  b12  b12 )")')
      WRITE(stdout,'(  5x,"( b12  b11  b12 )")')
      WRITE(stdout,'(  5x,"( b12  b12  b11 )")')

   CASE DEFAULT
      WRITE(stdout,'(/,5x,"( b11  b12  b13 )")')
      WRITE(stdout,'(  5x,"( b12  b22  b23 )")')
      WRITE(stdout,'(  5x,"( b13  b23  b33 )")')
END SELECT
   
END SUBROUTINE print_b_fact_shape

END MODULE nye
