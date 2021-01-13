!
! Copyright (C) 2013-14 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE bz_asy_mod
    !
    !   This module contains the additional variables needed to make
    !   the plot of the Brillouin zone. It defines the type bz_asy
    !   that contains an array of logical variables, one for each face,
    !   .TRUE. if the face is visible and .FALSE. if the face is not
    !   visible. Then for each possible label it contains its standard
    !   position close to the point where it is put.
    !   It offers the following subroutines:
    !
    !   allocate_bz_asy     allocate memory for a bz_asy type, taking
    !                       info from a bz_struc 
    !  
    !   deallocate_bz_asy   deallocate the memory of a bz_asy type 
    !
    !   allocate_2d_bz_asy  allocate the memory of a bz_asy type, taking 
    !                       info from a bz_2d_struc 
    !   init_bz_asy         initializes the variables of a bz_asy type using
    !                       an input bz structure
    !
    !   init_2d_bz_asy      initialize the variables of a bz_asy type using
    !                       an input bz_2d structure
    !
    !   find_letter_position receives a label and gives the position of this
    !                       label with respect to the 3d point
    !
    !   find_2d_letter_position receives a label and gives the position of this
    !                       label with respect to the 2d point
    !
    USE kinds, ONLY : DP
    USE bz_form, ONLY : bz
    USE bz_2d_form, ONLY : bz_2d
    IMPLICIT NONE
    PRIVATE
    SAVE

TYPE bz_asy
    LOGICAL, ALLOCATABLE :: visible(:) ! true if a face is visible in a plot
                                       ! of the bz
    CHARACTER(LEN=3), ALLOCATABLE :: letter_position(:) ! contains
                                       ! the position of each label in the
                                       ! plot
END TYPE
   
    PUBLIC bz_asy, allocate_bz_asy, init_bz_asy, find_letter_position, &
           deallocate_bz_asy, allocate_2d_bz_asy, init_2d_bz_asy,      &
           find_2d_letter_position

CONTAINS

!----------------------------------------------------------------------
SUBROUTINE allocate_bz_asy( bz_struc, bz_asy_struc )
!----------------------------------------------------------------------
!
IMPLICIT NONE
TYPE(bz), INTENT(INOUT) :: bz_struc
TYPE(bz_asy), INTENT(INOUT) :: bz_asy_struc
!
ALLOCATE(bz_asy_struc%visible(bz_struc%nfaces))
ALLOCATE(bz_asy_struc%letter_position(bz_struc%nlett))
!
RETURN
END SUBROUTINE allocate_bz_asy
!
!----------------------------------------------------------------------
SUBROUTINE deallocate_bz_asy( bz_asy_struc )
!----------------------------------------------------------------------
!
IMPLICIT NONE
TYPE(bz_asy), INTENT(INOUT) :: bz_asy_struc
!
IF (ALLOCATED(bz_asy_struc%visible)) DEALLOCATE(bz_asy_struc%visible)
IF (ALLOCATED(bz_asy_struc%letter_position)) DEALLOCATE(bz_asy_struc%&
                                                         letter_position)
!
RETURN
END SUBROUTINE deallocate_bz_asy

!------------------------------------------------------------------------
SUBROUTINE allocate_2d_bz_asy( bz_2d_struc, bz_asy_struc )
!------------------------------------------------------------------------
IMPLICIT NONE
TYPE(bz_2d), INTENT(INOUT) :: bz_2d_struc
TYPE(bz_asy), INTENT(INOUT) :: bz_asy_struc
!
ALLOCATE(bz_asy_struc%visible(bz_2d_struc%nvertices))
ALLOCATE(bz_asy_struc%letter_position(bz_2d_struc%nlett))
!
RETURN
END SUBROUTINE allocate_2d_bz_asy

!----------------------------------------------------------------------
SUBROUTINE init_bz_asy(bz_struc, bz_asy_struc)
!----------------------------------------------------------------------
!
USE asy, ONLY : asy_proj
USE bz_form, ONLY : adjust_orthorombic_vect
IMPLICIT NONE
TYPE(bz), INTENT(IN) :: bz_struc
TYPE(bz_asy), INTENT(INOUT) :: bz_asy_struc
INTEGER :: ibz, iface
REAL(DP) :: prod, proj(3)

bz_asy_struc%visible(:)=.TRUE.
bz_asy_struc%letter_position(:)='N'
ibz=bz_struc%ind
proj=asy_proj
!
! Determine which faces are visible in the plot, projecting on the 
! direction of the orthogonal projection.
!
! in the orthorhombic case the bz can be constructed in rotated axis, we
! need to rotate also the projection direction to find the visible faces
!
IF (ibz>7.AND.ibz<13) &
   CALL adjust_orthorombic_vect(bz_struc,proj)

DO iface=1, bz_struc%nfaces
   prod = proj(1) * bz_struc%normal(1,iface) + &
          proj(2) * bz_struc%normal(2,iface) + &
          proj(3) * bz_struc%normal(3,iface)  
   bz_asy_struc%visible(iface) = (prod >= 0.0_DP)
END DO 

IF ( ibz ==1) THEN
!
!  simple cubic bz
!
!
!  This tell the code where to write the letter in the plot
!
   bz_asy_struc%letter_position(1)='S'
   bz_asy_struc%letter_position(2)='SE'
   bz_asy_struc%letter_position(3)='NE'
   bz_asy_struc%letter_position(4)='N'

   IF (bz_struc%letter_type=='BI') THEN
      bz_asy_struc%letter_position(5)='NW'
   ENDIF

ELSEIF (ibz==2) THEN
!
!  fcc bz
!

   bz_asy_struc%letter_position(1)='S'
   bz_asy_struc%letter_position(2)='S'
   bz_asy_struc%letter_position(3)='SE'
   bz_asy_struc%letter_position(4)='NE'
   bz_asy_struc%letter_position(5)='NE'
   bz_asy_struc%letter_position(6)='N'

   IF (bz_struc%letter_type=='BI') THEN
      bz_asy_struc%letter_position(7) ='NW'
      bz_asy_struc%letter_position(8) ='SE'
      bz_asy_struc%letter_position(9) =' U1'
      bz_asy_struc%letter_position(10) ='NW'
      bz_asy_struc%letter_position(11) ='NW'
      bz_asy_struc%letter_position(12) ='SE'
      bz_asy_struc%letter_position(13) ='SE'
   ENDIF

ELSEIF (ibz==3) THEN
!
!  bcc bz
!

   bz_asy_struc%letter_position(2)='SE'
   bz_asy_struc%letter_position(3)='N'
   bz_asy_struc%letter_position(4)='NE'

   IF (bz_struc%letter_type=='BI') THEN
      bz_asy_struc%letter_position(5)='N'
   ENDIF

ELSEIF (ibz==4) THEN
!
!  simple tetragonal bz
!

   bz_asy_struc%letter_position(1)='S'
   bz_asy_struc%letter_position(2)='SE'
   bz_asy_struc%letter_position(3)='SE'
   bz_asy_struc%letter_position(4)='NE'
   bz_asy_struc%letter_position(5)='N'
   bz_asy_struc%letter_position(6)='N'

ELSEIF (ibz==5) THEN
!
!  centered tetragonal (c<a) bz
!
   bz_asy_struc%letter_position(1)='S'
   bz_asy_struc%letter_position(2)='SE'
   bz_asy_struc%letter_position(3)='SE'
   bz_asy_struc%letter_position(4)='NW'
   bz_asy_struc%letter_position(5)='NW'
   bz_asy_struc%letter_position(6)='NW'
   bz_asy_struc%letter_position(7)='NW'

ELSEIF (ibz==6) THEN
!
!  centered tetragonal (c>a) bz
!
   bz_asy_struc%letter_position(1)='S'
   bz_asy_struc%letter_position(2)='SE'
   bz_asy_struc%letter_position(3)='NW'
   bz_asy_struc%letter_position(4)='NW'
   bz_asy_struc%letter_position(5)='NW'
   bz_asy_struc%letter_position(6)='N'
   bz_asy_struc%letter_position(7)='NE'
   bz_asy_struc%letter_position(8)='SE'
   bz_asy_struc%letter_position(9)='S'
   IF (bz_struc%letter_type=='BI') THEN
      bz_asy_struc%letter_position(2)='NW'
      bz_asy_struc%letter_position(4)='N'
      bz_asy_struc%letter_position(5)='NE'
      bz_asy_struc%letter_position(6)='N'
      bz_asy_struc%letter_position(9)='N'
   ENDIF

   IF (bz_struc%letter_type=='BI') THEN
      bz_asy_struc%letter_position(10)='NW'
      bz_asy_struc%letter_position(11)='NW'
      bz_asy_struc%letter_position(12)='SE'
      bz_asy_struc%letter_position(13)='SE'
      bz_asy_struc%letter_position(14)='NE'
      bz_asy_struc%letter_position(15)='NW'
   ENDIF

ELSEIF (ibz==7) THEN
!
!  simple orthorhombic bz
!
   bz_asy_struc%letter_position(1)='S'
   bz_asy_struc%letter_position(2)='NW'
   bz_asy_struc%letter_position(3)='SE'
   bz_asy_struc%letter_position(4)='SE'
   bz_asy_struc%letter_position(5)='N'
   bz_asy_struc%letter_position(6)='N'
   bz_asy_struc%letter_position(7)='N'
   bz_asy_struc%letter_position(8)='NW'

ELSEIF (ibz==8) THEN
!
!  face centered orthorhombic (1/a^2 > 1/b^2 + 1/c^2) bz
!
   bz_asy_struc%letter_position(1)='S'
   bz_asy_struc%letter_position(2)='NW'
   bz_asy_struc%letter_position(3)='NW'
   bz_asy_struc%letter_position(4)='NW'
   bz_asy_struc%letter_position(5)='NE'
   bz_asy_struc%letter_position(6)='SW'
   bz_asy_struc%letter_position(7)='SE'
   bz_asy_struc%letter_position(8)='NE'
   bz_asy_struc%letter_position(9)='NE'

   IF (bz_struc%switch_b_c) THEN
      bz_asy_struc%letter_position(2)='SE'
      bz_asy_struc%letter_position(4)='NE'
      bz_asy_struc%letter_position(6)='NW'
      bz_asy_struc%letter_position(7)='NW'
      bz_asy_struc%letter_position(9)='NW'
   ENDIF
   IF ( bz_struc%switch_a_b ) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_asy_struc%letter_position(2)='NE'
         bz_asy_struc%letter_position(3)='NW'
         bz_asy_struc%letter_position(4)='NE'
         bz_asy_struc%letter_position(7)='NW'
         bz_asy_struc%letter_position(8)='SE'
         bz_asy_struc%letter_position(9)='NW'
      ELSE
         bz_asy_struc%letter_position(2)='NE'
         bz_asy_struc%letter_position(3)='NE'
         bz_asy_struc%letter_position(8)='NW'
         bz_asy_struc%letter_position(9)='NW'
      END IF
   ELSEIF (bz_struc%rotate_a_b_c) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_asy_struc%letter_position(2)='NW'
         bz_asy_struc%letter_position(3)='NW'
         bz_asy_struc%letter_position(4)='NW'
         bz_asy_struc%letter_position(7)='NE'
         bz_asy_struc%letter_position(8)='SE'
         bz_asy_struc%letter_position(9)='NE'
      ELSE
         bz_asy_struc%letter_position(2)='SE'
         bz_asy_struc%letter_position(3)='NE'
         bz_asy_struc%letter_position(4)='NW'
         bz_asy_struc%letter_position(7)='NE'
         bz_asy_struc%letter_position(8)='NW'
         bz_asy_struc%letter_position(9)='NW'
      END IF
   ENDIF

ELSEIF (ibz==9) THEN
!
!  face centered orthorhombic (1/a^2 < 1/b^2 + 1/c^2) bz
!
!
   bz_asy_struc%letter_position(1)='S'
   bz_asy_struc%letter_position(2)='NW'
   bz_asy_struc%letter_position(3)='SE'
   bz_asy_struc%letter_position(4)='SE'
   bz_asy_struc%letter_position(5)='NE'
   bz_asy_struc%letter_position(6)='N'
   bz_asy_struc%letter_position(7)='N'
   bz_asy_struc%letter_position(8)='NE'
   bz_asy_struc%letter_position(9)='N'
   bz_asy_struc%letter_position(10)='N'
   bz_asy_struc%letter_position(11)='NW'

   IF (bz_struc%switch_b_c) THEN
      bz_asy_struc%letter_position(3)='N'
      bz_asy_struc%letter_position(4)='N'
      bz_asy_struc%letter_position(5)='NW'
      bz_asy_struc%letter_position(6)='SE'
      bz_asy_struc%letter_position(8)='N'
      bz_asy_struc%letter_position(9)='SE'
      bz_asy_struc%letter_position(10)='NE'
      bz_asy_struc%letter_position(11)='NE'
   ENDIF

   IF (bz_struc%switch_a_b) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_asy_struc%letter_position(2)='NE'
         bz_asy_struc%letter_position(3)='NW'
         bz_asy_struc%letter_position(4)='NW'
         bz_asy_struc%letter_position(5)='NW'
         bz_asy_struc%letter_position(6)='NE'
         bz_asy_struc%letter_position(8)='SE'
         bz_asy_struc%letter_position(9)='NE'
         bz_asy_struc%letter_position(10)='SE'
         bz_asy_struc%letter_position(11)='NE'
      ELSE
         bz_asy_struc%letter_position(2)='NE'
         bz_asy_struc%letter_position(3)='SE'
         bz_asy_struc%letter_position(4)='SE'
         bz_asy_struc%letter_position(5)='NW'
         bz_asy_struc%letter_position(6)='NE'
         bz_asy_struc%letter_position(8)='NW'
         bz_asy_struc%letter_position(9)='NE'
         bz_asy_struc%letter_position(10)='NW'
      ENDIF
   ELSEIF (bz_struc%rotate_a_b_c) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_asy_struc%letter_position(2)='NW'
         bz_asy_struc%letter_position(3)='NE'
         bz_asy_struc%letter_position(4)='NE'
         bz_asy_struc%letter_position(5)='NE'
         bz_asy_struc%letter_position(6)='NW'
         bz_asy_struc%letter_position(8)='SE'
         bz_asy_struc%letter_position(9)='NW'
         bz_asy_struc%letter_position(10)='SE'
         bz_asy_struc%letter_position(11)='NW'
      ELSE
         bz_asy_struc%letter_position(2)='NE'
         bz_asy_struc%letter_position(3)='NE'
         bz_asy_struc%letter_position(4)='NE'
         bz_asy_struc%letter_position(6)='SE'
         bz_asy_struc%letter_position(5)='NW'
         bz_asy_struc%letter_position(8)='NW'
         bz_asy_struc%letter_position(9)='SE'
         bz_asy_struc%letter_position(10)='NW'
         bz_asy_struc%letter_position(11)='NW'
      ENDIF
   ENDIF

ELSEIF (ibz==10) THEN
!
!  face centered orthorhombic (1/a^2 = 1/b^2 + 1/c^2) bz
!

   bz_asy_struc%letter_position(1)='S'
   bz_asy_struc%letter_position(2)='NW'
   bz_asy_struc%letter_position(3)='NW'
   bz_asy_struc%letter_position(4)='N'
   bz_asy_struc%letter_position(5)='NW'
   bz_asy_struc%letter_position(6)='NE'
   bz_asy_struc%letter_position(7)='SE'
   bz_asy_struc%letter_position(8)='NE'

   IF (bz_struc%switch_b_c) THEN
      bz_asy_struc%letter_position(3)='SE'
      bz_asy_struc%letter_position(5)='NE'
      bz_asy_struc%letter_position(7)='NW'
      bz_asy_struc%letter_position(8)='NW'
   ENDIF

   IF (bz_struc%switch_a_b) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_asy_struc%letter_position(2)='NW'
         bz_asy_struc%letter_position(3)='NE'
         bz_asy_struc%letter_position(5)='NE'
         bz_asy_struc%letter_position(6)='SE'
         bz_asy_struc%letter_position(7)='NW'
         bz_asy_struc%letter_position(8)='NW'
      ELSE
         bz_asy_struc%letter_position(3)='NE'
         bz_asy_struc%letter_position(6)='NW'
         bz_asy_struc%letter_position(8)='NW'
         bz_asy_struc%letter_position(2)='NE'
      ENDIF
   ELSEIF (bz_struc%rotate_a_b_c) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_asy_struc%letter_position(2)='NW'
         bz_asy_struc%letter_position(3)='NW'
         bz_asy_struc%letter_position(5)='NW'
         bz_asy_struc%letter_position(6)='SE'
         bz_asy_struc%letter_position(7)='NE'
         bz_asy_struc%letter_position(8)='NE'
      ELSE
         bz_asy_struc%letter_position(2)='NE'
         bz_asy_struc%letter_position(3)='SE'
         bz_asy_struc%letter_position(5)='NW'
         bz_asy_struc%letter_position(6)='NW'
         bz_asy_struc%letter_position(7)='NE'
         bz_asy_struc%letter_position(8)='NW'
      ENDIF
   ENDIF

ELSEIF (ibz==11) THEN
!
!  body centered orthorhombic bz
!
!   bz_asy_struc%visible(2)=.FALSE.
!   bz_asy_struc%visible(3)=.FALSE.
!   bz_asy_struc%visible(7)=.FALSE.
!   bz_asy_struc%visible(8)=.FALSE.
!   bz_asy_struc%visible(11)=.FALSE.
!   bz_asy_struc%visible(12)=.FALSE.
!   bz_asy_struc%visible(14)=.FALSE.

!   IF (bz_struc%switch_b_c) THEN
!      bz_asy_struc%visible=.TRUE.
!      bz_asy_struc%visible(2)=.FALSE.
!      bz_asy_struc%visible(3)=.FALSE.
!      bz_asy_struc%visible(7)=.FALSE.
!      bz_asy_struc%visible(10)=.FALSE.
!      bz_asy_struc%visible(11)=.FALSE.
!      bz_asy_struc%visible(12)=.FALSE.
!      bz_asy_struc%visible(14)=.FALSE.
!   ENDIF
!   IF (bz_struc%switch_a_b) THEN
!      IF (bz_struc%switch_b_c) THEN
!         bz_asy_struc%visible=.TRUE.
!         bz_asy_struc%visible(3)=.FALSE.
!         bz_asy_struc%visible(4)=.FALSE.
!         bz_asy_struc%visible(8)=.FALSE.
!         bz_asy_struc%visible(9)=.FALSE.
!         bz_asy_struc%visible(11)=.FALSE.
!         bz_asy_struc%visible(12)=.FALSE.
!         bz_asy_struc%visible(14)=.FALSE.
!      ELSE
!         bz_asy_struc%visible=.TRUE.
!         bz_asy_struc%visible(3)=.FALSE.
!         bz_asy_struc%visible(4)=.FALSE.
!         bz_asy_struc%visible(7)=.FALSE.
!         bz_asy_struc%visible(8)=.FALSE.
!         bz_asy_struc%visible(11)=.FALSE.
!         bz_asy_struc%visible(12)=.FALSE.
!         bz_asy_struc%visible(14)=.FALSE.
!      ENDIF
!   ELSEIF (bz_struc%rotate_a_b_c) THEN
!      IF (bz_struc%switch_b_c) THEN
!         bz_asy_struc%visible=.TRUE.
!         bz_asy_struc%visible(3)=.FALSE.
!         bz_asy_struc%visible(4)=.FALSE.
!         bz_asy_struc%visible(9)=.FALSE.
!         bz_asy_struc%visible(10)=.FALSE.
!         bz_asy_struc%visible(11)=.FALSE.
!         bz_asy_struc%visible(12)=.FALSE.
!         bz_asy_struc%visible(14)=.FALSE.
!      ELSE
!         bz_asy_struc%visible=.TRUE.
!         bz_asy_struc%visible(2)=.FALSE.
!         bz_asy_struc%visible(3)=.FALSE.
!         bz_asy_struc%visible(9)=.FALSE.
!         bz_asy_struc%visible(10)=.FALSE.
!         bz_asy_struc%visible(11)=.FALSE.
!         bz_asy_struc%visible(12)=.FALSE.
!         bz_asy_struc%visible(14)=.FALSE.
!      ENDIF
!   ENDIF


   bz_asy_struc%letter_position(1)='S'
   bz_asy_struc%letter_position(2)='NW'
   bz_asy_struc%letter_position(3)='SW'
   bz_asy_struc%letter_position(4)='SE'
   bz_asy_struc%letter_position(5)='SE'
   bz_asy_struc%letter_position(6)='NE'
   bz_asy_struc%letter_position(7)='NW'
   bz_asy_struc%letter_position(8)='N'
   bz_asy_struc%letter_position(9)='NE'
   bz_asy_struc%letter_position(10)='NW'
   bz_asy_struc%letter_position(11)='E'
   bz_asy_struc%letter_position(12)='NE'
   bz_asy_struc%letter_position(13)='NW'

   IF (bz_struc%letter_type=='BI') THEN
      bz_asy_struc%letter_position(2)='gS0'
      bz_asy_struc%letter_position(4)=' R '
      bz_asy_struc%letter_position(6)='gL0'
      bz_asy_struc%letter_position(7)=' T '
      bz_asy_struc%letter_position(10)=' F0'
      bz_asy_struc%letter_position(12)=' G0'
      bz_asy_struc%letter_position(13)=' G '
   ENDIF

   IF (bz_struc%switch_b_c) THEN
      bz_asy_struc%letter_position(3)='W'
      bz_asy_struc%letter_position(4)='NW'
      bz_asy_struc%letter_position(5)='NW'
      bz_asy_struc%letter_position(7)='SE'
      bz_asy_struc%letter_position(8)='N'
      bz_asy_struc%letter_position(9)='NE'
      bz_asy_struc%letter_position(10)='SE'
      bz_asy_struc%letter_position(11)='N'
      bz_asy_struc%letter_position(12)='NE'
      bz_asy_struc%letter_position(6)='NW' 
      bz_asy_struc%letter_position(13)='SE'
      IF (bz_struc%letter_type=='BI') THEN
         bz_asy_struc%letter_position(6)='gL0'
         bz_asy_struc%letter_position(13)=' G '
      ENDIF
   ENDIF

   IF (bz_struc%switch_a_b) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_asy_struc%letter_position(2)='NE'
         bz_asy_struc%letter_position(3)='N'
         bz_asy_struc%letter_position(4)='NW'
         bz_asy_struc%letter_position(5)='SW'
         bz_asy_struc%letter_position(6)='NW'
         bz_asy_struc%letter_position(7)='NE'
         bz_asy_struc%letter_position(9)='SE'
         bz_asy_struc%letter_position(10)='NE'
         bz_asy_struc%letter_position(11)='SE'
         bz_asy_struc%letter_position(12)='SE'
         bz_asy_struc%letter_position(13)='NE'
         IF (bz_struc%letter_type=='BI') THEN
            bz_asy_struc%letter_position(2)='gL0'
            bz_asy_struc%letter_position(6)='gS0'
            bz_asy_struc%letter_position(13)=' G '
         ENDIF
      ELSE
         bz_asy_struc%letter_position(2)='NE'
         bz_asy_struc%letter_position(3)='SE'
         bz_asy_struc%letter_position(5)='SW'
         bz_asy_struc%letter_position(6)='NW'
         bz_asy_struc%letter_position(13)='NW'
         bz_asy_struc%letter_position(7)='NE'
         bz_asy_struc%letter_position(9)='NW'
         bz_asy_struc%letter_position(10)='NE'
         bz_asy_struc%letter_position(11)='N'
         bz_asy_struc%letter_position(12)='NW'
         IF (bz_struc%letter_type=='BI') THEN
            bz_asy_struc%letter_position(2)=' G '
            bz_asy_struc%letter_position(6)='gS0'
            bz_asy_struc%letter_position(13)='gL0'
         ENDIF
      END IF
   ELSEIF (bz_struc%rotate_a_b_c) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_asy_struc%letter_position(3)='NE'
         bz_asy_struc%letter_position(4)='NE'
         bz_asy_struc%letter_position(5)='NE'
         bz_asy_struc%letter_position(7)='NW'
         bz_asy_struc%letter_position(8)='N'
         bz_asy_struc%letter_position(9)='SE'
         bz_asy_struc%letter_position(10)='NW'
         bz_asy_struc%letter_position(11)='NW'
         bz_asy_struc%letter_position(12)='SE'
         bz_asy_struc%letter_position(2)='NW'
         bz_asy_struc%letter_position(6)='NE'
         bz_asy_struc%letter_position(13)='NW'
         IF (bz_struc%letter_type=='BI') THEN
            bz_asy_struc%letter_position(2)='gL0'
            bz_asy_struc%letter_position(6)=' G '
            bz_asy_struc%letter_position(13)='gS0'
         ENDIF
      ELSE
         bz_asy_struc%letter_position(2)='SE'
         bz_asy_struc%letter_position(3)='NE'
         bz_asy_struc%letter_position(4)='NE'
         bz_asy_struc%letter_position(5)='NE'
         bz_asy_struc%letter_position(6)='NW'
         bz_asy_struc%letter_position(7)='SE'
         bz_asy_struc%letter_position(9)='NW'
         bz_asy_struc%letter_position(10)='SE'
         bz_asy_struc%letter_position(11)='NW'
         bz_asy_struc%letter_position(12)='NW'
         bz_asy_struc%letter_position(13)='NW'
         IF (bz_struc%letter_type=='BI') THEN
            bz_asy_struc%letter_position(2)=' G '
            bz_asy_struc%letter_position(6)='gL0'
            bz_asy_struc%letter_position(13)='gS0'
         ENDIF
      ENDIF
   END IF

ELSEIF (ibz==12) THEN
!
!  one face centered orthorhombic bz
!
!   bz_asy_struc%visible(3)=.FALSE.
!   bz_asy_struc%visible(4)=.FALSE.
!   bz_asy_struc%visible(5)=.FALSE.
!   bz_asy_struc%visible(8)=.FALSE.
!   IF (bz_struc%switch_a_b) THEN
!      bz_asy_struc%visible=.TRUE.
!      bz_asy_struc%visible(4)=.FALSE.
!      bz_asy_struc%visible(5)=.FALSE.
!      bz_asy_struc%visible(6)=.FALSE.
!      bz_asy_struc%visible(8)=.FALSE.
!   ENDIF


   bz_asy_struc%letter_position(1)='S'
   bz_asy_struc%letter_position(2)='NW'
   bz_asy_struc%letter_position(3)='SE'
   bz_asy_struc%letter_position(4)='SE'
   bz_asy_struc%letter_position(5)='NE'
   bz_asy_struc%letter_position(6)='NE'
   bz_asy_struc%letter_position(7)='NW'
   bz_asy_struc%letter_position(8)='NE'
   bz_asy_struc%letter_position(9)='SW'
   bz_asy_struc%letter_position(10)='NW'

   IF (bz_struc%switch_a_b) THEN
      bz_asy_struc%letter_position(2)='NE'
      bz_asy_struc%letter_position(3)='SE'
      bz_asy_struc%letter_position(4)='SE'
      bz_asy_struc%letter_position(5)='NW'
      bz_asy_struc%letter_position(6)='NW'
      bz_asy_struc%letter_position(7)='N'
      bz_asy_struc%letter_position(8)='NE'
      bz_asy_struc%letter_position(9)='NE'
      bz_asy_struc%letter_position(10)='NW'
   ENDIF

   IF (bz_struc%letter_type=='BI') THEN
      bz_asy_struc%letter_position(2)='gD0'
      bz_asy_struc%letter_position(9)=' B0'
   ENDIF

ELSEIF (ibz==13) THEN
!
!  hexagonal
!
!  bz_asy_struc%visible(3)=.FALSE.
!   bz_asy_struc%visible(4)=.FALSE.
!   bz_asy_struc%visible(5)=.FALSE.
!   bz_asy_struc%visible(8)=.FALSE.

   bz_asy_struc%letter_position(1)='S'
   bz_asy_struc%letter_position(2)='SE'
   bz_asy_struc%letter_position(3)='SE'
   bz_asy_struc%letter_position(4)='NW'
   bz_asy_struc%letter_position(5)='NE'
   bz_asy_struc%letter_position(6)='NW'

!
ELSEIF (ibz==14) THEN
!
!  trigonal alpha < 90 bz
!

!   bz_asy_struc%visible(3)=.FALSE.
!   bz_asy_struc%visible(4)=.FALSE.
!   bz_asy_struc%visible(5)=.FALSE.
!   bz_asy_struc%visible(9)=.FALSE.
!   bz_asy_struc%visible(10)=.FALSE.
!   bz_asy_struc%visible(11)=.FALSE.
!   bz_asy_struc%visible(14)=.FALSE.

   bz_asy_struc%letter_position(1)='S'
   bz_asy_struc%letter_position(2)='SE'
   bz_asy_struc%letter_position(3)='NW'
   bz_asy_struc%letter_position(4)='S'
   bz_asy_struc%letter_position(5)='E'
   bz_asy_struc%letter_position(6)='E'
   bz_asy_struc%letter_position(7)='NW'
   bz_asy_struc%letter_position(8)='SW'
   bz_asy_struc%letter_position(9)='NE'
   bz_asy_struc%letter_position(10)='NW'
   bz_asy_struc%letter_position(11)='SE'
   bz_asy_struc%letter_position(12)='SE'

ELSEIF (ibz==15) THEN
!
!  trigonal alpha > 90 bz
!
!   bz_asy_struc%visible(3)=.FALSE.
!   bz_asy_struc%visible(4)=.FALSE.
!   bz_asy_struc%visible(5)=.FALSE.
!   bz_asy_struc%visible(9)=.FALSE.
!   bz_asy_struc%visible(11)=.FALSE.
!   bz_asy_struc%visible(12)=.FALSE.

   bz_asy_struc%letter_position(1)='S'
   bz_asy_struc%letter_position(2)='SE'
   bz_asy_struc%letter_position(3)='NE'
   bz_asy_struc%letter_position(4)='SW'
   bz_asy_struc%letter_position(5)='NW'
   bz_asy_struc%letter_position(6)='NE'
   bz_asy_struc%letter_position(7)='NE'
   bz_asy_struc%letter_position(8)='NW'

   IF (bz_struc%letter_type=='BI') THEN
      bz_asy_struc%letter_position(3)=' R0'
      bz_asy_struc%letter_position(5)=' T '
      bz_asy_struc%letter_position(6)=' FA'
      bz_asy_struc%letter_position(8)=' P2'
   ENDIF

ELSEIF (ibz==16) THEN
!
!  simple monoclinic bz
!
   bz_asy_struc%letter_position(1)='SE'
   bz_asy_struc%letter_position(2)='NW'
   bz_asy_struc%letter_position(3)='NE'
   bz_asy_struc%letter_position(4)='NE'
   bz_asy_struc%letter_position(5)='NE'
   bz_asy_struc%letter_position(6)='NE'
ELSE
   CALL errore('init_2d_bz_asy', &
                       ' 2d Brillouin zone type not available',1)
ENDIF

RETURN
END SUBROUTINE init_bz_asy

!----------------------------------------------------------------------
SUBROUTINE init_2d_bz_asy(bz_2d_struc, bz_asy_struc, ibrav_2d, celldm_2d)
!----------------------------------------------------------------------
!
IMPLICIT NONE
TYPE(bz_2d), INTENT(IN) :: bz_2d_struc
TYPE(bz_asy), INTENT(INOUT) :: bz_asy_struc
INTEGER, INTENT(IN) :: ibrav_2d
REAL(DP), INTENT(IN) :: celldm_2d(3)

bz_asy_struc%visible(:)=.TRUE.
bz_asy_struc%letter_position(:)='N'
!
!
IF (ibrav_2d==1) THEN
!
!  oblique lattice
!
   bz_asy_struc%letter_position(1)='S '
   bz_asy_struc%letter_position(2)='  '
   bz_asy_struc%letter_position(3)='  '
!
ELSEIF (ibrav_2d==2) THEN
!
!  rectangular lattice
!
   bz_asy_struc%letter_position(1)='S '
   bz_asy_struc%letter_position(2)='SE'
   bz_asy_struc%letter_position(3)='NE'
   bz_asy_struc%letter_position(4)='N'
!
ELSEIF (ibrav_2d==3) THEN
!
!  centered rectangular lattice
!
   bz_asy_struc%letter_position(1)='S '
   IF (celldm_2d(2) < 1.0_DP) THEN
      bz_asy_struc%letter_position(2)='SE'
      bz_asy_struc%letter_position(3)='SE'
      bz_asy_struc%letter_position(4)='NE'
      bz_asy_struc%letter_position(5)='NE'
   ELSE
      bz_asy_struc%letter_position(2)='N'
      bz_asy_struc%letter_position(3)='SE'
      bz_asy_struc%letter_position(4)='NE'
      bz_asy_struc%letter_position(5)='NE'
   END IF

ELSEIF (ibrav_2d==4) THEN
!
!  square lattice
!
   bz_asy_struc%letter_position(1)='S '
   bz_asy_struc%letter_position(2)='NE'
   bz_asy_struc%letter_position(3)='SE'

ELSEIF (ibrav_2d==5) THEN
!
!  hexagonal lattice
!
   bz_asy_struc%letter_position(1)='S '
   bz_asy_struc%letter_position(2)='SE'
   bz_asy_struc%letter_position(3)='NE'

ELSE
   CALL errore('init_2d_bz_asy', &
                       ' 2d Brillouin zone type not available',1)
ENDIF

RETURN
END SUBROUTINE init_2d_bz_asy

!----------------------------------------------------------------------
SUBROUTINE find_letter_position(bz_struc, bz_asy_struc, letter, let_pos)
!----------------------------------------------------------------------
!
!  This routine checks if among the labels of special points defined
!  for each BZ there is the label letter and in that case it 
!  returns the position of that label on the asy plot. It stops
!  if the letter is not recognized.
!
IMPLICIT NONE
TYPE(bz), INTENT(IN) :: bz_struc
TYPE(bz_asy), INTENT(IN) :: bz_asy_struc
CHARACTER(LEN=3), INTENT(OUT) :: let_pos
CHARACTER(LEN=3), INTENT(IN) :: letter

INTEGER :: i

DO i=1, bz_struc%nlett
   IF ((letter(1:2) == bz_struc%letter_list(i)(2:3) .AND. &
         bz_struc%letter_list(i)(1:1)/='g') .OR. &
        (letter(1:3) == bz_struc%letter_list(i)(1:3) )) THEN
         let_pos = bz_asy_struc%letter_position(i)
         RETURN
   ENDIF
ENDDO

CALL errore('find_letter_position','Letter not recognized '//TRIM(letter),1)

RETURN
END SUBROUTINE find_letter_position

!----------------------------------------------------------------------
SUBROUTINE find_2d_letter_position(bz_2d_struc, bz_asy_struc, letter, let_pos)
!----------------------------------------------------------------------
!
!  This routine checks if among the labels of special points defined
!  for each BZ there is the label letter and in that case it 
!  returns the position of that label on the asy plot. It stops
!  if the letter is not recognized.
!
IMPLICIT NONE
TYPE(bz_2d), INTENT(IN) :: bz_2d_struc
TYPE(bz_asy), INTENT(IN) :: bz_asy_struc
CHARACTER(LEN=3), INTENT(OUT) :: let_pos
CHARACTER(LEN=3), INTENT(IN) :: letter

INTEGER :: i

DO i=1, bz_2d_struc%nlett
   IF ((letter(1:2) == bz_2d_struc%letter_list(i)(2:3) .AND. &
         bz_2d_struc%letter_list(i)(1:1)/='g') .OR. &
        (letter(1:3) == bz_2d_struc%letter_list(i)(1:3) )) THEN
         let_pos = bz_asy_struc%letter_position(i)
         RETURN
   ENDIF
ENDDO

CALL errore('find_2d_letter_position','Letter not recognized '//TRIM(letter),1)

RETURN
END SUBROUTINE find_2d_letter_position

END MODULE bz_asy_mod
