!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE bz_2d_form
    !
    USE kinds, ONLY : DP
    USE io_global, ONLY : stdout
    IMPLICIT NONE
    PUBLIC
    SAVE
    

TYPE bz_2d
    INTEGER :: ind      ! number of the bz
    INTEGER :: nvertices ! The number of vertices
    REAL(DP), ALLOCATABLE :: normal(:,:)  ! The G vector normal to each edge
                                          ! in unit 2 pi / celldm(1)
    REAL(DP), ALLOCATABLE :: vertex_coord(:,:) ! coordinates of each vertex
                                           !(carthesian units 2 pi / celldm(1))
    INTEGER, ALLOCATABLE :: ivertex(:,:) ! for each vertex which edges define it
    INTEGER, ALLOCATABLE :: indsur(:,:) ! for each edge the vertex that
                                        ! define it
    INTEGER :: xaxis, yaxis          ! the indices of the edges that
                                     ! intersect the x, and y axis
    REAL(DP) :: xi(3), yi(3)         ! the actual coordinates of intersection
    INTEGER :: nlett                 ! number of letters for which the position
                                     ! in the BZ is known
    CHARACTER(LEN=3), ALLOCATABLE :: letter_list(:) ! list of each letter
    REAL(DP), ALLOCATABLE :: letter_coord(:,:) ! coordinates of each letter

    INTEGER :: npx = 8

    INTEGER :: ibrav
    REAL(DP) :: at(3,3), bg(3,3)
    REAL(DP) :: celldm_2d(3)
 
END TYPE

CONTAINS

!----------------------------------------------------------------------
SUBROUTINE allocate_2d_bz(ibz, bz_struc, celldm, at, bg, celldm_2d)
!----------------------------------------------------------------------
!
!  This routine identifies the 2d Brillouin zone (ibz) starting from the at
!  and celldm and sets celldm_2d and the bz_struc data for the following
!  routines
!
IMPLICIT NONE
INTEGER, INTENT(OUT) :: ibz
TYPE(bz_2d), INTENT(INOUT) :: bz_struc
REAL(DP), INTENT(IN) :: celldm(6), at(3,3), bg(3,3)
REAL(DP), INTENT(OUT) ::  celldm_2d(3)
INTEGER :: ibrav_2d

CALL find_ibrav_2d(at,ibrav_2d,celldm,celldm_2d)

ibz=ibrav_2d
bz_struc%ind=ibz
bz_struc%ibrav=ibrav_2d
bz_struc%celldm_2d=celldm_2d
bz_struc%at=at
bz_struc%bg=bg

IF ( ibz ==1) THEN
!
!  oblique bz
!
   bz_struc%nvertices=6
   bz_struc%nlett=6
ELSEIF (ibz==2) THEN
!
!  rectangular bz
!
   bz_struc%nvertices=4
   bz_struc%nlett=4
ELSEIF (ibz==3) THEN
!
!  centered rectangular bz
!
   bz_struc%nvertices=6
   bz_struc%nlett=5
ELSEIF (ibz==4) THEN
!
!  square bz
!
   bz_struc%nvertices=4
   bz_struc%nlett=3
ELSEIF (ibz==5) THEN
!
!  hexagonal
!
   bz_struc%nvertices=6
   bz_struc%nlett=3
ELSE
   CALL errore('allocate_bz','Brillouin zone type not available',1)
ENDIF

ALLOCATE(bz_struc%normal(3,bz_struc%nvertices))
ALLOCATE(bz_struc%ivertex(2,bz_struc%nvertices))
ALLOCATE(bz_struc%vertex_coord(3,bz_struc%nvertices))
ALLOCATE(bz_struc%indsur(3,bz_struc%nvertices))
ALLOCATE(bz_struc%letter_list(bz_struc%nlett))
ALLOCATE(bz_struc%letter_coord(3,bz_struc%nlett))

RETURN
END SUBROUTINE allocate_2d_bz

!----------------------------------------------------------------------
SUBROUTINE deallocate_2d_bz(bz_struc)
!----------------------------------------------------------------------
!
IMPLICIT NONE
TYPE(bz_2d), INTENT(INOUT) :: bz_struc

DEALLOCATE(bz_struc%normal)
DEALLOCATE(bz_struc%ivertex)
DEALLOCATE(bz_struc%vertex_coord)
DEALLOCATE(bz_struc%indsur)
DEALLOCATE(bz_struc%letter_list)
DEALLOCATE(bz_struc%letter_coord)

RETURN
END SUBROUTINE deallocate_2d_bz

SUBROUTINE init_2d_bz(bz_struc,celldm_2d)
IMPLICIT NONE
TYPE(bz_2d), INTENT(INOUT) :: bz_struc
REAL(DP), INTENT(IN) :: celldm_2d(3)
INTEGER :: ibz, ivert

bz_struc%letter_list(1)='gG '
bz_struc%letter_coord(:,1)=0.0_DP
ibz=bz_struc%ind


WRITE(stdout,'(/,5x,"2D direct lattice vectors for ibz:",i5)')  ibz
write(stdout,'(5x,"(",2(f16.8,","),f16.8,")")') bz_struc%at(:,1)
write(stdout,'(5x,"(",2(f16.8,","),f16.8,")")') bz_struc%at(:,2)
write(stdout,'(5x,"(",2(f16.8,","),f16.8,")")') bz_struc%at(:,3)

WRITE(stdout,'(/,5x,"2D reciprocal lattice vectors")')  
write(stdout,'(5x,"(",2(f16.8,","),f16.8,")")') bz_struc%bg(:,1)
write(stdout,'(5x,"(",2(f16.8,","),f16.8,")")') bz_struc%bg(:,2)
write(stdout,'(5x,"(",2(f16.8,","),f16.8,")")') bz_struc%bg(:,3)

IF (ibz ==1) THEN
!
!  oblique
!
!
!  G vector normal to each surface
!
!   bz_struc%normal(:,1) =
!   bz_struc%normal(:,2) =
!   bz_struc%normal(:,3) =
!   bz_struc%normal(:,4) =
!   bz_struc%normal(:,5) =
!   bz_struc%normal(:,6) =
!
!  The number of vertice of each surface and its number
!
   bz_struc%indsur(:,1) = (/ 2, 6, 1 /)
   bz_struc%indsur(:,2) = (/ 2, 1, 2 /)
   bz_struc%indsur(:,3) = (/ 2, 2, 3 /)
   bz_struc%indsur(:,4) = (/ 2, 3, 4 /)
   bz_struc%indsur(:,5) = (/ 2, 4, 5 /)
   bz_struc%indsur(:,6) = (/ 2, 5, 6 /)

!   CALL find_2d_vertices(bz_struc)
!   CALL compute_vertices(bz_struc) 

!   bz_struc%letter_list(2)='   '
!   bz_struc%letter_list(3)='   '
!   bz_struc%letter_list(4)='   '

!   bz_struc%letter_coord(:,2)=
!   bz_struc%letter_coord(:,3)=
!   bz_struc%letter_coord(:,4)=
                             
   CALL errore('init_2d_bz','Oblique BZ not implemented yet',1)

!   CALL find_axis_coordinates(bz_struc)

ELSEIF (ibz==2) THEN
!
!  rectangular bz
!
   bz_struc%normal(:,1) =   bz_struc%bg(:,1)  
   bz_struc%normal(:,2) =   bz_struc%bg(:,2) 
   bz_struc%normal(:,3) = - bz_struc%bg(:,1) 
   bz_struc%normal(:,4) = - bz_struc%bg(:,2)

   bz_struc%indsur(:,1) = (/ 2, 4, 1 /)
   bz_struc%indsur(:,2) = (/ 2, 1, 2 /)
   bz_struc%indsur(:,3) = (/ 2, 2, 3 /)
   bz_struc%indsur(:,4) = (/ 2, 3, 4 /)

   CALL find_2d_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' X '
   bz_struc%letter_list(3)=' S '
   bz_struc%letter_list(4)=' Y '
!
!   standard bg
!
   bz_struc%letter_coord(:,2) = 0.5_DP * bz_struc%bg(:,1) 
   bz_struc%letter_coord(:,3) = bz_struc%vertex_coord(:,1)
   bz_struc%letter_coord(:,4) = 0.5_DP * bz_struc%bg(:,2)

!   CALL find_axis_coordinates(bz_struc)

ELSEIF (ibz==3) THEN
!
! centered rectangular
!
   bz_struc%normal(:,2) =   bz_struc%bg(:,1) 
   bz_struc%normal(:,5) = - bz_struc%bg(:,1)   
   IF ( celldm_2d(2) < 1.0_DP ) THEN
      bz_struc%normal(:,1) =   bz_struc%bg(:,1) - bz_struc%bg(:,2) 
      bz_struc%normal(:,3) =   bz_struc%bg(:,2)
      bz_struc%normal(:,4) = - bz_struc%bg(:,1) + bz_struc%bg(:,2)
      bz_struc%normal(:,6) = - bz_struc%bg(:,2)
   ELSE
      bz_struc%normal(:,1) = - bz_struc%bg(:,2) 
      bz_struc%normal(:,3) =   bz_struc%bg(:,1) + bz_struc%bg(:,2)
      bz_struc%normal(:,4) =   bz_struc%bg(:,2)
      bz_struc%normal(:,6) = - bz_struc%bg(:,1) - bz_struc%bg(:,2)
   END IF

   bz_struc%indsur(:,1) = (/ 2, 6, 1 /)
   bz_struc%indsur(:,2) = (/ 2, 1, 2 /)
   bz_struc%indsur(:,3) = (/ 2, 2, 3 /)
   bz_struc%indsur(:,4) = (/ 2, 3, 4 /)
   bz_struc%indsur(:,5) = (/ 2, 4, 5 /)
   bz_struc%indsur(:,6) = (/ 2, 5, 6 /)

   CALL find_2d_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(3)=' K '
   bz_struc%letter_list(4)=' M '
   bz_struc%letter_list(5)=' K1'
   bz_struc%letter_coord(:,3) = bz_struc%vertex_coord(:,1)
   bz_struc%letter_coord(:,4) = 0.5_DP*bz_struc%bg(:,1)
   bz_struc%letter_coord(:,5) = bz_struc%vertex_coord(:,2)
   IF (celldm_2d(2) < 1.0_DP) THEN
      bz_struc%letter_list(2)=' X '
      bz_struc%letter_coord(:,2)=0.5_DP*(bz_struc%bg(:,1)-bz_struc%bg(:,2))
   ELSE
      bz_struc%letter_list(2)=' Y '
      bz_struc%letter_coord(:,2)=0.5_DP*(bz_struc%bg(:,1)+bz_struc%bg(:,2))
   END IF

!   CALL find_axis_coordinates(bz_struc)

ELSEIF (ibz==4) THEN
!
!  square bz
!
   bz_struc%normal(:,1)=bz_struc%bg(:,1)
   bz_struc%normal(:,2)=bz_struc%bg(:,2)
   bz_struc%normal(:,3)=-bz_struc%bg(:,1)
   bz_struc%normal(:,4)=-bz_struc%bg(:,2)

   bz_struc%indsur(:,1) = (/ 2, 4, 1 /)
   bz_struc%indsur(:,2) = (/ 2, 1, 2 /)
   bz_struc%indsur(:,3) = (/ 2, 2, 3 /)
   bz_struc%indsur(:,4) = (/ 2, 3, 4 /)

   CALL find_2d_vertices(bz_struc) 

   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' M '
   bz_struc%letter_list(3)=' X '

   bz_struc%letter_coord(:,2) = 0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,2))
   bz_struc%letter_coord(:,3) = 0.5_DP*bz_struc%bg(:,1)

!   CALL find_axis_coordinates(bz_struc)

ELSEIF (ibz==5) THEN
!
!  hexagonal bz
!
   bz_struc%normal(:,1) = bz_struc%bg(:,1) - bz_struc%bg(:,2)
   bz_struc%normal(:,2) = bz_struc%bg(:,1)
   bz_struc%normal(:,3) = bz_struc%bg(:,2)
   bz_struc%normal(:,4) = bz_struc%bg(:,2) - bz_struc%bg(:,1)
   bz_struc%normal(:,5) = -bz_struc%bg(:,1) 
   bz_struc%normal(:,6) = -bz_struc%bg(:,2)

   bz_struc%indsur(:,1) = (/ 2, 6, 1 /)
   bz_struc%indsur(:,2) = (/ 2, 1, 2 /)
   bz_struc%indsur(:,3) = (/ 2, 2, 3 /)
   bz_struc%indsur(:,4) = (/ 2, 3, 4 /)
   bz_struc%indsur(:,5) = (/ 2, 4, 5 /)
   bz_struc%indsur(:,6) = (/ 2, 5, 6 /)

   CALL find_2d_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' K '
   bz_struc%letter_list(3)=' M '

   bz_struc%letter_coord(:,2) = bz_struc%vertex_coord(:,1) 
   bz_struc%letter_coord(:,3) = 0.5_DP* bz_struc%bg(:,1)

!   CALL find_axis_coordinates(bz_struc)

ELSE
   CALL errore('init_2d_bz','Brillouin zone type not available',1)
ENDIF

RETURN
END SUBROUTINE init_2d_bz

!----------------------------------------------------------------------
SUBROUTINE compute_vertices(bz_struc)
!----------------------------------------------------------------------
!
!  This routine finds the coordinates of the vertex of the BZ, given
!  the index of the three planes that define each vertex.
!
IMPLICIT NONE
TYPE(bz_2d), INTENT(INOUT) :: bz_struc
REAL(DP) :: xk(3)
INTEGER :: i

DO i = 1, bz_struc%nvertices
   CALL find_2d_intersection( bz_struc%ivertex(:,i), bz_struc%normal, &
                           bz_struc%nvertices, xk)
   bz_struc%vertex_coord(:,i)=xk(:)
ENDDO
RETURN
END SUBROUTINE compute_vertices

SUBROUTINE find_2d_letter_coordinate(bz_struc, letter, xk_let)
!
!  This routine checks if among the labels of special points defined
!  for each BZ there is the label letter and in that case it 
!  returns the coordinates of the point with that label. It stops
!  if the letter is not recognized.
!
IMPLICIT NONE
TYPE(bz_2d), INTENT(IN) :: bz_struc
REAL(DP), INTENT(OUT) :: xk_let(3)
CHARACTER(LEN=3), INTENT(IN) :: letter

INTEGER :: i

DO i=1, bz_struc%nlett
   IF ((letter(1:2) == bz_struc%letter_list(i)(2:3) .AND. &
         bz_struc%letter_list(i)(1:1)/='g') .OR. &
        (letter(1:3) == bz_struc%letter_list(i)(1:3) )) THEN
         xk_let(:) = bz_struc%letter_coord(:,i)
         RETURN
   ENDIF
ENDDO

CALL errore('find_letter_coordinate','Letter not recognized '//TRIM(letter),1)

RETURN
END SUBROUTINE find_2d_letter_coordinate

!----------------------------------------------------------------------
SUBROUTINE find_2d_intersection( ivertex, normal, nfaces, outputk) 
!----------------------------------------------------------------------
!
! This routine receives as input the number of the two lines that define
! a vertex of the BZ, the reciprocal vectors perpendicular to all the 
! lines and gives as output the intersection point.
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nfaces, ivertex(2)
REAL(DP), INTENT(IN) :: normal(3,nfaces)
REAL(DP), INTENT(OUT) :: outputk(3)

REAL(DP) :: a(2,2), g1mod, g2mod, den
INTEGER :: ipol, jpol

outputk(3)=0.0_DP

DO ipol=1,2
   DO jpol=1,2
      a(jpol,ipol) = normal(jpol,ivertex(ipol))
   ENDDO
END DO

g1mod = a(1,1) ** 2 + a(2,1) ** 2
g2mod = a(1,2) ** 2 + a(2,2) ** 2 

IF (ABS(a(1,1)) > 1.d-9) THEN
   den = ( a(2,2) * a(1,1) - a(2,1) * a(1,2) ) / a(1,1)
   outputk(2) = ( g2mod * 0.5_DP - a(1,2) / a(1,1) * g1mod * 0.5_DP ) / den
   outputk(1) = - a(2,1) * outputk(2) / a(1,1) + g1mod * 0.5_DP / a(1,1) 
ELSE
   outputk(2) = g1mod * 0.5_DP / a(2,1)
   IF ( ABS(a(1,2)) < 1.d-9 ) &
      CALL errore('find_2d_intersection','problem with solution',1)
   outputk(1) = - a(2,2) * outputk(2) / a(1,2) + g2mod * 0.5_DP / a(1,2)
END IF

RETURN
END SUBROUTINE find_2d_intersection

!----------------------------------------------------------------------
SUBROUTINE find_2d_bz_type(ibrav, ibz)
!----------------------------------------------------------------------
!
!  This routine identifies the bz type that corresponds to the given
!  bravais lattice and structural parameters. 
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: ibrav
INTEGER, INTENT(OUT) :: ibz
REAL(DP) :: value

IF (ibrav>0.AND.ibrav<6) THEN
   ibz=ibrav
ELSE
   CALL errore('find_2d_bz_type','Wrong ibrav',1)
ENDIF

RETURN
END SUBROUTINE find_2d_bz_type

!----------------------------------------------------------------------
SUBROUTINE find_2d_vertices(bz_struc) 
!----------------------------------------------------------------------
!
!  This routine uses the definition of the vertices of each edge to 
!  identify, for each vertex, the two edges that define it.
!
IMPLICIT NONE
TYPE(bz_2d), INTENT(INOUT) :: bz_struc

INTEGER :: ivert, iface, i, iv

DO ivert = 1, bz_struc%nvertices
   iv=1
   DO iface=1, bz_struc%nvertices
      DO i=2, bz_struc%indsur(1,iface)+1
         IF (bz_struc%indsur(i,iface) == ivert) THEN
             bz_struc%ivertex(iv,ivert)= iface
             iv=iv+1
             IF (iv > 2) GOTO 100
             EXIT
         ENDIF
      ENDDO
   ENDDO
   CALL errore('find_2d_vertices','face not found',ivert)
100 CONTINUE   
ENDDO

RETURN
END SUBROUTINE find_2d_vertices

!----------------------------------------------------------------------
SUBROUTINE find_ibrav_2d(at, ibrav, celldm, celldm_2d)
!----------------------------------------------------------------------
IMPLICIT NONE
REAL(DP), INTENT(IN) :: at(3,3), celldm(6)
REAL(DP) :: eps=1.d-12, prod1, prod2, moda1, moda2
INTEGER, INTENT(OUT) :: ibrav
REAL(DP), INTENT(OUT) :: celldm_2d(3)
!
!  check that at(:,3) is perpendicular to the other two vectors
!
celldm_2d=0.0_DP
prod1 = at(1,1) * at(1,3) + at(2,1) * at(2,3) + at(3,1) * at(3,3) 
prod2 = at(1,2) * at(1,3) + at(2,2) * at(2,3) + at(3,2) * at(3,3) 

IF (ABS(prod1) > eps .OR. ABS(prod2) > eps) &
     CALL errore('find_ibrav_2d','a3 must be perpendicular to a1 and a2',1)

moda1 = SQRT(at(1,1) ** 2 + at(2,1) **2 + at(3,1) ** 2)
moda2 = SQRT(at(1,2) ** 2 + at(2,2) **2 + at(3,2) ** 2)

prod1 = (at(1,1)*at(1,2)+at(2,1)*at(2,2)+at(3,1)*at(3,2)) / moda1 / moda2


IF (ABS(prod1) < eps) THEN
   IF (ABS(moda1-moda2) < eps) THEN
!
!  cubic
!
      ibrav = 4
      celldm_2d(1) = celldm(1)
   ELSE
!
!  rectangular
!
      ibrav = 2
      celldm_2d(1) = celldm(1)
      celldm_2d(2) = celldm(2)
   ENDIF
ELSEIF (ABS(prod1+0.5_DP) < eps) THEN
!
!  hexagonal
!
   ibrav = 5
   celldm_2d(1) = celldm(1)
ELSEIF (ABS(moda1-moda2)< eps ) THEN
!
!  centered rectangular
!
   ibrav = 3
   celldm_2d(1) = celldm(1)
   celldm_2d(2) = celldm(2)
ELSE
   ibrav=1
   celldm_2d(1) = celldm(1)
   celldm_2d(2) = celldm(2)
   celldm_2d(3) = celldm(4)
END IF

RETURN
END SUBROUTINE find_ibrav_2d

!----------------------------------------------------------------------
SUBROUTINE transform_2d_label_coord(ibrav, celldm, xk, letter, label_list, &
                                 npk_label, nks, k_points )
!----------------------------------------------------------------------
!
!  This routine transforms the labels in the array letter into k points
!  coordinates that are put in the array xk in the position indicated
!  by label_list. If k_point='crystal' the coordinates are tranformed
!  in the basis of the crystal. npk_label is the size of the array 
!  letter and label_list,  while nks is the size of the array xk.
!  This routine used the labels of the 2d paths.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: npk_label
INTEGER, INTENT(IN) :: nks
INTEGER, INTENT(IN) :: ibrav
INTEGER, INTENT(IN) :: label_list(npk_label)
REAL(DP), INTENT(IN) :: celldm(6)
REAL(DP), INTENT(INOUT) :: xk(3, nks)
CHARACTER(LEN=3), INTENT(IN) :: letter(npk_label)
CHARACTER(LEN=*), INTENT(IN) :: k_points
INTEGER :: bzt, i
REAL(DP) :: omega, at(3,3), bg(3,3), xk_buffer(3), celldm_2d(3)
TYPE(bz_2d) :: bz_2d_struc
!
!    generate direct lattice vectors 
!
   CALL latgen(ibrav,celldm,at(:,1),at(:,2),at(:,3),omega)
!
!   we use at in units of celldm(1)
!
   at=at/celldm(1)
!
! generate reciprocal lattice vectors
!
   CALL recips( at(:,1), at(:,2), at(:,3), bg(:,1), bg(:,2), bg(:,3) )
!
! load the information on the Brillouin zone in bz_2d_struc and identify 
! the brillouin zone type (bzt) and the celldm_2d
!
   CALL allocate_2d_bz(bzt, bz_2d_struc, celldm, at, bg, celldm_2d)
   CALL init_2d_bz(bz_2d_struc, celldm_2d)
!
! find for each label the corresponding coordinates and save them
! on the k point list
!
   DO i=1, npk_label
      CALL find_2d_letter_coordinate(bz_2d_struc, letter(i), xk_buffer)
!
!  The output of this routine is in cartesian coordinates. If the other
!  k points are in crystal coordinates we transform xk_buffer to the bg
!  base.
!
      IF (TRIM(k_points)=='crystal') CALL cryst_to_cart(1,xk_buffer,at,-1) 
      xk(:,label_list(i))=xk_buffer(:)
   ENDDO
  CALL deallocate_2d_bz(bz_2d_struc)
RETURN
END SUBROUTINE transform_2d_label_coord

END MODULE bz_2d_form
