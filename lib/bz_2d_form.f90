!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE bz_2d_form
!
!  This module provides a type bz_2d that contains all the variables
!  that describe a 2d Brillouin zone and subroutines to allocate and
!  deallocate memory for these variables and to initialize them.
!  It offers the following subroutines:
!
!  allocate_2d_bz    identifies the 2d bz type and sets the number of 
!                    vertices and the number of high symmetry points
!                    in this lattice. It allocates all the arrays of
!                    a bz_2d structure.
!
!  deallocate_2d_bz  deallocate the arrays of a bz_2d structure.
!
!  init_2d_bz        for each edge sets the vertices, sets the normal
!                    for each edge, sets for each vertex the coordinates
!                    and the edges that define it. Then sets the possible
!                    high symmetry points and their labels.
!
!  find_2d_letter_coordinate receives the label of a high symmetry point
!                    and gives its coordinates.
!
!  transform_2d_label_coord receives a set of labels and transform them
!                   into a set of coordinates.
!
!  find_ibrav_2d    finds given the at and celldm of a 3d lattice, the
!                   2d ibrav and the 2d celldm that define the 2d lattice.
!
!  find_2d_bz_type  finds the 2d bz type given the 2d bravais lattice.
!
!  find_2d_axis_coordinates finds the intersections of the k_x and k_y axis
!                   with the 2d bz.
!
!  The module has also a set of privite routines that
!
!  compute_vertices  computes the coordinates of all the vertices of a
!                    2d bz.
!
!  find_2d_intersection receives the numbers of two edges and the
!                    normal to all edges. Find the vertex defined by 
!                    the two edges.
!  
!  find_2d_vertices  sets the array ivertex. For each vertex says which are
!                    the two edges that define it.
! 
!  inter_line_line   this routine find the intersection of two lines,
!                    one containing an edge and the other given with
!                    a point x0 and a direction vect.
!

    USE kinds, ONLY : DP
    USE io_global, ONLY : stdout
    IMPLICIT NONE
    PRIVATE
    SAVE

TYPE bz_2d
    INTEGER :: ind                        ! number of the bz type
    INTEGER :: nvertices                  ! The number of vertices
    REAL(DP), ALLOCATABLE :: normal(:,:)  ! The G vector normal to each edge
                                          ! in unit 2 pi / celldm(1)
    REAL(DP), ALLOCATABLE :: vertex_coord(:,:) ! coordinates of each vertex
                                          !(carthesian units 2 pi / celldm(1))
    INTEGER, ALLOCATABLE :: ivertex(:,:) ! for each vertex which edges define it
    INTEGER, ALLOCATABLE :: indsur(:,:)   ! for each edge the vertex that
                                          ! define it
    INTEGER :: xaxis, yaxis          ! the indices of the edges that
                                     ! intersect the x, and y axis
    REAL(DP) :: xi(3), yi(3)         ! the actual coordinates of intersection
    INTEGER :: nlett                 ! number of letters for which the position
                                     ! in the BZ is known
    CHARACTER(LEN=3), ALLOCATABLE :: letter_list(:) ! list of each letter
    REAL(DP), ALLOCATABLE :: letter_coord(:,:) ! coordinates of each letter

    INTEGER :: npx = 4               ! this parameter is used in the oblique
                                     ! lattice to find the reciprocal lattice
                                     ! vectors perpendicular to the edges of
                                     ! the bz
    INTEGER :: ibrav                 ! the bravais lattice index of the
                                     ! lattice that define this bz
    REAL(DP) :: at(3,3), bg(3,3)     ! direct and reciprocal primitive vectors
    REAL(DP) :: celldm_2d(3)         ! sizes of the 2d lattice
 
END TYPE

    PUBLIC bz_2d, allocate_2d_bz, deallocate_2d_bz, init_2d_bz, &
           find_2d_letter_coordinate, transform_2d_label_coord, &
           find_2d_bz_type, find_ibrav_2d, find_2d_axis_coordinates

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
   bz_struc%nlett=3
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

!--------------------------------------------------------------------------
SUBROUTINE init_2d_bz(bz_struc,celldm_2d)
!--------------------------------------------------------------------------
IMPLICIT NONE
TYPE(bz_2d), INTENT(INOUT) :: bz_struc
REAL(DP), INTENT(IN) :: celldm_2d(3)
INTEGER :: ibz, ivert, i, ind1, ind2, n1(6), n2(6)

bz_struc%letter_list(1)='gG '
bz_struc%letter_coord(:,1)=0.0_DP
ibz=bz_struc%ind


WRITE(stdout,'(/,5x,"2D direct lattice vectors for ibz:",i5)')  ibz
WRITE(stdout,'(5x,"(",2(f16.8,","),f16.8,")")') bz_struc%at(:,1)
WRITE(stdout,'(5x,"(",2(f16.8,","),f16.8,")")') bz_struc%at(:,2)
WRITE(stdout,'(5x,"(",2(f16.8,","),f16.8,")")') bz_struc%at(:,3)

WRITE(stdout,'(/,5x,"2D reciprocal lattice vectors")')  
WRITE(stdout,'(5x,"(",2(f16.8,","),f16.8,")")') bz_struc%bg(:,1)
WRITE(stdout,'(5x,"(",2(f16.8,","),f16.8,")")') bz_struc%bg(:,2)
WRITE(stdout,'(5x,"(",2(f16.8,","),f16.8,")")') bz_struc%bg(:,3)

IF (ibz ==1) THEN
!
!  oblique
!
!  first find the indices of the six G vectors that define the bz
!
   CALL find_2d_n1n2_monoclinic(n1, n2, bz_struc)
!
!  G vectors normal to each edge built from n1 and n2. If bg(:,1) and
!  bg(:,2) are within these vectors we put X and Y at one half of these, 
!  otherwise we put them at one half of two of the G vectors that define
!  the edges.
!
   ind1=1
   ind2=3
   DO i=1,6
      bz_struc%normal(:,i) = n1(i)*bz_struc%bg(:,1)+n2(i)*bz_struc%bg(:,2)
      IF (n1(i)==1.AND.n2(i)==0) ind1=i
      IF (n1(i)==0.AND.n2(i)==1) ind2=i
   ENDDO
   IF (ind1==ind2) ind2=MOD(ind1+2,6)+1
!
!  The number of vertices of each surface and its number
!
   bz_struc%indsur(:,1) = (/ 2, 6, 1 /)
   bz_struc%indsur(:,2) = (/ 2, 1, 2 /)
   bz_struc%indsur(:,3) = (/ 2, 2, 3 /)
   bz_struc%indsur(:,4) = (/ 2, 3, 4 /)
   bz_struc%indsur(:,5) = (/ 2, 4, 5 /)
   bz_struc%indsur(:,6) = (/ 2, 5, 6 /)

   CALL find_2d_vertices(bz_struc)
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' X '
   bz_struc%letter_list(3)=' Y '

   bz_struc%letter_coord(:,2) = 0.5_DP * bz_struc%normal(:,ind1) 
   bz_struc%letter_coord(:,3) = 0.5_DP * bz_struc%normal(:,ind2) 

   CALL find_2d_axis_coordinates(bz_struc)

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

   CALL find_2d_axis_coordinates(bz_struc)

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

   CALL find_2d_axis_coordinates(bz_struc)

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

   CALL find_2d_axis_coordinates(bz_struc)

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

   CALL find_2d_axis_coordinates(bz_struc)

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

!------------------------------------------------------------------------
SUBROUTINE find_2d_letter_coordinate(bz_struc, letter, xk_let)
!------------------------------------------------------------------------
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

CALL errore('find_2d_letter_coordinate','Letter not recognized '&
                                                        //TRIM(letter),1)
RETURN
END SUBROUTINE find_2d_letter_coordinate

!----------------------------------------------------------------------
SUBROUTINE find_2d_intersection( ivertex, normal, nedges, outputk) 
!----------------------------------------------------------------------
!
! This routine receives as input the number of the two edges that define
! a vertex of the BZ, the reciprocal vectors perpendicular to all the 
! edges and gives as output the coordinates of the vertex. This vertex 
! is at the intersection of two lines. Each line pass througt the middle 
! point of a reciprocal lattice vector and is perpendicular to it.
! The two reciprocal lattice vectors are those perpendicular to the two
! edges that define the vertex.
!
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nedges, ivertex(2)
REAL(DP), INTENT(IN) :: normal(3,nedges)
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
   outputk(2) = ( g2mod - a(1,2) / a(1,1) * g1mod ) * 0.5_DP / den
   outputk(1) = - a(2,1) * outputk(2) / a(1,1) + g1mod * 0.5_DP / a(1,1) 
ELSE
   outputk(2) = g1mod * 0.5_DP / a(2,1)
   IF ( ABS(a(1,2)) < 1.d-9 ) &
      CALL errore('find_2d_intersection','problem with solution',1)
   outputk(1) = ( - a(2,2) * outputk(2) + g2mod * 0.5_DP )/ a(1,2)
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

INTEGER :: ivert, iedge, i, iv

DO ivert = 1, bz_struc%nvertices
   iv=1
   DO iedge=1, bz_struc%nvertices   ! the number of edges is equal to the
                                    ! number of vertices
      DO i=2, bz_struc%indsur(1,iedge)+1
         IF (bz_struc%indsur(i,iedge) == ivert) THEN
             bz_struc%ivertex(iv,ivert)= iedge
             iv=iv+1
             IF (iv > 2) GOTO 100
             EXIT
         ENDIF
      ENDDO
   ENDDO
   CALL errore('find_2d_vertices','edge not found',ivert)
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
ELSEIF (ABS(moda1-moda2)< eps.AND.ABS(celldm(4)) < eps ) THEN
!
!  centered rectangular
!
   ibrav = 3
   celldm_2d(1) = celldm(1)
   celldm_2d(2) = celldm(2)
ELSE
!
!  oblique
!
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
!  This routine uses the labels of the 2d paths.
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

!--------------------------------------------------------------------------
SUBROUTINE inter_line_line(x0, vect, bedge, xk)
!--------------------------------------------------------------------------
!
!  This routine finds the intersection between the line passing through
!  x0 and parallel to vect, and the line passing through bedge/2 and
!  perpendicular to bedge.
!  Although x0, vect, and bedge are 3d vectors, there are supposed to be
!  all in the xy plane.
!
IMPLICIT NONE

REAL(DP), INTENT(IN) :: x0(3), vect(3), bedge(3)
REAL(DP), INTENT(OUT) :: xk(3)

REAL(DP) :: sprod, sprod1, lambda, gmod2

IF (ABS(x0(3)) > 1.D-9.OR.ABS(vect(3)) > 1.D-9 .OR. ABS(bedge(3))>1D-9) &
   CALL errore('inter_line_line','input vectors not in the xy plane',1)
xk(3)=0.0_DP

sprod = bedge(1) * vect(1) + bedge(2) * vect(2)
sprod1 = x0(1) * bedge(1) + x0(2) * bedge(2)
gmod2 = bedge(1) * bedge(1) + bedge(2) * bedge(2)

IF (ABS(sprod) < 1.D-9) &
   CALL errore('inter_line_line','input vectors are perpendicular',1)

lambda = (gmod2 * 0.5_DP - sprod1) / sprod 

xk(:)=x0(:) + lambda * vect(:)

RETURN
END SUBROUTINE inter_line_line
!
!---------------------------------------------------------------------
SUBROUTINE find_2d_axis_coordinates(bz_struc)
!---------------------------------------------------------------------
!
!  This routine checks all the edges with a normal not perpendicular 
!  to the axis and find the intersection between the axis and the line
!  that contains the edge. The edge with the shortest distance is the
!  one with the intersection between the axis and the brillouin zone 
!
IMPLICIT NONE
TYPE(bz_2d), INTENT(INOUT) :: bz_struc
REAL(DP) :: x0(3), vect(3), xi(3), xmin
INTEGER :: iedges
!
!   find the interception with the xaxis
!
x0 = 0.0_DP
vect=0.0_DP
vect(1)=1.0_DP
xmin=1.D20
DO iedges=1, bz_struc%nvertices
   IF (ABS(bz_struc%normal(1,iedges)) > 1.d-9) THEN
      CALL inter_line_line(x0, vect, bz_struc%normal(:,iedges), xi)
      IF (xi(1) > 0.0_DP .AND. xi(1) < xmin) THEN
         bz_struc%xi=xi
         bz_struc%xaxis=iedges
         xmin=xi(1)
      ENDIF
   ENDIF
ENDDO
!
!   find the intersection with the y axis
!
x0 = 0.0_DP
vect=0.0_DP
vect(2)=1.0_DP
xmin=1.D20
DO iedges=1, bz_struc%nvertices
   IF (ABS(bz_struc%normal(2,iedges)) > 1.d-9) THEN
      CALL inter_line_line(x0, vect, bz_struc%normal(:,iedges), xi)
      IF (xi(2) > 0.0_DP .AND. xi(2) < xmin) THEN
         bz_struc%yi=xi
         bz_struc%xaxis=iedges
         xmin=xi(2)
      ENDIF
   ENDIF
ENDDO

RETURN
END SUBROUTINE find_2d_axis_coordinates

!-----------------------------------------------------------------------
SUBROUTINE find_2d_n1n2_monoclinic(n1, n2, bz_struc)
!-----------------------------------------------------------------------
!
!   This routine finds the six reciprocal lattice vectors whose mid point
!   is closest to the origin and order them in order of increasing angle 
!   with the x-axis. These are the six normals to the faces of the 
!   oblique Brillouin zone.
!
USE constants, ONLY : pi
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: n1(6), n2(6)
TYPE(bz_2d), INTENT(INOUT) :: bz_struc
INTEGER :: in1, in2, ind(6), inaux(6), nfound, ivect, isub
REAL(DP) :: min_mod, max_mod, modul, vect(3), save_mod(6), save_angle(6), angle
INTEGER :: npx_
LOGICAL :: done
!
!  Search among (2*npx+1)**2 vectors. Not all cases are covered by this
!  search, but the largest part should be.
!
npx_=bz_struc%npx
min_mod=0.0_DP
nfound=0
DO in1=-npx_,npx_
   DO in2=-npx_,npx_
      IF ( in1 == 0 .AND. in2 == 0 ) CYCLE
      vect(1:3)= in1 * bz_struc%bg(1:3,1) + in2 * bz_struc%bg(1:3,2)
      modul = SQRT( vect(1) ** 2 + vect(2) ** 2 + vect(3) ** 2 ) 
      angle = ACOS( vect(1) / modul )
      IF (vect(2) < 0.0_DP ) angle= 2.0_DP * pi - angle
!
!    Check if new vector is a multiple of one already found
!
      done=.FALSE.
      DO ivect=1, nfound
         done=done.OR.(ABS(angle-save_angle(ivect))< 1.D-7)
      ENDDO
      IF ( done ) THEN
!
!    In this case the new vector is a multiple of an already found vector.
!    Substitute its modulus is lower
!
         DO ivect=1, nfound
            IF (ABS(angle-save_angle(ivect))< 1.D-7 .AND. &
                                          modul < save_mod(ivect)) THEN
               n1(ivect)=in1
               n2(ivect)=in2
               save_mod(ivect)=modul
               save_angle(ivect)=angle
            END IF
         ENDDO
!
!   ricompute the maximum modulus of the found vectors
!
         min_mod=0.0_DP
         DO ivect=1, nfound
            IF (save_mod(ivect) > min_mod) min_mod=save_mod(ivect)
         ENDDO
      ELSE
         IF ( nfound < 6 .OR. modul < min_mod ) THEN
            IF (nfound < 6 ) THEN
!
!    In this case we have not yet found six vectors. 
!
               nfound = nfound + 1
               n1(nfound)=in1
               n2(nfound)=in2
               save_mod(nfound) = modul
               save_angle(nfound) = angle
               IF (modul > min_mod) min_mod = modul
            ELSE
!
!    In this case we substitute the vector with maximum modulus with this one
!
               max_mod=0.0_DP
               isub=0
               DO ivect=1,6
                  IF (save_mod(ivect) > max_mod) THEN
                     isub = ivect
                     max_mod=save_mod(ivect)
                  ENDIF
               END DO
               IF (isub==0) CALL errore('find_2d_n1n2_monoclinic',&
                                                       'Problem with isub',1)
               n1(isub)=in1
               n2(isub)=in2
               save_mod(isub)=modul
               save_angle(isub)=angle
               min_mod=0.0_DP
               DO ivect=1, 6
                  IF (save_mod(ivect) > min_mod) min_mod=save_mod(ivect)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
   ENDDO
ENDDO
IF (nfound /= 6) CALL errore('find_2d_n1n2_monoclinic','Problem with nfound',1)
!
!  now order the six vectors, according to the angle they form with the x axis
!
!
!  If n1 or n2 is at the limit of the checked cell, tell the user to
!  double the parameter npx
!
DO ivect=1,6
   IF (n1(ivect)==npx_ .OR. n2(ivect)==npx_) &
      CALL errore('find_2d_n1n2_monoclinic','Difficult monoclinic cell, &
                                                  &double npx',1)
ENDDO
!
!  now order the six vectors, according to the angle they form with the x axis
!
ind(1)=0
CALL hpsort(6, save_angle, ind)
!
!  order n1 and n2
!
inaux=n1
n1(:) =inaux(ind(:))
inaux=n2
n2(:) =inaux(ind(:))

RETURN
END SUBROUTINE find_2d_n1n2_monoclinic
!
END MODULE bz_2d_form
