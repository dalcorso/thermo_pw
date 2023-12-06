!
! Copyright (C) 2015 - 2017 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE lattices
!
!  This module contains variables and routines to deal with Bravais
!  lattices. 
!
!  Presently it has the routines: 
!
!  compute_conventional  : this routine receives the direct lattice
!                          vectors of a centered cell and computes the
!                          vectors of the corresponding conventional cell.
!                          atc in output are in the same units of the at
!                          in input. They are in cartesian coordinates.
!                          This routine must be matched to the routine
!                          latgen because it depends on the choice of the
!                          direct lattice vectors.
!
!  conventional_ibrav   : given the ibrav, sets the ibrav of the conventional
!                         lattice
!
!  find_ibrav_code        : this routine receives the direct lattice
!                          vectors and finds the code number (1-14) of the
!                          bravais lattice. The code correspond to the
!                          ibrav code of QE (see the routine Modules/latgen)
!
!
!  find_combination      : if two set of primitive Bravais lattices 
!                          describe the same lattice express the first
!                          as a linear combination of the second.
!
!  is_bravais_lattice    : is a function that returns .TRUE. if a vector
!                          is a Bravais lattice vector.
!
!  same_lattice          : is a logical function that receives two sets of
!                          primitive vectors and gives .TRUE. if the two sets 
!                          give the same Bravais lattice. It generalizes find
!                          combination because the two Bravais lattices
!                          are considered the same even if they are related
!                          by a global rotation. The routine provides the
!                          global rotation and the integer matrix that
!                          links the rotated vectors of the second set
!                          to those of the first. The two object separately
!                          are not unique and the routine gives one possible
!                          rotation and decomposition.
!
!  lattice_point_group   : is a routine that finds the code (1-32) of the
!                          point group of the Bravais lattice given by three
!                          principal vectors. It is more general than the
!                          routine provided by QE because it allow any
!                          arbitrary orientation of the rotation axis.
!
!  bravais_dir :  is a function that receive the primitive lattice vectors
!                 and a direction and gives .TRUE. if there is a bravais
!                 lattice vector parallel to that direction within the
!                 shell defined by -3 < n_i < 3.. In output
!                 it gives also the lattice vector of shortest modulus and
!                 its modulus.
!
!  compute_omega: receives three primitive lattice vectors and gives as
!                 output the volume of the unit cell
!
!  is_centered : receives ibrav and gives true if the lattice is centered
!
!  lattice_name : receives the bravais lattice and gives the lattice name
!
!  zone_border : receives a vector in direct or reciprocal space, 
!                the direct and reciprocal lattice vectors and gives 
!                .true. if the input vector is at zone border
!
!  is_compatible_group_ibrav : a function that check if a point group 
!                is compatible with a given Bravais lattice and gives also 
!                all the list of compatible lattices
!  
!  crystal_parameters : receives the Bravais lattice index and give the
!                       number of required independent crystal parameters
!
!  same_star : receives two k vectors and the list of point group symmetries.
!              It gives .TRUE. if the two k vectors belong to the same star.
!
!  print_bravais_description : receives the celldm of a given lattice 
!              and writes a small description of the crystal parameters.
!
!  compress_celldm : receives celldm and ibrav and produces a vector with
!                  the variable crystal parameters for that ibrav
!
!  expand_celldm : receives a vector with the variable crystal parameters
!                  for a given ibrav and sets the celldm vector
!
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  SAVE

  CHARACTER(LEN=7) :: &
           celldm_name(6)

  DATA celldm_name / 'a=     ', 'b/a=   ', 'c/a=   ', &
                     'cos(a)=', 'cos(b)=', 'cos(c)=' /

  CHARACTER(LEN=16) :: &
           celldm_gnuplot_name(6)

  DATA celldm_gnuplot_name / 'a (a.u.)        ', 'b/a             ',   &
                             'c/a             ',                       &
                             'cos({/Symbol a})', 'cos({/Symbol b})',   &
                             'cos({/Symbol c})' /
                     
  PUBLIC compute_conventional, find_ibrav_code, find_combination,      &
         is_bravais_lattice, same_lattice, lattice_point_group,        &
         compute_omega, conventional_ibrav, is_centered, lattice_name, &
         zone_border, same_star, is_compatible_group_ibrav,            &
         bravais_dir, print_bravais_description, crystal_parameters,   &
         compress_celldm, expand_celldm, compress_int_vect, needed_celldm, &
         celldm_name, celldm_gnuplot_name

CONTAINS

!----------------------------------------------------------------------
  SUBROUTINE compute_conventional(at, atp, ibrav)
!----------------------------------------------------------------------
!
!   This routine write the conventional vectors of a given centered
!   bravais lattice. It assumes that the primitive vectors are defined
!   as in the latgen routine of QE.
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ibrav
  REAL(DP), INTENT(IN) :: at(3,3)
  REAL(DP), INTENT(INOUT) :: atp(3,3)

  SELECT CASE (ibrav)
     CASE(2)
         atp(:,1)= -at(:,1) + at(:,2) - at(:,3)
         atp(:,2)= -at(:,1) + at(:,2) + at(:,3)
         atp(:,3)=  at(:,1) + at(:,2) - at(:,3)
     CASE(3)
         atp(:,1)=  at(:,1) - at(:,2) 
         atp(:,2)=  at(:,2) - at(:,3)
         atp(:,3)=  at(:,1) + at(:,3)
     CASE(5)
         atp(:,1)= at(:,1) - at(:,3)
         atp(:,2)= - at(:,1) + at(:,2)
         atp(:,3)= -at(:,1) - at(:,2) - at(:,3)
     CASE(7)
         atp(:,1)=  at(:,1) - at(:,3) 
         atp(:,2)= -at(:,1) + at(:,2)
         atp(:,3)=  at(:,2) + at(:,3)
     CASE(9)
         atp(:,1)=  at(:,1) - at(:,2) 
         atp(:,2)=  at(:,1) + at(:,2)
         atp(:,3)=  at(:,3)
     CASE(91)
         atp(:,1)=  at(:,1) 
         atp(:,2)=  at(:,2) + at(:,3)
         atp(:,3)= -at(:,2) + at(:,3)
     CASE(10)
         atp(:,1)=  at(:,1) + at(:,2) - at(:,3)
         atp(:,2)= -at(:,1) + at(:,2) + at(:,3)
         atp(:,3)=  at(:,1) - at(:,2) + at(:,3)
     CASE(11)
         atp(:,1)=  at(:,1) - at(:,2) 
         atp(:,2)=  at(:,2) - at(:,3)
         atp(:,3)=  at(:,1) + at(:,3)
     CASE(13)
         atp(:,1)=  at(:,1) + at(:,3) 
         atp(:,2)=  at(:,2)
         atp(:,3)= -at(:,1) + at(:,3)
     CASE(-13)
         atp(:,1)=  at(:,1) + at(:,2) 
         atp(:,2)= -at(:,1) + at(:,2)
         atp(:,3)=  at(:,3)
     CASE DEFAULT
!
!  If this is not a known centered lattice we simply copy the at in atp
!
         atp=at
  END SELECT   
  RETURN
  END SUBROUTINE compute_conventional

!----------------------------------------------------------------------
  SUBROUTINE conventional_ibrav(ibrav, cibrav)
!----------------------------------------------------------------------
!
! Given the ibrav of a given Bravais lattice gives the ibrav of the
! conventional lattice.
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ibrav
  INTEGER, INTENT(OUT) :: cibrav
  INTEGER :: c_ibrav(14)
  
  DATA c_ibrav / 1, 1, 1, 4, 4, 6, 6, 8, 8, 8, 8, 12, 12, 14 /

  IF (ibrav>0.AND.ibrav<15) THEN
     cibrav= c_ibrav(ibrav)
  ELSEIF (ibrav==91) THEN
     cibrav=8
  ELSEIF (ibrav==-12.OR.ibrav==-13) THEN
     cibrav=-12
  ENDIF

  RETURN
  END SUBROUTINE conventional_ibrav

!----------------------------------------------------------------------
  LOGICAL FUNCTION is_centered(ibrav)
!----------------------------------------------------------------------
!
!  Gives true if the input bravais lattice is centered
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ibrav

  is_centered = (ibrav==2.OR.ibrav==3.OR.ibrav==5.OR.ibrav==7.OR.ibrav==9&
                 .OR.ibrav==91.OR.ibrav==10.OR.ibrav==11.OR.ibrav==13.OR. &
                 ibrav==-13)

  RETURN
  END FUNCTION is_centered


!----------------------------------------------------------------------
  SUBROUTINE find_ibrav_code(a1, a2, a3, ibrav, celldm, code_group_ext, ur, &
                             global_s, verbosity)
!----------------------------------------------------------------------
!
!   This routine is closely related to the routine latgen.f90 of the
!   QE distribution. Given ibrav and celldm the routine latgen
!   provides three primitive vectors at2 of the Bravais lattice. 
!   The present routine receives three primitive vectors of a Bravais lattice
!   a1, a2, a3 and finds the Bravais lattice type (ibrav) and 
!   the dimensions of 
!   the lattice celldm that give the same Bravais lattice when given
!   as input to latgen. In general the orientation of the a1, a2, a3 can
!   be different from that given by latgen and also it might happen that
!   the a1, a2, a3 are linear combinations of the vectors given by latgen.
!   This routine gives as output two matrices, ur and global_s. These
!   two matrices allows to transform the at2 given by latgen into the 
!   a1, a2, a3 given in input.
!   Note also that in some cases the set at1, at2, at3 and the matrices
!   ur and global_s are not uniquely defined. In particular 
!   for the ortorhombic systems any exchange of at1, at2, at3 can be
!   recovered with an appropriate rotation. In this case we make the following
!   choice. If the input conventional cell is parallel to x, y, and z we 
!   keep the order of the input a1, a2, a3, so that the global rotation is
!   the identity. If the input conventional cell is rotated, we choose
!   the shorter celldm along x, then along y, and finally along z so that
!   a < b < c. These rules might be broken for ibrav=9 since latgen gives
!   a lattice in which the C face is centered, while the input lattice 
!   might have a different face centered. In this case the a, b, c are 
!   chosen so that the C face is centered.
!
!   a1, a2, a3 -> at
!   ibrav, celldm -> at2
!   at1_i = \sum_j ur_ij at2_j
!   at_i  = global_s at1_i
!   at, at1, at2 and global_s are in cartesian coordinates.
!   ur is a matrix of integer numbers.
!
!   Given the a1, a2, a3 the routine finds the lenghts and the angles 
!   between the three vectors.
!   Then it finds the point group of the resulting Bravais lattice.
!   From the point group symmetry elements we can find a global rotation
!   that brings the main symmetry elements in the standard orientation.
!   At that point we can find the conventional cell and therefore
!   the ibrav and celldm. These variables are used to generate the at2
!
!   If the variable code_group_ext is set to zero the routine finds the
!   code of the point group, otherwise it will search the Bravais lattice
!   only among those compatible with the input code_group_ext, even if 
!   the actual at might be invariant for a larger point group. In any case 
!   code_group_ext must be a subgroup of the actual point group of the at.
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE rotate, ONLY : is_rotation, find_rotation
  USE point_group, ONLY : find_group_info_ext, group_index_from_ext, &
                          is_subgroup, point_group_bravais

  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a1(3), a2(3), a3(3)
  LOGICAL, INTENT(IN)  :: verbosity
  INTEGER, INTENT(OUT) :: ibrav
  INTEGER, INTENT(INOUT) :: code_group_ext
  REAL(DP), INTENT(OUT) :: celldm(6), global_s(3,3), ur(3,3)

  REAL(DP), PARAMETER :: eps1=1.D-7, eps2=5.D-7
  REAL(DP) :: tmod(3), cangle(3), sp(3), e(3,3), emod(3), at(3,3), sr(3,3,48),&
              ax(3,3), cc(3,3), at1(3,3), at2(3,3), bg(3,3), ang, &
              modulus, omega, angle_rot, ps, e1(3,3), rot_mat(3,3), lm, &
              celldm0(6)
  INTEGER :: ipol, jpol, ivec, code_group, nsym, isym, ts, tipo_sym, iperp, &
             ipar1, ipar2, enne(3), ind(3), ncount, code_group_, &
             which_elem(48), group_desc(48), code_group_ext_, mm(3,3), idet, &
             bgenne(3)
  LOGICAL :: is_axis, parallel_axis, bunique
  CHARACTER(LEN=11) :: gname, gname_at, group_name
!
! find the point group
!
  at(:,1)=a1(:)
  at(:,2)=a2(:)
  at(:,3)=a3(:)
  CALL lattice_point_group(at,gname_at,code_group,nsym,sr)

  IF (code_group_ext==0) THEN
     code_group_ext=0
     gname=gname_at
  ELSE
     CALL find_group_info_ext(nsym, sr, code_group, code_group_ext_, &
                                                which_elem, group_desc)
     IF (.NOT. is_subgroup(code_group_ext_, code_group_ext)) &
        CALL errore('find_ibrav_code','input point group code wrong',1)
     CALL point_group_bravais(group_index_from_ext(code_group_ext), code_group)
     gname=group_name(code_group)
  ENDIF

  IF (verbosity) THEN
     WRITE(stdout,'(/,5x," find ibrav for a1, a2, a3",/)') 

     WRITE(stdout,'(3f20.12)') a1(:)
     WRITE(stdout,'(3f20.12)') a2(:)
     WRITE(stdout,'(3f20.12)') a3(:)
  END IF

  WRITE(stdout,'(/,5x,"Point group of the at ",a11,/)') gname_at
  IF (gname /= gname_at) &
     WRITE(stdout,'(5x,"Searching the ibrav among the lattices &
                   &compatible with ",a11)') gname
!
!  set the three cartesian directions
!
  cc=0.0_DP
  cc(1,1)=1.0_DP
  cc(2,2)=1.0_DP
  cc(3,3)=1.0_DP
!
!  Find the modulus and the angles of the input vectors
!
  tmod(1) = SQRT ( a1(1)**2 + a1(2)**2 + a1(3)**2 )
  tmod(2) = SQRT ( a2(1)**2 + a2(2)**2 + a2(3)**2 )
  tmod(3) = SQRT ( a3(1)**2 + a3(2)**2 + a3(3)**2 )
  IF ( tmod(1) < eps1 .OR. tmod(2) < eps1 .OR. tmod(3) < eps1 ) &
     CALL errore('find_ibrav_code','a primitive vector has zero module',1)

  sp(1) = a1(1)*a2(1) + a1(2)*a2(2) + a1(3)*a2(3)
  sp(2) = a1(1)*a3(1) + a1(2)*a3(2) + a1(3)*a3(3)  
  sp(3) = a2(1)*a3(1) + a2(2)*a3(2) + a2(3)*a3(3)  
 
  cangle(1) = sp(1) / tmod(1) / tmod(2)
  cangle(2) = sp(2) / tmod(1) / tmod(3)
  cangle(3) = sp(3) / tmod(2) / tmod(3)

!
!  In order to find a possible global rotation of the input vector with
!  respect to the standard orientation we use the rotation axis of the 
!  point group operations.
!
  celldm=0.0_DP
  SELECT CASE(code_group)
     CASE(2)
!
!  triclinic lattice  C_i
! 
       ibrav=14
       celldm(1)=tmod(1)
       celldm(2)=tmod(2)/tmod(1)
       celldm(3)=tmod(3)/tmod(1)
       celldm(4)=cangle(3)
       celldm(5)=cangle(2)
       celldm(6)=cangle(1)

       CALL latgen(ibrav,celldm,at2(1,1),at2(1,2),at2(1,3),omega)

       IF (.NOT. same_lattice(at, at2, ur, at1, global_s)) &
          CALL errore('find_ibrav_code','problem with triclinic lattice',1)

     CASE(16)
!
!  monoclinic lattice  C_2h
!
!  the two-fold rotation axis or the perpendicular to the mirror will be the 
!  z axis. One of them is among the first three symmetry operations.
!
        DO isym=1,3
           ts=tipo_sym(sr(1,1,isym))
           IF (ts==4) THEN
              CALL versor(sr(1,1,isym), ax(1,3))
              IF (.NOT. bravais_dir(at, ax(1,3), enne, e(1,3), emod(3) ) ) &
                 CALL errore('find_ibrav_code','problem with monoclinic dir',1)
              EXIT
           ELSEIF (ts==5) THEN
              CALL mirror_axis(sr(1,1,isym), ax(1,3))
              IF (.NOT. bravais_dir(at, ax(1,3), enne, e(1,3), emod(3) ) ) &
                 CALL errore('find_ibrav_code','problem with monoclinic dir',1)
              EXIT
           ENDIF
        ENDDO
        bunique=is_axis(ax(1,3),2)
!
!   mm contains the matrix needed to write our vectors in terms of
!   the input at vectors
!
        mm(:,3)=enne(:)

        CALL recips(at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))
!
!   find the shortest G vector parallel to the symmetry axis
!
        IF (.NOT. bravais_dir(bg, ax(:,3), bgenne, e(1,1), emod(1) ) ) &
               CALL errore('find_ibrav_code','No g vector found ',1)
!
!   and find the two bravais lattice vectors closest to the origin
!   in the plane perpendicular to the G vector. Together with the Bravais
!   lattice vector parallel to the axis they form our first guess of the 
!   monoclinic primitive vectors.
!
        CALL find_closest_at( bgenne, at, mm, 2, e )

        idet = mm(1,1) * ( mm(2,2) * mm(3,3) - mm(3,2) * mm(2,3) )-   &
               mm(1,2) * ( mm(2,1) * mm(3,3) - mm(3,1) * mm(2,3) )+   &
               mm(1,3) * ( mm(2,1) * mm(3,2) - mm(3,1) * mm(2,2) )
!
!   If the determinant is negative our guessed vector have left hand
!   orientation. Exchange e(:,1) and e(:,2).
!
        IF (idet < 0) THEN
           e1(:,1)=e(:,1)
           e(:,1)=e(:,2)
           e(:,2)=e1(:,1)
           idet=-idet
        ENDIF
!
!    If idet=1 the lattice is simple monoclinic and we have found it
!
        IF (idet==1) THEN
           ibrav=12
           e1(:,:)=e(:,:)
!
!   If idet==2 the lattice is centered and we must see if it is
!   B, C or I. In all cases we transform it to B.
!
        ELSEIF (idet==2) THEN
           ibrav=13
!
!   first check if it is B
!
           e1(:,1)=0.5_DP*(e(:,1) - e(:,3))
           e1(:,2)=e(:,2) 
           e1(:,3)=0.5_DP*(e(:,1) + e(:,3))
           IF (same_lattice(e1, at, ur, at1, global_s)) GOTO 100
!
!   otherwise check if it is C
!
           e1(:,1)=0.5_DP*(e(:,2) + e(:,3))
           e1(:,2)=e(:,1) 
           e1(:,3)=0.5_DP*(e(:,2) - e(:,3))
           IF (same_lattice(e1, at, ur, at1, global_s)) GOTO 100
!
!   If the codes arrive here it must be I
!
           e1(:,1)=0.5_DP*(e(:,1) +e(:,2) + e(:,3))
           e1(:,2)=e(:,1) 
           e1(:,3)=0.5_DP*(e(:,1) + e(:,2) - e(:,3))

           IF (same_lattice(e1, at, ur, at1, global_s)) GOTO 100
!
!   If the code arrives here, there is some problem. It should not happen.
!
           CALL errore('find_ibrav_code','problem1 with monoclinic',1)
100        CONTINUE
!
!    return to the conventional vectors to extract the size of the unit cell
!
           ax(:,1)=e1(:,1)
           e1(:,1)=e1(:,1) + e1(:,3)
           e1(:,3)=ax(:,1) - e1(:,3)
        ENDIF
!
!   Find the moduli of the vectors
!
        emod(1)=SQRT(e1(1,1)**2+e1(2,1)**2+e1(3,1)**2)
        emod(2)=SQRT(e1(1,2)**2+e1(2,2)**2+e1(3,2)**2)
        emod(3)=SQRT(e1(1,3)**2+e1(2,3)**2+e1(3,3)**2)
!
!  And set the dimensions of the unit cell. b-unique cells are
!  used when the symmetry axis of the original lattice was parallel 
!  to the y direction. In all other cases the lattice is cunique
!
        IF (bunique) THEN
           ibrav=-ibrav
           celldm(1)=emod(1)
           celldm(2)=emod(3)/celldm(1)
           celldm(3)=emod(2)/celldm(1)
           celldm(5)= (e1(1,1)*e1(1,2)+e1(2,1)*e1(2,2)+e1(3,1)*e1(3,2)) / &
                                            emod(1) / emod(2)
        ELSE
           celldm(1)=emod(1)
           celldm(2)=emod(2)/celldm(1)
           celldm(3)=emod(3)/celldm(1)
           celldm(4)= (e1(1,1)*e1(1,2)+e1(2,1)*e1(2,2)+e1(3,1)*e1(3,2)) / &
                                            emod(1) / emod(2)
        ENDIF
!
!   Now check that we found the correct lattice.
!
        IF (ibrav==0) &
           CALL errore('find_ibrav_code','No lattice found',1)
        CALL latgen(ibrav,celldm,at2(1,1),at2(1,2),at2(1,3),omega)

        IF (.NOT. same_lattice(at, at2, ur, at1, global_s)) &
           CALL errore('find_ibrav_code','problem with monoclinic lattice',1)

     CASE(20)
!
!  orthorhombic case D_2h
!
        ipol=0
        DO isym=1,nsym
           ts=tipo_sym(sr(1,1,isym))
           IF (ts==4) THEN
              ipol=ipol+1
              CALL versor(sr(1,1,isym), ax(1,ipol))
              IF (.NOT. bravais_dir(at, ax(1,ipol), enne, e(1,ipol), &
                         emod(ipol) ) ) &
                 CALL errore('find_ibrav_code','problem with orthorhombic dir',1)
           ENDIF
        ENDDO
!
!   if possible take the axis in the order x, y, z
!
        parallel_axis=.FALSE. 
        ncount=0
        DO ipol=1,3
           ind(ipol)=ipol
           DO jpol=1,3
              IF (is_axis(ax(1,jpol),ipol)) THEN
                 ind(ipol)=jpol
                 ncount=ncount+1
              ENDIF
           ENDDO
           IF (ncount==3) parallel_axis=.TRUE.
        ENDDO
        IF (.NOT. parallel_axis) THEN
!
!  choose an arbitrary orientation:  a < b < c 
!
           ind(1)=0
           CALL hpsort(3, emod, ind)
!
!    the emod are already ordered
!
           DO ipol=1,3
              ind(ipol)=ipol
           ENDDO
        ENDIF
        
        celldm(1)=emod(ind(1))
        celldm(2)=emod(ind(2))/emod(ind(1))
        celldm(3)=emod(ind(3))/emod(ind(1))

        ibrav=8
        CALL latgen(ibrav,celldm,at2(1,1),at2(1,2),at2(1,3),omega)
        IF (.NOT.same_lattice(at, at2, ur, at1, global_s)) THEN
           ibrav=9
           CALL latgen(ibrav,celldm,at2(1,1),at2(1,2),at2(1,3),omega)
           IF (.NOT.same_lattice(at, at2, ur, at1, global_s)) THEN
              ibrav=10
              CALL latgen(ibrav,celldm,at2(1,1),at2(1,2),at2(1,3),omega)
              IF (.NOT.same_lattice(at, at2, ur, at1, global_s)) THEN
                 ibrav=11
                 CALL latgen(ibrav,celldm,at2(1,1),at2(1,2),at2(1,3),omega)
                 IF (.NOT.same_lattice(at, at2, ur, at1, global_s)) THEN
!
!    we have still the possibility that the centered faces do not coincide
!    with the previous choice of the axis and we must check the other 
!    two possibilities.
!
                    ibrav=9
                    celldm0(:)=celldm(:)
                    celldm(1)=celldm0(2)*celldm0(1)
                    celldm(2)=celldm0(1)/celldm(1)
                    celldm(3)=celldm0(3)*celldm0(1)/celldm(1)
                    CALL latgen(ibrav,celldm,at2(1,1),at2(1,2),at2(1,3),omega)
                    IF (.NOT.same_lattice(at, at2, ur, at1, global_s)) THEN
                       celldm(1)=celldm0(3)*celldm0(1)
                       celldm(2)=celldm0(2)*celldm0(1)/celldm(1)
                       celldm(3)=celldm0(1)/celldm(1)
                       CALL latgen(ibrav,celldm,at2(1,1),at2(1,2),at2(1,3),&
                                                                        omega)
                       IF (.NOT.same_lattice(at, at2, ur, at1, global_s)) &
                          CALL errore('find_ibrav_code','orthorhombic lattice &
                                                    &not found',1)
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
         
     CASE(22)
!
!  tetragonal case D_4h
!
!  first find the direction of the z axis 
!
        DO isym=1,nsym
           ts=tipo_sym(sr(1,1,isym))
           IF (ts==3) THEN
              IF (ABS(angle_rot(sr(1,1,isym))-90.0_DP)<eps1 ) THEN
                 CALL versor(sr(1,1,isym), ax(1,3))
                 IF (.NOT.bravais_dir(at, ax(1,3), enne, e(1,3), emod(3) )) &
                    CALL errore('find_ibrav_code','problem with tetragonal dir',1)
                 EXIT
              ENDIF
           ENDIF
        ENDDO
!
!  Find the two perpendicular twofold axis with the shorter Bravais lattice 
!  vectors
!
        modulus=1.D8
        DO isym=1, nsym
           ts=tipo_sym(sr(1,1,isym))
           IF (ts==4) THEN
              CALL versor(sr(1,1,isym), ax(1,1))
              ps=ax(1,1)*ax(1,3)+ax(2,1)*ax(2,3)+ax(3,1)*ax(3,3)
!
!   neglect the two-fold rotation with axis parallel to the four-fold axis
!
              IF (ABS(ps)>eps1) CYCLE
              IF (.NOT. bravais_dir(at, ax, enne, e(1,1), emod(1) ) ) &
                 CALL errore('find_ibrav_code','problem with two-fold axis',1)
              IF (emod(1) < modulus) modulus=emod(1)
           ENDIF
        ENDDO

        celldm(1) = modulus
        celldm(3) = emod(3)/modulus

        ibrav=6
        CALL latgen(ibrav,celldm,at2(1,1),at2(1,2),at2(1,3),omega)
        IF (.NOT.same_lattice(at, at2, ur, at1, global_s)) THEN
           ibrav=7
           CALL latgen(ibrav,celldm,at2(1,1),at2(1,2),at2(1,3),omega)
           IF (.NOT.same_lattice(at, at2, ur, at1, global_s)) THEN
              CALL errore('find_ibrav_code','tetragonal lattice not found',1)
           END IF
        END IF
!
     CASE(23,25)
!
!  hexagonal case  D_6h or trigonal D_3d
!
        IF (ABS(tmod(1)-tmod(2))>eps1.OR.ABS(tmod(1)-tmod(3))>eps1) THEN
           DO isym=1,nsym
              ts=tipo_sym(sr(1,1,isym))
              IF (ts==3) THEN
                 IF (ABS(angle_rot(sr(1,1,isym))-60.0_DP)<eps1 ) THEN
                    CALL versor(sr(1,1,isym), ax(1,3))
                 ENDIF
              ENDIF
           ENDDO
!
!   find an at paralell and one perpendicular to ax(1,3)
!
!
!   hexagonal lattice
!       
           ibrav=4
           celldm(1)=1.D20
           DO ipol=1,3
              ps=ABS(at(1,ipol)*ax(1,3)+at(2,ipol)*ax(2,3)+at(3,ipol)*ax(3,3))
              IF (ABS(ps-tmod(ipol))<eps1) celldm(3)=tmod(ipol)
              IF (ABS(ps)<eps1.AND.tmod(ipol)<celldm(1)) celldm(1)=tmod(ipol)
           ENDDO
           celldm(3)=celldm(3)/celldm(1)

           CALL latgen(ibrav,celldm,at2(1,1),at2(1,2),at2(1,3),omega)
           IF (.NOT.same_lattice(at, at2, ur, at1, global_s)) &
              CALL errore('find_ibrav_code','hexagonal lattice not found',1)
        ELSE
!
!
!  trigonal case  D_3d
!
!   set the standard orientation of the three axis
!
           ibrav=5
           celldm(1)=tmod(1)
           celldm(4)=cangle(1)

           CALL latgen(ibrav,celldm,at2(1,1),at2(1,2),at2(1,3),omega)
           IF (.NOT.same_lattice(at, at2, ur, at1, global_s)) &
              CALL errore('find_ibrav_code','trigonal lattice not found',1)
        ENDIF
     CASE(32)
!
!   cubic case O_h. The three fourfold axis give the global orientation
!   of the lattice
!
        ipol=0
        DO isym=1,nsym
           ts=tipo_sym(sr(1,1,isym))
           IF (ts==3) THEN
              ang=angle_rot(sr(1,1,isym))
              IF (ABS(ang-90.0_DP)< eps1) THEN
                 ipol=ipol+1
                 CALL versor(sr(1,1,isym), ax(1,ipol))
                 IF (.NOT.bravais_dir(at, ax(1,ipol), enne, e(1,ipol), &
                                                   emod(ipol) ) ) &
                    CALL errore('find_ibrav_code','problem with cubic dir',1)
                    
              ENDIF
           ENDIF
        ENDDO

!        IF (ABS(emod(1)-emod(2))> eps1.OR. ABS(emod(1)-emod(3)) > eps1) &
!           CALL errore('find_ibrav_code','wrong cubic lattice',1)

        celldm(1)=emod(1)

        ibrav=1
        CALL latgen(ibrav,celldm,at2(1,1),at2(1,2),at2(1,3),omega)
        IF (.NOT.same_lattice(at, at2, ur, at1, global_s)) THEN
           ibrav=2
           CALL latgen(ibrav,celldm,at2(1,1),at2(1,2),at2(1,3),omega)
           IF (.NOT.same_lattice(at, at2, ur, at1, global_s)) THEN
              ibrav=3
              CALL latgen(ibrav,celldm,at2(1,1),at2(1,2),at2(1,3),omega)
              IF (.NOT.same_lattice(at, at2, ur, at1, global_s)) &
                 CALL errore('find_ibrav_code','unknown cubic lattice',1)
           ENDIF
        ENDIF

  CASE DEFAULT
     CALL errore('find_ibrav_code','wrong code_group',1)
  END SELECT

  RETURN
  END SUBROUTINE find_ibrav_code

!----------------------------------------------------------------------
  SUBROUTINE toint(b)
!----------------------------------------------------------------------
!
!  This routine receives a vector with three components. It divides
!  all the components for the one with smallest nonzero modulus.
!
  USE kinds, ONLY : DP
  IMPLICIT NONE 

  REAL(DP), INTENT(INOUT) :: b(3)

  REAL(DP), PARAMETER :: eps1=1.D-8
  REAL(DP) :: minb
  INTEGER :: ipol

  minb=1.D8
  DO ipol=1,3
     IF (ABS(b(ipol))>eps1.AND. ABS(b(ipol)) < minb ) minb= ABS(b(ipol))
  ENDDO
  b=b/minb

  RETURN
  END SUBROUTINE toint

!----------------------------------------------------------------------
SUBROUTINE find_combination(at, at1, ur, is_combination)
!----------------------------------------------------------------------
!
!   This routine receives six vectors at and at1. The three vectors
!   at are obtained as a linear combination with integer coefficients
!   of the at1. 
!   The routine gives as output the matrix that gives the linear combination.
!   at, at1 are in cartesian coordinates. ur is a matrix of integers 
!   even if we store it in a real matrix. 
!   If the at are not a linear combination with integer coefficient,
!   the variable is_combination is set to .FALSE.
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP), INTENT(IN) :: at(3,3), at1(3,3)
REAL(DP), INTENT(OUT) :: ur(3,3)
LOGICAL, INTENT(OUT) :: is_combination

REAL(DP), PARAMETER :: eps1=1.D-8
REAL(DP) :: bg(3,3)
INTEGER :: ipol, jpol

CALL recips(at1(1,1), at1(1,2), at1(1,3), bg(1,1), bg(1,2), bg(1,3)) 

is_combination=.TRUE.
DO ipol=1,3
   DO jpol=1,3
      ur(ipol,jpol) = at(1,ipol) * bg(1,jpol) +  &
                      at(2,ipol) * bg(2,jpol) +  &
                      at(3,ipol) * bg(3,jpol) 
      IF (ABS(ur(ipol,jpol)-NINT(ur(ipol,jpol))) > eps1 ) is_combination=.FALSE.
   END DO
END DO

RETURN
END SUBROUTINE find_combination

!----------------------------------------------------------------------
LOGICAL FUNCTION is_bravais_lattice(at, erre, setn)
!----------------------------------------------------------------------
!
!   This routine receives three primitive vectors at and a 
!   vector erre. It becomes .TRUE. if erre is a Bravais lattice vector
!   of the Bravais lattice defined by the at. In any case it gives
!   setn, the components of erre in the at basis.
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP), INTENT(IN) :: at(3,3), erre(3)
REAL(DP), INTENT(OUT) :: setn(3)

REAL(DP), PARAMETER :: eps1=1.D-8
REAL(DP) :: bg(3,3)
INTEGER :: ipol, jpol

CALL recips(at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3)) 

is_bravais_lattice=.TRUE.
DO ipol=1,3
   setn(ipol) = erre(1) * bg(1,ipol) +  &
                erre(2) * bg(2,ipol) +  &
                erre(3) * bg(3,ipol) 
   is_bravais_lattice = is_bravais_lattice .AND. (ABS(setn(ipol)-&
                                                 NINT(setn(ipol))) < eps1) 
END DO

RETURN
END FUNCTION is_bravais_lattice

!----------------------------------------------------------------------
LOGICAL FUNCTION same_lattice(at, at2, ur, at1, sr)
!----------------------------------------------------------------------
!
!   This routine receives six vectors at and at2 and returns .TRUE.
!   when these two sets of primitive vectors generate the same Bravais
!   lattice. Three vectors at1 are obtained as a linear combination 
!   with integer coefficients of the at2 and the at are obtained 
!   rotating the at1.
!   This function espresses at in terms of at2 using ur the integer
!   coefficients that links at1 and at2 and the rotation matrix sr that
!   links at to at1:
!
!   at1_{\alpha,i} = \sum_j ur_ij at2_{\alpha,j}
!   at_{\alpha,i}  = \sum_\beta sr_{\alpha,\beta} at1_{\beta,i}
!   at, at1, at2 and sr are in cartesian coordinates.
!
!   ur_ij and sr_{\alpha,\beta} as well as the at1 are output of this
!   function when it returns .TRUE..
!
USE kinds, ONLY : DP
USE rotate, ONLY : find_rotation, is_rotation

IMPLICIT NONE

REAL(DP), INTENT(IN) :: at(3,3), at2(3,3)
REAL(DP), INTENT(OUT) :: ur(3,3), at1(3,3), sr(3,3)

REAL(DP) :: amodulus(3), max_modulus, bg2(3,3), tvec(3), tmod, omega, omega0
REAL(DP), PARAMETER :: eps1=1.D-6, eps2=1.D-4
LOGICAL :: is_combination
INTEGER, PARAMETER :: max_comb=96
INTEGER :: ipol, ivec, nx1, nx2, nx3, comb(3,3,max_comb), nc(3), i, j, k
!
!   check that at and at2 have the same volume
!
CALL compute_omega(at,omega0)
CALL compute_omega(at2,omega)
!
!  if they have not, the two lattices are different.  
!
same_lattice=.FALSE.
IF (ABS(omega-omega0)> eps2) RETURN
!
!  first check if a linear combination works
!
sr=0.0_DP
DO ipol=1,3
   sr(ipol,ipol)=1.0_DP
END DO

CALL find_combination(at, at2, ur, is_combination)
IF (is_combination) THEN
   at1=at
   same_lattice=.TRUE.
   RETURN
END IF
!
!  then try with a simple rotation
!
ur=0.0_DP
DO ipol=1,3
   ur(ipol,ipol)=1.0_DP
END DO
CALL find_rotation(at,at2,sr)
IF (is_rotation(sr)) THEN
   at1=at2
   same_lattice=.TRUE.
   RETURN
END IF
!
!  if the routine arrives here, none of the above worked.
!  We check if a rotation and a linear combination work. 
!  First find the maximum modulus of at
!
DO ivec=1,3
   amodulus(ivec)=at(1,ivec)**2 + at(2,ivec)**2 + at(3,ivec)**2
ENDDO
max_modulus=MAX(amodulus(1), amodulus(2), amodulus(3))
!
!  Now find the size of the mesh in which we must check for the modulus
!
CALL recips(at2(1,1), at2(1,2), at2(1,3), bg2(1,1), bg2(1,2), bg2(1,3))
nx1=INT(SQRT(max_modulus)*SQRT(bg2(1,1)**2+bg2(2,1)**2+bg2(3,1)**2))+1
nx2=INT(SQRT(max_modulus)*SQRT(bg2(1,2)**2+bg2(2,2)**2+bg2(3,2)**2))+1
nx3=INT(SQRT(max_modulus)*SQRT(bg2(1,3)**2+bg2(2,3)**2+bg2(3,3)**2))+1

!WRITE(stdout,*) 'same lattice using grid', nx1, nx2, nx3
!
!  save the combinations that have the correct modulus
!
nc=0
DO i=nx1, -nx1, -1
   DO j=nx2, -nx2, -1
      DO k=nx3, -nx3, -1
         tvec(:) = i * at2(:,1) + j * at2(:,2) + k * at2(:,3)
         tmod=tvec(1)**2 + tvec(2)**2 + tvec(3)**2
         DO ipol=1,3
            IF (ABS(tmod - amodulus(ipol))<eps1) THEN
               nc(ipol)=nc(ipol)+1
               IF (nc(ipol)> max_comb) &
                  CALL errore('same_lattice','increase max_comb',1)
               comb(ipol,1,nc(ipol)) = i   
               comb(ipol,2,nc(ipol)) = j   
               comb(ipol,3,nc(ipol)) = k   
            END IF
         END DO
      END DO
   END DO
END DO
!
!   And now try all combinations with the correct modulus until one is found
!   such that at and at1 are related by a pure rotation. Exclude the case
!   in which the at are a supercell
!
DO i=1,nc(1)
   DO j=1,nc(2)
      DO k=1,nc(3)
         at1(:,1)=comb(1,1,i) * at2(:,1) +  &
                  comb(1,2,i) * at2(:,2) +  &
                  comb(1,3,i) * at2(:,3) 
         at1(:,2)=comb(2,1,j) * at2(:,1) +  &
                  comb(2,2,j) * at2(:,2) +  &
                  comb(2,3,j) * at2(:,3) 
         at1(:,3)=comb(3,1,k) * at2(:,1) +  &
                  comb(3,2,k) * at2(:,2) +  &
                  comb(3,3,k) * at2(:,3) 
         CALL find_rotation(at,at1,sr)
         IF (is_rotation(sr)) THEN
            CALL compute_omega(at1,omega)
            IF (ABS(omega-omega0)<eps1) THEN
               same_lattice=.TRUE.
               DO ipol=1,3
                  ur(1,ipol)=DBLE(comb(1,ipol,i))
                  ur(2,ipol)=DBLE(comb(2,ipol,j))
                  ur(3,ipol)=DBLE(comb(3,ipol,k))
               ENDDO
               RETURN
            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
!  if the routine arrives here at and at2 do not generate the
!  same Bravais lattice.
!
   same_lattice=.FALSE.
RETURN
END FUNCTION same_lattice

!----------------------------------------------------------------------
SUBROUTINE lattice_point_group(at,gname,code_group,nsym,sr)
!----------------------------------------------------------------------
!
!   This routine receives three primitive vectors, at, of a Bravais lattice
!   and finds the number of symmetries, the rotation matrices, and the point 
!   group code (1-32) of the Bravais lattice.
!
USE kinds,  ONLY : DP
USE rotate, ONLY : find_rotation, is_rotation
IMPLICIT NONE
REAL(DP), INTENT(IN) :: at(3,3)
INTEGER, INTENT(OUT) :: code_group, nsym
REAL(DP), INTENT(OUT) :: sr(3,3,48)
CHARACTER (LEN=11), INTENT(OUT) :: gname

REAL(DP) :: sr_test(3,3), bg(3,3), amodulus(3), max_modulus, &
            tvec(3), tmod, at1(3,3)
REAL(DP), PARAMETER :: eps1=1.D-7
INTEGER, PARAMETER :: max_comb=96
INTEGER :: ivec, nx1, nx2, nx3, i, j, k, nc(3), comb(3,3, max_comb), ipol
!
!  Start finding the moduli of the input vectors 
!
DO ivec=1,3
   amodulus(ivec)=at(1,ivec)**2 + at(2,ivec)**2 + at(3,ivec)**2
ENDDO
max_modulus=MAX(amodulus(1), amodulus(2), amodulus(3))

!
!  Now find the size of the mesh in which we must check for the modulus
!
CALL recips(at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))
nx1=INT(SQRT(max_modulus)*SQRT(bg(1,1)**2+bg(2,1)**2+bg(3,1)**2))+1
nx2=INT(SQRT(max_modulus)*SQRT(bg(1,2)**2+bg(2,2)**2+bg(3,2)**2))+1
nx3=INT(SQRT(max_modulus)*SQRT(bg(1,3)**2+bg(2,3)**2+bg(3,3)**2))+1

!WRITE(stdout,*) 'lattice point group using grid', nx1, nx2, nx3

nc=0
DO i=nx1, -nx1, -1
   DO j=nx2, -nx2, -1
      DO k=nx3, -nx3, -1
         tvec(:) = i * at(:,1) + j * at(:,2) + k * at(:,3)
         tmod = tvec(1)**2 + tvec(2)**2 + tvec(3)**2
         DO ipol=1,3
            IF (ABS(tmod - amodulus(ipol)) < eps1*tmod) THEN
               nc(ipol)=nc(ipol)+1
               IF (nc(ipol)> max_comb) &
                  CALL errore('lattice_point_group','increase max_comb',1)
               comb(ipol,1,nc(ipol)) = i
               comb(ipol,2,nc(ipol)) = j
               comb(ipol,3,nc(ipol)) = k
            END IF
         END DO
      END DO
   END DO
END DO
!
!   And now try all combinations with the correct modulus until one is found
!   such that at and at1 are related by a pure rotation
!
nsym=0
DO i=1,nc(1)
   DO j=1,nc(2)
      DO k=1,nc(3)
         at1(:,1)=comb(1,1,i) * at(:,1) +  &
                  comb(1,2,i) * at(:,2) +  &
                  comb(1,3,i) * at(:,3)
         at1(:,2)=comb(2,1,j) * at(:,1) +  &
                  comb(2,2,j) * at(:,2) +  &
                  comb(2,3,j) * at(:,3)
         at1(:,3)=comb(3,1,k) * at(:,1) +  &
                  comb(3,2,k) * at(:,2) +  &
                  comb(3,3,k) * at(:,3)
         CALL find_rotation(at,at1,sr_test)
         IF (is_rotation(sr_test)) THEN
            nsym=nsym+1
            IF (nsym > 48) &
                  CALL errore('lattice_point_group','too many symmetries',1)
            sr(:,:,nsym)=sr_test(:,:)
         ENDIF
      ENDDO
   ENDDO
ENDDO

CALL find_group(nsym, sr, gname, code_group)

RETURN
END SUBROUTINE lattice_point_group

!----------------------------------------------------------------------
LOGICAL FUNCTION bravais_dir(at, ax, enne, rvec_out, rvecmod_out )
!----------------------------------------------------------------------
!
!   This routine receives three primitive lattices at and a direction ax
!   and finds the shorter Bravais lattice vector in the direction of 
!   the versor ux. The function gives .FALSE. if none of the computed 
!   vectors is in the direction of ax or if this vector has n1, n2, n3 
!   outside the checked range.
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP), INTENT(IN)  :: at(3,3), ax(3)
REAL(DP), INTENT(OUT) :: rvec_out(3), rvecmod_out
INTEGER, INTENT(OUT)  :: enne(3)

INTEGER, PARAMETER  :: npx=150
REAL(DP), PARAMETER :: eps1=1.d-8
INTEGER  :: ipol, i1, i2, i3, n1, n2, n3
REAL(DP) :: bg(3,3), rvec(3), rvecmod, ps, min_mod, proj(3), rn2, rn3

CALL recips(at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))

DO ipol=1,3
   proj(ipol)= bg(1, ipol) * ax(1) + bg(2, ipol) * ax(2) + bg(3, ipol) * ax(3)
END DO

IF (ABS(proj(1))>eps1) THEN
   i1=1
   i2=2
   i3=3
ELSEIF (ABS(proj(2))>eps1) THEN
   i1=2
   i2=3
   i3=1
ELSEIF (ABS(proj(3))>eps1) THEN
   i1=3
   i2=1
   i3=2
ELSE
   CALL errore('bravais_dir',' versor is zero',1)
ENDIF

min_mod=1.D8
bravais_dir=.FALSE.
DO n1=npx, -npx, -1
   rn2= n1 * proj(i2) / proj(i1) 
   rn3= n1 * proj(i3) / proj(i1) 
   IF (ABS(rn2 - NINT(rn2))< eps1 .AND. ABS(rn3 - NINT(rn3))< eps1) THEN
      n2 = NINT(rn2)
      n3 = NINT(rn3)
      rvec(:) = n1 * at(:,i1) + n2 * at(:,i2) + n3 * at(:,i3)
      ps= SQRT(rvec(1)**2 + rvec(2)**2 + rvec(3)**2)
      IF (ps > eps1 .AND. ps<min_mod) THEN
         enne(i1)=n1
         enne(i2)=n2
         enne(i3)=n3
         min_mod=ps
         rvec_out(:)=rvec(:)
         rvecmod_out=ps
         bravais_dir=.TRUE.
      ENDIF
   ENDIF
ENDDO

RETURN
END FUNCTION bravais_dir

!----------------------------------------------------------------------
SUBROUTINE compute_omega(at, omega)
!----------------------------------------------------------------------
!
!  Compute the volume of the unit cell given by the at
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP), INTENT(IN) :: at(3,3)
REAL(DP), INTENT(OUT) :: omega

omega= at(1,1) * (at(2,2)*at(3,3) - at(2,3)*at(3,2) ) - &
       at(2,1) * (at(1,2)*at(3,3) - at(1,3)*at(3,2) ) + &
       at(3,1) * (at(1,2)*at(2,3) - at(1,3)*at(2,2) )  

omega=ABS(omega)

RETURN
END SUBROUTINE compute_omega

!----------------------------------------------------------------------
SUBROUTINE lattice_name(ibrav, latt_name)
!----------------------------------------------------------------------
!
!  this subroutine receives as input the Bravais lattice vector index
!  and gives as output the lattice name
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav
CHARACTER(LEN=40), INTENT(OUT) :: latt_name

SELECT CASE (ibrav)
    CASE(0)
        latt_name='free lattice'
    CASE(1)
        latt_name='simple cubic'
    CASE(2)
        latt_name='face centered cubic'
    CASE(3)
        latt_name='body centered cubic'
    CASE(4)
        latt_name='hexagonal'
    CASE(5,-5)
        latt_name='trigonal'
    CASE(6)
        latt_name='tetragonal'
    CASE(7)
        latt_name='centered tetragonal'
    CASE(8)
        latt_name='simple orthorhombic'
    CASE(9, -9)
        latt_name='one face centered orthorhombic (C)'
    CASE(91)
        latt_name='one face centered orthorhombic (A)'
    CASE(10)
        latt_name='face centered orthorhombic'
    CASE(11)
        latt_name='body centered orthorhombic'
    CASE(12)
        latt_name='monoclinic (c unique)'
    CASE(-12)
        latt_name='monoclinic (b unique)'
    CASE(13)
        latt_name='base centered monoclinic (c unique)'
    CASE(-13)
        latt_name='base centered monoclinic (b unique)'
    CASE(14)
        latt_name='triclinic'
CASE DEFAULT
     CALL errore('lattice_name','ibrav not known',1)
END SELECT

RETURN
END SUBROUTINE lattice_name

!----------------------------------------------------------------------
LOGICAL FUNCTION zone_border(vect,at,bg,iflag)
!----------------------------------------------------------------------
!
!   This function receives a vector vect, the direct and reciprocal lattice
!   vectors at and bg and a iflag=+-1.
!   The vector vect is in cartesian coordinates. 
!   When iflag is +1 it is supposed to be in real space and the
!   function returns .TRUE. if more than one Bravais lattice point is
!   at the same distance from vect.
!   When iflag is -1 the vectors vect is supposed to be in reciprocal
!   space and the function return .TRUE. if more than one reciprocal
!   lattice vector is at the same distance from vect.
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP), INTENT(IN) :: vect(3), at(3,3), bg(3,3)
INTEGER, INTENT(IN) :: iflag

REAL(DP) :: max_modulus, gvec(3), tvec(3), tmod
REAL(DP), PARAMETER :: eps1=1.D-10
INTEGER :: nc, i, j, k
INTEGER :: nx1, nx2, nx3

max_modulus=vect(1)**2 + vect(2)**2 + vect(3)**2

IF (iflag==1) THEN
   nx1=INT(SQRT(4.0_DP*max_modulus)*SQRT(bg(1,1)**2+bg(2,1)**2+bg(3,1)**2))+1
   nx2=INT(SQRT(4.0_DP*max_modulus)*SQRT(bg(1,2)**2+bg(2,2)**2+bg(3,2)**2))+1
   nx3=INT(SQRT(4.0_DP*max_modulus)*SQRT(bg(1,3)**2+bg(2,3)**2+bg(3,3)**2))+1
ELSE
   nx1=INT(SQRT(4.0_DP*max_modulus)*SQRT(at(1,1)**2+at(2,1)**2+at(3,1)**2))+1
   nx2=INT(SQRT(4.0_DP*max_modulus)*SQRT(at(1,2)**2+at(2,2)**2+at(3,2)**2))+1
   nx3=INT(SQRT(4.0_DP*max_modulus)*SQRT(at(1,3)**2+at(2,3)**2+at(3,3)**2))+1
ENDIF

nc=0
DO i=nx1, -nx1, -1
   DO j=nx2, -nx2, -1
      DO k=nx3, -nx3, -1
         IF (iflag==1) THEN
            gvec(:) = i * at(:,1) + j * at(:,2) + k * at(:,3)
         ELSE
            gvec(:) = i * bg(:,1) + j * bg(:,2) + k * bg(:,3)
         ENDIF
         tvec(:)=vect(:) - gvec(:)
         tmod = tvec(1)**2 + tvec(2)**2 + tvec(3)**2
         IF ((tmod - max_modulus)<-eps1) THEN
            nc=1
            max_modulus=tmod
         ELSEIF (ABS(tmod - max_modulus)<eps1) THEN
            nc=nc+1
         ENDIF
      END DO
   END DO
END DO

zone_border = nc>1

RETURN
END FUNCTION zone_border

!----------------------------------------------------------------------
LOGICAL FUNCTION same_star(nsym, s, xk1, xk2, at)
!----------------------------------------------------------------------
!
!  This routine receives a point group of order nsym and its rotations
!  matrices s in the crystal basis. It uses these symmetries to see if
!  xk1 and xk2 belong to the same star of k points. It returns .TRUE.
!  if this is the case. It return .TRUE. even if xk2 differ by a reciprocal
!  lattice vector from a vector of the star of xk1
!  xk1 and xk2 are in cartesian coordinates.
!  at are the direct lattice vectors.
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nsym, s(3,3,nsym)

REAL(DP) :: xk1(3), xk2(3), at(3,3), xkr(3), xk1c(3), xk2c(3), xkrmod
INTEGER :: isym, ipol 
!
!   Bring xk1 in reciprocal space
!
xk1c=xk1
CALL cryst_to_cart(1, xk1c, at, -1)
xk2c=xk2
CALL cryst_to_cart(1, xk2c, at, -1)

same_star=.FALSE.
DO isym=1,nsym
   DO ipol=1,3
      xkr(ipol) = s(ipol,1,isym) * xk1c(1) &
                + s(ipol,2,isym) * xk1c(2) &
                + s(ipol,3,isym) * xk1c(3)
      xkr(ipol) = xk2c(ipol) - xkr(ipol) - nint( xk2c(ipol) - xkr(ipol) )
   ENDDO
   xkrmod= xkr(1)**2 + xkr(2)**2 + xkr(3)**2
   IF (xkrmod < 1.D-9) THEN
      same_star=.TRUE.
      RETURN
   ENDIF
ENDDO   

RETURN
END FUNCTION same_star

!----------------------------------------------------------------------
FUNCTION is_compatible_group_ibrav(code_group, ibrav, is_compatible, ncomp)
!----------------------------------------------------------------------
!
!  This function returns .TRUE. is the Bravais lattice and the
!  point group are compatible.
!
USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: code_group, ibrav
INTEGER, INTENT(OUT) :: is_compatible(32,6), ncomp(32)
LOGICAL :: is_compatible_group_ibrav

INTEGER :: i

is_compatible=0
!
!   C_1, C_i triclinic
!
is_compatible(1,1)=14
ncomp(1)=1
is_compatible(2,1)=14
ncomp(2)=1
!
!   C_s, C_2, monoclinic
!
is_compatible(3,1)=12
is_compatible(3,2)=13
is_compatible(3,3)=-12
is_compatible(3,4)=-13
ncomp(3)=4
is_compatible(4,1)=12
is_compatible(4,2)=13
is_compatible(4,3)=-12
is_compatible(4,4)=-13
ncomp(4)=4
!
!   C_3, trigonal or hexagonal
!
is_compatible(5,1)=4
is_compatible(5,2)=5
is_compatible(5,3)=-5
ncomp(5)=3
!
!   C_4, tetragonal
!
is_compatible(6,1)=6
is_compatible(6,2)=7
ncomp(6)=2
!
!   C_6, hexagonal
!
is_compatible(7,1)=4
ncomp(7)=1
!
!   D_2, orthorhombic
!
is_compatible(8,1)=8
is_compatible(8,2)=9
is_compatible(8,3)=10
is_compatible(8,4)=11
is_compatible(8,5)=-9
is_compatible(8,6)=91
ncomp(8)=6
!
!   D_3 trigonal or hexagonal
!
is_compatible(9,1)=4
is_compatible(9,2)=5
is_compatible(9,3)=-5
ncomp(9)=3
!
!   D_4 tetragonal
!
is_compatible(10,1)=6
is_compatible(10,2)=7
ncomp(10)=2
!
!   D_6 hexagonal
!
is_compatible(11,1)=4
ncomp(11)=1
!
!   C_2v orthorhombic
!
is_compatible(12,1)=8
is_compatible(12,2)=9
is_compatible(12,3)=10
is_compatible(12,4)=11
is_compatible(12,5)=-9
is_compatible(12,6)=91
ncomp(12)=6
!
!   C_3v hexagonal or trigonal
!
is_compatible(13,1)=4
is_compatible(13,2)=5
is_compatible(13,3)=-5
ncomp(13)=3
!
!   C_4v tetragonal
!
is_compatible(14,1)=6
is_compatible(14,2)=7
ncomp(14)=2
!
!   C_6v hexagonal
!
is_compatible(15,1)=4
ncomp(15)=1
!
!   C_2h monoclinic
!
is_compatible(16,1)=12
is_compatible(16,2)=13
is_compatible(16,3)=-12
is_compatible(16,4)=-13
ncomp(16)=4
!
!  C_3h hexagonal
!
is_compatible(17,1)=4
ncomp(17)=1
!
!  C_4h tetragonal
!
is_compatible(18,1)=6
is_compatible(18,2)=7
ncomp(18)=2
!
!  C_6h hexagonal
!
is_compatible(19,1)=4
ncomp(19)=1
!
!  D_2h orthorhombic
!
is_compatible(20,1)=8
is_compatible(20,2)=9
is_compatible(20,3)=10
is_compatible(20,4)=11
is_compatible(20,5)=-9
is_compatible(20,6)=91
ncomp(20)=6
!
!  D_3h hexagonal
!
is_compatible(21,1)=4
ncomp(21)=1
!
!  D_4h tetragonal
!
is_compatible(22,1)=6
is_compatible(22,2)=7
ncomp(22)=2
!
!  D_6h hexagonal
!
is_compatible(23,1)=4
ncomp(23)=1
!
!  D_2d tetragonal
!
is_compatible(24,1)=6
is_compatible(24,2)=7
ncomp(24)=1
!
!   D_3d hexagonal or trigonal
!
is_compatible(25,1)=4
is_compatible(25,2)=5
is_compatible(25,3)=-5
ncomp(25)=3
!
!   S_4 tetragonal
!
is_compatible(26,1)=6
is_compatible(26,2)=7
ncomp(26)=2
!
!   S_6 hexagonal or trigonal
!
is_compatible(27,1)=4
is_compatible(27,2)=5
is_compatible(27,3)=-5
ncomp(27)=3
!
!   T cubic
!
is_compatible(28,1)=1
is_compatible(28,2)=2
is_compatible(28,3)=3
ncomp(28)=3
!
!   T_h cubic
!
is_compatible(29,1)=1
is_compatible(29,2)=2
is_compatible(29,3)=3
ncomp(29)=3
!
!   T_d cubic
!
is_compatible(30,1)=1
is_compatible(30,2)=2
is_compatible(30,3)=3
ncomp(30)=3
!
!   O cubic
!
is_compatible(31,1)=1
is_compatible(31,2)=2
is_compatible(31,3)=3
ncomp(31)=3
!
!   O_h cubic
!
is_compatible(32,1)=1
is_compatible(32,2)=2
is_compatible(32,3)=3
ncomp(32)=3

is_compatible_group_ibrav=.FALSE.

DO i=1, ncomp(code_group)
   IF (is_compatible(code_group,i)==ibrav) is_compatible_group_ibrav=.TRUE.
ENDDO

RETURN
END FUNCTION is_compatible_group_ibrav

!----------------------------------------------------------------------
SUBROUTINE print_bravais_description(ibrav,celldm)
!----------------------------------------------------------------------
!
USE kinds, ONLY : DP
USE constants, ONLY : pi
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav
REAL(DP), INTENT(IN) :: celldm(6)

CHARACTER(LEN=40) :: latt_name

CALL lattice_name(ibrav, latt_name)
WRITE(stdout,'(/,5x, "ibrav=",i3,": ", a)') ibrav, TRIM(latt_name)
WRITE(stdout,'(5x,"Cell parameters:")')

SELECT CASE (ibrav)
   CASE(1)  
      WRITE(stdout,'(/,5x,"alat=",f11.6," a.u.",/)') celldm(1)
   CASE(2)  
      WRITE(stdout,'(/,5x,"alat=",f11.6," a.u.",/)') celldm(1)
   CASE(3)  
      WRITE(stdout,'(/,5x,"alat=",f11.6," a.u.",/)') celldm(1)
   CASE(4)  
      WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., c/a=",f11.6,/)') celldm(1), &
                                                            celldm(3)
   CASE(5)  
      WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., cos(alpha)=",f11.6,/)') &
                                        celldm(1), celldm(4)
   CASE(6)  
      WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., c/a=",f11.6,/)') &
                                  celldm(1), celldm(3)
   CASE(7)  
      WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., c/a=",f11.6,/)') celldm(1), &
                                                            celldm(3)
   CASE(8,9,-9,91,10,11)  
      WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., b/a=",f11.6,&
                          &", c/a=",f11.6,/)') celldm(1), celldm(2), celldm(3)
   CASE(12,13)  
      WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., b/a=",f11.6,&
                   &", c/a=",f11.6)') celldm(1), celldm(2), celldm(3)
      WRITE(stdout,'(/,5x,"cos(gamma)=",f11.6,/)') celldm(4)
      WRITE(stdout,'(5x,"gamma=",f11.6,/)') ACOS(celldm(4))*180.0_DP/pi
   CASE(-12,-13)
      WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., b/a=",f11.6,&
                  &", c/a=",f11.6)') celldm(1), celldm(2), celldm(3)
      WRITE(stdout,'(/,5x,"cos(beta)=",f11.6,/)') celldm(5)
      WRITE(stdout,'(5x,"beta=",f11.6,/)') ACOS(celldm(5))*180.0_DP/pi
   CASE(14)  
      WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., b/a=",f11.6,&
                  &", c/a=",f11.6)') celldm(1), celldm(2), celldm(3)
      WRITE(stdout,'(/,5x,"cos(alpha)=",f11.6, &
                     &", cos(beta)=",f11.6,", cos(gamma)=",f11.6)') &
                      celldm(4), celldm(5), celldm(6)
      WRITE(stdout,'(5x,"alpha=",f11.6," deg,",2x,"beta=",f11.6,&
                                      &" deg,",2x,"gamma=",f11.6," deg",/)') &
                      ACOS(celldm(4))*180.0_DP/pi,  &
                      ACOS(celldm(5))*180.0_DP/pi,  &
                      ACOS(celldm(6))*180.0_DP/pi
  CASE(0)
  CASE DEFAULT
     CALL errore('print_bravais_description','ibrav not programmed',1)
  END SELECT

RETURN
END SUBROUTINE print_bravais_description

!----------------------------------------------------------------------
INTEGER FUNCTION crystal_parameters(ibrav)
!----------------------------------------------------------------------
!
!   This function returs the number of independent crystal parameters
!   for each Bravais lattice index
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav

INTEGER :: nvar

SELECT CASE (ibrav)
   CASE(1,2,3)
      nvar=1
   CASE(4,5,6,7)
      nvar=2
   CASE(8,9,91,10,11)
      nvar=3
   CASE(12,-12,13,-13)
      nvar=4
   CASE DEFAULT
      nvar=6
END SELECT

crystal_parameters=nvar

RETURN
END FUNCTION crystal_parameters

!----------------------------------------------------------------------
SUBROUTINE find_closest_at( gvec, at, enne, m, rmu )
!----------------------------------------------------------------------
!
!  This subroutine finds the m vectors of the Bravais lattice defined
!  by the primitive vectors at, closest to the origin and belonging to 
!  the lattice plane perpendicular to gvec (given in the basis of bg)
!  Note that only a subset of the plane point are tested. Increase npx
!  if this routine gives problems.
!
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: gvec(3), m
REAL(DP), INTENT(IN) :: at(3,3)

REAL(DP), INTENT(INOUT) :: rmu(3,m)
INTEGER, INTENT(INOUT) :: enne(3,m)

INTEGER, PARAMETER :: npx=10
INTEGER :: nn, ip, iq, ir, i1, i2, i3, i, icount, nvecx, fact, ip1, iq1, ir1
REAL(DP), ALLOCATABLE :: r(:,:), rmod(:)
INTEGER,  ALLOCATABLE :: ind(:), renne(:,:)

IF (gvec(1) /= 0) THEN
   i1=1
   i2=2
   i3=3
ELSEIF (gvec(2) /= 0) THEN
   i1=2
   i2=3
   i3=1
ELSE
   i1=3
   i2=1
   i3=2
ENDIF

icount=0
nn=ABS(gvec(i1))
nvecx = (2 * npx * nn + 1) **2 
ALLOCATE(r(3,nvecx))
ALLOCATE(rmod(nvecx))
ALLOCATE(renne(3,nvecx))
ALLOCATE(ind(nvecx))
DO iq=-npx*nn, npx*nn
   DO ir=-npx*nn, npx*nn
      ip= - ( iq * gvec(i2) + ir * gvec(i3) ) 
!
!  The condition ip>=0 removes one half of the points so as to avoid to have
!  R and -R in the list. 
!
      IF ((MOD(ip, nn) == 0 .AND. ip > 0) .OR. (ip==0.AND.iq > 0) .OR. &
              (ip==0.AND.iq==0.AND.ir>0) ) THEN
!
!   Avoid to take parallel vectors. Only when ip/nn, iq, and ir have no
!   common factor the point is inserted in the list
!
         ip1=ip/nn
         iq1=iq
         ir1=ir
         CALL remove_common_factors(ip1,iq1,ir1,fact)
         IF (fact /= 1) CYCLE
         icount=icount+1
         r(:,icount)= ip * at(:,i1) / nn + iq * at(:,i2) + ir * at(:,i3)
         renne(i1,icount)=ip / nn
         renne(i2,icount)=iq
         renne(i3,icount)=ir
         rmod(icount)=SQRT(r(1,icount)**2+r(2,icount)**2+r(3,icount)**2)
         ind(icount)=icount
      ENDIF
   ENDDO
ENDDO

IF (icount < m) CALL errore('find_closest_at','too few vectors increase npx',1)

CALL hpsort(icount, rmod, ind)

!DO i = 1, icount
!   WRITE(6,'(5i8,f15.7)') i, ind(i), renne(1, ind(i)), renne(2, ind(i)), &
!               renne(3, ind(i)), rmod(i)
!END DO

DO icount = 1, m
   rmu(:,icount)=r(:,ind(icount))
   enne(:,icount)=renne(:,ind(icount))
ENDDO

DEALLOCATE (r)
DEALLOCATE (rmod)
DEALLOCATE (renne)
DEALLOCATE (ind)

RETURN
END SUBROUTINE find_closest_at

!----------------------------------------------------------------------
SUBROUTINE remove_common_factors(m,n,o,fact)
!----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(INOUT) :: m, n, o
INTEGER, INTENT(OUT) :: fact

INTEGER :: ind

fact=1
DO ind=2,10
   IF (MOD(m,ind)==0.AND.MOD(n,ind)==0.AND.MOD(o,ind)==0) THEN
      m = m/ind
      n = n/ind
      o = o/ind
      fact=fact*ind
   END IF
ENDDO

RETURN
END SUBROUTINE remove_common_factors

!----------------------------------------------------------------------
SUBROUTINE compress_celldm(cm,x,nvar,ibrav)
!----------------------------------------------------------------------
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ibrav
REAL(DP), INTENT(IN) :: cm(6)
REAL(DP), INTENT(INOUT) :: x(nvar)

SELECT CASE (ibrav)
   CASE(1,2,3)
      x(1) = cm(1)
   CASE(4,5,6,7)
      x(1) = cm(1)
      x(2) = cm(3)
      IF (ibrav==5) x(2) = ACOS(cm(4))
   CASE(8,9,91,10,11)
      x(1) = cm(1)
      x(2) = cm(2)
      x(3) = cm(3)
   CASE(12,-12,13,-13)
      x(1) = cm(1)
      x(2) = cm(2)
      x(3) = cm(3)
      IF (ibrav>0) THEN
!
!   c unique
!
         x(4) = ACOS(cm(4))
      ELSE
!
!   b unique
!
         x(4) = ACOS(cm(5))
      ENDIF
   CASE DEFAULT
      x(1) = cm(1)
      x(2) = cm(2)
      x(3) = cm(3)
      x(4) = ACOS( cm(4) )
      x(5) = ACOS( cm(5) )
      x(6) = ACOS( cm(6) )
END SELECT

RETURN
END SUBROUTINE compress_celldm

!----------------------------------------------------------------------
SUBROUTINE compress_int_vect(inv,x,n,ibrav)
!----------------------------------------------------------------------

USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, ibrav
INTEGER, INTENT(IN) :: inv(n)
INTEGER, INTENT(INOUT) :: x(n)

INTEGER :: i

x=inv
SELECT CASE (ibrav)
   CASE(4,5,6,7)
      DO i=1,n
         IF (inv(i)==3) x(i)=2
         IF (inv(i)==4) x(i)=2 
      ENDDO
   CASE(12,-12,13,-13)
      DO i=1,n
         IF (inv(i)==5) x(i)=4
      ENDDO
   CASE DEFAULT
END SELECT

RETURN
END SUBROUTINE compress_int_vect

!----------------------------------------------------------------------
SUBROUTINE expand_celldm(cm, x, nvar, ibrav)
!----------------------------------------------------------------------
!
!  This routine receives a set of compressed crystallographic parameters 
!  in the array x(nvar) and transforms them in the celldm array
!
USE kinds, ONLY : DP

IMPLICIT NONE
INTEGER,  INTENT(IN)    :: ibrav, nvar
REAL(DP), INTENT(IN)    :: x(nvar)
REAL(DP), INTENT(INOUT) :: cm(6)

cm=0.0_DP
SELECT CASE (ibrav)
   CASE(1,2,3) 
      cm(1)=x(1)
   CASE(4,5,6,7)
      cm(1)=x(1)
      IF (ibrav==5) THEN
         cm(4)=COS(x(2))
      ELSE
         cm(3)= x(2)
      ENDIF
   CASE(8,9,91,10,11)
      cm(1)=x(1)
      cm(2)=x(2)
      cm(3)=x(3)
   CASE(12,-12,13,-13) 
      cm(1)=x(1)
      cm(2)=x(2)
      cm(3)=x(3)
      IF (ibrav>0) THEN
!
!   c unique
!
         cm(4)=COS(x(4))
      ELSE
!
!   b unique
!
         cm(5)=COS(x(4))
      ENDIF
   CASE DEFAULT
      cm(1)=x(1)
      cm(2)=x(2)
      cm(3)=x(3)
      cm(4)=COS(x(4))
      cm(5)=COS(x(5))
      cm(6)=COS(x(6))
END SELECT

RETURN
END SUBROUTINE expand_celldm

SUBROUTINE needed_celldm(ibrav,celldm_in_use)

IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav
LOGICAL, INTENT(OUT) :: celldm_in_use(6)

celldm_in_use=.FALSE.
celldm_in_use(1)=.TRUE.

SELECT CASE (ibrav)
   CASE(1,2,3)
   CASE(4,6,7) 
      celldm_in_use(3)=.TRUE.
   CASE(5) 
      celldm_in_use(4)=.TRUE.
   CASE(8,9,91,10,11) 
      celldm_in_use(2)=.TRUE.
      celldm_in_use(3)=.TRUE.
   CASE(12,13) 
      celldm_in_use(2)=.TRUE.
      celldm_in_use(3)=.TRUE.
      celldm_in_use(4)=.TRUE.
   CASE(-12,-13) 
      celldm_in_use(2)=.TRUE.
      celldm_in_use(3)=.TRUE.
      celldm_in_use(5)=.TRUE.
   CASE DEFAULT
      celldm_in_use=.TRUE.
END SELECT
RETURN
END SUBROUTINE needed_celldm

END MODULE lattices
