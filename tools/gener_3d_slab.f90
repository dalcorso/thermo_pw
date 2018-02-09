!
! Copyright (C) 2014-2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM gener_3d_slab
!
!  This program generates all the points of a tridimensional solid
!  defined by a three dimensional Bravais lattice and nat_3d atoms in
!  the unit cell. Calling a1, a2, and a3 the three primitive vectors of the
!  Bravais lattice, it prints on output the information necessary, 
!  using periodic boundary conditions, to simulate a slab whose surface
!  is perpendicular to a G vector defined as G = m b1 + n b2 + o b3. 
!  It is supposed that m, n, and o, have no common factors. 
!  The user specifies the number of layers nlayers and 
!  the size of the unit cell on the surface. The size of the unit cell on
!  the surface is specified as a 2x2 matrix of integers t11, t12, t21, t22  
!  that define a multiple of the minimal unit cell chosen by the code. 
!  If N_sur is the number of atoms per surface plane defined by the above 
!  integers, the program writes on output the coordinates 
!  of N_sur x nlayers atoms.
!  The distance between repeated copies of the slab can be indicated as 
!  input (vacuum) in a.u. and the program can put that distance
!  (ldist_vacuum=.TRUE.) or a vacuum space equivalent to an integer 
!  number of empty layers sufficient to contain at least a vacuum distance.
!  The atomic coordinates are written on output and periodic boundary
!  conditions are assumed in the plane while repeated slabs are obtained
!  along z.
!  List of input variables:
!
!  ibrav_3d     The Bravais lattice of the bulk, defined as in QE. This
!               code use the same latgen routine.
!  
!  celldm_3d(6) The size of the unit cell.
!
!  lcryst       1 all crystal coordinates as in QE input, 
!               2 only inequivalent atoms crystal coordinates (as in space_sg)
!               3 cartesian coordinates
!
!  nat_3d       The number of atoms inside the unit cell
!
!  atm_3d, tau_3d(3) ! label and coordinates of the atoms inside the unit cell.
!
!  m,n,o        m, n, and o of the G vector
!
!  nlayers      number of rows 
!               
!  t11, t12, t21, t22  where tij are integer numbers. The size of the unit cell 
!               is defined by the two vectors 
!               A'_1 = t11 A_1 + t12 A_2, 
!               A'_2 = t21 A_1 + t22 A_2 
!               where A_1 and A_2 are minimal direct lattice vectors 
!               on the surface plane. 
!               The two dimensional Bravais lattice on the surface plane 
!               is determined by the code.
!
!  ldist_vacuum if .TRUE. the distance between repeated slab is vacuum,
!               otherwise it is the smaller integer number of rows larger
!               than vacuum
!
!  vacuum       the vacuum distance
!
!  origin_shift Number of the layer that contains the origin.
!               When origin_shift is zero, the first atom of the 
!               central layer is in the origin. You might want to 
!               shift the origin to a different layer. origin_shift
!               is an integer that specifies on which layer to put
!               the origin. In any case the layer with the origin will
!               be the central layer. Note that when the number of layers 
!               is even the origin is in the middle of two layers. The 
!               layer indicated by origin_shift has negative z and is 
!               the closest to the origin. If you do not understand
!               the previous sentence or do not know which value to use
!               set origin_shift=0.
!
!  filename     the name of a file were the output coordinates and
!               unit cell are written
!              
USE kinds, ONLY : DP
USE constants, ONLY : pi
USE wyckoff,  ONLY : nattot, tautot, ityptot, sup_spacegroup, clean_spacegroup
USE io_global,   ONLY : ionode, stdout
USE mp_global,   ONLY : mp_startup, mp_global_end
USE environment, ONLY : environment_start, environment_end

IMPLICIT NONE

INTEGER, PARAMETER :: nmax=10000
REAL(DP) :: at(3,3), bg(3,3), g(3), t(3), c1(3), c2(3), d(3), bc1(3), bc2(3), &
            bt(3), beff(3,3), d1(3), d2(3), bd1(3), bd2(3), vprod(3), dz(3)
REAL(DP) :: alat, c1mod, c2mod, tmod, prod1, prod2, prod3, vacuum, &
            gmod, dist, dist1, tau(3), det, omega, prod, area, super_area, &
            alpha, c_alat, shiftz, scalar
INTEGER  :: m, n, h, o, p, q, s, t11, t12, t21, t22, nlayers, nat, ibrav_3d, &
            nat_3d, p0, q0, s0, p01, q01, s01, p1, q1, s1, ia, iat, j, itry, &
            jtry, na, found, nspace, sign1, sign2, ia_first, m1, n1, o1, fact, &
            copies, min_j, lcryst, space_group_number, origin_choice, nt, &
            ntyp, ibrav_sur, origin_shift, i1

REAL(DP) :: y(3,nmax), celldm(6), celldm_sur(6), gc(3)
REAL(DP), PARAMETER :: eps=1.D-6
INTEGER :: ps(nmax), qs(nmax), ss(nmax), ityp_all(nmax)
INTEGER, ALLOCATABLE :: ityp(:), if_pos(:,:)
REAL(DP), ALLOCATABLE :: tau_3d(:,:), tau_sur(:,:), extfor(:,:)
CHARACTER ( LEN=3 ), ALLOCATABLE :: atm_3d(:)
CHARACTER ( LEN=3 ) :: atm(nmax), atm_typ(100)
LOGICAL :: ldist_vacuum, three_indices, uniqueb, rhombohedral
INTEGER :: iuout
INTEGER :: find_free_unit
CHARACTER( LEN=256 ) :: filename
CHARACTER(LEN=9) :: code='3D_SLAB'

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

WRITE(stdout,'("ibrav_3d? ")')
READ(5,*) ibrav_3d
WRITE(stdout,'(i5)') ibrav_3d
WRITE(stdout,'("celldm ")')
READ(5,*) celldm(1), celldm(2), celldm(3), celldm(4), celldm(5), celldm(6)
WRITE(stdout,'(6f12.6)') celldm
alat=celldm(1)
WRITE(stdout,'("Crystal coordinates? (1 yes all, 2 yes only irreducible, 3 no)")')
READ(5,*) lcryst
WRITE(stdout,'(i5)') lcryst

IF (lcryst==2) THEN
   WRITE(stdout,'("Space group number?")')
   READ(5,*) space_group_number
   WRITE(stdout,'(i5)') space_group_number
   WRITE(stdout,'("Unique axis b? (.TRUE. or .FALSE.")') 
   READ(5,*) uniqueb
   WRITE(stdout,'(l5)') uniqueb
   WRITE(stdout,'("Rombohedral? (.TRUE. or .FALSE.")') 
   READ(5,*) rhombohedral
   WRITE(stdout,'(l5)') rhombohedral
   WRITE(stdout,'("Origin choice? (1 or 2)")') 
   READ(5,*) origin_choice
   WRITE(stdout,'(i5)') origin_choice
END IF

WRITE(stdout,'("Number of atoms in the bulk unit cell?")') 
READ(5,*) nat_3d
WRITE(stdout,'(i5)') nat_3d

ALLOCATE(tau_3d(3,nat_3d))
ALLOCATE(ityp(nat_3d))
ALLOCATE(atm_3d(nat_3d))
DO ia=1,nat_3d
   READ(5,*) atm_3d(ia), tau_3d(1,ia), tau_3d(2,ia), tau_3d(3,ia)
   WRITE(stdout,'(a3, 3f18.10)') atm_3d(ia), tau_3d(1,ia), tau_3d(2,ia), tau_3d(3,ia)
ENDDO
!
!  Count how many types we have and how they are called
!
ntyp=1
atm_typ(1)=atm_3d(1)
ityp(1)=1
DO ia=2, nat_3d
   found=0
   DO nt=1, ntyp 
      IF ( TRIM( atm_3d(ia) ) == TRIM( atm_typ(nt) ) ) THEN
         ityp(ia)=nt
         found=1
      END IF
   ENDDO
   IF (found==0) THEN
      ntyp = ntyp + 1
      atm_typ(ntyp) = atm_3d(ia)
      ityp(ia)=ntyp
   END IF
END DO

IF (lcryst==2) THEN
!
!   find all the equivalent positions
!
   ALLOCATE(extfor(3,nat_3d))
   ALLOCATE(if_pos(3,nat_3d))
   extfor=0.0_DP
   if_pos=0
   CALL sup_spacegroup(tau_3d,ityp,extfor,if_pos,space_group_number,nat_3d, &
              uniqueb,rhombohedral,origin_choice,ibrav_3d)
   nat_3d=nattot
   DEALLOCATE(tau_3d)
   DEALLOCATE(ityp)
   DEALLOCATE(atm_3d)
   ALLOCATE(tau_3d(3,nat_3d))
   ALLOCATE(ityp(nat_3d))
   ALLOCATE(atm_3d(nat_3d))
   tau_3d(:,:)=tautot(:,:)
   ityp(:) = ityptot(:)
   DO ia=1,nat_3d
      atm_3d(ia) = atm_typ(ityp(ia))
   END DO
   CALL clean_spacegroup()
   DEALLOCATE(extfor)
   DEALLOCATE(if_pos)
ENDIF



IF (ibrav_3d == 4) THEN
   WRITE(stdout,'("Three (.TRUE.) or four (.FALSE.) indices  ?")') 
   READ(5,*) three_indices
   WRITE(stdout,'(l5)') three_indices
   IF (three_indices) THEN
       WRITE(stdout,'("Crystal coordinates of the &
                             &G vector m b1 + n b2 + o b3 (m,n,o)?")')
       READ(5,*) m, n, o
       WRITE(stdout,'(3i5)') m,n,o
   ELSE
       WRITE(stdout,'("Crystal coordinates of the &
                             &G vector m b1 + n b2 + o b3 (m,n,h,o) h=-m-n?")')
       READ(5,*) m, n, h, o
       WRITE(stdout,'(4i5)') m,n,h,o
       IF (h /= -m-n) CALL errore('gener_3d_slab','h must be equal to -m-n', 1)
   ENDIF
ELSE
   WRITE(stdout,'("Crystal coordinates of the &
                             &G vector m b1 + n b2 + o b3 (m,n,o)?")')
   READ(5,*) m, n, o
   WRITE(stdout,'(3i5)') m,n,o
ENDIF

WRITE(stdout,'("Number of layers ?")')
READ(5,*) nlayers
WRITE(stdout,'(i5)') nlayers

WRITE(stdout,'("Transformation matrix ? (t11, t12, t21, t22) ")')
READ(5,*) t11, t12, t21, t22
WRITE(stdout,'(2i5)') t11, t12
WRITE(stdout,'(2i5)') t21, t22

WRITE(stdout,'("Exact vacuum (.TRUE.) or row distance multiples (.FALSE.)?")')
READ(5,*) ldist_vacuum
WRITE(stdout,*) ldist_vacuum

WRITE(stdout,'("Vacuum space in a.u. ?")')
READ(5,*) vacuum
WRITE(stdout,'(f15.6)') vacuum 

WRITE(stdout,'("In which layer do you want to put the origin ?")')
READ(5,*) origin_shift
WRITE(stdout,'(i5)') origin_shift

WRITE(stdout,'("Output file name?")')
READ(5,*) filename
WRITE(stdout,'(a)') TRIM(filename)

CALL latgen(ibrav_3d, celldm, at(1,1), at(1,2), at(1,3), omega)

at=at/alat

WRITE(stdout,'("Direct lattice vectors")')
WRITE(stdout,'("(",2(f15.6,","),f15.6,")")') at(:,1)
WRITE(stdout,'("(",2(f15.6,","),f15.6,")")') at(:,2)
WRITE(stdout,'("(",2(f15.6,","),f15.6,")")') at(:,3)

CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )

WRITE(stdout,'("Reciprocal lattice vectors")')
WRITE(stdout,'("(",2(f15.6,","),f15.6,")")') bg(:,1)
WRITE(stdout,'("(",2(f15.6,","),f15.6,")")') bg(:,2)
WRITE(stdout,'("(",2(f15.6,","),f15.6,")")') bg(:,3)

IF (lcryst==1.OR.lcryst==2) THEN
!
!  atomic coordinates in crystal axis, bring them to cartesian axis
!
   CALL cryst_to_cart( nat_3d, tau_3d, at, 1 )
END IF

IF (ibrav_3d == 2)  THEN
!
!  With a face centered cubic lattice m, n, and o refer to the cubic
!  reciprocal lattice vectors
!
   beff=0.0_DP
   beff(1,1) = 2.0_DP
   beff(2,2) = 2.0_DP
   beff(3,3) = 2.0_DP
   m1 = m
   n1 = n
   o1 = o
   g(:) = m1 * beff(:,1) + n1 * beff(:,2) + o1 * beff(:,3)
!
!  here we return to the m and n of the face centered lattice
!
   m = -m1 + o1
   n =  n1 + o1
   o = -m1 + n1 

   CALL remove_common_factors(m,n,o,fact)
   g(:)=g(:)/fact
!
ELSEIF (ibrav_3d == 3)  THEN
!
!  With a body centered cubic lattice m, n, and o refer to the cubic
!  reciprocal lattice vectors
!
   beff=0.0_DP
   beff(1,1) = 2.0_DP
   beff(2,2) = 2.0_DP
   beff(3,3) = 2.0_DP
   m1 = m
   n1 = n
   o1 = o
   g(:) = m1 * beff(:,1) + n1 * beff(:,2) + o1 * beff(:,3)
!
!  here we return to the m, n, and o of the face centered lattice
!
   m =  m1 + n1 + o1
   n = -m1 + n1 + o1
   o = -m1 - n1 + o1 

   CALL remove_common_factors(m,n,o,fact)
   g(:)=g(:)/fact
ELSEIF (ibrav_3d == 7)  THEN
!
!  With a centered tetragonal lattice m, n, and o refer to the tetragonal
!  reciprocal lattice vectors
!
   beff=0.0_DP
   beff(1,1) = 2.0_DP
   beff(2,2) = 2.0_DP
   beff(3,3) = 2.0_DP / celldm(3)
   m1 = m
   n1 = n
   o1 = o
   g(:) = m1 * beff(:,1) + n1 * beff(:,2) + o1 * beff(:,3)
!
!  here we return to the m, n, and o of the centered tetragonal lattice
!
   m =  m1 - n1 + o1
   n =  m1 + n1 + o1
   o = -m1 - n1 + o1 

   CALL remove_common_factors(m,n,o,fact)
   g(:)=g(:)/fact
ELSEIF (ibrav_3d == 9)  THEN
!
!  With one-base (C) centered orthorhombic lattice m, n, and o refer to 
!  the primitive orthorhombic reciprocal lattice vectors
!
   beff=0.0_DP
   beff(1,1) = 2.0_DP
   beff(2,2) = 2.0_DP / celldm(2)
   beff(3,3) = 1.0_DP / celldm(3)
   m1 = m
   n1 = n
   o1 = o
   g(:) = m1 * beff(:,1) + n1 * beff(:,2) + o1 * beff(:,3)
!
!  here we return to the m, n, and o of the centered tetragonal lattice
!
   m =  m1 + n1 
   n =  m1 - n1 
   o =  o1 

   CALL remove_common_factors(m,n,o,fact)
   g(:)=g(:)/fact
ELSEIF (ibrav_3d == 91)  THEN
!
!  With one-base (A) centered orthorhombic lattice m, n, and o refer to 
!  the primitive orthorhombic reciprocal lattice vectors
!
   beff=0.0_DP
   beff(1,1) = 1.0_DP
   beff(2,2) = 2.0_DP / celldm(2)
   beff(3,3) = 2.0_DP / celldm(3)
   m1 = m
   n1 = n
   o1 = o
   g(:) = m1 * beff(:,1) + n1 * beff(:,2) + o1 * beff(:,3)
!
!  here we return to the m, n, and o of the centered tetragonal lattice
!
   m =  m1 
   n =  n1 - o1 
   o =  n1 + o1
   CALL remove_common_factors(m,n,o,fact)
   g(:)=g(:)/fact

ELSEIF (ibrav_3d == 10)  THEN
!
!  With face centered orthorhombic lattice m, n, and o refer to 
!  the primitive orthorhombic reciprocal lattice vectors
!
   beff=0.0_DP
   beff(1,1) = 2.0_DP
   beff(2,2) = 2.0_DP / celldm(2)
   beff(3,3) = 2.0_DP / celldm(3)
   m1 = m
   n1 = n
   o1 = o
   g(:) = m1 * beff(:,1) + n1 * beff(:,2) + o1 * beff(:,3)
!
!  here we return to the m, n, and o of the centered orthorhombic lattice
!
   m = m1 + o1  
   n = m1 + n1  
   o = n1 + o1 

   CALL remove_common_factors(m,n,o,fact)
   g(:)=g(:)/fact
ELSEIF (ibrav_3d == 11)  THEN
!
!  With body centered orthorhombic lattice m, n, and o refer to 
!  the primitive orthorhombic reciprocal lattice vectors
!
   beff=0.0_DP
   beff(1,1) = 2.0_DP
   beff(2,2) = 2.0_DP / celldm(2)
   beff(3,3) = 2.0_DP / celldm(3)
   m1 = m
   n1 = n
   o1 = o
   g(:) = m1 * beff(:,1) + n1 * beff(:,2) + o1 * beff(:,3)
!
!  here we return to the m, n, and o of the body centered orthorhombic lattice
!
   m =  m1 + n1 + o1  
   n = -m1 + n1 + o1 
   o = -m1 - n1 + o1 

   CALL remove_common_factors(m,n,o,fact)
   g(:)=g(:)/fact
ELSEIF (ibrav_3d == 13)  THEN
!
!  With one base centered monoclinic lattice m, n, and o refer to 
!  the primitive monoclinic reciprocal lattice vectors
!
   beff=0.0_DP
   beff(:,1) =   bg(:,1) + bg(:,3)
   beff(:,2) =   bg(:,2)
   beff(:,3) = - bg(:,1) + bg(:,3)
   m1 = m
   n1 = n
   o1 = o
   g(:) = m1 * beff(:,1) + n1 * beff(:,2) + o1 * beff(:,3)
!
!  here we return to the m, n, and o of the body centered orthorhombic lattice
!
   m = m1 - o1  
   n = n1  
   o = m1 + o1 

   CALL remove_common_factors(m,n,o,fact)
   g(:)=g(:)/fact
ELSEIF (ibrav_3d == -13)  THEN
!
!  With one base centered monoclinic lattice m, n, and o refer to 
!  the primitive monoclinic reciprocal lattice vectors
!
   beff=0.0_DP
   beff(:,1) =   bg(:,1) + bg(:,2)
   beff(:,2) = - bg(:,1) + bg(:,2)
   beff(:,3) =   bg(:,3)
   m1 = m
   n1 = n
   o1 = o
   g(:) = m1 * beff(:,1) + n1 * beff(:,2) + o1 * beff(:,3)
!
!  here we return to the m, n, and o of the body centered orthorhombic lattice
!
   m = m1 - n1  
   n = m1 + n1  
   o = o1 

   CALL remove_common_factors(m,n,o,fact)
   g(:)=g(:)/fact
ELSE
   g(:) = m * bg(:,1) + n * bg(:,2) + o * bg(:,3)
ENDIF
gmod = SQRT( g(1)**2 + g(2)**2 + g(3)**2 )

WRITE(stdout,'("G vector cartesian", 3f18.7)') g(1), g(2), g(3) 
gc=g
CALL cryst_to_cart(1,gc,at,-1)
WRITE(stdout,'("G vector cryst", 3f18.7)') gc(1), gc(2), gc(3) 

IF (.NOT.(ldist_vacuum)) THEN
!
!  in this case the vacuum is a multiple of the layer distance sufficient to
!  contain a vacuum space equal to vacuum
!
   nspace= INT(vacuum * gmod / alat) + 1
   vacuum= nspace * alat / gmod
END IF

IF ( MOD(nlayers,2) == 0 ) THEN
!
! Put the origin of the slabs with an even number of layers in the center of
! two layers
!
   dz(:) = - (0.5_DP + origin_shift) * g(:) / gmod ** 2
ELSE
   dz(:) = - origin_shift * g(:) / gmod ** 2
ENDIF

prod1 = g(1) * bg(1,1) + g(2) * bg(2,1) + g(3) * bg(3,1)  
prod2 = g(1) * bg(1,2) + g(2) * bg(2,2) + g(3) * bg(3,2)  
prod3 = g(1) * bg(1,3) + g(2) * bg(2,3) + g(3) * bg(3,3)  
!
!  first find two Bravais lattice vectors on the plane passing through 
!  the origin
!
IF ( o /= 0) THEN
   found=0
   o1=ABS(o)
   DO itry=0, o1
      DO jtry = 0, o1
         DO sign1=1,-1,-2
            DO sign2=1,-1,-2
               IF (MOD(-sign1*itry*m - sign2*jtry*n, o1)==0) THEN
                  found=found+1
                  ps(found) = sign1 * itry
                  qs(found) = sign2 * jtry
                  ss(found) = ( -sign1*itry*m - sign2*jtry*n ) / o
               END IF
            END DO
         END DO
      END DO
   END DO
ELSEIF ( n /= 0) THEN
   found=0
   n1=ABS(n)
   DO itry=0, n1
      DO jtry = 0, n1
         DO sign1=1,-1,-2
            DO sign2=1,-1,-2
               IF (MOD(-sign1*itry*m - sign2*jtry*o, n1)==0) THEN
                  found=found+1
                  ps(found) = sign1 * itry
                  qs(found) = ( -sign1*itry*m - sign2*jtry*o ) / n
                  ss(found) = sign2 * jtry
               END IF
            END DO
         END DO
      END DO
   END DO
ELSEIF ( m /= 0) THEN
   found=0
   m1 = ABS(m)
   DO itry=0, m1
      DO jtry = 0, m1
         DO sign1=1,-1,-2
            DO sign2=1,-1,-2
               IF (MOD( -sign1 * itry * n - sign2 * jtry * o, m1)==0) THEN
                  found=found+1
                  ps(found) = ( -sign1*itry*n - sign2*jtry*o ) / m
                  qs(found) = sign1 * itry
                  ss(found) = sign2 * jtry
               END IF
            END DO
         END DO
      END DO
   END DO
END IF

IF ( found == 0  ) THEN
   WRITE(stdout,*) 'p,q, and s not found'
   STOP
END IF

dist=1.d20
DO ia=1, found
   tau(:) = ps(ia) * at(:,1) + qs(ia) * at(:,2) + ss(ia) * at(:,3)
   dist1= SQRT( ( g(2) * tau(3) - g(3) * tau(2) ) ** 2 + &
                ( g(3) * tau(1) - g(1) * tau(3) ) ** 2 + &
                ( g(1) * tau(2) - g(2) * tau(1) ) ** 2 ) / gmod
   
   IF (dist1 < dist .AND. dist1>eps) THEN
      p0=ps(ia)
      q0=qs(ia)
      s0=ss(ia)
      dist=dist1
      ia_first=ia
   ENDIF
END DO

!WRITE(stdout,*) 'p0, q0, and s0', p0, q0, s0

c1(:) = p0 * at(:,1) + q0 * at(:,2) + s0 * at(:,3)

dist=1.d20
DO ia=1, found
   IF ( ia /= ia_first ) THEN
      tau(:) = ps(ia) * at(:,1) + qs(ia) * at(:,2) + ss(ia) * at(:,3)
      CALL vector_prod(g, tau, vprod)
      dist1 = SQRT( vprod(1) ** 2 + vprod(2) ** 2 + vprod(3) ** 2 ) / gmod
      CALL vector_prod(c1, tau, vprod)
      prod = vprod(1) ** 2 + vprod(2) ** 2 + vprod(3) ** 2
!
!   prod > 0 removes the vectors that are parallel to c1
!
      IF (dist1 < dist .AND. prod > eps) THEN
         p01=ps(ia)
         q01=qs(ia)
         s01=ss(ia)
         dist=dist1
      END IF
   END IF
END DO

!WRITE(stdout,*) 'p01, q01, and s01', p01, q01, s01

min_j = 1000
nat=0
DO j=-(nlayers-1)/2+origin_shift, nlayers/2 + origin_shift
   IF ( o /= 0 ) THEN
      p1 = NINT ( j * prod1 / gmod / gmod )
      q1 = NINT ( j * prod2 / gmod / gmod )
      found=0
      o1=ABS(o)
      DO itry=0, o1
         DO jtry = 0, o1
            DO sign1=-1,1,2
               DO sign2=-1,1,2
                  IF (MOD(-(p1+sign1*itry)*m-(q1+sign2*jtry)*n+j,o1)==0) THEN
                     found=found+1
                     ps(found) = p1 + sign1 * itry
                     qs(found) = q1 + sign2 * jtry
                     ss(found) = ( - ( p1 + sign1 * itry ) * m   &
                                   - ( q1 + sign2 * jtry ) * n + j ) / o 
                  END IF
               END DO
            END DO
         END DO
      END DO
   ELSE IF ( n /= 0) THEN
      p1 = NINT ( j * prod1 / gmod / gmod )
      s1 = NINT ( j * prod3 / gmod / gmod )
      found=0
      n1=ABS(n)
      DO itry=0, n1
         DO jtry = 0, n1
            DO sign1=-1,1,2
               DO sign2=-1,1,2
                  IF (MOD(-(p1+sign1*itry)*m-(s1+sign2*jtry)*o+j,n1)==0) THEN
                     found=found+1
                     ps(found) = p1 + sign1 * itry
                     qs(found) = ( - ( p1 + sign1 * itry ) * m  &
                                   - ( s1 + sign2 * jtry ) * o + j ) / n
                     ss(found) = s1 + sign2 * jtry
                  END IF
               END DO
            END DO
         END DO
      END DO
   ELSE IF ( m /= 0) THEN
      q1 = NINT ( j * prod2 / gmod / gmod )
      s1 = NINT ( j * prod3 / gmod / gmod )
      found=0
      m1=ABS(m)
      DO itry=0, m1
         DO jtry = 0, m1
            DO sign1=-1,1,2
               DO sign2=-1,1,2
                  IF (MOD(-(q1+sign1*itry)*n-(s1+sign2*jtry)*o+j,m1)==0) THEN
                     found=found+1
                     ps(found) = ( - ( q1 + sign1 * itry ) * n  &
                                   - ( s1 + sign2 * jtry ) * o + j ) / m
                     qs(found) = q1 + sign1 * itry
                     ss(found) = s1 + sign2 * jtry
                  END IF
               END DO
            END DO
         END DO
      END DO
   ELSE
      CALL errore('gener_3d_slab','m,n, and o all zero?',1)
   END IF

   IF ( found == 0  ) THEN
      WRITE(stdout,*) 'p,q, and s not found'
      STOP
   END IF
!
!   compute the distance for all candidates p,q and s and take the shortest
!
   dist=1.d20
   DO ia=1, found
      tau(:) = ps(ia) * at(:,1) + qs(ia) * at(:,2) + ss(ia) * at(:,3)
      CALL vector_prod(g,tau,vprod)
      dist1= SQRT( vprod(1) ** 2 + vprod(2) ** 2 + vprod(3) ** 2 ) / gmod
      IF (dist1 < dist) THEN
         p=ps(ia)
         q=qs(ia)
         s=ss(ia)
         dist=dist1
      ENDIF
   ENDDO

   IF (dist < eps) THEN
      IF (ABS(j) < min_j .AND. j /= 0) min_j=ABS(j)
   ENDIF

!   WRITE(stdout,*) 'p, q, and s', p, q, s

   DO ia=1,nat_3d
      nat = nat + 1
      y(:,nat) = p * at(:,1) + q * at(:,2) + s * at(:,3) + tau_3d(:,ia)
      atm(nat) = atm_3d(ia)
      ityp_all(nat) = ityp(ia)
   END DO
END DO

IF (min_j < 1000) THEN
   WRITE(stdout,'("In this direction the bulk has lattice planes of",i5,&
                                                  " types")') min_j
ELSE
   WRITE(stdout,'("Vanishing distance not found")')
   IF (ibrav_3d > 3) THEN
      WRITE(stdout,'("Either the bulk cannot be obtained in this direction")')
      WRITE(stdout,'("or it requires more than",i5," lattice planes")') (nlayers)/2
   ELSE
      WRITE(stdout,'("The bulk has more than",i5," lattice planes")') (nlayers)/2
   END IF 
END IF
!
!  c1 has been calculated above
!
c2(:) = p01 * at(:,1) + q01 * at(:,2) + s01 * at(:,3)
t(:) = ( nlayers * alat / gmod + vacuum ) * g(:) / gmod / alat

c1mod = SQRT ( c1(1)**2 + c1(2)**2 + c1(3)**2 )
c2mod = SQRT ( c2(1)**2 + c2(2)**2 + c2(3)**2 )
tmod  = SQRT ( t(1)**2  + t(2)**2  + t(3)**2 )
!
!  If c1, c2, and t do not have the same orientation as x, y, and z
!  change the sign of c2
!
det = c1(1) * ( c2(2) * t(3) - c2(3) * t(2) ) -   &
      c2(1) * ( c1(2) * t(3) - c1(3) * t(2) ) +   &
      t(1)  * ( c1(2) * c2(3) - c1(3) * c2(2) )


IF (det < 0.0_DP) c2(:) = - c2(:)
!
!  If the angle between c1 and c2 is 60 degrees and c1 and c2 have the same
!  modulus, take c1 and c2-c1 that have an angle of 120 as the hexagonal
!  lattice in latgen
!
prod1 = ( c1(1) * c2(1) + c1(2) * c2(2) + c1(3) * c2(3) ) / c1mod / c2mod
IF ( ABS( c1mod - c2mod ) < eps .AND. ABS( prod1 - 0.5_DP ) < eps ) THEN
   c2(:) = c2(:) - c1(:)
   c2mod = SQRT ( c2(1)**2 + c2(2)**2 + c2(3)**2 )
   IF (ABS( c1mod - c2mod ) > eps) CALL errore('gener_3d_slab',&
                                               'some problem with c2mod',1)
   prod1 = ( c1(1) * c2(1) + c1(2) * c2(2) + c1(3) * c2(3) ) / c1mod / c2mod
END IF
!
!  If c1 and c2 are not orthogonal, check if c1 and a linear combination of
!  c1 and c2 have the same area and are orthogonal. Only a small number of
!  combinations is tried. In that case take the linear combination as the
!  new c2
!
IF (ABS(prod1) > 1.D-9) THEN
   prod1 =   c1(1) * c2(1) + c1(2) * c2(2) + c1(3) * c2(3) 
   DO i1=-8,8
      scalar= prod1 + i1 * c1mod ** 2
      IF (ABS(scalar)< 1.D-9) THEN
         c2=i1 * c1 + c2 
         EXIT
      ENDIF 
   ENDDO
ENDIF
!
!  write all the periodic copies on the surface if the user asked for 
!  a supercell
!
d1(:) = t11 * c1(:) + t12 * c2(:)
d2(:) = t21 * c1(:) + t22 * c2(:)

CALL recips( d1, d2, t, bd1, bd2, bt )

iat=nat
copies = 1
DO m1 = MIN(0, t11, t21, t11+t21), MAX(0, t11, t21, t11+t21)
   DO n1 = MIN(0, t12, t22, t12+t22), MAX(0, t12, t22, t12+t22)
      IF ( m1==0 .AND. n1==0 ) CYCLE
      tau(:) = m1 * c1(:) + n1 * c2(:) 
      prod = tau(1) * bd1(1) + tau(2) * bd1(2) + tau(3) * bd1(3)
      prod1 = tau(1) * bd2(1) + tau(2) * bd2(2) + tau(3) * bd2(3)
      IF (prod >=-eps .AND. prod < 1.0_DP-eps .AND. &
                            prod1 >=-eps .AND. prod1<1.0_DP-eps) THEN
         DO ia=1,nat
            iat=iat+1
            y(:,iat) = y(:,ia) + tau(:)
            atm(iat) = atm(ia)
            ityp_all(iat)= ityp(ia)
         END DO
         copies=copies + 1
      END IF
   END DO
END DO
!
nat=iat
!
CALL vector_prod(c1,c2,vprod)
area = SQRT( vprod(1)**2 + vprod(2)**2 + vprod(3)**2 )
CALL vector_prod(d1,d2,vprod)
super_area = SQRT( vprod(1)**2 + vprod(2)**2 + vprod(3)**2 )

IF ( NINT( super_area / area ) - copies /= 0 ) &
   CALL errore('gener_3d_slab','Same problem with the copies',1)

c1 = d1 
c2 = d2

c1mod = SQRT ( c1(1)**2 + c1(2)**2 + c1(3)**2 )
c2mod = SQRT ( c2(1)**2 + c2(2)**2 + c2(3)**2 )
prod1 = ( c1(1) * c2(1) + c1(2) * c2(2) + c1(3) * c2(3) ) / c1mod / c2mod

CALL recips(c1,c2,t,bc1,bc2,bt)

IF (ionode) THEN
   iuout=find_free_unit()
   OPEN(unit=iuout, file=TRIM(filename), status='unknown', form='formatted')
ENDIF

celldm_sur=0.0_DP
IF (ABS(prod1) < eps) THEN
!
!  c1 and c2 form an angle of 90 degrees. If their modulus is the same
!  use a simple tetragonal cell, 
!
   IF (ABS(c1mod - c2mod) < eps) THEN
      c_alat = c1mod * alat
      celldm_sur(1) = c_alat
      celldm_sur(3) = tmod/c1mod
      ibrav_sur=6
      IF (ionode) THEN
         WRITE (iuout, '("ibrav=6")')
         WRITE (iuout, '("celldm(1)= ",f15.8)') celldm_sur(1)
         WRITE (iuout, '("celldm(3)= ",f15.8)') celldm_sur(3)
      ENDIF
   ELSE
!
!  c1 and c2 have different modulus, use a simple tetragonal cell   
!
      c_alat = c1mod * alat
      celldm_sur(1) = c_alat
      celldm_sur(2) = c2mod / c1mod
      celldm_sur(3) = tmod / c1mod
      ibrav_sur=8
      IF (ionode) THEN
         WRITE (iuout, '("ibrav=8")')
         WRITE (iuout, '("celldm(1)= ",f15.8)') c1mod * alat
         WRITE (iuout, '("celldm(2)= ",f15.8)') c2mod / c1mod
         WRITE (iuout, '("celldm(3)= ",f15.8)') tmod / c1mod
      END IF
   END IF
ELSEIF (ABS(c1mod-c2mod) < eps .AND. ABS(prod1 + 0.5_DP) < eps ) THEN
!
!  The angle is 120 degrees and the vectors are equal use the hexagonal lattice
!
   c_alat = c1mod * alat
   celldm_sur(1) = c_alat
   celldm_sur(3) = tmod / c1mod
   ibrav_sur=4
   IF (ionode) THEN
      WRITE (iuout, '("ibrav=4")')
      WRITE (iuout, '("celldm(1)= ",f15.8)') c1mod * alat
      WRITE (iuout, '("celldm(3)= ",f15.8)') tmod / c1mod
   ENDIF
ELSEIF (ABS(c1mod-c2mod) < eps ) THEN
!
!  The angle is neither 90 nor 120, but the vectors have equal length, 
!  use the one-base (C) centered lattice
!
   alpha=ACOS(prod1) * 0.5_DP
   c_alat= 2.0_DP * c1mod * SIN(alpha) * alat
   celldm_sur(1) = c_alat
   celldm_sur(2) = 2.0_DP*c1mod*COS(alpha)*alat / c_alat
   celldm_sur(3) = tmod * alat / c_alat
   ibrav_sur=9
   IF (ionode) THEN
      WRITE (iuout, '("ibrav=9")')
      WRITE (iuout, '("celldm(1)= ",f15.8)') c_alat
      WRITE (iuout, '("celldm(2)= ",f15.8)') 2.0_DP*c1mod*COS(alpha)*alat &
                                                                    / c_alat
      WRITE (iuout, '("celldm(3)= ",f15.8)') tmod * alat / c_alat
   ENDIF
ELSE
!
!  c1 and c2 form an angle different from 90 and 120 degrees or have different
!  lenghts, use a  monoclinic cell (c unique)
!
   c_alat = c1mod * alat
   celldm_sur(1) = c_alat
   celldm_sur(2) = c2mod / c1mod
   celldm_sur(3) = tmod / c1mod
   celldm_sur(4) = prod1
   ibrav_sur=12
   IF (ionode) THEN
      WRITE (iuout, '("ibrav=12")')
      WRITE (iuout, '("celldm(1)= ",f15.8)') c1mod * alat
      WRITE (iuout, '("celldm(2)= ",f15.8)') c2mod / c1mod
      WRITE (iuout, '("celldm(3)= ",f15.8)') tmod / c1mod
      WRITE (iuout, '("celldm(4)= ",f15.8)') prod1
   ENDIF
END IF

IF (ionode) THEN
   WRITE (iuout, '("nat= ",i5)') nat
   WRITE (iuout, '("ATOMIC_POSITIONS {crystal}")') 
ENDIF

!
!  project dz in the t direction
!
shiftz = dz(1) * bt(1) + dz(2) * bt(2) + dz(3) * bt(3)
!
ALLOCATE(tau_sur(3,nat))
DO na=1,nat
   prod1 = y(1,na) * bc1(1) + y(2,na) * bc1(2) + y(3,na) * bc1(3)
   prod2 = y(1,na) * bc2(1) + y(2,na) * bc2(2) + y(3,na) * bc2(3)
   prod3 = y(1,na) * bt(1) + y(2,na) * bt(2) + y(3,na) * bt(3)
   tau_sur(1,na) = prod1
   tau_sur(2,na) = prod2
   tau_sur(3,na) = prod3 + shiftz
   IF (ionode) &
      WRITE (iuout,'(a,3f20.13,"  0  0  0")') atm(na), prod1, prod2, &
                                              prod3 + shiftz
ENDDO
IF (ionode) CLOSE(iuout)
!
!  write the coordinate for xcrydens with the original orientation 
!  of the surface
!
IF (ionode) THEN
   OPEN(unit=iuout, file=TRIM(filename)//'.or.xsf', status='unknown', &
                                              form='formatted')
   at(:,1)=c1(:)
   at(:,2)=c2(:)
   at(:,3)=t(:)
   CALL xsf_struct (alat, at, nat, y, atm_typ, ityp_all, iuout)

   CLOSE(iuout)
ENDIF
!
!  rotate the coordinates so that the normal to the surface is along the
!  z axis. 
!
IF (ionode) &
OPEN(unit=iuout, file=TRIM(filename)//'.xsf', status='unknown', &
                                              form='formatted')

CALL latgen(ibrav_sur, celldm_sur, at(1,1), at(1,2), at(1,3), omega)
at=at/c_alat
CALL cryst_to_cart( nat, tau_sur, at, 1 )

CALL xsf_struct (c_alat, at, nat, tau_sur, atm_typ, ityp_all, iuout)

IF (ionode) CLOSE(iuout)

CALL environment_end( code )
CALL mp_global_end ()

END PROGRAM gener_3d_slab

SUBROUTINE remove_common_factors(m,n,o,fact)

IMPLICIT NONE
INTEGER, INTENT(INOUT) :: m, n, o
INTEGER, INTENT(OUT) :: fact

fact=1
IF (MOD(m,2)==0.AND.MOD(n,2)==0.AND.MOD(o,2)==0) THEN
   m = m/2
   n = n/2
   o = o/2
   fact=fact*2
END IF

IF (MOD(m,3)==0.AND.MOD(n,3)==0.AND.MOD(o,3)==0) THEN
   m = m/3
   n = n/3
   o = o/3
   fact=fact*3
END IF

IF (MOD(m,4)==0.AND.MOD(n,4)==0.AND.MOD(o,4)==0) THEN
   m = m/4
   n = n/4
   o = o/4
   fact=fact*4
END IF

IF (MOD(m,5)==0.AND.MOD(n,5)==0.AND.MOD(o,5)==0) THEN
   m = m/5
   n = n/5
   o = o/5
   fact=fact*5
END IF

IF (MOD(m,7)==0.AND.MOD(n,7)==0.AND.MOD(o,7)==0) THEN
   m = m/7
   n = n/7
   o = o/7
   fact=fact*7
END IF

IF (MOD(m,9)==0.AND.MOD(n,9)==0.AND.MOD(o,9)==0) THEN
   m = m/9
   n = n/9
   o = o/9
   fact=fact*9
END IF

RETURN
END SUBROUTINE remove_common_factors

SUBROUTINE vector_prod(a,b,c)

USE kinds, ONLY : DP

IMPLICIT NONE
REAL(DP), INTENT(IN) :: a(3), b(3)
REAL(DP), INTENT(OUT) :: c(3)

c(1) = a(2) * b(3) - a(3) * b(2)
c(2) = a(3) * b(1) - a(1) * b(3)
c(3) = a(1) * b(2) - a(2) * b(1)

RETURN
END SUBROUTINE vector_prod
