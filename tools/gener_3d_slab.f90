!
! Copyright (C) 2014-2018 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM gener_3d_slab
!
!  This program computes all the atomic positions of a slab of a bulk solid
!  defined by a three dimensional (3D) Bravais lattice and nat_3d atoms in
!  the unit cell. Calling a1, a2, and a3 the three primitive vectors of the
!  Bravais lattice of the bulk, it prints on output the information necessary, 
!  using periodic boundary conditions, to simulate a slab with surfaces
!  perpendicular to a G vector: G = m b1 + n b2 + o b3 where m, n, and o, 
!  should have no common factors. 
!  For centered lattices b1, b2, and b3 are the reciprocal lattice vectors
!  of the conventional lattice.
!  The number of layers of the slab, nlayers, is given in input
!  and the size of the unit cell on the surface is given by a 2x2 matrix 
!  of integers t11, t12, t21, t22 that define a super-cell of the smallest
!  unit cell of the slab. The latter is found by the code. 
!  If N_sur is the number of atoms per surface plane defined by the above 
!  integers, the program writes on output the coordinates of N_sur x nlayers 
!  atoms.
!  The distance between repeated copies of the slab can be indicated in 
!  input (vacuum) in a.u.. The program puts that distance if
!  ldist_vacuum=.TRUE. or a vacuum space equivalent to the minimum integer 
!  number of empty layers needed to contain the vacuum distance
!  if ldist_vacuum=.FALSE..
!  The atomic coordinates are written on output and periodic boundary
!  conditions are assumed in the xy plane while repeated slabs are obtained
!  along z. The output file, whose name (filename) is given in input, contains 
!  the size of the slab unit cell and the atomic coordinates in a format 
!  suited to write the pw.x input.
!  The program produces also two .xsf files for plotting the slab using
!  the code XCrysDen. The file filename.xsf plots the slab with the 
!  orientation that it will have in the pw.x simulation (in this reference 
!  system  the input G vector is parallel to z), while filename.or.xsf 
!  plots the slab exactly as cut in the solid, so the normal to the surface 
!  will be parallel to the input G vector.
!
!  List of input variables:
!
!  ibrav_3d     The Bravais lattice of the bulk, defined as in the pw.x input. 
!               The latgen routine of QE is used to generate the 3D primitive
!               lattice vectors.
!  
!  celldm_3d(6) The size of the unit cell as in pw.x input.
!
!  lcryst       1 atomic positions in the crystal basis. 
!               2 atomic positions of inequivalent atoms in crystal 
!                 coordinates as in pw.x input when the space group number
!                 and the crystal_sg option are used.
!               3 atomic positions in cartesian coordinates.
!  IF lcryst==2
!    
!     space_group_number  The number of the space group as in the ITA tables
!
!     uniqueb         !
!     rhombohedral    ! See the pw.x input for an explanation of these 
!                     ! variables
!     origin_choice   !
!  ENDIF
!
!  nat_3d       The number of atoms inside the unit cell
!
!  atm_3d, tau_3d(3)   label and coordinates of the atoms inside the unit cell.
!
!  IF (ibrav_3d==4) THEN
!
!     three_indices    if true the input G vector has only three indices
!     IF (three indices) THEN
!        m,n,o        m, n, and o of the G vector.
!     ELSE
!         m, n, -m-n, o
!     ENDIF
!  ELSE
!     m, n, o         m, n, and o of the G vector.
!  ENDIF
!
!  nlayers      number of layers.
!               
!  t11, t12, t21, t22  where tij are integer numbers. The size of the 
!               surface unit cell is defined by the two primitive vectors 
!               A'_1 = t11 A_1 + t12 A_2, 
!               A'_2 = t21 A_1 + t22 A_2 
!               where A_1 and A_2 are minimal direct lattice vectors 
!               on the surface plane, identified by the code.  
!               Find them by using 1, 0, 0, 1 for t11, t12, t21, and t22.
!
!  ldist_vacuum if .TRUE. the distance between repeated slab is vacuum,
!               otherwise it is the minimum integer number of empty layers 
!               that contains vacuum.
!
!  vacuum       the vacuum distance. This distance is calculated along z from
!               the boundary of the last layer of the slab to the first
!               periodically repeated 3D Bravais lattice point of the slab
!               
!               --o----   first point of the slab due to periodiciy 
!                   |
!                   |      vacuum  
!                   |
!               -------   boundary of the last layer of the slab
!               --o----   last 3D Bravais lattice point of the slab
!               ---o---   second-last 3D Bravais lattice point of slab
!
!               If the slab has the correct number of layers, setting 
!               vacuum to 0.0 reproduces a bulk with a slab perpendicular 
!               to G. If the number of layer is not correct setting 
!               the vacuum to 0.0 might lead to very close or overlapping 
!               atoms and a cell that cannot be simulated.
!
!  origin_shift Number of the layer that contains the origin.
!               When origin_shift is zero, the first atom of the 
!               central layer is in the origin. You might want to 
!               shift the origin to a different layer. origin_shift
!               is an integer that specifies on which layer to put
!               the origin. In the resulting slab the layer with the 
!               origin will be the central layer. Note that when the 
!               number of layers is even the origin is in the middle 
!               of two layers. The layer indicated by origin_shift will
!               be the one with negative z closest to the origin.  
!               Set origin_shift=0 if you do not know what to do here.
!
!  do_shift     if .TRUE. apply a shift to slabs with an even number
!               of layers so that the origin is in the middle of
!               two slabs.
!
!  filename     the name of a file were the output coordinates and
!               unit cell are written. Note that other two files 
!               filename.xsf and filename.or.xsf will be written as 
!               explained above.
!              
USE kinds,       ONLY : DP
USE constants,   ONLY : pi, bohr_radius_si
USE wyckoff,     ONLY : nattot, tautot, ityptot, sup_spacegroup, &
                        clean_spacegroup
USE lattices,    ONLY : is_centered, compute_conventional
USE atomic_pos,  ONLY : find_ityp
USE io_global,   ONLY : ionode, stdout
USE mp_global,   ONLY : mp_startup, mp_global_end
USE environment, ONLY : environment_start, environment_end

IMPLICIT NONE

INTEGER, PARAMETER :: nmax=100000, ntypx=100
!
!  Input parameters
!
REAL(DP) :: celldm_3d(6)
REAL(DP), ALLOCATABLE :: tau_3d(:,:), extfor(:,:)
INTEGER  :: m, n, h, o, t11, t12, t21, t22, nlayers, ibrav_3d, nat_3d, &
            origin_shift, lcryst, space_group_number, origin_choice 
LOGICAL  :: ldist_vacuum, three_indices, uniqueb, rhombohedral, do_shift
CHARACTER ( LEN=3 ), ALLOCATABLE :: atm_3d(:)
CHARACTER( LEN=256 ) :: filename          
!
!  local variables
!
REAL(DP), PARAMETER :: eps=1.D-6, small=1.D-10

REAL(DP) :: at(3,3), bg(3,3),   &  ! direct and reciprocal 3d lattice
            atp(3,3),beff(3,3), &  ! direct and reciprocal conventional 3d 
            g(3), &                ! normal to the surface
            c1(3), c2(3), t(3),   &  ! primitive vectors of the slab
            bc1(3), bc2(3), bt(3),&! reciprocal vectors of the slab
            d1(3), d2(3), bd1(3), bd2(3), & ! direct and reciprocal vectors of 
                                   ! the super-cell slab
            vprod(3), dz(3),     & ! auxiliary
            tau(3),              & ! auxiliary to store Bravais lattices
            y(3,nmax),           & ! coordinates of all atoms
            celldm_sur(6),       & ! the dimension of the slab 
            gc(3)                  ! the G vector in crystal coordinates

INTEGER ::  ps(nmax), qs(nmax), ss(nmax),& ! coordinates of R in each plane
            ityp_all(nmax)                 ! the types of the slab atoms

REAL(DP) :: alat,                & ! all length are in units of alat
            c1mod, c2mod, tmod,  & ! moduli of c1, c2, t
            gmod,                & ! modulus of G (in units of 2 pi / alat)
            prod1, prod2, prod3, & ! scalar product between G and b1, b2, b3
            vacuum,              & ! the target vacuum length
            dist, dist1,         & ! auxiliary to save distances
            prod,                & ! auxiliary
            det,                 & ! auxiliary deteminant
            omega,               & ! volume of the bulk
            area, super_area,    & ! surface area of the slab and of the 
                                   ! supercell slab
            volume,              & ! volume of a slab cell (for checking)
            pr, qr, sr,          & ! auxiliary crystal coordinates of points
            alpha,               & ! angle between c1 and c2
            c_alat,              & ! the length unit for the slab 
            shiftz,              & ! shift of all z coordinates
            scalar                 ! auxiliary

INTEGER  :: p, q, s,       & ! the indices of the R
            p0, q0, s0,    & ! the indices of c1
            p01, q01, s01, & ! the indices of c2
            p1, q1, s1,    & ! the initial guess of the indices
            m1, n1, o1,    & ! auxiliary for m, n, o
            j,             & ! index of the layer
            min_j,         & ! index of the bulk layer 
            idet,          & ! the number of surface unit cells
            itry, jtry, ktry, & ! grid of trying indices
            sign1, sign2,  & ! extend the grid to negative values
            ibrav_sur,     & ! Bravais lattice index of the slab
            nat, na, iat,  & ! atoms of the slab, indeces on atoms
            ntyp, nt,      & ! number of atoms types, counter on types
            found,         & ! number of atoms found in the plane
            nspace,        & ! number of layers of vacuum
            na_first,      & ! atom of c1
            copies,        & ! number of slabs in the super-cell slab
            fact, ifact,   & ! auxiliary
            find_free_unit, iuout ! iunits auxiliary

INTEGER, ALLOCATABLE :: ityp(:),        & ! type of atoms 
                        if_pos(:,:)       ! if to be kept fixed (in pw.x input)
REAL(DP), ALLOCATABLE :: tau_sur(:,:)     ! crystal coordinates of the atoms of
                                          ! the slab in the basis of the slab at

LOGICAL :: lfound                         ! if true bulk exists for this G
                                          ! direction
CHARACTER ( LEN=3 ) :: atm(nmax),       & ! the atomic labels all atoms
                       atm_typ(ntypx)     ! atomic labels, all types
CHARACTER(LEN=9) :: code='3D_SLAB'        ! code name


CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )
!
!  Reading input
!
WRITE(stdout,'("ibrav_3d? ")')
READ(5,*) ibrav_3d
WRITE(stdout,'(i5)') ibrav_3d
WRITE(stdout,'("celldm ")')
READ(5,*) celldm_3d(1), celldm_3d(2), celldm_3d(3), &
          celldm_3d(4), celldm_3d(5), celldm_3d(6)
WRITE(stdout,'(6f12.6)') celldm_3d
alat=celldm_3d(1)
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
DO na=1,nat_3d
   READ(5,*) atm_3d(na), tau_3d(1,na), tau_3d(2,na), tau_3d(3,na)
   WRITE(stdout,'(a3, 3f18.10)') atm_3d(na), tau_3d(1,na), tau_3d(2,na), tau_3d(3,na)
ENDDO
!
!  Count how many types we have and how they are called
!
CALL find_ityp(nat_3d, atm_3d, ntyp, ityp, atm_typ, ntypx)

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
   DO na=1,nat_3d
      atm_3d(na) = atm_typ(ityp(na))
   END DO
   CALL clean_spacegroup()
   DEALLOCATE(extfor)
   DEALLOCATE(if_pos)
ENDIF
!
!  Read the G vector
!
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
IF (nlayers>nmax) CALL errore('gener_3d_slab','nlayer too large',nmax)

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

WRITE(stdout,'("In which layer do you want to put the origin?")')
READ(5,*) origin_shift
WRITE(stdout,'(i5)') origin_shift

WRITE(stdout,'("Do you want to put the origin within two slabs?")')
READ(5,*) do_shift
WRITE(stdout,'(l5)') do_shift

WRITE(stdout,'("Output file name?")')
READ(5,*) filename
WRITE(stdout,'(a)') TRIM(filename)
!
!  generate the 3D bulk primitive direct and reciprocal lattice vectors.
!
CALL latgen(ibrav_3d, celldm_3d, at(1,1), at(1,2), at(1,3), omega)
at=at/alat

WRITE(stdout,'("Direct lattice vectors")')
WRITE(stdout,'(/,5x,"Direct lattice vectors")')
WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') at(:,1)
WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') at(:,2)
WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') at(:,3)

CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )

WRITE(stdout,'(/,5x,"Reciprocal lattice vectors")')
WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') bg(:,1)
WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') bg(:,2)
WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') bg(:,3)
!
!  atomic coordinates in crystal axis, bring them to cartesian axis
!
IF (lcryst==1.OR.lcryst==2) THEN
   CALL cryst_to_cart( nat_3d, tau_3d, at, 1 )
END IF
!
!  For centered Bravais lattice the input n, m, o refer to the convential
!  cell, transform them into the crystal coordinates of the centered cell
!  Multiply by an integer factor if this G is not a G of the centered lattice
!
IF (is_centered(ibrav_3d)) THEN
   CALL compute_conventional(at, atp, ibrav_3d)

   WRITE(stdout,'(/,5x,"Conventional lattice vectors")')
   WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') atp(:,1)
   WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') atp(:,2)
   WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') atp(:,3)

   CALL recips(atp(1,1), atp(1,2), atp(1,3), beff(1,1), beff(1,2), beff(1,3))

   WRITE(stdout,'(/,5x,"Conventional reciprocal vectors")')
   WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') beff(:,1)
   WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') beff(:,2)
   WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') beff(:,3)
   g(:) = m * beff(:,1) + n * beff(:,2) + o * beff(:,3)
   gc=g
   CALL cryst_to_cart(1, gc, at, -1)
!
!  multiply by a factor if gc are not integer
!
   DO ifact=1,2
      m = NINT(ifact*gc(1))
      n = NINT(ifact*gc(2))
      o = NINT(ifact*gc(3))
      IF (ABS(m-ifact*gc(1))<1.D-10.AND.ABS(n-ifact*gc(2))<1.D-10.AND.&
                              ABS(o-ifact*gc(3))<1.D-10) THEN
         EXIT
      ELSEIF (ifact==2) THEN
         CALL errore('gener_3d_slab','No factor found',1)
      ENDIF
   ENDDO
   CALL remove_common_factors(m,n,o,fact)
   g(:)=g(:)*ifact/fact
ELSE
   g(:) = m * bg(:,1) + n * bg(:,2) + o * bg(:,3)
   beff = bg
ENDIF
gmod = SQRT( g(1)**2 + g(2)**2 + g(3)**2 )
gc=g
CALL cryst_to_cart(1, gc, at, -1)
!
!  Write information on the G vector
!
WRITE(stdout,'(/,5x,"G vector perpendicular to the surface:")')
WRITE(stdout,'(5x,"Cart. Coord. ",3f17.7, " (2 pi/a)")') g(1), g(2), g(3)
WRITE(stdout,'(5x,"Cryst. Coord.", 3f17.7)') gc(1), gc(2), gc(3)

WRITE(stdout,'(/,5x,"Lattice constant", f17.7, " a.u.")') alat
WRITE(stdout,'(5x,"Planar distance", f18.7, " a.u.")') alat / gmod
WRITE(stdout,'(5x,"Planar distance/alat", f13.7)') 1.0_DP / gmod

CALL search_vicinals(g, beff)
!
!  Now predict how many different layers there are in this direction.
!  This is done by finding the indeces of p,q,r of the shortest Bravais 
!  lattice vector R parallel to G and finding the j of the lattice plane 
!  that contains the vector R.
!
prod1 = g(1) * bg(1,1) + g(2) * bg(2,1) + g(3) * bg(3,1)  
prod2 = g(1) * bg(1,2) + g(2) * bg(2,2) + g(3) * bg(3,2)  
prod3 = g(1) * bg(1,3) + g(2) * bg(2,3) + g(3) * bg(3,3)  

lfound=.TRUE.
IF (ABS(prod1)>=ABS(prod2).AND.ABS(prod1)>=ABS(prod3)) THEN
   DO itry=1, 10000
      qr= itry * prod2 / prod1
      sr= itry * prod3 / prod1
      IF (ABS(qr - NINT(qr))<small .AND. ABS(sr - NINT(sr))< small) THEN
         p=itry
         q=NINT(qr)
         s=NINT(sr)
         EXIT
      ENDIF
      IF (itry==10000) THEN
         WRITE(stdout,'(5x,"Bulk plane not found")')
         lfound=.FALSE.
      ENDIF
   ENDDO
ELSEIF (ABS(prod2)>=ABS(prod1).AND.ABS(prod2)>=ABS(prod3)) THEN
   DO jtry=1, 10000
      pr= jtry * prod1 / prod2
      sr= jtry * prod3 / prod2
      IF (ABS(pr - NINT(pr))<small .AND. ABS(sr - NINT(sr))< small) THEN
         p=NINT(pr)
         q=jtry
         s=NINT(sr)
         EXIT
      ENDIF
      IF (jtry==10000) THEN
         WRITE(stdout,'(5x,"Bulk plane not found")')
         lfound=.FALSE.
      ENDIF
   ENDDO
ELSEIF (ABS(prod3)>=ABS(prod2).AND.ABS(prod3)>=ABS(prod1)) THEN
   DO ktry=1, 10000
      pr= ktry * prod1 / prod3
      qr= ktry * prod2 / prod3
      IF (ABS(pr - NINT(pr))<small .AND. ABS(qr - NINT(qr))< small) THEN
         p=NINT(pr)
         q=NINT(qr)
         s=ktry
         EXIT
      ENDIF
      IF (ktry==10000) THEN
         WRITE(stdout,'(5x,"Bulk plane not found")')
         lfound=.FALSE.
      ENDIF
   ENDDO
ENDIF

IF (lfound) THEN
   min_j= m * p + n * q + o * s
!   WRITE(stdout,'(/,5x,"Found R vector parallel to G",3i8)') p, q, s 
   WRITE(stdout,'(/,5x,"In this direction the bulk has lattice planes of",i8,&
                                                  " types")') ABS(min_j)
ENDIF
!
!  Now find two Bravais lattice vectors on the plane passing through 
!  the origin and perpendicular to G. First create a mesh of R vectors
!  in this plane
!
IF ( o /= 0) THEN
   found=0
   o1=ABS(o) * 5
   DO itry=0, o1
      DO jtry = 0, o1
         DO sign1=1,-1,-2
            DO sign2=1,-1,-2
               IF (MOD(-sign1*itry*m - sign2*jtry*n, ABS(o))==0) THEN
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
   n1=ABS(n) * 5
   DO itry=0, n1
      DO jtry = 0, n1
         DO sign1=1,-1,-2
            DO sign2=1,-1,-2
               IF (MOD(-sign1*itry*m - sign2*jtry*o, ABS(n))==0) THEN
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
   m1 = ABS(m) * 5
   DO itry=0, m1
      DO jtry = 0, m1
         DO sign1=1,-1,-2
            DO sign2=1,-1,-2
               IF (MOD( -sign1 * itry * n - sign2 * jtry * o, ABS(m))==0) THEN
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

IF ( found == 0  ) CALL errore('gener_3d_slab', 'p,q, and s not found',1)
!
!  Then find the shortest R vector. This is taken as c1. 
!
dist=1.d20
DO na=1, found
   tau(:) = ps(na) * at(:,1) + qs(na) * at(:,2) + ss(na) * at(:,3)
   dist1= SQRT( tau(1)** 2 + tau(2) ** 2 + tau(3)** 2 ) 
   IF (dist1 < dist .AND. dist1>eps) THEN
      p0=ps(na)
      q0=qs(na)
      s0=ss(na)
      dist=dist1
      na_first=na
   ENDIF
END DO
!WRITE(stdout,*) 'p0, q0, and s0', p0, q0, s0
c1(:) = p0 * at(:,1) + q0 * at(:,2) + s0 * at(:,3)
!
!  c2 is the vector closest to the origin not parallel to c1
!
dist=1.d20
DO na=1, found
   IF ( na /= na_first ) THEN
      tau(:) = ps(na) * at(:,1) + qs(na) * at(:,2) + ss(na) * at(:,3)
      CALL vector_prod(c1, tau, vprod)
      prod = vprod(1) ** 2 + vprod(2) ** 2 + vprod(3) ** 2
!
!   prod > 0 removes the vectors parallel to c1
!
      IF (prod > eps) THEN   
         dist1 = SQRT( tau(1) ** 2 + tau(2) ** 2 + tau(3) ** 2 ) 
         IF (dist1 < dist) THEN
            p01=ps(na)
            q01=qs(na)
            s01=ss(na)
            dist=dist1
         END IF
      END IF
   END IF
END DO
!WRITE(stdout,*) 'p01, q01, and s01', p01, q01, s01
c2(:) = p01 * at(:,1) + q01 * at(:,2) + s01 * at(:,3)
!
! Now we build the layers with the same technique. We take only one R vector
! in each layer, creating a mesh of R vectors in the plane of equation
! m p + n q + o r = j
! and taking the R vector closest to the line r=nu G, where nu is a parameter.
! Then in each layer we add to this R vector the coordinates of the nat_3d 
! atoms of the bulk unit cell. 
!
min_j = 100000000
nat=0
DO j=-(nlayers-1)/2+origin_shift, nlayers/2 + origin_shift
!
!  Generate a mesh of R vectors close to the point
!  where the perpendicular to the plane crosses the layer
!
   found=0
   IF ( o /= 0 ) THEN
      p1 = NINT ( j * prod1 / gmod / gmod )
      q1 = NINT ( j * prod2 / gmod / gmod )
      o1=ABS(o) * 5
      DO itry=0, o1
         DO jtry = 0, o1
            DO sign1=-1,1,2
               DO sign2=-1,1,2
                  IF (MOD(-(p1+sign1*itry)*m-(q1+sign2*jtry)*n+j,ABS(o))==0) &
                                                                        THEN
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
      n1=ABS(n) * 5
      DO itry=0, n1
         DO jtry = 0, n1
            DO sign1=-1,1,2
               DO sign2=-1,1,2
                  IF (MOD(-(p1+sign1*itry)*m-(s1+sign2*jtry)*o+j,ABS(n))==0) &
                                                                        THEN
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
      m1=ABS(m) * 5
      DO itry=0, m1
         DO jtry = 0, m1
            DO sign1=-1,1,2
               DO sign2=-1,1,2
                  IF (MOD(-(q1+sign1*itry)*n-(s1+sign2*jtry)*o+j,ABS(m))==0) &
                                                                        THEN
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

   IF ( found == 0  ) CALL errore('gener_3d_slab', 'p,q, and s not found',1)
!
!  compute the distance for all candidate R from the line r=nu G and take 
!  the shortest
!
   dist=1.d20
   DO na=1, found
      tau(:) = ps(na) * at(:,1) + qs(na) * at(:,2) + ss(na) * at(:,3)
      CALL vector_prod(g,tau,vprod)
      dist1= SQRT( vprod(1) ** 2 + vprod(2) ** 2 + vprod(3) ** 2 ) / gmod
      IF (dist1 < dist) THEN
         p=ps(na)
         q=qs(na)
         s=ss(na)
         dist=dist1
      ENDIF
   ENDDO
!
!  Zero distance means that this plane is identical to the plane
!  passing through the origin. If there is no closest plane with these
!  characteristics report that a bulk can be created by this number of
!  layers. Used to check consistency with the value predicted before.
!
   IF (dist < eps) THEN
      IF (ABS(j) < min_j .AND. j /= 0) min_j=ABS(j)
   ENDIF

!   WRITE(stdout,*) 'p, q, and s', p, q, s
!
!   Now to each R add the atoms of a 3D bulk unit cell
!
   DO na=1,nat_3d
      nat = nat + 1
      y(:,nat) = p * at(:,1) + q * at(:,2) + s * at(:,3) + tau_3d(:,na)
      atm(nat) = atm_3d(na)
      ityp_all(nat) = ityp(na)
   END DO
END DO

IF (min_j < 100000000) &
   WRITE(stdout,'(/,5x,"In this direction the bulk has lattice planes of",i8,&
                                                  " types")') min_j
!
!  Adjust vacuum if requested by the user.
!
IF (.NOT.(ldist_vacuum)) THEN
   nspace= INT(vacuum * gmod / alat) 
!
!  This instruction is needed to avoid to take one layer more than necessary
!  when vacuum is exactly an integer number of layers
!
   IF ((vacuum * gmod / alat - nspace)> small) nspace=nspace+1
   vacuum= nspace * alat / gmod
END IF
WRITE(stdout,'(/,5x,"Vacuum distance", f18.7, " a.u.")') vacuum
WRITE(stdout,'(5x,"Vacuum planes", f20.7)') vacuum * gmod / alat
!
!  The third primitive vector of the slab unit cell t is parallel to the 
!  versor g(:) / gmod and has length equal to the number of layers, 
!  multiplied by the layer distance alat/gmod plus the vacuum space
!
t(:) = ( nlayers * alat / gmod + vacuum ) * g(:) / gmod 
!
!  Put t in units of alat
!
t(:) = t(:) / alat
! 
!  Find the modulus of c1, c2 and t
!
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
!  modulus, take c1 and c2-c1 since these two vectors form an angle of 120 
!  degrees as the primitive vectors of the hexagonal lattice given by latgen.
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
!  Write information on the three vectors that define the slab
!
WRITE(stdout,'(/,5x,"The three vectors that define the slab are (alat units):")')

WRITE(stdout,'(/,5x,3f17.8)') c1(:)
WRITE(stdout,'(5x,3f17.8)') c2(:)
WRITE(stdout,'(5x,3f17.8)') t(:)
WRITE(stdout,'(/,5x,"The muduli of these vectors are:",3f14.7)') c1mod, &
                                                        c2mod, tmod
WRITE(stdout,'(5x,"In a.u. the moduli are:         ",3f14.7)') c1mod*alat, &
                                                        c2mod*alat, tmod*alat
WRITE(stdout,'(5x,"In Ang. the moduli are:         ",3f14.7)') &
                c1mod*alat*bohr_radius_si*1.D10, &
                c2mod*alat*bohr_radius_si*1.D10, &
                tmod*alat*bohr_radius_si*1.D10
!   
!   Now build a super-cell of slabs if the user requested it. First find the
!   new direct and reciprocal primitive lattice vectors.
!
d1(:) = t11 * c1(:) + t12 * c2(:)
d2(:) = t21 * c1(:) + t22 * c2(:)

idet = ABS(t11 * t22 - t21 * t12)

IF (idet == 0) CALL errore('gener_3d_slab','incorrect surface unit cell',1)
IF (idet /= 1) &
   WRITE(stdout,'(/,5x,"Surface cell has",i5," unit cells" )') idet 

CALL recips( d1, d2, t, bd1, bd2, bt )
!
!  Now take all the R vectors within the super-cell of slabs and for each one
!  of them add nat_3d atomic coordinates to the list of atoms. The cell with
!  m1=0, n1=0 is already in the atom list. 
!
iat=nat
copies = 1
DO m1 = MIN(t11, t21, t11+t21), MAX(t11, t21, t11+t21)
   DO n1 = MIN(t12, t22, t12+t22), MAX(t12, t22, t12+t22)
      IF ( m1==0 .AND. n1==0 ) CYCLE
      tau(:) = m1 * c1(:) + n1 * c2(:) 
      prod = tau(1) * bd1(1) + tau(2) * bd1(2) + tau(3) * bd1(3)
      prod1 = tau(1) * bd2(1) + tau(2) * bd2(2) + tau(3) * bd2(3)
      IF (prod >=-eps .AND. prod < 1.0_DP-eps .AND. &
                            prod1 >=-eps .AND. prod1<1.0_DP-eps) THEN
         DO na=1,nat
            iat=iat+1
            y(:,iat) = y(:,na) + tau(:)
            atm(iat) = atm(na)
            ityp_all(iat)= ityp_all(na)
         END DO
         copies=copies + 1
      END IF
   END DO
END DO

IF (copies /= idet) CALL errore('gener_3d_slab','error with the supercell',1)
!
!  Update the total number of atoms in the supercell.
!
nat=iat
!
!  Compute the area of the slab and of the super-cell of slabs, then
!  copy the lattice vectors of the super-cell of slabs into the lattice
!  vectors of the slab.
!
CALL vector_prod(c1,c2,vprod)
area = SQRT( vprod(1)**2 + vprod(2)**2 + vprod(3)**2 )
CALL vector_prod(d1,d2,vprod)
super_area = SQRT( vprod(1)**2 + vprod(2)**2 + vprod(3)**2 )
c1 = d1 
c2 = d2
!
!  Print information on the slab and super-cell of slabs
!
IF ( NINT( super_area / area ) - copies /= 0 ) &
   CALL errore('gener_3d_slab','Same problem with the copies',1)

WRITE(stdout,'(/,5x,"The area of the slab is:           ",f17.8," (a.u.)^2")') &
                                                      area*alat**2
WRITE(stdout,'(5x,"The area of the super-cell slab is:",f17.8," (a.u.)^2")') &
                                                      super_area*alat**2

volume=area * alat **3 / gmod
WRITE(stdout,'(/,5x,"Volume of one slab cell",f17.8," (a.u.)^3")') volume
WRITE(stdout,'(5x,"Volume of one bulk cell",f17.8," (a.u.)^3")') omega
!
!   Recalculate the modulus of c1 and c2 and the cosine of their angle, 
!   as well as the reciprocal lattice vectors of the slab. All vectors are
!   in units of alat. The reciprocal lattice vectors are in units of 
!   2 pi / alat.
!
c1mod = SQRT ( c1(1)**2 + c1(2)**2 + c1(3)**2 )
c2mod = SQRT ( c2(1)**2 + c2(2)**2 + c2(3)**2 )
prod1 = ( c1(1) * c2(1) + c1(2) * c2(2) + c1(3) * c2(3) ) / c1mod / c2mod

CALL recips(c1,c2,t,bc1,bc2,bt)
!
!  Open the output file and write the Bravais lattice index of
!  the slab, the lattice dimensions, and the atomic coordinates.
!
IF (ionode) THEN
   iuout=find_free_unit()
   OPEN(unit=iuout, file=TRIM(filename), status='unknown', form='formatted')
ENDIF
!
!  Identify the Bravais lattice of the slab
!
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
!  c1 and c2 have different moduli, use a simple orthorhombic cell.   
!
      c_alat = c1mod * alat
      celldm_sur(1) = c_alat
      celldm_sur(2) = c2mod / c1mod
      celldm_sur(3) = tmod / c1mod
      ibrav_sur=8
      IF (ionode) THEN
         WRITE (iuout, '("ibrav=8")')
         WRITE (iuout, '("celldm(1)= ",f15.8)') celldm_sur(1)
         WRITE (iuout, '("celldm(2)= ",f15.8)') celldm_sur(2)
         WRITE (iuout, '("celldm(3)= ",f15.8)') celldm_sur(3)
      END IF
   END IF
ELSEIF (ABS(c1mod-c2mod) < eps .AND. ABS(prod1 + 0.5_DP) < eps ) THEN
!
!  The angle between c1 and c2 is 120 degrees and the vectors are 
!  equal: use the hexagonal lattice.
!
   c_alat = c1mod * alat
   celldm_sur(1) = c_alat
   celldm_sur(3) = tmod / c1mod
   ibrav_sur=4
   IF (ionode) THEN
      WRITE (iuout, '("ibrav=4")')
      WRITE (iuout, '("celldm(1)= ",f15.8)') celldm_sur(1)
      WRITE (iuout, '("celldm(3)= ",f15.8)') celldm_sur(3)
   ENDIF
ELSEIF (ABS(c1mod-c2mod) < eps ) THEN
!
!  The angle between c1 and c2 is neither 90 nor 120 degrees, but the 
!  vectors have equal length, use the base-centered (C) orthorhombic lattice.
!
   alpha=ACOS(prod1) * 0.5_DP
   c_alat= 2.0_DP * c1mod * SIN(alpha) * alat
   celldm_sur(1) = c_alat
   celldm_sur(2) = 2.0_DP*c1mod*COS(alpha)*alat / c_alat
   celldm_sur(3) = tmod * alat / c_alat
   ibrav_sur=9
   IF (ionode) THEN
      WRITE (iuout, '("ibrav=9")')
      WRITE (iuout, '("celldm(1)= ",f15.8)') celldm_sur(1)
      WRITE (iuout, '("celldm(2)= ",f15.8)') celldm_sur(2)
      WRITE (iuout, '("celldm(3)= ",f15.8)') celldm_sur(3)
   ENDIF
ELSE
!
!  c1 and c2 form an angle different from 90 and 120 degrees or have different
!  lengths, use a monoclinic cell (c unique).
!
   c_alat = c1mod * alat
   celldm_sur(1) = c_alat
   celldm_sur(2) = c2mod / c1mod
   celldm_sur(3) = tmod / c1mod
   celldm_sur(4) = prod1
   ibrav_sur=12
   IF (ionode) THEN
      WRITE (iuout, '("ibrav=12")')
      WRITE (iuout, '("celldm(1)= ",f15.8)') celldm_sur(1)
      WRITE (iuout, '("celldm(2)= ",f15.8)') celldm_sur(2)
      WRITE (iuout, '("celldm(3)= ",f15.8)') celldm_sur(3)
      WRITE (iuout, '("celldm(4)= ",f15.8)') celldm_sur(4)
   ENDIF
END IF

IF (ionode) THEN
   WRITE (iuout, '("nat= ",i5)') nat
   WRITE (iuout, '("ATOMIC_POSITIONS {crystal}")') 
ENDIF
!
! For slabs with an even number of layers the origin can be moved
! in the center of two layers. Moreover in any case shift the z coordinates 
! of a distance equivalent to origin_shift layers.
!
dz(:) = - origin_shift * g(:) / gmod ** 2
IF ( MOD(nlayers,2) == 0 .AND. do_shift ) dz(:)=dz(:)-0.5_DP*g(:)/gmod**2
!
!  Put dz in crystal coordinates. Note that it has only the component parallel
!  to t.
!
shiftz = dz(1) * bt(1) + dz(2) * bt(2) + dz(3) * bt(3)
!
! Now compute the atomic positions in crystal coordinates of the slab lattice 
! vectors and add shiftz to the t component if necessary.
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
!  write the coordinates in a .xsf file readable by xcrysden with the 
!  original surface orientation
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
!  rotate the coordinates so that the normal to the surface is along z
!  and rewrite the .xsf file.
!
IF (ionode) &
OPEN(unit=iuout, file=TRIM(filename)//'.xsf', status='unknown', &
                                              form='formatted')

CALL latgen(ibrav_sur, celldm_sur, at(1,1), at(1,2), at(1,3), omega)
at=at/c_alat
CALL cryst_to_cart( nat, tau_sur, at, 1 )

IF (ionode) THEN
   CALL xsf_struct (c_alat, at, nat, tau_sur, atm_typ, ityp_all, iuout)
   CLOSE(iuout)
ENDIF

DEALLOCATE(tau_3d)
DEALLOCATE(ityp)
DEALLOCATE(atm_3d)
DEALLOCATE(tau_sur)

CALL environment_end( code )
CALL mp_global_end ()

END PROGRAM gener_3d_slab

!--------------------------------------------------------------------
SUBROUTINE remove_common_factors(m,n,o,fact)
!--------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(INOUT) :: m, n, o
INTEGER, INTENT(OUT) :: fact

INTEGER :: ind

fact=1
DO ind=2,100
   IF (MOD(m,ind)==0.AND.MOD(n,ind)==0.AND.MOD(o,ind)==0) THEN
      m = m/ind
      n = n/ind
      o = o/ind
      fact=fact*ind
   END IF
ENDDO

RETURN
END SUBROUTINE remove_common_factors

!--------------------------------------------------------------------
SUBROUTINE vector_prod(a,b,c)
!--------------------------------------------------------------------

USE kinds, ONLY : DP

IMPLICIT NONE
REAL(DP), INTENT(IN) :: a(3), b(3)
REAL(DP), INTENT(OUT) :: c(3)

c(1) = a(2) * b(3) - a(3) * b(2)
c(2) = a(3) * b(1) - a(1) * b(3)
c(3) = a(1) * b(2) - a(2) * b(1)

RETURN
END SUBROUTINE vector_prod

!--------------------------------------------------------------------
SUBROUTINE search_vicinals(g, bg)
!--------------------------------------------------------------------
!
!  This routine receives the G vector perpendicular to a surface and
!  the primitive vectors of the reciprocal lattice (of the conventional
!  cell for centered lattices) and computes the angle of several low
!  Miller index surfaces with the surface perpendicular to G. Writes on
!  output the angles of the three surfaces with the smallest miscut angles.
!
USE kinds,     ONLY : DP
USE constants, ONLY : pi
USE io_global, ONLY : stdout

IMPLICIT NONE
REAL(DP) :: g(3), bg(3,3)

INTEGER, PARAMETER  :: nvicinals=7
REAL(DP), PARAMETER :: small=1.D-10

INTEGER  :: gvic(3, nvicinals)     ! list of low Miller index surfaces
REAL(DP) :: angle(nvicinals),    & ! angle of our surface with the low index 
            gmod,                & ! modulus of the input g vector
            gv(3),               & ! G vector of the low index surface
            gmodv,               & ! modulus of the G vector of the low index
            alpha                  ! cosine of the angle
INTEGER  :: ivic,                & ! counter on vicinals
            ind(nvicinals)         ! index used to order the angles

DATA gvic / 1, 1, 1,   &           ! define the miller indeces of the low
            1, 1, 0,   &           ! index surfaces
            1, 0, 1,   &
            0, 1, 1,   &
            0, 0, 1,   &
            0, 1, 0,   &
            1, 0, 0    /

gmod=SQRT(g(1)**2+g(2)**2+g(3)**2)
DO ivic=1, nvicinals
   gv(:) = gvic(1,ivic) *bg(:,1) + gvic(2,ivic) * bg(:,2) + gvic(3,ivic)*bg(:,3)
   gmodv=SQRT(gv(1)**2 + gv(2)**2 + gv(3)**2)
   gv=gv/gmodv
   alpha=(gv(1)*g(1) + gv(2)*g(2) + gv(3)*g(3))/gmod 
!
!  Due to round-off errors the acos function could give a NaN. Set exactly to
!  1 or !  -1 the alpha that are close to 1 or -1.
!
   IF (ABS(alpha-1.0_DP)<small) alpha=1.0_DP
   IF (ABS(alpha+1.0_DP)<small) alpha=-1.0_DP
!
!  compute the angle and transform in degrees
!
   angle(ivic)=ACOS(alpha) * 180._DP / pi
   ind(ivic)=ivic
ENDDO
!
!  order the angles
!
CALL hpsort(nvicinals, angle, ind)
!
!  Write on output only the angles on the first three vicinals
!
WRITE(stdout,*)
DO ivic=1,3
   WRITE(stdout,'(5x,"The miscut angle with the (",3i1") surface is", f13.7, &
                 &" deg")') gvic(:,ind(ivic)), angle(ivic)
ENDDO

RETURN
END SUBROUTINE search_vicinals
