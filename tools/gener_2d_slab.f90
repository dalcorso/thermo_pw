!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM gener_2d_slab
!
!  This program generates all the points of a two dimensional solid
!  defined by a two dimensional Bravais lattice and nat_2d atoms in
!  the unit cell. Calling a1 and a2 the two primitive vectors of the
!  2D Bravais lattice, it prints on output the informations necessary, 
!  using periodic boundary conditions to simulate a ribbon whose edge 
!  is perpendicular to a 2d G vector defined as G = m b1 + n b2. It is 
!  supposed that m and n have no common factors. 
!  The user specifies the number of atomic rows nrows to be used and
!  the number of atoms per row nat_row and the program writes on output
!  the coordinates of nrows x nat_row atoms and a unit cell sufficient
!  to simulate the ribbon. The distance between repeated copies of the
!  ribbon can be indicated as input (vacuum) in a.u. and the program 
!  can put that distance (ldist_vacuum=.TRUE.) or a vacuum space equivalent 
!  to an integer number of empty rows.
!  The atomic coordinates are written on output and periodic boundary
!  conditions are assumed in the xy plane, while a unit cell that can be
!  used by pw.x is calculated in such a way that the distance between
!  ribbons along z is alat_box.
!  List of input variables:
!
!  ibrav_2d     The Bravais lattice type: 1 oblique, 2 rectangular, 
!               3 centered rectangular, 4 square, 5 hexagonal.
!  
!  celldm_2d(3) The size of the unit cell.
!
!  nat_2d       The number of atoms inside the unit cell
!
!  atm_2d, tau_2d(3) ! label and coordinates of the atoms inside the unit cell.
!               these atomic coordinated might have also a z component 
!               producing a buckled layer and a buckled wire, or adding a
!               molecule to a nanowire.
!  alat_box     the size of the box for pw.x calculations
!
!  m,n          m and n  of the G vector
!
!  nrows        number of rows (at the moment must be odd)
!               
!  nat_row      number of atoms per row
!
!  ldist_vacuum if .TRUE. the distance between repeated ribbons is vacuum,
!               otherwise it is the smaller integer number of rows larger
!               than vacuum
!
!  vacuum       the vacuum distance
!
!  filename     output file name
!
USE kinds,       ONLY : DP
USE constants,   ONLY : bohr_radius_si
USE atomic_pos,  ONLY : find_ityp
USE io_global,   ONLY : stdout, ionode
USE mp_global,   ONLY : mp_startup, mp_global_end
USE environment, ONLY : environment_start, environment_end

IMPLICIT NONE

REAL(DP) :: a1(3), a2(3), b1(3), b2(3), g(3), t(3), c(3), d(3), t1(3), c1(3)
REAL(DP) :: alat, cmod, tmod, pi, prod1, prod2, vacuum, gmod, dist, dist1, &
            tau(3), b1eff(3), b2eff(3)
INTEGER  :: m, n, m1, n1, p, q, nrows, nat_row, nat, ibrav_2d, nat_2d, q0, p0, &
            p1, q1, ia, iat, j, itry, na, nt, ntyp, found, nspace, iuout
INTEGER, PARAMETER :: ntypx=10
REAL(DP) :: alat_box, celldm_2d(3), celldm(6), omega, at(3,3)
INTEGER  :: direction, jbulk
INTEGER, ALLOCATABLE :: ityp(:), ps(:), qs(:)
INTEGER :: find_free_unit
INTEGER :: stdin
REAL(DP), ALLOCATABLE :: tau_2d(:,:), tau_ribbon(:,:), y(:,:)
CHARACTER ( LEN=3 ), ALLOCATABLE :: atm_2d(:)
CHARACTER ( LEN=3 ), ALLOCATABLE :: atm(:)
LOGICAL :: ldist_vacuum, invert
CHARACTER(LEN=256) :: filename, xsf_filename
CHARACTER(LEN=3) :: atm_typ(ntypx)
CHARACTER(LEN=9) :: code='2D_SLAB'

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

stdin=5
WRITE(stdout,'("ibrav_2d ")')
WRITE(stdout,'("1 -- oblique, give a, b, cos(gamma) ")')
WRITE(stdout,'("2 -- rectangular, give a and b ")')
WRITE(stdout,'("3 -- centered rectangular, give a and b ")')
WRITE(stdout,'("4 -- square, give a ")')
WRITE(stdout,'("5 -- hexagonal, give a ")')
WRITE(stdout,'("ibrav_2d?")')

READ(stdin,*) ibrav_2d
WRITE(stdout,'(i5)') ibrav_2d
WRITE(stdout,'("a, b/a, COS(gamma)? ")')

READ(stdin,*) celldm_2d(1), celldm_2d(2), celldm_2d(3)
WRITE(stdout,'(3f15.6)') celldm_2d
alat=celldm_2d(1)
WRITE(stdout,'("Number of atoms in the 2d unit cell?")') 

READ(stdin,*) nat_2d
WRITE(stdout,'(i5)') nat_2d
ALLOCATE(tau_2d(3,nat_2d))
ALLOCATE(atm_2d(nat_2d))
DO ia=1,nat_2d
   READ(stdin,*) atm_2d(ia), tau_2d(1,ia), tau_2d(2,ia), tau_2d(3,ia)
   WRITE(stdout,'(a3, 3f18.10)') atm_2d(ia), tau_2d(1,ia), tau_2d(2,ia), &
                                             tau_2d(3,ia)
ENDDO

WRITE(stdout,'("Dimension of the box?")')
READ(stdin,*) alat_box
WRITE(stdout,'(f15.6)') alat_box

WRITE(stdout,'("Crystal coordinates of the G vector m b1 + n b2 (m,n)?")')
READ(stdin,*) m, n
WRITE(stdout,'(2i5)') m,n

WRITE(stdout,'("Number of rows ?")')
READ(stdin,*) nrows
WRITE(stdout,'(i5)') nrows

WRITE(stdout,'("Number of atoms per row ?")')
READ(stdin,*) nat_row
WRITE(stdout,'(i5)') nat_row

WRITE(stdout,'("Exact vacuum (.TRUE.) or row distance multiples (.FALSE.)?")')
READ(stdin,*) ldist_vacuum
WRITE(stdout,*) ldist_vacuum

WRITE(stdout,'("Vacuum space in a.u. ?")')
READ(stdin,*) vacuum
WRITE(stdout,'(f15.6)') vacuum 

WRITE(stdout,'("Output file name?")')
READ(stdin,*) filename
WRITE(stdout,'(a)') TRIM(filename)

pi=4.0_DP * atan(1.0_DP)
CALL latgen_2d(ibrav_2d, celldm_2d, a1, a2)

WRITE(stdout,'("Direct lattice vectors")')
WRITE(stdout,'("(",f15.6,",",f15.6,")")') a1(1), a1(2)
WRITE(stdout,'("(",f15.6,",",f15.6,")")') a2(1), a2(2)

CALL recips_2d(a1,a2,b1,b2)

WRITE(stdout,'("Reciprocal lattice vectors")')
WRITE(stdout,'("(",f15.6,",",f15.6,")")') b1(1), b1(2)
WRITE(stdout,'("(",f15.6,",",f15.6,")")') b2(1), b2(2)

nat=nrows * nat_row * nat_2d  
ALLOCATE(y(3,nat))
ALLOCATE(atm(nat))

IF (ibrav_2d == 3)  THEN
!
!  With a centered rectangular lattice m and n refer to the rectangular
!  reciprocal lattice vectors
!
   b1eff(1) = 2.0_DP
   b1eff(2) = 0.0_DP
   b1eff(3) = 0.0_DP
   b2eff(1) = 0.0_DP 
   b2eff(2) = 2.0_DP / celldm_2d(2)
   b2eff(3) = 0.0_DP
   m1 = m
   n1 = n
   g(:) = m1 * b1eff(:) + n1 * b2eff(:)
!
!  here we return to the m and n of the centered tetragonal lattice
!
   m = m1 + n1
   n = n1 - m1 
ELSE
   g(:) = m * b1(:) + n * b2(:)
ENDIF
gmod = SQRT( g(1)**2 + g(2)**2 )

prod1= g(1) * a1(2) - g(2) * a1(1)
prod2= g(1) * a2(2) - g(2) * a2(1)

IF (.NOT.(ldist_vacuum)) THEN
!
!  in this case the vacuum is a multiple of the row distace sufficient to
!  contain a vacuum space equal to vacuum
!
   nspace= INT(vacuum * gmod / alat) + 1
   vacuum= nspace * alat / gmod
END IF

ALLOCATE(ps(-(nrows-1)/2:nrows/2))
ALLOCATE(qs(-(nrows-1)/2:nrows/2))

jbulk=10000000
IF (m == 0) THEN
   p0=1
   q0=0
   DO j=-(nrows-1)/2, nrows/2
      q = j
      p = - NINT ( j * prod2 / prod1 )   
      tau(:) = p * a1(:) + q * a2(:) 
      dist= ABS(g(2) * tau(1) - g(1) * tau(2)) / gmod
!      WRITE(stdout,*) 'p and q', p, q
      WRITE(stdout,'(5x,"Layer j",i8," distance of the lattice",f15.7)') j, dist
      IF (ABS(dist)<1.D-8.AND.ABS(j)<ABS(jbulk).AND.j/=0) jbulk=j
      ps(j)=p
      qs(j)=q
   END DO
ELSEIF (n == 0) THEN
   p0=0
   q0=1
   DO j=-(nrows-1)/2, nrows/2
      p = j
      q = - NINT ( j * prod1 / prod2 )   
      tau(:) = p * a1(:) + q * a2(:) 
      dist= ABS(g(2) * tau(1) - g(1) * tau(2)) / gmod
!      WRITE(stdout,*) 'p and q', p, q
      WRITE(stdout,'(5x,"Layer j",i8," distance of the lattice",f15.7)') j, dist
      IF (ABS(dist)<1.D-8.AND.ABS(j)<ABS(jbulk).AND.j/=0) jbulk=j
      ps(j)=p
      qs(j)=q
   END DO
ELSE
   DO itry=1,m
      IF ( MOD(itry * n, m) == 0 ) THEN
         q0 = itry
         p0 = -q0 * n / m
      END IF
   END DO
   DO j=-(nrows-1)/2, nrows/2     ! valid for even and for odd nrows
      q1 =  j / n 
      found=0
      DO itry=0,m-1
         IF ( MOD( -( q1 - itry ) * n + j, m ) == 0 ) THEN
            p= ( - ( q1 + itry ) * n + j ) / m
            q= q1 + itry
            found=1
         END IF
         IF (found==1) EXIT
      END DO
      IF ( found == 0  ) THEN
         WRITE(stdout,*) 'p and q not found'
         STOP
      END IF
!
!   compute the distance for the candidate p and q 
!
      tau(:) = p * a1(:) + q * a2(:) 
      dist= ABS(g(2) * tau(1) - g(1) * tau(2)) / gmod
!      WRITE(stdout,'("p and q ",2i5)') p, q
!      WRITE(stdout,'("tau ",2f15.8,"  dist=",f15.9)') tau(1:2), dist
!
!   now try in one direction and then in the opposite. When the distance
!   grows we have found the minimum distance from the line passing through
!   the origin and parallel to G/|G|
!
      direction=1
      dist1=0.0_DP
      invert=.TRUE.
      DO WHILE (dist1 < dist) 
         p1=p+p0*direction
         q1=q+q0*direction
         tau(:) = p1 * a1(:) + q1 * a2(:) 
         dist1= ABS(g(2) * tau(1) - g(1) * tau(2)) / gmod
!         WRITE(stdout,'("p1 and q1 ",2i5)') p1, q1
!         WRITE(stdout,'("tau ",2f15.8,"  dist1=",f15.9, i3)') tau(1:2), dist1, &

!             direction

         IF (dist1 < dist) THEN
            p=p1
            q=q1
            dist=dist1
            dist1=0.0_DP
         ELSEIF (invert) THEN
            direction=-1
            dist1=0.0_DP
         ENDIF
         invert=.FALSE.
      ENDDO
!      WRITE(stdout,*) 'p and q', p, q
      WRITE(stdout,'(5x,"Layer j",i8," distance of the lattice",f15.7)') j, dist
      IF (ABS(dist)<1.D-8.AND.ABS(j)<ABS(jbulk).AND.j/=0) jbulk=j 
      ps(j)=p
      qs(j)=q
   END DO
END IF
!
!  create the ribbon with one lattice point per row
!
nat=0
DO j=-(nrows-1)/2, nrows/2     ! valid for even and for odd nrows
   DO ia=1,nat_2d
      nat = nat + 1
      y(:,nat) = ps(j) * a1(:) + qs(j) * a2(:) + tau_2d(:,ia)
      atm(nat) = atm_2d(ia)
   END DO
ENDDO

IF (jbulk<10000000) THEN
   WRITE(stdout,'(/,5x,"The 2D sheet can be simulated with a &
                                   &minimum of",i6," rows")') ABS(jbulk)
ELSE
   WRITE(stdout,'(/,5x,"Unknown size of the 2D sheet. Try to increase nrows &
                          &to find it")') 
ENDIF
c(:) = p0 * a1(:) + q0 * a2(:)
t(:) = ( (nrows - 1) * alat / gmod + vacuum ) * g(:) / gmod / alat

!
!  If c and t do not have the same orientation as x and y 
!  change the sign of c
!
prod1 = c(1) * t(2) - t(1) * c(2)
IF (prod1 < 0.0_DP) c(:) = - c(:)
!
!  write all the periodic copies along a row if the user asked for more
!  atoms per row
!
IF (nat_row > 1) THEN
   iat=nat
   DO j = 1, nat_row-1
      d(:) = j * c(:)
      DO ia=1,nat
         iat=iat+1
         y(:,iat) = y(:,ia) + d(:)
         atm(iat) = atm(ia)
      ENDDO
   ENDDO
   nat = nat * nat_row
   c = c * nat_row
END IF

CALL recips_2d(c,t,c1,t1)
cmod = SQRT ( c(1)**2 + c(2)**2 )
tmod = SQRT ( t(1)**2 + t(2)**2 )
WRITE(stdout,'(/,5x,"The step periodicity is:",f15.8," a.u.", f15.8," A")') &
                                     cmod*alat, cmod*alat*bohr_radius_si*1.D10
!
!  we use an orthorhombic cell
!
celldm=0.0_DP
celldm(1)=cmod * alat
celldm(2)=tmod/cmod
celldm(3)=alat_box/cmod/alat
ALLOCATE(tau_ribbon(3,nat))
ALLOCATE(ityp(nat))
IF (ionode) THEN
   iuout=find_free_unit()
   OPEN(unit=iuout, file=TRIM(filename), status='unknown', form='formatted')
   WRITE (iuout, '("ibrav=8")')
   WRITE (iuout, '("celldm(1)= ",f15.8)') celldm(1)
   WRITE (iuout, '("celldm(2)= ",f15.8)') celldm(2)
   WRITE (iuout, '("celldm(3)= ",f15.8)') celldm(3)
   WRITE (iuout, '("nat= ",i5)') nat
   WRITE (iuout, '("ATOMIC_POSITIONS {crystal}")') 
   DO na=1,nat
      prod1 = y(1,na) * c1(1) + y(2,na) * c1(2)
      prod2 = y(1,na) * t1(1) + y(2,na) * t1(2)
      tau_ribbon(1,na)=prod1
      tau_ribbon(2,na)=prod2
      tau_ribbon(3,na)=y(3,na) * alat / alat_box
      WRITE (iuout,'(a,3f18.10)') atm(na), prod1, prod2, &
                                           y(3,na) * alat / alat_box
   ENDDO
   CLOSE(iuout)
END IF
!
!  Count how many types of atoms we have and how they are called
!
CALL find_ityp(nat, atm, ntyp, ityp, atm_typ, ntypx)

IF (ionode) THEN
   iuout=find_free_unit()
   xsf_filename=TRIM(filename)//'.xsf'
   OPEN(UNIT=iuout, FILE=TRIM(xsf_filename), STATUS='unknown', &
                                                             FORM='formatted')
   CALL latgen(8, celldm, at(1,1), at(1,2), at(1,3), omega)
   at=at/celldm(1) 
!
!   bring the atomic positions in cartesian coordinates in unit of alat
!
   CALL cryst_to_cart(nat,tau_ribbon,at,1)
   CALL xsf_struct (celldm(1), at, nat, tau_ribbon, atm_typ, ityp, iuout)
   CLOSE(iuout)
ENDIF

DEALLOCATE(y)
DEALLOCATE(atm)
DEALLOCATE(ps)
DEALLOCATE(qs)
DEALLOCATE(tau_ribbon)
DEALLOCATE(ityp)

CALL environment_end( code )
CALL mp_global_end ()

END PROGRAM gener_2d_slab

