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
!  can put that distance
!  (ldist_vacuum=.TRUE.) or a vacuum space equivalent to an integer 
!  number of empty rows.
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
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER, PARAMETER :: nmax=10000

REAL(DP) :: a1(3), a2(3), b1(3), b2(3), g(3), t(3), c(3), d(3), t1(3), c1(3)
REAL(DP) :: alat, cmod, tmod, pi, prod1, prod2, vacuum, gmod, dist, dist1, &
            tau(3), b1eff(3), b2eff(3)
INTEGER :: m, n, m1, n1, p, q, nrows, nat_row, nat, ibrav_2d, nat_2d, q0, p0, &
           q1, ia, iat, j, itry, na, found, nspace
REAL(DP) :: y(3,nmax), alat_box, celldm_2d(3)
INTEGER :: ps(nmax), qs(nmax)
REAL(DP), ALLOCATABLE :: tau_2d(:,:)
CHARACTER ( LEN=3 ), ALLOCATABLE :: atm_2d(:)
CHARACTER ( LEN=3 ) :: atm(nmax)
LOGICAL :: ldist_vacuum
INTEGER :: iuout
CHARACTER(LEN=256) :: filename

WRITE(6,'("ibrav_2d ")')
WRITE(6,'("1 -- oblique, give a, b, cos(gamma) ")')
WRITE(6,'("2 -- rectangular, give a and b ")')
WRITE(6,'("3 -- centered rectangular, give a and b ")')
WRITE(6,'("4 -- square, give a ")')
WRITE(6,'("5 -- hexagonal, give a ")')
WRITE(6,'("ibrav_2d?")')
READ(5,*) ibrav_2d
WRITE(6,'(i5)') ibrav_2d
WRITE(6,'("a, b/a, COS(gamma)? ")')
READ(5,*) celldm_2d(1), celldm_2d(2), celldm_2d(3)
WRITE(6,'(3f15.6)') celldm_2d
alat=celldm_2d(1)
WRITE(6,'("Number of atoms in the 2d unit cell?")') 
READ(5,*) nat_2d
WRITE(6,'(i5)') nat_2d

ALLOCATE(tau_2d(3,nat_2d))
ALLOCATE(atm_2d(nat_2d))
DO ia=1,nat_2d
   READ(5,*) atm_2d(ia), tau_2d(1,ia), tau_2d(2,ia), tau_2d(3,ia)
   WRITE(6,'(a3, 3f18.10)') atm_2d(ia), tau_2d(1,ia), tau_2d(2,ia), tau_2d(3,ia)
ENDDO

WRITE(6,'("Dimension of the box?")')
READ(5,*) alat_box
WRITE(6,'(f15.6)') alat_box

WRITE(6,'("Crystal coordinates of the G vector m b1 + n b2 (m,n)?")')
READ(5,*) m, n
WRITE(6,'(2i5)') m,n

WRITE(6,'("Number of rows ?")')
READ(5,*) nrows
IF ( MOD(nrows,2) == 0 ) THEN
   nrows=nrows+1
END IF
WRITE(6,'(i5)') nrows

WRITE(6,'("Number of atoms per row ?")')
READ(5,*) nat_row
WRITE(6,'(i5)') nat_row

WRITE(6,'("Exact vacuum (.TRUE.) or row distance multiples (.FALSE.)?")')
READ(5,*) ldist_vacuum
WRITE(6,*) ldist_vacuum

WRITE(6,'("Vacuum space in a.u. ?")')
READ(5,*) vacuum
WRITE(6,'(f15.6)') vacuum 

WRITE(6,'("Output file name?")')
READ(5,*) filename
WRITE(6,'(a)') TRIM(filename)

pi=4.0_DP * atan(1.0_DP)
CALL latgen_2d(ibrav_2d, celldm_2d, a1, a2)

WRITE(6,'("Direct lattice vectors")')
WRITE(6,'("(",f15.6,",",f15.6,")")') a1(1), a1(2)
WRITE(6,'("(",f15.6,",",f15.6,")")') a2(1), a2(2)

CALL recips_2d(a1,a2,b1,b2)

WRITE(6,'("Reciprocal lattice vectors")')
WRITE(6,'("(",f15.6,",",f15.6,")")') b1(1), b1(2)
WRITE(6,'("(",f15.6,",",f15.6,")")') b2(1), b2(2)

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

IF (m == 0) THEN
   nat=0
   DO j=-(nrows-1)/2, (nrows-1)/2
      q = j
      p = - NINT ( j * prod2 / prod1 )   
      DO ia=1,nat_2d
         nat = nat + 1
         y(:,nat) = p * a1(:) + q * a2(:) + tau_2d(:,ia)
         atm(nat) = atm_2d(ia)
      END DO
   END DO
   p0=nat_row
   q0=0
ELSE
   DO itry=1,m
      IF ( MOD(itry * n, m) == 0 ) THEN
         q0 = itry
         p0 = -q0 * n / m
      END IF
   END DO
   nat=0
   DO j=-(nrows-1)/2, (nrows-1)/2
      q1 = NINT (  j * (g(1) * b2(1) + g(2) * b2(2)) / gmod )
      found=0
      DO itry=0,m
         IF ( MOD( -( q1 + itry ) * n + j, m ) == 0 ) THEN
            found=found+1
            ps(found)= ( - ( q1 + itry ) * n + j ) / m
            qs(found)= q1 + itry
         END IF
         IF ( MOD( -( q1 - itry ) * n + j, m ) == 0 ) THEN
            found=found+1
            ps(found)= ( -( q1 - itry ) * n + j ) / m
            qs(found)= q1 - itry
         END IF
      END DO
      IF ( found == 0  ) THEN
         WRITE(6,*) 'p and q not found'
         STOP
      END IF
!
!   compute the distance for all candidates p and q and take the shortest
!
      dist=1.d20
      DO ia=1, found
         tau(:) = ps(ia) * a1(:) + qs(ia) * a2(:) 
         dist1= ABS(g(2) * tau(1) - g(1) * tau(2)) / gmod
         IF (dist1 < dist) THEN
            p=ps(ia)
            q=qs(ia)
            dist=dist1
         ENDIF
      ENDDO
      write(6,*) 'p and q', p, q
      DO ia=1,nat_2d
         nat = nat + 1
         y(:,nat) = p * a1(:) + q * a2(:) + tau_2d(:,ia)
         atm(nat) = atm_2d(ia)
      END DO
   END DO
END IF

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

iuout=28
OPEN(unit=iuout, file=TRIM(filename), status='unknown', form='formatted')
!
!  we use an orthorombic cell
!
WRITE (iuout, '("ibrav=8")')
WRITE (iuout, '("celldm(1)= ",f15.8)') cmod * alat
WRITE (iuout, '("celldm(2)= ",f15.8)') tmod/cmod
WRITE (iuout, '("celldm(3)= ",f15.8)') alat_box/cmod/alat
WRITE (iuout, '("nat= ",i5)') nat
WRITE (iuout, '("ATOMIC_POSITIONS {crystal}")') 
DO na=1,nat
   prod1 = y(1,na) * c1(1) + y(2,na) * c1(2)
   prod2 = y(1,na) * t1(1) + y(2,na) * t1(2)
   WRITE (iuout,'(a,3f18.10)') atm(na), prod1, prod2, y(3,na) * alat / alat_box
ENDDO

CLOSE(iuout)

END PROGRAM gener_2d_slab

