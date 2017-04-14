!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM gener_nanowire
!
!  This program generates all the points of a two dimensional solid
!  defined by a two dimensional Bravais lattice and nat_2d atoms in
!  the unit cell. Calling a1 and a2 the two primitive vectors of the
!  2D Bravais lattice, it prints on output all the atoms that are
!  contained in the surface defined by the four integer numbers
!  (m,n) and (p,q). The four numbers define two vectors
!  C = m a1 + n a2
!  T = p a1 + q a2 
!  and the area is defined by r = x C + y T, with 0<= x <1 and 0<= y < 1.
!  The two numbers m and n are given in input, while p and q can be given
!  in input if lrect=.FALSE. or if the lattice is not contained in a 
!  rectangle that has C as one edge. Otherwise if lrect=.TRUE. the code
!  finds, for the lattices that allow it the two integers p and q such
!  that T is orthogonal to C. 
!  A multiplicative factor nz can be given in input to take nz p, and nz q,
!  instead of p and q.
!  The size of the unit cell is given as input and defined by the vector 
!  celldm_2d(3) that contains |a|, |b| and cos(gamma)
!  The atomic coordinates are written on output and periodic boundary
!  conditions are assumed in the plane, while a unit cell that can be
!  used by pw.x is calculated in such a way that the distance between
!  planes is alat_box.
!  Finally the flag lwire=.true. allows to obtain as output the coordinates
!  of the atoms in the given area wrapped around a cilider obtained
!  obtained by transforming the vector C in a circle and building
!  a cilinder of height |T| cos (alpha) where alpha is the angle between
!  C and T. The wire is put inside a tetragonal box. The a size of this box
!  is alat_box, while the size c is the height of the wire.
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
!  lwire        .true. if the output coordinates are those of a wire
!
!  lrect        .true. if p and q are chosen by the program
!
!  m or m,n     m and n 
!
!  nz           unit of repetition along z.
!
USE kinds,       ONLY : DP
USE io_global,   ONLY : stdout
USE mp_global,   ONLY : mp_startup, mp_global_end
USE environment, ONLY : environment_start, environment_end
IMPLICIT NONE

INTEGER, PARAMETER :: nmax=10000

REAL(DP), PARAMETER :: eps=1.d-12
REAL(DP) :: a1(3), a2(3), c(3), t(3), c1(3), t1(3), tau(3)
REAL(DP) :: alat, fact, r, cmod, tmod, pi, phi, z, prod, prod1, xq
INTEGER :: m, n, p, q, ip, nz, n1, m1, ia, na, nat, ibrav_2d, nat_2d
REAL(DP) :: y(3,nmax), nano(3,nmax), alat_box, celldm_2d(3)
REAL(DP), ALLOCATABLE :: tau_2d(:,:)
CHARACTER(LEN=3), ALLOCATABLE :: atm_2d(:)
CHARACTER(LEN=3) :: atm(nmax)
LOGICAL :: lwire, lrect
INTEGER :: iuout
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=9) :: code='NANOWIRE'

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )


WRITE(stdout,'("ibrav_2d ")')
WRITE(stdout,'("1 -- oblique, give a, b, cos(gamma) ")')
WRITE(stdout,'("2 -- rectangular, give a and b ")')
WRITE(stdout,'("3 -- centered rectangular, give a and b ")')
WRITE(stdout,'("4 -- square, give a ")')
WRITE(stdout,'("5 -- hexagonal, give a ")')
WRITE(stdout,'("ibrav_2d?")')
READ(5,*) ibrav_2d
WRITE(stdout,'(i5)') ibrav_2d
WRITE(stdout,'("a, b/a, COS(gamma)? ")')
READ(5,*) celldm_2d(1), celldm_2d(2), celldm_2d(3)
WRITE(stdout,'(3f15.6)') celldm_2d
alat=celldm_2d(1)
WRITE(stdout,'("Number of atoms in the 2d unit cell?")') 
READ(5,*) nat_2d
WRITE(stdout,'(i5)') nat_2d

ALLOCATE(tau_2d(3,nat_2d))
ALLOCATE(atm_2d(nat_2d))
DO ia=1,nat_2d
   READ(5,*) atm_2d(ia), tau_2d(1,ia), tau_2d(2,ia), tau_2d(3,ia)
   WRITE(stdout,'(a3, 3f18.10)') atm_2d(ia), tau_2d(1,ia), tau_2d(2,ia), &
                                                            tau_2d(3,ia)
ENDDO

WRITE(stdout,'("Dimension of the box?")')
READ(5,*) alat_box
WRITE(stdout,'(f15.6)') alat_box

WRITE(stdout,'("Two dimensional sheet (.FALSE.) or wire (.TRUE.)?")')
READ(5,*) lwire
WRITE(stdout,*) lwire

SELECT CASE (ibrav_2d) 
   CASE(1)
!
!  rhombus
!
      WRITE(stdout,'("General lattice, give m,n,p,q")')
      READ(5,*) m, n, p, q
   CASE(2)
!
!  rectangular
!
      WRITE(stdout,'("Only nanowires of type (m,0) or (0,n) have rectangular cell")')
      WRITE(stdout,'("Rectangular? If .TRUE. then give m and n otherwise &
                 &give m n p q")')
      READ(5,*) lrect
      WRITE(stdout,*) lrect
      IF (lrect) THEN
         READ(5,*) m, n
         IF (m > 0) THEN
            p=0
            q=1
         ELSE
            p=1
            q=0
         ENDIF
      ELSE
         READ(5,*) m, n, p, q
      ENDIF
   CASE(3)
!
!  centered rectangular
!
      WRITE(stdout,'("Only nanowires of type (m,m) are rectangular.")')
      WRITE(stdout,'("Rectangular? If .TRUE. then give m otherwise &
                 &give m n p q")')
      READ(5,*) lrect
      WRITE(stdout,*) lrect
      IF (lrect) THEN
         READ(5,*) m
         n=m
         p=1
         q=-p
      ELSE
         READ(5,*) m, n, p, q
      ENDIF
   CASE(4)
!
!  square
!
      WRITE(stdout,'("Nanowires of all types (m,n) are rectangular.")')
      WRITE(stdout,'("Rectangular? If .TRUE. then give m and n otherwise &
                 &give m n p q")')
      READ(5,*) lrect
      WRITE(stdout,*) lrect
      IF (lrect) THEN
         READ(5,*) m,n
         IF (n>0) THEN
            fact=DBLE(m) / DBLE(n)
            DO ip=1, ABS(n)
               xq = - ip * fact
               IF (ABS(xq - NINT(xq)) < 1.d-10 ) THEN
                  p=ip
                  q=NINT(xq)
                  EXIT
               ENDIF
            ENDDO
         ELSE
           q=1
           p=0
         END IF
      ELSE
         READ(5,*) m, n, p, q
      ENDIF
   CASE(5)
!
!  hexagonal
!
      WRITE(stdout,'("Nanowires of all types (m,n) are rectangular.")')
      WRITE(stdout,'("Rectangular? If .TRUE. then give m and n otherwise &
                 &give m n p q")')
      READ(5,*) lrect
      WRITE(stdout,*) lrect
      IF (lrect) THEN
         READ(5,*) m,n
         IF (2*n-m /= 0) THEN
            fact= DBLE( 2 * m - n ) / DBLE( 2 * n - m )
            p=0
            q=0
            DO ip=1, ABS(2*n-m)
               xq = - ip * fact
               WRITE(6,*) ip, ABS(xq - NINT(xq)), NINT(xq)
               IF (ABS(xq - NINT(xq)) < 1.d-10 ) THEN
                  p=ip
                  q=NINT(xq)
                  EXIT
               ENDIF
            ENDDO
         ELSE
           q=1
           p=0
         END IF
      ELSE
         READ(5,*) m, n, p, q
      ENDIF
END SELECT

IF (p==0.AND.q==0) THEN
   WRITE(stdout,*) 'Unable to find p and q'
   STOP 1
ENDIF

WRITE(stdout,'("repeated units along z?")')
READ(5,*) nz

IF (nz > 1) THEN
   p=p*nz
   q=q*nz
ENDIF

WRITE(stdout,'("Output file name?")')
READ(5,*) filename

pi=4.0_DP * atan(1.0_DP)
CALL latgen_2d(ibrav_2d, celldm_2d, a1, a2)

WRITE(stdout,'("Direct lattice vectors")')
WRITE(stdout,'("(",f15.6,",",f15.6,")")') a1(1), a1(2)
WRITE(stdout,'("(",f15.6,",",f15.6,")")') a2(1), a2(2)

WRITE(stdout,'("(m, n), and (p, q)", 4i6)') m, n, p, q

c(:) = m * a1(:) + n * a2(:)
t(:) = p * a1(:) + q * a2(:)

CALL recips_2d(c,t,c1,t1)

cmod=SQRT(c(1)**2 + c(2)**2)
tmod=SQRT(t(1)**2 + t(2)**2)
nat=0
DO m1 = MIN(-1,m-1,p-1,m+p-1), MAX(1,m+1,p+1,m+p+1)
   DO n1 = MIN(-1,n-1,q-1,n+q-1), MAX(1,n+1,q+1,n+q+1)
      DO ia=1, nat_2d
         tau(:) = m1 * a1(:) + n1 * a2(:) + tau_2d(:,ia)
         prod = tau(1)*c1(1)+tau(2)*c1(2) 
         prod1 = tau(1)*t1(1)+tau(2)*t1(2) 
         IF (prod >=-eps.AND.prod < 1.0_DP-eps.AND.&
                               prod1 >=-eps.AND.prod1<1.0_DP-eps)THEN
            nat=nat+1
            IF ( nat > nmax) THEN
               WRITE(stdout,'("Too many atoms")')
               STOP 1
            END IF
            y(1,nat)=prod
            y(2,nat)=prod1
            y(3,nat)= tau(3)
            atm(nat) = atm_2d(ia)
         ENDIF
      END DO
   END DO
END DO

iuout=28
OPEN(unit=iuout, file=TRIM(filename), status='unknown', form='formatted')

IF (.NOT.lwire) THEN
   prod=c(1)*t(1)+c(2)*t(2)
   IF (ABS(prod)<1.d-8) THEN
!
!  vectors are orthogonal, a wire can be built
!
      WRITE (iuout, '("ibrav=8")')
      WRITE (iuout, '("celldm(1)= ",f15.8)') cmod * alat
      WRITE (iuout, '("celldm(2)= ",f15.8)') tmod/cmod
      WRITE (iuout, '("celldm(3)= ",f15.8)') alat_box/cmod/alat
      WRITE (iuout, '("nat= ",i5)') nat
   ELSE
      WRITE (iuout, '("ibrav=12")')
      WRITE (iuout, '("celldm(1)= ",f15.8)') cmod * alat
      WRITE (iuout, '("celldm(2)= ",f15.8)') tmod/cmod
      WRITE (iuout, '("celldm(3)= ",f15.8)') alat_box/cmod/alat
      WRITE (iuout, '("celldm(4)= ",f15.8)') prod/cmod/tmod
   ENDIF
   WRITE (iuout, '("nat= ",i5)') nat
   WRITE (iuout, '("ATOMIC_POSITIONS {crystal}")') 
   DO na=1,nat
      WRITE (iuout,'(a,3f18.10)') atm(na), y(1,na), y(2,na), y(3,na)
   ENDDO
ELSE
   r=cmod / 2.0_DP / pi 

   WRITE (iuout, '("wire radius r= ",f15.8, "  a.u.")') r * alat
   WRITE (iuout, '("ibrav=6")')
   WRITE (iuout, '("celldm(1)= ",f15.8)') alat_box
   WRITE (iuout, '("celldm(3)= ",f15.8)') tmod * alat / alat_box
   WRITE (iuout, '("nat= ",i5)') nat
   WRITE (iuout, '("ATOMIC_POSITIONS {alat}")') 
   DO na=1, nat
      phi = 2.0_DP * pi * y(1,na) 
      z = y(2,na) 
      WRITE(iuout, '(a, 3f18.10)' )  atm(na), &
                                 (r + y(3,na))*COS(phi)*alat/alat_box, &
                                 (r + y(3,na))*SIN(phi) * alat / alat_box, &
                                 z * tmod * alat / alat_box
   END DO
ENDIF

CLOSE(iuout)
CALL environment_end( code )
CALL mp_global_end ()


END PROGRAM gener_nanowire

