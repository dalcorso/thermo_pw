!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM translate
!
!  This program reads a translation vector and a set of atomic positions 
!  It applies the translation to the atomic positions and writes them on
!  output. Both the translation vectors and the atomic positions must
!  be in the same units (alat or crystal), but the units are irrelevant.
!  The output is in the same units as the input.
!
USE kinds, ONLY : DP

IMPLICIT NONE
INTEGER :: nat
REAL(DP), ALLOCATABLE :: tau(:,:)
REAL(DP) :: a(3)
INTEGER :: na 
CHARACTER(LEN=3), ALLOCATABLE :: label(:)

WRITE(6,'(5x," Translation vector? ")')
READ(5,*) a(1), a(2), a(3)

WRITE(6,'(5x," Number of atoms and atomic coordinates")')
READ(5,*) nat
WRITE(6,'(i5)') nat 

ALLOCATE(tau(3,nat))
ALLOCATE(label(nat))

DO na=1, nat
   READ(5,*) label(na), tau(1,na), tau(2,na), tau(3,na)
   WRITE(6,'(a3,3f20.12)') label(na), tau(1,na) + a(1), tau(2,na) + a(2), &
                                     tau(3,na) + a(3)
ENDDO

DEALLOCATE(tau)
DEALLOCATE(label)
END PROGRAM translate
