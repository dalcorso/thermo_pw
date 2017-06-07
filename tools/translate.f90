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
!  and optionally a rotation matrix.
!  It applies the (roto-)translation to the atomic positions and writes them on
!  output. Both the translation vectors and the atomic positions must
!  be in the same units (alat or crystal), but the units are irrelevant.
!  The rotation matrix must be in the same units.
!  The output is in the same units as the input.
!  The input is the following:
!  a1,a2,a3 : the three components of the translation
!  0-1      : 0 if no rotation is given, 1 otherwise.
!  R_11, R_12, R_13  : the first row of the rotation matrix
!  R_21, R_22, R_23  : the second row of the rotation matrix
!  R_31, R_32, R_33  : the third row of the rotation matrix
!  nat               : the number of atoms
!  label, taux, tauy, tauz  : for each atom the name and the coordinates
!
!
USE kinds, ONLY : DP
USE rotate, ONLY : is_rotation, rotate_vect

IMPLICIT NONE
INTEGER :: nat
REAL(DP), ALLOCATABLE :: tau(:,:), tau0(:,:)
REAL(DP) :: a(3), rot(3,3), original_units, final_units, fact
INTEGER :: na, ipol, jpol, input_rot
CHARACTER(LEN=3), ALLOCATABLE :: label(:)

WRITE(6,'(5x," Translation vector? ")')
READ(5,*) a(1), a(2), a(3)
WRITE(6,'(3f15.8)') a(1), a(2), a(3)

WRITE(6,'(5x," Rotation matrix? (0 to skip) ")')
READ(5,*) input_rot
WRITE(6,'(i5)') input_rot

IF (input_rot == 0) THEN
   rot=0.0_DP
   rot(1,1)=1.0_DP
   rot(2,2)=1.0_DP
   rot(3,3)=1.0_DP
ELSE
   DO ipol=1,3
      READ(5,*) (rot(ipol,jpol), jpol=1,3)
      WRITE(6,'(3f15.8)') (rot(ipol,jpol), jpol=1,3)
   ENDDO

   IF (.NOT.is_rotation(rot)) THEN
      WRITE(6,'(/,5x,"WARNING: input matrix not a rotation")')
   END IF
END IF

WRITE(6,'(5x," Number of atoms and atomic coordinates")')
READ(5,*) nat
WRITE(6,'(i5)') nat 

ALLOCATE(tau0(3,nat))
ALLOCATE(tau(3,nat))
ALLOCATE(label(nat))
DO na=1, nat
   READ(5,*) label(na), tau0(1,na), tau0(2,na), tau0(3,na)
   WRITE(6,'(a3,3f20.12)') label(na), tau0(1,na), tau0(2,na), tau0(3,na)
ENDDO

WRITE(6,'(/,5x,"(Roto-)translated coordinates")') 
WRITE(6,*) nat
CALL rotate_vect(rot,nat, tau0, tau,1) 
DO na=1, nat
   WRITE(6,'(a3,3f20.12)') label(na), tau(1,na) + a(1), tau(2,na) + a(2), &
                                     tau(3,na) + a(3)
ENDDO

WRITE(6,'(5x," change of units? (original units - final units, &
                                & 0.0, 0.0 to skip)")')
READ(5,*) original_units, final_units
WRITE(6,'(2f16.7)') original_units, final_units

IF (original_units>0.0_DP .AND. final_units>0.0_DP) THEN
   fact=original_units/final_units
   DO na=1, nat
      WRITE(6,'(a3,3f20.12)') label(na), (tau(1,na) + a(1))*fact, &
                                         (tau(2,na) + a(2))*fact, &
                                         (tau(3,na) + a(3))*fact
   ENDDO
ENDIF

DEALLOCATE(tau)
DEALLOCATE(tau0)
DEALLOCATE(label)
END PROGRAM translate
