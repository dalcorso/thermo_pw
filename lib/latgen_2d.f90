!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE latgen_2d(ibrav_2d, celldm_2d, a1, a2)
!----------------------------------------------------------------------
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav_2d
REAL(DP), INTENT(IN) :: celldm_2d(3)
REAL(DP), INTENT(OUT) :: a1(3), a2(3)

a1=0.0_DP
a2=0.0_DP

SELECT CASE(ibrav_2d) 
    CASE(1)
!
!   rhombic
!
       a1(1)=1.0_DP
       a2(1)=celldm_2d(2) * celldm_2d(3)
       a2(2)=celldm_2d(2) * SQRT( 1.0_DP - celldm_2d(3)**2 )
    CASE(2)
!
!  rectangular
!
       a1(1)=1.0_DP
       a2(2)=celldm_2d(2)
    CASE(3)
!
!  centered rectangular
!
       a1(1)=0.5_DP
       a1(2)=0.5_DP * celldm_2d(2)
       a2(1)=-0.5_DP
       a2(2)=0.5_DP * celldm_2d(2)

    CASE(4)
!
!  square lattice
!
       a1(1)=1.0_DP
       a2(2)=1.0_DP
    CASE(5)
!
!  hexagonal lattice
!
       a1(1)=1.0_DP
       a2(1)=-0.5_DP
       a2(2)=SQRT(3.0_DP) * 0.5_DP
END SELECT

RETURN
END SUBROUTINE latgen_2d

!----------------------------------------------------------------------
SUBROUTINE recips_2d(a1,a2,b1,b2)
!----------------------------------------------------------------------
!
!  This routine receives as input the direct lattice vectors of a 2d 
!  Bravais lattice and finds the reciprocal lattice vectors b1 and b2
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP),  INTENT(IN) :: a1(3), a2(3)
REAL(DP),  INTENT(OUT) :: b1(3), b2(3)

REAL(DP) :: det

det = a1(1) * a2(2) - a1(2) * a2(1)

b1(1) =   a2(2) / det
b1(2) = - a2(1) / det

b2(1) = - a1(2) / det
b2(2) =   a1(1) / det

RETURN
END SUBROUTINE recips_2d
