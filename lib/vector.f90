!
! Copyright (C) 2019 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE vector_mod
!
!   This module contains the routines for dealing with vectors.
!   Presently it has only a routine to print the vector components.
!
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: write_vector

CONTAINS

!-----------------------------------------------------------------------
SUBROUTINE write_vector(nvar, x)
!-----------------------------------------------------------------------
!
!   This routine writes on output a vector of dimension nvar
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
REAL(DP) :: x(nvar)

INTEGER :: i

DO i=1, nvar
   WRITE(stdout,'(23x,"x",i1,"=",f16.9)') i, x(i)
END DO

RETURN
END SUBROUTINE write_vector

END MODULE vector_mod
