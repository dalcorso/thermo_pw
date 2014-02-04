!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_temperature()

USE kinds,          ONLY : DP
USE temperature,    ONLY : tmin, tmax, deltat, ntemp, temp

IMPLICIT NONE

INTEGER  :: itemp
!
!  Allocate thermodynamic quantities
!
IF (deltat <= 0.0_8) CALL errore('print_thermo','Negative deltat',1)
ntemp=1+NINT((tmax-tmin)/deltat)

IF (.NOT.ALLOCATED(temp)) ALLOCATE(temp(ntemp))
!
!  Compute temperature
!
DO itemp = 1, ntemp
   temp(itemp) = tmin + (itemp-1) * deltat
END DO

RETURN
END SUBROUTINE set_temperature

