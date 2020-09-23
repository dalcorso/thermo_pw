!
! Copyright (C) 2014-2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------
SUBROUTINE set_temperature()
!------------------------------------------------------------------

USE kinds,          ONLY : DP
USE temperature,    ONLY : tmin, tmax, deltat, ntemp, temp

IMPLICIT NONE

INTEGER  :: itemp
!
!  Allocate thermodynamic quantities
!
IF (deltat <= 0.0_DP) CALL errore('set_temperature','Negative deltat',1)
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

!------------------------------------------------------------------
SUBROUTINE set_pressure()
!------------------------------------------------------------------
!
!  Set the pressure in Ry/(a.u.)^2 units. Input pressure is given in kbar.
!  This routine allocates also celldm0 and initialize it to the input geometry
!  for each pressure
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE control_pressure, ONLY : pressure, pressure_kb
IMPLICIT NONE
!
! save the pressure in kbar
!
pressure_kb=pressure
!
! convert the pressure in kbar
!
pressure=pressure_kb/ry_kbar

RETURN
END SUBROUTINE set_pressure

!------------------------------------------------------------------
SUBROUTINE add_pressure(filename)
!------------------------------------------------------------------
!
USE kinds, ONLY : DP
USE control_pressure, ONLY : pressure_kb

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=8) :: float_to_char

IF (pressure_kb /= 0.0_DP) &
   filename=TRIM(filename) //'.'//TRIM(float_to_char(pressure_kb,1))

RETURN
END SUBROUTINE add_pressure
