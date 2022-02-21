!
! Copyright (C) 2014-2021 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------
SUBROUTINE set_temperature()
!------------------------------------------------------------------

USE kinds,          ONLY : DP
USE constants,      ONLY : eps6
USE temperature,    ONLY : tmin, tmax, deltat, ntemp, temp, &
                           ntemp_plot, itemp_plot, temp_plot, itemp300
IMPLICIT NONE
INTEGER  :: itemp, itempp
REAL(DP) :: dt, dtmin
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
!
!  Check and set all the temperatures to plot
!
DO itempp=1, ntemp_plot
   itemp_plot(itempp)=ntemp+1
   dtmin=1.D10
   DO itemp=2, ntemp-1
      dt=ABS(temp(itemp) - temp_plot(itempp)) 
      IF (dt<dtmin) THEN
         itemp_plot(itempp)=itemp
         dtmin=dt
      ENDIF
   ENDDO
   IF (itemp_plot(itempp)==(ntemp+1)) CALL errore('set_temperature', &
                                              'temp_plot not found', itempp)
   temp_plot(itempp)=temp(itemp_plot(itempp))
ENDDO
!
!  Find the temperature closest to room temperature (300 K)
!  Some quantities are normalized to the value at 300 K.
!
itemp300=ntemp
dtmin=1.D10
DO itemp=1, ntemp
   dt=ABS(temp(itemp) - 300.D0)
   IF (dt<dtmin) THEN
      itemp300=itemp
      dtmin=dt
   ENDIF
ENDDO
IF (itemp300==ntemp) itemp300=0

RETURN
END SUBROUTINE set_temperature
!
!------------------------------------------------------------------
SUBROUTINE set_pressure()
!------------------------------------------------------------------
!
USE kinds,          ONLY : DP
USE constants,      ONLY : eps6
USE control_pressure, ONLY : pmin, pmax, deltap, npress, press, &
                             npress_plot, press_plot, ipress_plot

IMPLICIT NONE
INTEGER  :: ipress, ipressp
REAL(DP) :: dpr, dpmin
!
!  Allocate thermodynamic quantities
!
IF (deltap <= 0.0_DP) CALL errore('set_pressure','Negative deltap',1)
npress=1+NINT((pmax-pmin)/deltap)

IF (.NOT.ALLOCATED(press)) ALLOCATE(press(npress))
!
!  Compute pressure
!
DO ipress = 1, npress
   press(ipress) = pmin + (ipress-1) * deltap
END DO
!
!  Check and sets all the pressures to plot
!
DO ipressp=1, npress_plot
   ipress_plot(ipressp)=npress+1
   dpmin=1.D10
   DO ipress=1, npress
      dpr=ABS(press(ipress) - press_plot(ipressp))
      IF (dpr<dpmin) THEN
         ipress_plot(ipressp)=ipress
         dpmin=dpr
      ENDIF
   ENDDO
   IF (ipress_plot(ipressp)==(npress+1)) CALL errore('set_pressure', &
                                              'press_plot not found', ipressp)
   press_plot(ipressp)=press(ipress_plot(ipressp))
ENDDO

RETURN
END SUBROUTINE set_pressure
!
!------------------------------------------------------------------
SUBROUTINE set_pressure_kb()
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
END SUBROUTINE set_pressure_kb

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

!------------------------------------------------------------------
SUBROUTINE add_value(filename, press)
!------------------------------------------------------------------
!
USE kinds, ONLY : DP

IMPLICIT NONE
REAL(DP) :: press
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=8) :: float_to_char

IF (ABS(press) > 1.D-8) &
   filename=TRIM(filename) //'.'//TRIM(float_to_char(press,1))

RETURN
END SUBROUTINE add_value
