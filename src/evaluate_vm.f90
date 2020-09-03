!
! Copyright (C) 2015-2018 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------
SUBROUTINE evaluate_vm(temp_ph, vm, ntemp, temp, vminf_t)
!----------------------------------------------------------
!
!   This routine receives the volume as a function of temperature
!   in vminf_t, the set of temperature temp and their number ntemp
!   a given temperature temp_ph and finds the volume vm at temp_ph
!   via a linear interpolation of the two points closest to temp_ph.
!
USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), vminf_t(ntemp)
REAL(DP), INTENT(INOUT) :: temp_ph, vm

INTEGER :: itemp0, itemp1, itemp

itemp0=1
DO itemp=1,ntemp
   IF (temp(itemp) < temp_ph) itemp0=itemp
ENDDO

IF (itemp0 == ntemp) THEN
   WRITE(stdout,'(/,5x,"temp_ph too large setting to",f15.4," K")') &
                                                               temp(ntemp-1)
   temp_ph=temp(ntemp-1)
   vm=vminf_t(ntemp-1)
   RETURN
ENDIF

IF (itemp0 == 1) THEN
   WRITE(stdout,'(/,5x,"temp_ph too small setting to",f15.4," K")') temp(2)
   temp_ph=temp(2)
   vm=vminf_t(2)
   RETURN
ENDIF

itemp1=itemp0+1

vm = vminf_t(itemp0) + (temp_ph - temp(itemp0)) *          &
                       (vminf_t(itemp1)-vminf_t(itemp0)) / &
                       (temp(itemp1)-temp(itemp0))

RETURN
END SUBROUTINE evaluate_vm
