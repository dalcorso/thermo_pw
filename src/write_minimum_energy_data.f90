!
! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE write_minimum_energy_data()
  !-----------------------------------------------------------------------
  !
  !  This routine writes on output the information on the (Gibbs) energy 
  !  minimization.
  !
  !
  USE kinds, ONLY : DP
  USE control_mur,      ONLY : b0, b01, emin, celldm0, lmurn
  USE control_pressure, ONLY : pressure_kb
  USE cell_base,        ONLY : ibrav
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  !
  WRITE(stdout,'(/,2x,76("-"))')
  IF (pressure_kb /= 0.0_DP) THEN
     WRITE(stdout,'(5x,"At pressure ",f15.6," kbar")') pressure_kb
  ELSE
     WRITE(stdout,*) 
  ENDIF
  IF (ibrav==1.OR.ibrav==2.OR.ibrav==3) THEN
     WRITE(stdout,'(5x, "The equilibrium lattice constant is ",9x,f12.4,&
                                 &" a.u.")') celldm0(1)
  ELSE
     IF (pressure_kb /= 0.0_DP) THEN
        WRITE(stdout,'(5x,"The minimum enthalpy is obtained for celldm")')
     ELSE
        WRITE(stdout,'(5x,"The minimum energy is obtained for celldm")')
     ENDIF
     WRITE(stdout,'(5x,6f12.5)') celldm0(:)
  ENDIF
 
  IF (lmurn) THEN
     WRITE(stdout,'(5x, "The bulk modulus is ",24x,f12.3,"  kbar")') b0
     WRITE(stdout,'(5x, "The pressure derivative of the bulk modulus is ",&
                                  &f9.3)')  b01
  END IF
  IF (pressure_kb /= 0.0_DP) THEN
     WRITE(stdout,'(5x,"The enthalpy at the minimum is    ",6x,f20.9," Ry")') &
                                                               emin
  ELSE
     WRITE(stdout,'(5x,"The total energy at the minimum is",6x,f20.9," Ry")') &
                                                 emin
  END IF

  WRITE(stdout,'(2x,76("-"),/)')
  RETURN
END SUBROUTINE write_minimum_energy_data
