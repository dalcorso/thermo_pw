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
  USE kinds,            ONLY : DP
  USE control_mur,      ONLY : b0, b01, emin, lmurn
  USE equilibrium_conf, ONLY : celldm0
  USE control_pressure, ONLY : pressure_kb
  USE io_global,        ONLY : stdout
  IMPLICIT NONE
  CHARACTER(LEN=20) :: quantity
  !
  WRITE(stdout,'(/,2x,76("-"))')
  IF (pressure_kb /= 0.0_DP) THEN
     WRITE(stdout,'(5x,"At pressure ",f15.6," kbar")') pressure_kb
     quantity='enthalpy'
  ELSE
     WRITE(stdout,*) 
     quantity='total energy'
  ENDIF
  IF (lmurn) THEN
     WRITE(stdout,'(5x, "The equilibrium lattice constant is ",9x,f12.4,&
                                 &" a.u.")') celldm0(1)
     WRITE(stdout,'(5x, "The bulk modulus is ",24x,f12.3,"  kbar")') b0
     WRITE(stdout,'(5x, "The pressure derivative of the bulk modulus is ",&
                                  &f9.3)')  b01
  ELSE
     WRITE(stdout,'(5x,"The equilibrium celldm is:")')
     WRITE(stdout,'(5x,6f12.5)') celldm0(:)
  ENDIF
 
  WRITE(stdout,'(5x,"The ",a," at the minimum is:   ",6x,f20.9," Ry")') &
                                                   TRIM(quantity), emin
  WRITE(stdout,'(2x,76("-"),/)')
  RETURN
END SUBROUTINE write_minimum_energy_data
