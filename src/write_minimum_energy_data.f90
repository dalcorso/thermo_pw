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
  USE constants,        ONLY : bohr_radius_si
  USE control_ev,       ONLY : ieos
  USE control_mur,      ONLY : b0, b01, b02, emin, lmurn
  USE equilibrium_conf, ONLY : celldm0
  USE lattices,         ONLY : celldm_name, needed_celldm
  USE initial_conf,     ONLY : ibrav_save
  USE control_pressure, ONLY : pressure_kb
  USE io_global,        ONLY : stdout
  IMPLICIT NONE
  CHARACTER(LEN=20) :: quantity
  LOGICAL :: celldm_in_use(6)
  REAL(DP) :: omega0, compute_omega_geo
  INTEGER :: i
  !
  WRITE(stdout,'(/,2x,76("-"))')
  IF (ieos==1) THEN
     WRITE(stdout,'(/,5x,"Birch-Murnaghan 3 order equation of state")')
  ELSEIF (ieos==2) THEN
     WRITE(stdout,'(/,5x,"Birch-Murnaghan 4 order equation of state")')
  ELSEIF (ieos==4) THEN
     WRITE(stdout,'(/,5x,"Murnaghan equation of state")')
  ELSE
     CALL errore("write_minimum_energy_data","wrong ieos",1)
  ENDIF
  IF (pressure_kb /= 0.0_DP) THEN
     WRITE(stdout,'(5x,"At pressure ",f15.6," kbar")') pressure_kb
     quantity='enthalpy'
  ELSE
     WRITE(stdout,*) 
     quantity='total energy'
  ENDIF
  IF (lmurn) THEN
     CALL needed_celldm(ibrav_save, celldm_in_use)
     WRITE(stdout,'(5x,"The equilibrium celldm is:")')
     WRITE(stdout,'(44x,a7,f12.5," a.u.")') celldm_name(1), celldm0(1)
     DO i=2,6
        IF (celldm_in_use(i)) &
           WRITE(stdout,'(44x,a7,f12.5)') celldm_name(i), celldm0(i)
     ENDDO
!     WRITE(stdout,'(5x, "The equilibrium lattice constant is ",9x,f12.4,&
!                                 &" a.u.")') celldm0(1)
     WRITE(stdout,'(5x, "The bulk modulus is:",24x,f12.3,"  kbar")') b0
     WRITE(stdout,'(5x, "The pressure derivative of the bulk modulus is:",&
                                  &f9.3)')  b01
     IF (ieos==2) &
        WRITE(stdout,'(5x, "The second derivative of the bulk &
                           &modulus is:",f13.5," 1/kbar")')  b02
  ELSE
     WRITE(stdout,'(5x,"The equilibrium celldm is:")')
     WRITE(stdout,'(5x,6f12.5)') celldm0(:)
  ENDIF
 
  WRITE(stdout,'(5x,"The ",a," at the minimum is:   ",6x,f20.9," Ry")') &
                                                   TRIM(quantity), emin
  omega0=compute_omega_geo(ibrav_save,celldm0)
  WRITE(stdout,'(5x, "The volume is ",9x,f13.5," (a.u.)^3",&
          f15.5," (A)^3")') omega0, omega0*bohr_radius_si**3/1.D-30
  WRITE(stdout,'(2x,76("-"),/)')
  RETURN
END SUBROUTINE write_minimum_energy_data
