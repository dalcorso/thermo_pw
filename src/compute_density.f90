!
! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE compute_density(omega,density)
!
!  This routine receives the volume of a unit cell in (a.u.)^3,
!  the number of atoms inside the unit cell and the mass (in a.m.u.) of
!  each atom and gives as output the density in Kg/m^3 of the solid.
!
USE kinds, ONLY : DP
USE constants, ONLY : amu_si, bohr_radius_si
USE ions_base, ONLY : nat, ityp, atm, amass
USE io_global, ONLY : stdout

IMPLICIT NONE
REAL(DP), INTENT(IN) :: omega
REAL(DP), INTENT(OUT) :: density 

REAL(DP) :: total_mass, total_expected_mass, expected_mass
REAL(DP) :: current_mass, fact
REAL(DP) :: atom_weight
INTEGER  :: atomic_number
INTEGER :: ia, it
!
total_mass=0.0_DP
total_expected_mass=0.0_DP
DO ia=1,nat
   it=ityp(ia)
   expected_mass=atom_weight(atomic_number(TRIM(atm(it))))
   IF (amass(it)==0.0_DP) THEN
      current_mass=expected_mass
      total_mass=total_mass+current_mass
      total_expected_mass=total_expected_mass+current_mass
   ELSE
      current_mass=amass(it)
      IF (ABS(current_mass - expected_mass) > 1.0_DP) THEN
         IF (ia==1) WRITE(stdout,*)
         WRITE(stdout,'(5x,"Warning the mass of atom ",i5, f9.3,&
                          &" a.m.u. does not match its name ",a2)') &
                                    ia, amass(it), atm(it)
      ENDIF
      total_mass = total_mass + current_mass
      total_expected_mass=total_expected_mass+expected_mass
   ENDIF
ENDDO

fact = amu_si / (bohr_radius_si)**3 

IF (ABS(total_mass - total_expected_mass) > 1.0_DP) THEN
   WRITE(stdout,'(/,5x,"Total mass of this unit cell ",3x,f14.4," a.m.u.")') &
                                                   total_mass  
   WRITE(stdout,'(5x, "Expected mass of this unit cell ",f14.4," a.m.u.")') &
                                                   total_expected_mass  
   WRITE(stdout,'(5x, "Density of this solid ",9x,f15.2," kg/m^3",&
                         &f13.4," g/cm^3")') total_mass * fact / omega, &
                            total_mass * fact / omega / 1000._DP 
   WRITE(stdout,'(5x, "Expected density of this solid ", f15.2," kg/m^3",&
                   &f13.4," g/cm^3")') total_expected_mass *fact / omega, &
                         total_expected_mass * fact / omega / 1000._DP

ELSE
   WRITE(stdout,'(/,5x,"Total mass of this unit cell ",f15.4," a.m.u.")') &
                                    total_mass  
   WRITE(stdout,'(5x,"Density of this solid ",7x,f15.2," kg/m^3",&
                       &f15.4," g/cm^3")') total_mass * fact / omega, &
                                total_mass * fact / omega /1000._DP
ENDIF

density = total_mass * fact / omega

RETURN
END SUBROUTINE compute_density
