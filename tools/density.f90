!
! Copyright (C) 2023 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM density
!
!  This program reads in input the volume of a solid in (a.u.)^3
!  the mass of one unit cell in a.m.u. and computes the density 
!  in g/cm^3 and in kg/m^3. Alternatively, it reads the density
!  in g/cm^2 and the mass of one unit cell in a.m.u. and prints
!  the volume of a unit cell in (a.u.)^3 and in A^3.
!
!  The input variables are:
!
!  working_mode     1 input volume -> output density
!                   2 input density -> output volume
!  When working_mode is 1:
!       volume : the volume of a unit cell in (a.u.)^3
!       mass : the mass of the unit cell in a.m.u.
!  When working_mode is 2:
!       density : the density in g/cm^3
!       mass : the mass of the unit cell in a.m.u.
!
!
USE kinds,       ONLY : DP
USE constants,   ONLY : amu_si, bohr_radius_si
USE mp_global,   ONLY : mp_startup, mp_global_end
USE environment, ONLY : environment_start, environment_end
USE io_global,   ONLY : stdout

IMPLICIT NONE

INTEGER  :: iwork
REAL(DP) :: volume, mass, fact, dens
INTEGER  :: stdin
CHARACTER(LEN=9) :: code='density'

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

stdin=5
WRITE(stdout,'(5x,"Choose what to do: ")')
WRITE(stdout,'(5x,"1: Volume -> Density")')
WRITE(stdout,'(5x,"2: Density -> Volume")')
READ(stdin,*) iwork 
IF (iwork==1) THEN
   WRITE(stdout,'(5x,"Volume (a.u.)^3")') 
   READ(stdin,*) volume
ELSE
   WRITE(stdout,'(5x,"Density (g/cm^3)")') 
   READ(stdin,*) dens
ENDIF
WRITE(stdout,'(5x,"Mass of the unit cell in a.m.u.")')
READ(stdin,*) mass

IF (iwork==1) THEN
   fact = amu_si / (bohr_radius_si)**3
   WRITE(stdout,'(/,5x,"The density is",f15.6," Kg/m^3")') fact*mass/volume 
   WRITE(stdout,'(5x,"The density is",f15.6," g/cm^3")') fact * mass / &
                                                           volume / 1.D3
ELSE
   fact = amu_si / (bohr_radius_si)**3/1.D3
   WRITE(stdout,'(/,5x,"The volume is",f15.6," (a.u.)^3")') fact*mass/dens 
   fact = amu_si * 1.D30 / 1.D3
   WRITE(stdout,'(5x,"The volume is",f15.6," A^3")') fact*mass/dens 
ENDIF

CALL environment_end( code )
CALL mp_global_end ()
END PROGRAM density
