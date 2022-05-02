PROGRAM average_grun
USE kinds, ONLY : DP
USE constants, ONLY : avogadro

IMPLICIT NONE
REAL(DP) :: alpha, bulk, volume, cv, temp

WRITE(6,'(5x,"alpha x 1.E6 (1/K) ?")') 
READ(5,*) alpha
WRITE(6,'(5x,"Isothermal bulk modulus (kbar) ?")') 
READ(5,*) bulk
WRITE(6,'(5x,"Volume (A^3) ?")') 
READ(5,*) volume
WRITE(6,'(5x,"Heat capacity (J/mol/K) ?")') 
READ(5,*) cv

WRITE(6,'(/,5x,"Average Gruneisen parameter",f15.3)')  &
                              alpha*1.D-6*bulk*1.D8*volume* &
                              1.D-30*avogadro/cv

WRITE(6,'(5x,"Temperature ?")') 
READ(5,*) temp

WRITE(6,'(/,5x,"Isobaric heat capacity",f20.3," J/mol/K")')  &
                              cv + (alpha*1.D-6)**2*bulk*1.D8*volume* &
                              1.D-30*temp*avogadro

WRITE(6,'(/,5x,"Adiabatic bulk modulus",f20.1," kbar")')  &
                              bulk + (alpha*1.D-6*bulk*1.D8)**2*volume* &
                              1.D-30*avogadro*temp/cv/1.D8

END PROGRAM average_grun
