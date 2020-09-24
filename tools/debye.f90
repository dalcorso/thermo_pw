!
! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
PROGRAM debye
!----------------------------------------------------------------------

USE kinds,        ONLY : DP
USE debye_module, ONLY : debye_e0, debye_vib_energy, debye_free_energy, &
                         debye_entropy, debye_cv, debye_b_factor
USE io_global,    ONLY : ionode, stdout
USE mp_global,    ONLY : mp_startup, mp_global_end
USE environment,  ONLY : environment_start, environment_end

IMPLICIT NONE
CHARACTER(LEN=9)   :: code='Debye'
CHARACTER(LEN=256) :: filename

INTEGER  :: itemp, nat, nsp
INTEGER  :: ntemp

REAL(DP), ALLOCATABLE :: temp(:), deb_cv(:), deb_entropy(:), &
                         deb_energy(:), deb_free_energy(:), deb_bfact(:)

REAL(DP) :: debye_t, deb_e0, amass(1), tmin, tmax, deltat
!
CALL mp_startup (start_images=.TRUE.)
CALL environment_start(code)

WRITE(stdout,'(5x,"Debye temperature (K)? ")')
READ(5,*) debye_t

WRITE(stdout,'(5x,"Number of atoms? ")')
READ(5,*) nat

WRITE(stdout,'(5x,"Number of atom types? ")')
READ(5,*) nsp

WRITE(stdout,'(5x,"Atomic mass (a.m.u.)? ")')
READ(5,*) amass(1)

WRITE(stdout,'(5x,"Minimum temperature (K)? ")')
READ(5,*) tmin

WRITE(stdout,'(5x,"Maximum temperature (K)? ")')
READ(5,*) tmax

WRITE(stdout,'(5x,"Delta T (K)? ")')
READ(5,*) deltat

WRITE(stdout,'(5x,"Output file? ")')
READ(5,*) filename

WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(5x,"Computing the thermodynamic properties from elastic &
                                                              &constants")')
WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(filename)
WRITE(stdout,'(2x,76("+"),/)')

ntemp= (tmax-tmin)/deltat

ALLOCATE( temp(ntemp) )
ALLOCATE( deb_energy(ntemp) )
ALLOCATE( deb_free_energy(ntemp) )
ALLOCATE( deb_entropy(ntemp) ) 
ALLOCATE( deb_cv(ntemp) )
ALLOCATE( deb_bfact(ntemp) )

DO itemp=1, ntemp
   temp(itemp)= tmin + (itemp-1) * deltat
ENDDO
!
!  Compute thermodynamic quantities
!
CALL debye_e0 (debye_t, nat, deb_e0)
CALL debye_cv (debye_t, temp, ntemp, nat, deb_cv)
CALL debye_vib_energy (debye_t, temp, ntemp, nat, deb_energy)
CALL debye_free_energy (debye_t, temp, ntemp, nat, deb_free_energy)
CALL debye_entropy (debye_t, temp, ntemp, nat, deb_entropy)
IF (nsp==1) CALL debye_b_factor(debye_t, temp, ntemp, amass(1), deb_bfact)
!
!  Add the zero point energy
!
DO itemp=1,ntemp
   deb_energy(itemp) = deb_energy(itemp) + deb_e0
   deb_free_energy(itemp) = deb_free_energy(itemp) + deb_e0
END DO
!
!  Write on file
!
IF (ionode) THEN
   CALL write_thermo_info(deb_e0, debye_t, ntemp, temp, deb_energy, &
               deb_free_energy, deb_entropy, deb_cv, 3, filename)
   IF (nsp==1) CALL write_dw_debye(ntemp, temp, deb_bfact, filename)
END IF

DEALLOCATE( temp )
DEALLOCATE( deb_cv )
DEALLOCATE( deb_entropy ) 
DEALLOCATE( deb_free_energy )
DEALLOCATE( deb_energy )
DEALLOCATE( deb_bfact )

CALL environment_end(code)
CALL mp_global_end ()

END PROGRAM debye

!----------------------------------------------------------------------
SUBROUTINE write_thermo_info(e0, tot_states, ntemp, temp, energy, &
                                 free_energy, entropy, cv, iflag, filename)
!----------------------------------------------------------------------
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, iflag
REAL(DP), INTENT(IN) :: e0, tot_states, temp(ntemp), energy(ntemp),  &
                        free_energy(ntemp), entropy(ntemp), cv(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: iu_therm, itemp
INTEGER :: find_free_unit

iu_therm=find_free_unit()
OPEN (UNIT=iu_therm, FILE=TRIM(filename), STATUS='unknown',&
                                                     FORM='formatted')
WRITE(iu_therm,'("# Zero point energy:", f8.5, " Ry/cell,", f9.5, &
                 &" kJ/(N mol),", f9.5, " kcal/(N mol)")') e0, &
                    e0 * 1313.313_DP, e0 * 313.7545_DP 
IF (iflag==3) THEN
   WRITE(iu_therm,'("# Temperature T in K, Debye temperature=",f12.3, &
                                          &" K,")') tot_states
ELSE
   WRITE(iu_therm,'("# Temperature T in K, ")')
ENDIF
IF (iflag==1) THEN
   WRITE(iu_therm,'("# Total number of states is:", f15.5,",")') tot_states
ELSE
   WRITE(iu_therm,'("# ")') 
ENDIF
WRITE(iu_therm,'("# Energy and free energy in Ry/cell,")')
WRITE(iu_therm,'("# Entropy in Ry/cell/K,")')
WRITE(iu_therm,'("# Heat capacity Cv in Ry/cell/K.")')
WRITE(iu_therm,'("# Multiply by 13.6058 to have energies in &
                       &eV/cell etc..")')
WRITE(iu_therm,'("# Multiply by 13.6058 x 23060.35 = 313 754.5 to have &
                  &energies in cal/(N mol).")')
WRITE(iu_therm,'("# Multiply by 13.6058 x 96526.0 = 1 313 313 to &
                  &have energies in J/(N mol).")')
WRITE(iu_therm,'("# N is the number of formula units per cell.")')
WRITE(iu_therm,'("# For instance in silicon N=2. Divide by N to have &
                &energies in cal/mol etc. ")')
WRITE(iu_therm,'("#",5x,"   T  ", 10x, " energy ", 9x, "  free energy ",&
               & 9x, " entropy ", 12x, " Cv ")')

DO itemp = 1, ntemp
   WRITE(iu_therm, '(e16.8,4e20.12)') temp(itemp), energy(itemp), &
                 free_energy(itemp), entropy(itemp), cv(itemp)
END DO
CLOSE(UNIT=iu_therm, STATUS='KEEP')

RETURN
END SUBROUTINE write_thermo_info
!
!  Copyright (C) 2018 Cristiano Malica
!
SUBROUTINE write_dw_debye(ntemp, temp, deb_bfact, filename)
USE kinds, ONLY : DP

IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), deb_bfact(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: iu_therm, itemp, na, ipol, jpol
INTEGER :: find_free_unit

iu_therm=find_free_unit()

OPEN (UNIT=iu_therm, FILE=TRIM(filename)//'.dw', STATUS='unknown', &
                                                 FORM='formatted')
WRITE(iu_therm,'("#",6x,"T",14x,"B_33")')

DO itemp = 1, ntemp
   WRITE(iu_therm, '(e16.8,6e18.8)') temp(itemp), deb_bfact(itemp) 
END DO
CLOSE(UNIT=iu_therm, STATUS='KEEP')

RETURN
END SUBROUTINE write_dw_debye


