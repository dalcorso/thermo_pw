!
! Copyright (C) 2018 Cristiano Malica
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE write_anharm_bfact(temp, bfact_t, ntemp, filename)
!------------------------------------------------------------------------

USE ions_base, ONLY : nat
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ntemp 
REAL(DP), INTENT(IN) :: temp(ntemp), bfact_t(6,nat,ntemp)
CHARACTER(LEN=256) :: filename

INTEGER :: na, ijpol
INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

CHARACTER(LEN=6) :: int_to_char

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   DO na=1, nat
      OPEN(UNIT=iu_therm, FILE=TRIM(filename)//'.'//&
           TRIM(int_to_char(na))//'.dw', STATUS='UNKNOWN', FORM='FORMATTED')
      WRITE(iu_therm,'("# A^2")')

      WRITE(iu_therm,'("#",5x,"T (K)",10x,"B_11",14x,"B_12",14x,"B_13",14x,"B_22",&
                           &14x,"B_23",14x,"B_33")')

      DO itemp = 1, ntemp
      WRITE(iu_therm, '(e16.8,6e18.8)') temp(itemp), (bfact_t(ijpol,na,itemp),&
                                        ijpol=1,6)
      ENDDO
      CLOSE(UNIT=iu_therm, STATUS='KEEP')
   ENDDO
ENDIF
RETURN
END SUBROUTINE write_anharm_bfact

!------------------------------------------------------------------------
SUBROUTINE write_thermo_anharm(temp, ntemp, e0, ener_t, free_e_min_t, &
                                                 entropy_t, cv_t, filename)
!------------------------------------------------------------------------
USE kinds,     ONLY : DP
USE constants, ONLY : rytoev, electronvolt_si, rydberg_si, avogadro
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), PARAMETER :: caltoj=4.184_DP
REAL(DP), INTENT(IN) :: temp(ntemp), ener_t(ntemp), free_e_min_t(ntemp), &
                        cv_t(ntemp), entropy_t(ntemp), e0
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_therm,'("# Zero point energy:", f8.5, " Ry/cell,", f9.5, &
                 &" kJ/(N mol),", f9.5, " kcal/(N mol)")') e0, &
           e0 * rydberg_si*avogadro/1.D3, e0*rydberg_si*avogadro/caltoj/1.D3
   WRITE(iu_therm,'("# Temperature T in K, ")')
   WRITE(iu_therm,'("# Energy and free energy in Ry/cell,")')
   WRITE(iu_therm,'("# Entropy in Ry/cell/K,")')
   WRITE(iu_therm,'("# Heat capacity Cv in Ry/cell/K.")')
   WRITE(iu_therm,'("# Multiply by ",f7.4," to have energies in &
                       &eV/cell etc..")') rytoev
   WRITE(iu_therm,'("# Multiply by ",f7.4," x ",f8.2," = ",f9.1," to have &
                  &energies in cal/(N mol).")') rytoev, electronvolt_si &
                         * avogadro / caltoj, rydberg_si*avogadro/caltoj
   WRITE(iu_therm,'("# Multiply by ",f7.4," x ",f8.2," = ",f9.1," to &
                  &have energies in J/(N mol).")') rytoev, electronvolt_si&
                         * avogadro, rydberg_si*avogadro
   WRITE(iu_therm,'("# N is the number of formula units per cell.")')
   WRITE(iu_therm,'("# For instance in silicon N=2. Divide by N to have &
                &energies in cal/mol etc. ")')
   WRITE(iu_therm,'("#",5x,"   T  ", 10x, " energy ", 10x, "  free energy ",&
               & 12x, " entropy ", 17x, " Cv ")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,e23.13,e23.13,e23.13,e23.13)') temp(itemp), &
             ener_t(itemp), free_e_min_t(itemp), entropy_t(itemp), cv_t(itemp)
   END DO

   CLOSE(iu_therm)
ENDIF
RETURN
END SUBROUTINE write_thermo_anharm
!
! Copyright (C) 2022 Andrea Dal Corso
!
!------------------------------------------------------------------------
SUBROUTINE write_thermo_anharm_p(press, npress, e0_p, ener_p, free_e_min_p, &
                              entropy_p, cv_p, temp, ntemp, itemp, filename)

!------------------------------------------------------------------------
USE kinds,     ONLY : DP
USE constants, ONLY : rytoev, electronvolt_si, rydberg_si, avogadro
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER, INTENT(IN) :: npress, ntemp, itemp
REAL(DP), INTENT(IN) :: temp(ntemp)
REAL(DP), INTENT(IN) :: press(npress), ener_p(npress), free_e_min_p(npress), &
                        cv_p(npress), entropy_p(npress), e0_p(npress)
CHARACTER(LEN=*) :: filename
REAL(DP), PARAMETER :: caltoj=4.184_DP
INTEGER :: ipress, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_therm,'("# Temperature T= ",f15.7," K")') temp(itemp)
   WRITE(iu_therm,'("# Pressure p in kbar, ")')
   WRITE(iu_therm,'("# Energy and free energy in Ry/cell,")')
   WRITE(iu_therm,'("# Entropy in Ry/cell/K,")')
   WRITE(iu_therm,'("# Heat capacity Cv in Ry/cell/K.")')
   WRITE(iu_therm,'("# Multiply by ",f7.4," to have energies in &
                       &eV/cell etc..")') rytoev
   WRITE(iu_therm,'("# Multiply by ",f7.4," x ",f8.2," = ",f9.1," to have &
                  &energies in cal/(N mol).")') rytoev, electronvolt_si &
                         * avogadro / caltoj, rydberg_si*avogadro/caltoj
   WRITE(iu_therm,'("# Multiply by ",f7.4," x ",f8.2," = ",f9.1," to &
                  &have energies in J/(N mol).")') rytoev, electronvolt_si&
                         * avogadro, rydberg_si*avogadro
   WRITE(iu_therm,'("# N is the number of formula units per cell.")')
   WRITE(iu_therm,'("# For instance in silicon N=2. Divide by N to have &
                &energies in cal/mol etc. ")')
   WRITE(iu_therm,'("#",5x,"p(kbar)", 10x, " energy ", 10x, "  free energy ",&
               & 12x, " entropy ", 17x, " Cv ", 15x, "E_zpe")')

   DO ipress = 1, npress
      WRITE(iu_therm, '(e12.5,5e23.13)') press(ipress), ener_p(ipress), &
          free_e_min_p(ipress), entropy_p(ipress), cv_p(ipress), e0_p(ipress)
   END DO

   CLOSE(iu_therm)
ENDIF
RETURN
END SUBROUTINE write_thermo_anharm_p
!
