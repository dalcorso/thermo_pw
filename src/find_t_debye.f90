!
! Copyright (C) 2022 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE find_t_debye(free_energy, temp, ntemp, thetad)
!----------------------------------------------------------------------
!
!   This routine receives the Helmholtz free energy in a mesh of ntemp 
!   temperatures temp. For each temperature it gives on output the 
!   Debye temperature (thetad) that used in a Debye model gives 
!   that free energy at temperature temp.
!   Debye temperature is determined with four digits. Increase temp_err
!   if you need more accurate values.
!
USE kinds,        ONLY : DP
USE constants,    ONLY : k_boltzmann_ry
USE ions_base,    ONLY : nat
USE debye_module, ONLY : debye_free_energy_0d

IMPLICIT NONE
INTEGER, INTENT(IN)     :: ntemp
REAL(DP), INTENT(IN)    :: free_energy(ntemp), temp(ntemp)
REAL(DP), INTENT(INOUT) :: thetad(ntemp)

INTEGER, PARAMETER  :: maxiter=500
INTEGER             :: itemp, iter
REAL(DP), PARAMETER :: temp_err=1.D-4
REAL(DP)            :: tdl, tdr, tdc, fel, fer, fec, fet, tt

DO itemp=1, ntemp
   fet=free_energy(itemp)
   tt=temp(itemp)
   tdl=0.1_DP
   CALL debye_free_energy_0d (tdl, tt, nat, fel)
   tdr=8.D3
   CALL debye_free_energy_0d (tdr, tt, nat, fer)
   IF ((fel-fet)*(fer-fet) > 0.0_DP) &
                CALL errore('find_t_debye','Problem with bisection',1)
   DO iter=1, maxiter
      tdc=(tdl+tdr)*0.5_DP
      CALL debye_free_energy_0d (tdc, tt, nat, fec)
!      WRITE(6,'(i5,2f15.8,3e15.6)') iter, tdl, tdc, (fel-fet)*(fec-fet),&
!                  (fer-fet)*(fec-fet), fet
      IF ((fel-fet)*(fec-fet)>0.0_DP) THEN
         tdl=tdc
         fel=fec
      ELSEIF ((fer-fet)*(fec-fet)>0.0_DP) THEN
         tdr=tdc
         fer=fec
      ELSE
         CALL errore('find_t_debye','Problem1 with bisection',1)
      ENDIF
      IF (ABS(tdr-tdl)< temp_err) GOTO 100
   ENDDO
   CALL errore('find_t_debye','Maximum iteration number exceeded',1)
100 thetad(itemp)=tdl
ENDDO

RETURN
END SUBROUTINE find_t_debye
!
!----------------------------------------------------------------------
SUBROUTINE find_t_debye_cv(cv, temp, ntemp, thetad)
!----------------------------------------------------------------------
!
!   This routine receives the heat capacity in a mesh of ntemp 
!   temperatures temp. For each temperature it gives on output the 
!   Debye temperature (thetad) that used in a Debye model gives 
!   that heat capacity at temperature temp. 
!   Debye temperature is determined with four digits. Increase 
!   temp_err if you need more accurate values.
!
USE kinds,        ONLY : DP
USE constants,    ONLY : k_boltzmann_ry
USE ions_base,    ONLY : nat
USE debye_module, ONLY : debye_cv_0d

IMPLICIT NONE
INTEGER, INTENT(IN)     :: ntemp
REAL(DP), INTENT(IN)    :: cv(ntemp), temp(ntemp)
REAL(DP), INTENT(INOUT) :: thetad(ntemp)

INTEGER, PARAMETER  :: maxiter=500
INTEGER             :: itemp, iter
REAL(DP), PARAMETER :: temp_err=1.D-4
REAL(DP)            :: tdl, tdr, tdc, cvl, cvr, cvc, cvt, tt

DO itemp=1, ntemp
   cvt=cv(itemp)
   tt=temp(itemp)
   tdl=0.1_DP
   CALL debye_cv_0d (tdl, tt, nat, cvl)
   tdr=8.D3
   CALL debye_cv_0d (tdr, tt, nat, cvr)
   IF ((cvl-cvt)*(cvr-cvt) > 0.0_DP) CALL errore('find_t_debye_cv',&
                                             'Problem with bisection',1)
   DO iter=1, maxiter
      tdc=(tdl+tdr)*0.5_DP
      CALL debye_cv_0d (tdc, tt, nat, cvc)
!      WRITE(6,'(i5,2f15.8,3e15.6)') iter, tdl, tdc, (cvl-cvt)*(cvc-cvt),&
!                  (cvr-cvt)*(cvc-cvt), cvt
      IF ((cvl-cvt)*(cvc-cvt)>0.0_DP) THEN
         tdl=tdc
         cvl=cvc
      ELSEIF ((cvr-cvt)*(cvc-cvt)>0.0_DP) THEN
         tdr=tdc
         cvr=cvc
      ELSE
         CALL errore('find_t_debye_cv','Problem1 with bisection',1)
      ENDIF
      IF (ABS(tdr-tdl)< temp_err) GOTO 100
   ENDDO
   CALL errore('find_t_debye_cv','Maximum iteration number exceeded',1)
100 thetad(itemp)=tdl 
ENDDO

RETURN
END SUBROUTINE find_t_debye_cv
!
!----------------------------------------------------------------------
SUBROUTINE find_t_debye_ene(energy, temp, ntemp, thetad)
!----------------------------------------------------------------------
!
!   This routine receives the vibrational energy in a mesh of ntemp 
!   temperatures temp. For each temperature it gives on output the 
!   Debye temperature (thetad) that used in a Debye model gives that 
!   energy at temperature temp.
!   Debye temperature is determined with four digits. Increase temp_err
!   if you need more accurate values.
!
USE kinds,        ONLY : DP
USE constants,    ONLY : k_boltzmann_ry
USE ions_base,    ONLY : nat
USE debye_module, ONLY : debye_energy_0d

IMPLICIT NONE
INTEGER, INTENT(IN)     :: ntemp
REAL(DP), INTENT(IN)    :: energy(ntemp), temp(ntemp)
REAL(DP), INTENT(INOUT) :: thetad(ntemp)

INTEGER, PARAMETER  :: maxiter=500
INTEGER             :: itemp, iter
REAL(DP), PARAMETER :: temp_err=1.D-4
REAL(DP)            :: tdl, tdr, tdc, enl, enr, enc, ent, tt

DO itemp=1, ntemp
   ent=energy(itemp)
   tt=temp(itemp)
   tdl=0.1_DP
   CALL debye_energy_0d (tdl, tt, nat, enl)
   tdr=8.D3
   CALL debye_energy_0d (tdr, tt, nat, enr)
   IF ((enl-ent)*(enr-ent) > 0.0_DP) &
                CALL errore('find_t_debye_ene','Problem with bisection',1)
   DO iter=1, maxiter
      tdc=(tdl+tdr)*0.5_DP
      CALL debye_energy_0d (tdc, tt, nat, enc)
!      WRITE(6,'(i5,2f15.8,3e15.6)') iter, tdl, tdc, (enl-ent)*(enc-ent),&
!                  (enr-ent)*(enc-ent), ent
      IF ((enl-ent)*(enc-ent)>0.0_DP) THEN
         tdl=tdc
         enl=enc
      ELSEIF ((enr-ent)*(enc-ent)>0.0_DP) THEN
         tdr=tdc
         enr=enc
      ELSE
         CALL errore('find_t_debye_ene','Problem1 with bisection',1)
      ENDIF
      IF (ABS(tdr-tdl)< temp_err) GOTO 100
   ENDDO
   CALL errore('find_t_debye_ene','Maximum iteration number exceeded',1)
100 thetad(itemp)=tdl
ENDDO

RETURN
END SUBROUTINE find_t_debye_ene
