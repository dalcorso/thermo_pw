!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE ph_freq_module
!
!  This module provide methods to calculate their vibrational contribution 
!  to the free energy, energy, entropy and constant volume specific heat.
!  It defines a type ph_freq_type that contains the phonon frequencies on 
!  a uniform mesh of k points.
!
!
USE kinds, ONLY : DP
USE mp_images, ONLY : intra_image_comm
USE mp, ONLY : mp_sum
IMPLICIT NONE
SAVE
PRIVATE
REAL(DP), PARAMETER :: kb=8.6173324d-5/13.6058d0 ! Boltzmann constant in Ry/K
REAL(DP), PARAMETER :: ry_to_cmm1= 8065.5d0 * 13.6058d0
REAL(DP), PARAMETER :: kb1=1.0d0/8065.5d0/8.6173324d-5 ! inverse Boltzmann 
                                                       ! constant in cm^{-1}/K

TYPE ph_freq_type
   INTEGER :: number_of_points        ! the total number of frequencies
   INTEGER :: nqtot                   ! the total number of q points
   INTEGER :: nq                      ! the points reduced by symmetry
   INTEGER :: nq1, nq2, nq3           ! the number of points in the three
                                      ! reciprocal lattice directions
   INTEGER :: nat                     ! the number of atoms
   REAL(DP), ALLOCATABLE :: nu(:,:)   ! the frequencies (cm-1)
   REAL(DP), ALLOCATABLE :: wg(:)     ! the weight of each point
   INTEGER, ALLOCATABLE :: rap(:,:)   ! the number of the representation
END TYPE ph_freq_type


PUBLIC :: ph_freq_type, zero_point_energy_ph, free_energy_ph, vib_energy_ph, &
          vib_entropy_ph, specific_heat_cv_ph, init_ph_freq, destroy_ph_freq, &
          init_ph_rap, thermal_expansion_ph

CONTAINS

SUBROUTINE init_ph_freq(ph_freq, nat, nq1, nq2, nq3, nq)
IMPLICIT NONE

TYPE(ph_freq_type), INTENT(INOUT) :: ph_freq
INTEGER, INTENT(IN) :: nq1, nq2, nq3, nat, nq
INTEGER :: ndiv, nqtot

nqtot = nq1 * nq2 * nq3
ph_freq%nqtot=nqtot
ph_freq%nq=nq
ph_freq%number_of_points = 3 * nat * nqtot
ph_freq%nq1 = nq1
ph_freq%nq2 = nq2
ph_freq%nq3 = nq3
ph_freq%nat = nat
ndiv = ph_freq%number_of_points 
ALLOCATE(ph_freq%nu(3*nat, nq))
ALLOCATE(ph_freq%wg(nq))

RETURN
END SUBROUTINE init_ph_freq

SUBROUTINE init_ph_rap(ph_freq)
IMPLICIT NONE
TYPE(ph_freq_type), INTENT(INOUT) :: ph_freq
INTEGER :: nq, nat

nq = ph_freq%nq 
nat = ph_freq%nat
ALLOCATE(ph_freq%rap(3*nat,nq))

RETURN
END SUBROUTINE init_ph_rap

SUBROUTINE destroy_ph_freq(ph_freq)
IMPLICIT NONE
TYPE(ph_freq_type), INTENT(INOUT) :: ph_freq

IF (ALLOCATED(ph_freq%nu)) DEALLOCATE(ph_freq%nu)
IF (ALLOCATED(ph_freq%wg)) DEALLOCATE(ph_freq%wg)
IF (ALLOCATED(ph_freq%rap)) DEALLOCATE(ph_freq%rap)

RETURN
END SUBROUTINE destroy_ph_freq

SUBROUTINE zero_point_energy_ph(ph_freq, ener)
!
!  This subroutine receives as input a set of phonon frequencies and computes 
!  the zero point energy that corresponds to that frequencies. The output 
!  energy is in Ry.
!
TYPE(ph_freq_type), INTENT(IN) :: ph_freq
REAL(DP), INTENT(OUT) :: ener
REAL(DP) :: wg
INTEGER :: nq, iq, imode, nat

nq=ph_freq%nq
nat=ph_freq%nat

ener=0.0_DP
DO iq=1, nq
   wg= ph_freq%wg(iq)
   DO imode=1, 3*nat
      ener = ener + 0.5_DP * ph_freq%nu(imode, iq) * wg
   ENDDO
ENDDO
! result is in cm^{-1}, bring it to Ry
ener = ener / ry_to_cmm1 

RETURN
END SUBROUTINE zero_point_energy_ph

SUBROUTINE free_energy_ph(ph_freq, temp, free_ener)
!
!  This routine receives as input a phdos and a temperature and gives as 
!  output the vibrational free energy at that temperature. ener contains
!  only the vibrational contribution WITHOUT the zero point energy. 
!  
!
TYPE(ph_freq_type), INTENT(IN) :: ph_freq
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(OUT) :: free_ener

INTEGER :: nq, iq, nat, imode, startq, lastq
REAL(DP) :: nu, arg, wg, temp1, onesixth

free_ener=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
nq=ph_freq%nq
nat=ph_freq%nat
onesixth=1.0_DP / 6.0_DP
CALL divide (intra_image_comm, nq, startq, lastq)
DO iq=startq,lastq
   wg=ph_freq%wg(iq)
   DO imode=1,3*nat
      nu=ph_freq%nu(imode,iq)
      arg= kb1 * nu * temp1
      IF (arg > 1.d-5) THEN
         IF (arg < 650_DP) &
            free_ener = free_ener +kb*temp*wg*(-arg + LOG(EXP(arg)-1.0_DP))
      ELSEIF (nu > 1.D-3) THEN
         free_ener = free_ener+kb*temp*wg*(-arg + LOG( arg + arg**2*0.5_DP + &
                                                         arg**3 * onesixth )) 
         write(6,*) 'free_energy, using asymptote', temp, nu, arg
      ENDIF
   ENDDO
ENDDO
CALL mp_sum(free_ener, intra_image_comm)

RETURN
END SUBROUTINE free_energy_ph

SUBROUTINE vib_energy_ph(ph_freq, temp, ener)
!
!  This routine receives as input a set of phonon frequencies and a 
!  temperature and gives as 
!  output the vibrational energy at that temperature. ener contains
!  the energy WITHOUT the zero point energy. 
!  
!
TYPE(ph_freq_type), INTENT(IN) :: ph_freq
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(OUT) :: ener

INTEGER :: nq, iq, imode, nat, startq, lastq
REAL(DP) :: nu, temp1, arg, wg, onesixth

ener=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
onesixth=1.0_DP / 6.0_DP
nq=ph_freq%nq
nat=ph_freq%nat
CALL divide(intra_image_comm, nq, startq, lastq)
DO iq=startq,lastq
   wg=ph_freq%wg(iq)
   DO imode=1, 3*nat
      nu=ph_freq%nu(imode,iq)
      arg= kb1 * nu * temp1
      IF (arg > 1.d-5) THEN
        IF (arg < 650._DP) ener = ener +  wg * nu / ( EXP( arg ) - 1.0_DP ) 
      ELSEIF ( nu > 1.D-3 ) THEN
         ener = ener +  wg * nu / ( arg + arg**2*0.5_DP + arg**3 * onesixth )
         write(6,*) 'energy using asymptote', temp, nu, arg
      ENDIF
   ENDDO
ENDDO
ener = ener / ry_to_cmm1
CALL mp_sum(ener, intra_image_comm)

RETURN
END SUBROUTINE vib_energy_ph

SUBROUTINE vib_entropy_ph(ph_freq, temp, entr)
!
!  This routine receives as input a phdos and a temperature and gives as 
!  output the vibrational entropy at that temperature. 
!  
!
TYPE(ph_freq_type), INTENT(IN) :: ph_freq
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(OUT) :: entr
REAL(DP) :: ener, free_ener

CALL free_energy_ph(ph_freq, temp, free_ener)
CALL vib_energy_ph(ph_freq, temp, ener)

IF (temp > 0.0_DP) THEN
   entr = ( ener - free_ener ) / temp
ELSE
   entr = 0.0_DP
ENDIF

RETURN
END SUBROUTINE vib_entropy_ph

SUBROUTINE specific_heat_cv_ph(ph_freq, temp, cv)
!
!  This routine receives as input a set of phonon frequencies and a 
!  temperature and gives as output the constant volume specific heat at 
!  that temperature. 
!  The output cv is in Ry / K.
!
TYPE(ph_freq_type), INTENT(IN) :: ph_freq
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(OUT) :: cv

INTEGER :: nq, iq, imode, nat, startq, lastq, icount
REAL(DP) :: nu, temp1, arg, wg, onesixth

cv=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
onesixth=1.0_DP / 6.0_DP
temp1 = 1.0_DP / temp
nq=ph_freq%nq
nat=ph_freq%nat
CALL divide(intra_image_comm, nq, startq, lastq)
icount=0
DO iq=startq,lastq
   wg=ph_freq%wg(iq)
   DO imode=1,3*nat
      nu=ph_freq%nu(imode,iq)
      arg= kb1 * nu * temp1
      IF (arg > 1.D-6 ) THEN
         IF (arg < 650._DP) THEN
              cv = cv + wg * EXP(arg) * &
                                     ( arg / ( EXP( arg ) - 1.0_DP ) ) ** 2 
         icount=icount+1
         ENDIF
      ELSEIF (nu >1.D-3) THEN
           write(6,*) 'cv using expansion ', temp, nu, arg
             cv = cv + wg * EXP(arg) *  &
                  (arg /( arg + arg**2*0.5_DP + arg**3 * onesixth ))**2
         icount=icount+1
      ENDIF
   ENDDO
ENDDO
cv = cv * kb
CALL mp_sum(cv, intra_image_comm)
CALL mp_sum(icount, intra_image_comm)
WRITE(6,*) icount, ' phonons out of ', nq*3*nat, ' for specific heat at T', temp

RETURN
END SUBROUTINE specific_heat_cv_ph

SUBROUTINE thermal_expansion_ph(ph_freq, ph_grun, temp, betab)
!
!  This routine receives as input a set of phonon frequencies
!  on a uniform q mesh, the gruneisen parameters divided by the cell volume
!  on the same mesh and a temperature and gives as output the volume 
!  thermal expansion multiplied by the bulk modulus at that temperature. 
!  The output is in Ry / (a.u.)^3 / K. To obtain the volume thermal expansion 
!  divide by the bulk modulus in Ry / (a.u.)^3.
!
TYPE(ph_freq_type), INTENT(IN) :: ph_freq, ph_grun
REAL(DP), INTENT(IN)  :: temp
REAL(DP), INTENT(OUT) :: betab

INTEGER :: nq, iq, imode, nat, startq, lastq, icount
REAL(DP) :: nu, temp1, arg, gamman, wg, onesixth

betab=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
onesixth=1.0_DP/6.0_DP
temp1 = 1.0_DP / temp
nq=ph_freq%nq
nat=ph_freq%nat
CALL divide(intra_image_comm, nq, startq, lastq)
icount=0
DO iq=startq,lastq
   wg=ph_freq%wg(iq)
   DO imode=1, 3*nat
      nu=ph_freq%nu(imode,iq)
      gamman=ph_grun%nu(imode,iq) 
      arg= kb1 * nu * temp1
      IF (arg > 1.D-6 ) THEN
         IF (arg < 650._DP) THEN
            betab = betab + wg * gamman &
                            * EXP(arg) * ( arg / ( EXP( arg ) - 1.0_DP )) ** 2 
            icount=icount+1
         ENDIF
      ELSEIF (nu > 1.D-3) THEN
         write(6,*) 'te using expansion ', temp, nu, arg
         betab = betab + wg * gamman * EXP(arg) *  &
                       (arg /( arg + arg**2*0.5_DP + arg**3 * onesixth ))**2
         icount=icount+1
      ENDIF
   ENDDO
ENDDO
betab = betab * kb
CALL mp_sum(betab, intra_image_comm)
CALL mp_sum(icount, intra_image_comm)
WRITE(6,*) icount, ' phonons out of ', nq*3*nat, ' for specific beta at T', temp

RETURN
END SUBROUTINE thermal_expansion_ph

END MODULE ph_freq_module
