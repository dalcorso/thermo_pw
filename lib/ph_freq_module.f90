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
USE constants, ONLY : k_boltzmann_ry, ry_to_cmm1
USE io_global, ONLY : stdout
USE mp_images, ONLY : intra_image_comm
USE mp, ONLY : mp_sum
IMPLICIT NONE
SAVE
PRIVATE

REAL(DP), PARAMETER :: kb=k_boltzmann_ry ! Boltzmann constant in Ry/K
REAL(DP), PARAMETER :: kb1=1.0_DP/kb/ry_to_cmm1 ! inverse Boltzmann 
                                                ! constant in cm^{-1}/K

REAL(DP), PARAMETER :: thr_ph=1.D-3   ! a phonon with frequency smaller than
                                      ! this is considered of zero frequency.
!
! NB: we need this threshold to decide which are the three acoustic 
! frequencies. For numerical reason they are never exactly 0 and if 
! introduced in the formulas of this module they give large errors 
! in the average Gruneisen parameters at small temperature.
! All phonon frequencies smaller than thr_ph are removed from the calculation
! A warning is written on output if they are more or less than 3. 
!
REAL(DP), PARAMETER :: thr_taylor=1.D-3   ! When the argument of the functions
                                          ! is smaller than this threshold use
                                          ! a Taylor expansion expression 
                                          ! to avoid the numerical error due
                                          ! to the calculation of 
                                          ! 1.0-(1.0-epsilon) 

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
          init_ph_rap, thermal_expansion_ph, read_ph_freq_data, &
          write_ph_freq_data

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

SUBROUTINE read_ph_freq_data(ph_freq, filename)
!
!  This subroutine reads the ph_freq data from a file. It allocates space
!  for the frequencies and their representation
!
USE io_global, ONLY : ionode, ionode_id
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast
IMPLICIT NONE
TYPE(ph_freq_type), INTENT(INOUT) :: ph_freq
CHARACTER(LEN=256), INTENT(IN) :: filename
INTEGER :: iunit, ios
REAL(DP), ALLOCATABLE :: nu(:), dos(:)
INTEGER :: nat, nq1, nq2, nq3, nq
INTEGER :: iq, imode, ndiv, nqtot

iunit=65
IF (ionode) OPEN(file=TRIM(filename), unit=iunit, status='old', &
     form='formatted', err=100, iostat=ios)
100  CALL mp_bcast(ios, ionode_id, intra_image_comm)
    IF (ios /= 0) CALL errore('read_ph_freq_data','opening file',ios)

IF (ionode) READ(iunit, *) nat, nq1, nq2, nq3, nq

CALL mp_bcast(nat, ionode_id, intra_image_comm)
CALL mp_bcast(nq1, ionode_id, intra_image_comm)
CALL mp_bcast(nq2, ionode_id, intra_image_comm)
CALL mp_bcast(nq3, ionode_id, intra_image_comm)
CALL mp_bcast(nq, ionode_id, intra_image_comm)

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
ALLOCATE(ph_freq%rap(3*nat,nq))

IF (ionode) THEN
   DO iq=1,nq
      READ(iunit, *, END=20, ERR=10, IOSTAT=ios) ph_freq%wg(iq)
      DO imode=1,3*nat
         ! nu(i) = frequencies (cm^{-1})
         READ(iunit, *, END=20, ERR=10, IOSTAT=ios) ph_freq%nu(imode, iq), &
                                                    ph_freq%rap(imode, iq)
      END DO
   END DO
ENDIF
20 CONTINUE
   ios=0
10 CALL mp_bcast( ios, ionode_id, intra_image_comm )
   IF (ios /= 0 ) CALL errore('read_phdos_data', 'problem reading phdos', 1)
   CALL mp_bcast( ph_freq%wg, ionode_id, intra_image_comm )
   CALL mp_bcast( ph_freq%nu, ionode_id, intra_image_comm )
   CALL mp_bcast( ph_freq%rap, ionode_id, intra_image_comm )

   IF (ionode) CLOSE(iunit)

RETURN
END SUBROUTINE read_ph_freq_data

SUBROUTINE write_ph_freq_data(ph_freq, filename)
!
!  This subroutine writes the ph_freq data from a file. It allocates space
!  for the frequencies
!
USE io_global, ONLY : ionode, ionode_id
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast
IMPLICIT NONE
TYPE(ph_freq_type), INTENT(INOUT) :: ph_freq
CHARACTER(LEN=256), INTENT(IN) :: filename
INTEGER :: iunit, ios
REAL(DP), ALLOCATABLE :: nu(:), dos(:)
INTEGER :: nat, nq
INTEGER :: iq, imode

iunit=65
IF (ionode) OPEN(FILE=TRIM(filename), UNIT=iunit, STATUS='unknown', &
     FORM='formatted', ERR=100, IOSTAT=ios)
100  CALL mp_bcast(ios, ionode_id, intra_image_comm)
    IF (ios /= 0) CALL errore('write_ph_freq_data','opening file',ios)

IF (ionode) WRITE(iunit, '(5i10)') ph_freq%nat, ph_freq%nq1, ph_freq%nq2, &
                                   ph_freq%nq3, ph_freq%nq
nq = ph_freq%nq
nat = ph_freq%nat
IF (ionode) THEN
   DO iq=1,nq
      WRITE(iunit, '(E30.15)', ERR=10, IOSTAT=ios) ph_freq%wg(iq)
      DO imode=1,3*nat
         ! nu(i) = frequencies (cm^{-1}) 
         WRITE(iunit,'(E30.15,i30)',ERR=10,IOSTAT=ios) ph_freq%nu(imode, iq), &
                                                       ph_freq%rap(imode, iq)
      END DO
   END DO
END IF
10 CALL mp_bcast( ios, ionode_id, intra_image_comm )
IF (ios /= 0 ) CALL errore('write_phdos_data', 'problem reading phdos', 1)

IF (ionode) CLOSE(iunit)

RETURN
END SUBROUTINE write_ph_freq_data

SUBROUTINE zero_point_energy_ph(ph_freq, ener)
!
!  This subroutine receives as input a set of phonon frequencies and computes 
!  the zero point energy that corresponds to that frequencies. The output 
!  energy is in Ry.
!
TYPE(ph_freq_type), INTENT(IN) :: ph_freq
REAL(DP), INTENT(OUT) :: ener
REAL(DP) :: wg, nu
INTEGER :: nq, iq, imode, nat

nq=ph_freq%nq
nat=ph_freq%nat

ener=0.0_DP
DO iq=1, nq
   wg= ph_freq%wg(iq)
   DO imode=1, 3*nat
      nu=ph_freq%nu(imode, iq)
      IF (nu > thr_ph) ener = ener + 0.5_DP * nu * wg
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

INTEGER :: nq, iq, nat, imode, startq, lastq, counter
REAL(DP) :: nu, arg, wg, temp1, onesixth, one24

free_ener=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
nq=ph_freq%nq
nat=ph_freq%nat
onesixth=1.0_DP / 6.0_DP
one24=1.0_DP /24.0_DP
CALL divide (intra_image_comm, nq, startq, lastq)
counter=0
DO iq=startq,lastq
   wg=ph_freq%wg(iq)
   DO imode=1,3*nat
      nu=ph_freq%nu(imode,iq)
      arg= kb1 * nu * temp1
      IF (arg > thr_taylor) THEN
         free_ener = free_ener +kb*temp*wg*(LOG(1.0_DP-EXP(-arg)))
      ELSEIF (nu > thr_ph) THEN
         free_ener = free_ener+kb*temp*wg*(LOG( arg - arg**2*0.5_DP + &
                                          arg**3 * onesixth - arg**4 * one24 )) 
      ELSE 
         counter=counter+1
      ENDIF
   ENDDO
ENDDO
CALL mp_sum(free_ener, intra_image_comm)
CALL mp_sum(counter, intra_image_comm)
IF (counter > 3) THEN
   WRITE(stdout,'(5x,"WARNING: Too many acoustic modes")')
ELSEIF (counter<3) THEN
   WRITE(stdout,'(5x,"WARNING: Too few acoustic modes")')
ENDIF

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
REAL(DP) :: nu, temp1, arg, wg, onesixth, one24, earg

ener=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
onesixth=1.0_DP / 6.0_DP
one24=1.0_DP /24.0_DP
nq=ph_freq%nq
nat=ph_freq%nat
CALL divide(intra_image_comm, nq, startq, lastq)
DO iq=startq,lastq
   wg=ph_freq%wg(iq)
   DO imode=1, 3*nat
      nu=ph_freq%nu(imode,iq)
      arg= kb1 * nu * temp1
      earg = EXP(-arg)
      IF (arg > thr_taylor) THEN
        ener = ener +  wg*nu*earg / ( 1.0_DP - earg ) 
      ELSEIF ( nu > thr_ph ) THEN
         ener = ener +  wg*nu*earg/(arg-arg**2*0.5_DP+arg**3*onesixth &
                                                    - arg**4 * one24)
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

INTEGER :: nq, iq, imode, nat, startq, lastq
REAL(DP) :: nu, temp1, arg, wg, onesixth, one24, earg

cv=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
onesixth=1.0_DP / 6.0_DP
one24 = 1.0_DP / 24.0_DP
temp1 = 1.0_DP / temp
nq=ph_freq%nq
nat=ph_freq%nat
CALL divide(intra_image_comm, nq, startq, lastq)
DO iq=startq,lastq
   wg=ph_freq%wg(iq)
   DO imode=1,3*nat
      nu=ph_freq%nu(imode,iq)
      arg= kb1 * nu * temp1
      earg=EXP(-arg)
      IF (arg > thr_taylor ) THEN
          cv = cv + wg * earg * ( arg / ( 1.0_DP - earg  ) ) ** 2 
      ELSEIF (nu > thr_ph) THEN
             cv = cv + wg * earg *  &
                  (arg /( arg - arg**2*0.5_DP + arg**3 * onesixth &
                             - arg**4 * one24))**2
      ENDIF
   ENDDO
ENDDO
cv = cv * kb
CALL mp_sum(cv, intra_image_comm)

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

INTEGER :: nq, iq, imode, nat, startq, lastq
REAL(DP) :: nu, temp1, arg, gamman, wg, onesixth, one24, earg

betab=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
onesixth=1.0_DP/6.0_DP
one24=1.0_DP/24.0_DP
temp1 = 1.0_DP / temp
nq=ph_freq%nq
nat=ph_freq%nat
CALL divide(intra_image_comm, nq, startq, lastq)
DO iq=startq,lastq
   wg=ph_freq%wg(iq)
   DO imode=1, 3*nat
      nu=ph_freq%nu(imode,iq)
      gamman=ph_grun%nu(imode,iq) 
      arg= kb1 * nu * temp1
      earg=EXP(-arg)
      IF (arg > thr_taylor ) THEN
          betab = betab + wg * gamman &
                            * earg * ( arg / ( 1.0_DP - earg )) ** 2 
      ELSEIF (nu > thr_ph) THEN
         betab = betab + wg * gamman * earg *  &
                       (arg /( arg - arg**2*0.5_DP + arg**3*onesixth - &
                        arg**4 * one24))**2
      ENDIF
   ENDDO
ENDDO
betab = betab * kb
CALL mp_sum(betab, intra_image_comm)

RETURN
END SUBROUTINE thermal_expansion_ph

END MODULE ph_freq_module
