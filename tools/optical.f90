!
! Copyright (C) 2020 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM optical
!-----------------------------------------------------------------------
!
!  This program computes several quantities related to the optical 
!  properties of solids.
!  It can read the complex dielectric constant at a given frequency
!  and computes the complex index of refraction, the absorption coefficient,
!  and the reflectivity for normal incidence in a cubic solid.
!  It can also read the complex index of refraction and compute the complex
!  dielectric constant.
!  Finally it can convert the photon energy in eV into the frequency of
!  the electromagnetic wave in Hz or its wavelength in vacuum in nm.
!
!
USE kinds,          ONLY : DP
USE constants,      ONLY : rytoev
USE optical_module, ONLY : compute_refractive_index, compute_alpha, &
                           compute_energy, compute_frequency, compute_lambda, &
                           compute_complex_epsilon, compute_reflectivity
USE mp_global,      ONLY : mp_startup, mp_global_end
USE environment,    ONLY : environment_start, environment_end
USE io_global,      ONLY : stdout
IMPLICIT NONE

REAL(DP) :: epsilon1, epsilon2, enne, kappa, ref
REAL(DP) :: omega, omega_hz, lambda, lambda_in, alpha
INTEGER  :: work_choice, stdin
CHARACTER(LEN=9) :: code='Optical'

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

stdin=5

WRITE(stdout,'(/,5x,"Choose what to do")')
WRITE(stdout,'(5x,"1) Compute the complex refractive index n+ik and &
                                                   &the reflectivity R")')
WRITE(stdout,'(5x,"2) Compute the complex dielectric constant e1+ie2 and R")')
WRITE(stdout,'(5x,"3) Compute the absorption coefficient in cm^(-1) &
                                                        &given n and k")')
WRITE(stdout,'(5x,"4) Compute the absorption coefficient in cm^(-1) &
                                                        &given e1 and e2")')
WRITE(stdout,'(5x,"5) Enter Energy (eV or Ry), frequency (Hz) or &
                           &wavelength (nm), get the others")')

READ(stdin,*) work_choice

IF (work_choice == 1) THEN
   WRITE(stdout,'(/,5x,"Enter e1 and e2")')
   READ(stdin,*) epsilon1, epsilon2
   CALL compute_refractive_index(epsilon1, epsilon2, enne, kappa)
   CALL compute_reflectivity(enne,kappa,ref)
   CALL write_output(epsilon1, epsilon2, enne, kappa, ref)
ELSEIF (work_choice == 2) THEN
   WRITE(stdout,'(/,5x,"Enter n and k")')
   READ(stdin,*) enne, kappa
   CALL compute_complex_epsilon(epsilon1, epsilon2, enne, kappa)
   CALL compute_refractive_index(epsilon1, epsilon2, enne, kappa)
   CALL compute_reflectivity(enne,kappa,ref)
   CALL write_output(epsilon1, epsilon2, enne, kappa, ref)
ELSEIF (work_choice == 3) THEN
   WRITE(stdout,'(/,5x,"Enter n and k")')
   READ(stdin,*) enne, kappa
   CALL read_energy(omega)
   CALL compute_complex_epsilon(epsilon1, epsilon2, enne, kappa)
   CALL compute_refractive_index(epsilon1, epsilon2, enne, kappa)
   CALL compute_reflectivity(enne,kappa,ref)
   CALL compute_alpha(omega, kappa, alpha)
   CALL write_output(epsilon1, epsilon2, enne, kappa, ref)
   WRITE(stdout,'(/,5x,"alpha     =",e15.6, " cm^(-1)")') alpha
ELSEIF (work_choice == 4) THEN
   WRITE(stdout,'(/,5x,"Enter e1 and e2")')
   READ(stdin,*) epsilon1, epsilon2
   CALL read_energy(omega)
   CALL compute_refractive_index(epsilon1, epsilon2, enne, kappa)
   CALL compute_reflectivity(enne,kappa,ref)
   CALL compute_alpha(omega, kappa, alpha)
   CALL write_output(epsilon1, epsilon2, enne, kappa, ref)
   WRITE(stdout,'(/,5x,"alpha     =",e15.6, " cm^(-1)")') alpha
ELSEIF (work_choice == 5) THEN
   CALL read_energy(omega)
   CALL compute_frequency(omega, omega_hz)
   CALL compute_lambda(omega, lambda)
   WRITE(stdout,'(/,5x,"energy    =",e15.6," Ry", e15.6," eV")') omega, &
                                                omega * rytoev
   WRITE(stdout,'(5x,  "frequency =",e15.6," Hz")') omega_hz
   WRITE(stdout,'(5x,  "wavelength=",e15.6," nm")') lambda
ELSE
   WRITE(stdout,'(/,5x,"Work choice not recognized. Stoppingi ...")')
   STOP 1
ENDIF

CALL environment_end( code )
CALL mp_global_end ()

END PROGRAM optical
!
!------------------------------------------------------------
SUBROUTINE read_energy(omega) 
!------------------------------------------------------------
!
!  This routine reads the energy of a photon, the frequency of
!  the electromagnetic wave or its wavelength. In all cases it
!  gives as output the energy in eV.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : rytoev
USE optical_module, ONLY : compute_energy, compute_energy_hz
USE io_global,      ONLY : stdout

IMPLICIT NONE
REAL(DP), INTENT(OUT) :: omega

INTEGER  :: units, stdin
REAL(DP) :: lambda_in, freq_hz

stdin=5
WRITE(stdout,'(/,5x,"Enter units of photon energy")')
WRITE(stdout,'(5x,"1) eV")')
WRITE(stdout,'(5x,"2) Ry")')
WRITE(stdout,'(5x,"3) Give the frequency in Hz")')
WRITE(stdout,'(5x,"4) Give the wavelength in nm")')
READ(stdin,*) units
IF (units==1.OR.units==2) THEN
   WRITE(stdout,'(/,5x,"Enter photon energy")')
   READ(stdin,*) omega
   IF (units==1) omega = omega / rytoev
ELSEIF (units==3) THEN
   WRITE(stdout,'(/,5x,"Enter the frequency in Hz")')
   READ(stdin,*) freq_hz
   CALL compute_energy_hz(freq_hz, omega)
ELSEIF (units==4) THEN
   WRITE(stdout,'(/,5x,"Enter photon wavelength in nm")')
   READ(stdin,*) lambda_in
   CALL compute_energy(omega, lambda_in)
ELSE
   WRITE(stdout,'(5x,"Energy units not recognized. Stopping ...")')
   STOP 1
ENDIF

RETURN
END SUBROUTINE read_energy
!
!------------------------------------------------------------------
SUBROUTINE write_output(epsilon1, epsilon2, enne, kappa, ref)
!------------------------------------------------------------------
!
USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE
REAL(DP), INTENT(IN) :: epsilon1, epsilon2, enne, kappa, ref

WRITE(stdout,'(/,5x,"e1 + i e2=",e15.6," + i",e15.6)') epsilon1, epsilon2
WRITE(stdout,'(5x,  "n  + i k =",e15.6," + i",e15.6)') enne, kappa
WRITE(stdout,'(5x,  "R        =",e15.6)') ref
RETURN
END SUBROUTINE write_output
