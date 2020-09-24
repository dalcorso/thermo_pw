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
INTEGER  :: work_choice
CHARACTER(LEN=9) :: code='Optical'

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

WRITE(stdout,'(/,5x,"Choose what to do")')
WRITE(stdout,'(5x,"1) Compute the complex refractive index n+ik and &
                                                   &the reflectivity R")')
WRITE(stdout,'(5x,"2) Compute the complex dielectric constant e1+ie2 and R")')
WRITE(stdout,'(5x,"3) Compute the absorption coefficient in cm^(-1) &
                                                        &given n and k")')
WRITE(stdout,'(5x,"4) Compute the absorption coefficient in cm^(-1) &
                                                        &given e1 and e2")')
WRITE(stdout,'(5x,"5) Convert energy -> frequency (Hz), wavelength (nm)")')

READ(5,*) work_choice

IF (work_choice == 1) THEN
   WRITE(stdout,'(/,5x,"Enter e1 and e2")')
   READ(5,*) epsilon1, epsilon2
   CALL compute_refractive_index(epsilon1, epsilon2, enne, kappa)
   CALL compute_reflectivity(enne,kappa,ref)
   CALL write_output(epsilon1, epsilon2, enne, kappa, ref)
ELSEIF (work_choice == 2) THEN
   WRITE(stdout,'(/,5x,"Enter n and k")')
   READ(5,*) enne, kappa
   CALL compute_complex_epsilon(epsilon1, epsilon2, enne, kappa)
   CALL compute_refractive_index(epsilon1, epsilon2, enne, kappa)
   CALL compute_reflectivity(enne,kappa,ref)
   CALL write_output(epsilon1, epsilon2, enne, kappa, ref)
ELSEIF (work_choice == 3) THEN
   WRITE(stdout,'(/,5x,"Enter n and k")')
   READ(5,*) enne, kappa
   CALL read_energy(omega)
   CALL compute_complex_epsilon(epsilon1, epsilon2, enne, kappa)
   CALL compute_refractive_index(epsilon1, epsilon2, enne, kappa)
   CALL compute_reflectivity(enne,kappa,ref)
   CALL compute_alpha(omega, kappa, alpha)
   CALL write_output(epsilon1, epsilon2, enne, kappa, ref)
   WRITE(stdout,'(/,5x,"alpha     =",e15.6, " cm^(-1)")') alpha
ELSEIF (work_choice == 4) THEN
   WRITE(stdout,'(/,5x,"Enter e1 and e2")')
   READ(5,*) epsilon1, epsilon2
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
USE kinds,          ONLY : DP
USE constants,      ONLY : rytoev
USE optical_module, ONLY : compute_energy 
USE io_global,      ONLY : stdout

IMPLICIT NONE
REAL(DP), INTENT(OUT) :: omega

INTEGER  :: units
REAL(DP) :: lambda_in

WRITE(stdout,'(/,5x,"Enter units of photon energy")')
WRITE(stdout,'(5x,"1) eV")')
WRITE(stdout,'(5x,"2) Ry")')
WRITE(stdout,'(5x,"3) give wavelength in nm")')
READ(5,*) units
IF (units==1.OR.units==2) THEN
   WRITE(stdout,'(/,5x,"Enter photon energy")')
   READ(5,*) omega
   IF (units==1) omega = omega / rytoev
ELSEIF (units==3) THEN
   WRITE(stdout,'(/,5x,"Enter photon wavelength in nm")')
   READ(5,*) lambda_in
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
