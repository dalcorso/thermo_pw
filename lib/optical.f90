!
! Copyright (C) 2020 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE optical_module
!
!  This module provides subroutines to deal with optical properties.
!  It offers the following functions:
!
!  compute_refractive_index that receives the real and imaginary part
!                           of a complex dielectric constant and computes
!                           the complex index of refraction.
!  compute_complex_epsilon  that receives the real and imaginary part of the
!                           complex refractive index and gives the 
!                           real and imaginary part of the dielectric constant
!  compute_lambda     that receives the energy of the photon in Ry and
!                     gives the wavelength in nanometers
!  compute_omega      that receives the wavelength of the photon in nanometers
!                     and gives the energy in Ry
!  compute_frequency  that receives the energy of the photon in eV and
!                     computes its frequency in Hz
!  compute_alpha      receives the real part of the frequency in Ry, the 
!                     imaginary part of the refractive index and gives the
!                     absorption coefficient in cm^-1
!  compute_reflectivity receives the real and imaginary part of the complex 
!                     index of refraction and gives the reflectivity for
!                     normal incidence in a cubic solid.
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC compute_refractive_index, compute_complex_epsilon, compute_lambda, &
         compute_frequency, compute_energy, compute_reflectivity, compute_alpha

CONTAINS
!
!----------------------------------------------------------------------------
SUBROUTINE compute_refractive_index(epsilon1, epsilon2, enne, kappa)
!----------------------------------------------------------------------------
!
! This routine reveives as input the real and imaginary part of the
! complex dielectric constant and gives as output the real and 
! imaginary part of the complex index of reflection. 
!
USE kinds, ONLY : DP

IMPLICIT NONE
REAL(DP), INTENT(IN)  :: epsilon1, epsilon2
REAL(DP), INTENT(OUT) :: enne, kappa

REAL(DP) :: aux

aux=SQRT(epsilon1**2 + epsilon2**2)
enne = SQRT(0.5_DP*(epsilon1+aux))
kappa = SQRT(0.5_DP*(-epsilon1+aux))

RETURN
END SUBROUTINE compute_refractive_index
!
!----------------------------------------------------------------------------
SUBROUTINE compute_complex_epsilon(epsilon1, epsilon2, enne, kappa)
!----------------------------------------------------------------------------
!
! This routine reveives as input the real and imaginary part of the
! complex refractive index and gives as output the real and 
! imaginary part of the complex dielectric constant. 
!
USE kinds, ONLY : DP

IMPLICIT NONE
REAL(DP), INTENT(IN) :: enne, kappa
REAL(DP), INTENT(OUT)  :: epsilon1, epsilon2

epsilon1=enne**2 - kappa**2
epsilon2=2.0_DP * enne * kappa

RETURN
END SUBROUTINE compute_complex_epsilon
!
!--------------------------------------------------------------------
SUBROUTINE compute_lambda(freq_in, lambda)
!--------------------------------------------------------------------
!
!  Receives the energy of the photon in Rydberg and gives the wavelength 
!  of light in vacuum in nm.
!  freq_in multiplied rydberg_si is the h freq_in in Joule. Dividing by 
!  h_planck_si we obtain the freq_in in Hz. Given the speed of light
!  in m/s divided by the freq_in in Hz gives the wavelenght in m.
!  Multiplication by 1.D9 gives it in nm.
!
USE constants, ONLY : c_si, h_planck_si, rydberg_si
IMPLICIT NONE

REAL(DP), INTENT(IN)  :: freq_in
REAL(DP), INTENT(OUT) :: lambda

lambda = c_si * h_planck_si * 1.D9 / (freq_in * rydberg_si)

RETURN
END SUBROUTINE compute_lambda
!
!--------------------------------------------------------------------
SUBROUTINE compute_energy(omega, lambda)
!--------------------------------------------------------------------
!
!  Receives the wavelength of the light in vacuum in nm and gives
!  its energy in Rydberg.
!  The formula is the inverse of that used in compute_lambda
!
USE constants, ONLY : c_si, h_planck_si, rydberg_si
IMPLICIT NONE

REAL(DP), INTENT(IN) :: lambda
REAL(DP), INTENT(OUT)  :: omega

omega = c_si * h_planck_si * 1.D9 / (lambda * rydberg_si)

RETURN
END SUBROUTINE compute_energy

!------------------------------------------------------------
SUBROUTINE compute_reflectivity(enne, kappa, ref)
!------------------------------------------------------------
!
!   This routine computes the reflectivity for normal incidence of a 
!   cubic solid, given the complex index of refraction.
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP), INTENT(IN)  :: enne, kappa
REAL(DP), INTENT(OUT) :: ref

ref = (enne-1.0_DP)**2+kappa**2
ref = ref / ((enne+1.0_DP)**2+kappa**2)

RETURN
END SUBROUTINE compute_reflectivity

!--------------------------------------------------------------------
SUBROUTINE compute_alpha(freq_in, kappa, alpha)
!--------------------------------------------------------------------
!
!  Receives the frequency in Rydberg, kappa and gives the absorption
!  coefficient in cm^-1.
!  freq_in multiplied rydberg_si is the h freq_in in Joule. Dividing by 
!  h_planck_si we obtain the freq_in in Hz.
!  Multiplication by 2 pi gives the angular frequency w in Hz. The expression
!  for the absorption coefficient is:
!  alpha = 2 w k / c    
!  where c is the speed of light in m/s and k the imaginary part of the 
!  complex refractive index.
!
USE constants, ONLY : pi, c_si, h_planck_si, rydberg_si
IMPLICIT NONE

REAL(DP), INTENT(IN)  :: freq_in, kappa
REAL(DP), INTENT(OUT) :: alpha

alpha = 4.0_DP * pi * freq_in * rydberg_si * kappa / h_planck_si / c_si / 1.D2

RETURN
END SUBROUTINE compute_alpha

!--------------------------------------------------------------------
SUBROUTINE compute_frequency(freq_in, freq_out)
!--------------------------------------------------------------------
!
!  Receives the frequency in Rydberg and, as output, gives the frequency
!  in Hz.
!  freq_in multiplied rydberg_si is the h freq_in in Joule. Dividing by 
!  h_planck_si we obtain the freq_out in Hz. 
!
USE constants, ONLY : h_planck_si, rydberg_si
IMPLICIT NONE

REAL(DP), INTENT(IN)  :: freq_in
REAL(DP), INTENT(OUT) :: freq_out

freq_out = freq_in * rydberg_si / h_planck_si

RETURN
END SUBROUTINE compute_frequency

END MODULE optical_module
