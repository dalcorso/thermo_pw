! Copyright (C) 2022 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE model_free_energy
!
!   This module contains routines to add to the free 
!   energy contributions non computed by ab-initio, but derived
!   from model fits to experiment. The parameters of the free energy
!   must be given in input by the user. The module contains a few
!   of the commonly used functions.
!
!
  USE kinds, ONLY : DP
  USE constants, ONLY : k_boltzmann_ry, rytoev, bohr_radius_angs
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC  emp_quadratic_t, emp_anhar_am
          

CONTAINS
!-------------------------------------------------------------------
SUBROUTINE emp_quadratic_t(t, vol, alpha1, alpha2, free_ener, ener, &
                                                           entr, cv) 
!-------------------------------------------------------------------
!
!  In input t in (K), vol in (a.u.)*3, alpha0 in (eV/K^2), in 
!  alpha1 (eV/A^3/K^2)
!  The free_energy is in Ry.
!  This term is described for instance in Phys. Rev. B 81, 014301 (2010).
!  Note that the parameters in this paper refer to the conventional cell,
!  so the alpha1 to pass to this routine is 1/2 the one reported in the paper
!  alpha2 is the same.
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: t, vol, alpha1, alpha2
REAL(DP), INTENT(OUT) :: free_ener, ener, entr, cv

free_ener= (alpha1 + alpha2 *bohr_radius_angs**3 * vol ) * t**2 / rytoev
ener= -free_ener
entr= -2.0_DP * (alpha1 + alpha2 *bohr_radius_angs**3 * vol) * t / rytoev
cv= entr

RETURN
END SUBROUTINE emp_quadratic_t

!-------------------------------------------------------------------
SUBROUTINE emp_anhar_am(t, vol, alpha, emme, v0, nat, free_ener, &
                                                      ener, entr, cv) 
!-------------------------------------------------------------------
!
!  In input t in (K), vol and v0 in the same units, emme adimensional,
!  alpha in 1/K.
!  On output the free_energy is in Ry.
!  This term is described for instance in J. Appl. Phys.113, 093507 (2013).
!
IMPLICIT NONE
INTEGER :: nat
REAL(DP), INTENT(IN) :: t, vol, v0, alpha, emme
REAL(DP), INTENT(OUT) :: free_ener, ener, entr, cv

REAL(DP) :: x

x=vol/v0

free_ener= -1.5_DP * nat * k_boltzmann_ry * alpha * x**emme * t**2
ener=-free_ener
entr= 3.0_DP * nat * k_boltzmann_ry * alpha * x**emme * t
cv=entr

RETURN
END SUBROUTINE emp_anhar_am

END MODULE model_free_energy
