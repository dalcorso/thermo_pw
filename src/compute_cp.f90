!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE compute_cp()
!
!  This subroutine receives the constant volume heat capacity at
!  several volumes, the equilibrium volume as a function of temperature,
!  the bulk modulus as a function of the temperature and the 
!  compressibility as a function of temperature and computes the
!  isobaric heat capacity as a function of temperature and the
!  isoentropic bulk modulus as a function of temperature
!  cv_t and cp_t in Ry / K / cell
!  b0_t and b0_t in kbar
!
USE anharmonic,     ONLY : cp_t, cv_t, beta_t, b0_s
USE thermodynamics, ONLY : temp, ph_cv, ntemp, ngeo, omegav
USE thermo_mod,     ONLY : vmin_t, b0_t
USE constants,      ONLY : ry_kbar

IMPLICIT NONE
INTEGER :: itemp, igeo, igeo1, igeo2


IF ( .NOT. ALLOCATED (cp_t) ) ALLOCATE( cp_t (ntemp) )
IF ( .NOT. ALLOCATED (cv_t) ) ALLOCATE( cv_t (ntemp) )
IF ( .NOT. ALLOCATED (b0_s) ) ALLOCATE( b0_s (ntemp) )

DO itemp=1,ntemp
!
!  find the two volumes closer to omega_t
!
   igeo1=1
   igeo2=2
   DO igeo = 1, ngeo-1
      IF ( vmin_t(itemp) > omegav(igeo) ) THEN
         igeo1=igeo
         igeo2=igeo1+1
      END IF
   ENDDO
!
!  and interpolate linearly the specific heat
!
   cv_t(itemp) = ph_cv(itemp,igeo1) + ( ph_cv(itemp,igeo2) -              &
                 ph_cv(itemp,igeo1) ) * (vmin_t(itemp)-omegav(igeo1)) /  &
                                        (omegav(igeo2)-omegav(igeo1))

   cp_t(itemp) = cv_t(itemp) + temp(itemp,1) * vmin_t(itemp) * &
                           beta_t(itemp)**2 * b0_t(itemp) / ry_kbar

   b0_s(itemp) = b0_t(itemp) + temp(itemp,1) * vmin_t(itemp) * &
                           beta_t(itemp)**2 * b0_t(itemp)**2 / cp_t(itemp) &
                           / ry_kbar
END DO

RETURN
END SUBROUTINE compute_cp
