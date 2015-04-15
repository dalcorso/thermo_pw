!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE compute_cp(beta_t, vmin_t, b0_t, ph_cv, cv_t, cp_t, b0_s, gamma_t)
!
!  This subroutine receives the constant volume heat capacity at
!  several volumes, the equilibrium volume as a function of temperature,
!  the bulk modulus as a function of the temperature and the 
!  thermal expansion as a function of temperature and computes the
!  isobaric heat capacity as a function of temperature, the
!  isoentropic bulk modulus as a function of temperature, 
!  the average gruneisen parameter as a function of temperature
!  cv_t and cp_t in Ry / K / cell
!  b0_t and b0_s in kbar
!  gamma_t adimensional
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : tot_ngeo, omega_geo
USE temperature,    ONLY : temp, ntemp
USE constants,      ONLY : ry_kbar

IMPLICIT NONE
REAL(DP), INTENT(IN)  :: beta_t(ntemp), vmin_t(ntemp), b0_t(ntemp), &
                         ph_cv(ntemp, tot_ngeo)
REAL(DP), INTENT(OUT) :: cv_t(ntemp), cp_t(ntemp), b0_s(ntemp), gamma_t(ntemp)
INTEGER :: itemp, igeo, igeo1, igeo2

DO itemp=1,ntemp
!
!  find the two volumes closer to omega_t
!
   igeo1=1
   igeo2=2
   DO igeo = 1, tot_ngeo-1
      IF ( vmin_t(itemp) > omega_geo(igeo) ) THEN
         igeo1=igeo
         igeo2=igeo1+1
      END IF
   ENDDO
!
!  and interpolate linearly the specific heat
!
   cv_t(itemp) = ph_cv(itemp,igeo1) + ( ph_cv(itemp,igeo2) -              &
                 ph_cv(itemp,igeo1) ) * (vmin_t(itemp)-omega_geo(igeo1)) /  &
                                        (omega_geo(igeo2)-omega_geo(igeo1))
ENDDO

DO itemp=2,ntemp-1

   cp_t(itemp) = cv_t(itemp) + temp(itemp) * vmin_t(itemp) * &
                           beta_t(itemp)**2 * b0_t(itemp) / ry_kbar

   b0_s(itemp) = b0_t(itemp) + temp(itemp) * vmin_t(itemp) * &
                           beta_t(itemp)**2 * b0_t(itemp)**2 / cp_t(itemp) &
                           / ry_kbar

   gamma_t(itemp) =  beta_t(itemp) * b0_t(itemp) * vmin_t(itemp) / ry_kbar / &
                     cv_t(itemp) 

END DO

RETURN
END SUBROUTINE compute_cp
