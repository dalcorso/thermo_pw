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
USE thermo_mod,     ONLY : ngeo, omega_geo
USE anharmonic,     ONLY : cp_t, cv_t, beta_t, gamma_t, b0_s, vmin_t, b0_t
USE temperature,    ONLY : temp, ntemp
USE thermodynamics, ONLY : ph_cv
USE constants,      ONLY : ry_kbar

IMPLICIT NONE
INTEGER :: itemp, igeo, igeo1, igeo2

DO itemp=1,ntemp
!
!  find the two volumes closer to omega_t
!
   igeo1=1
   igeo2=2
   DO igeo = 1, ngeo-1
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

SUBROUTINE compute_ph_freq_cp()
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
USE thermo_mod,     ONLY : ngeo, omega_geo
USE ph_freq_anharmonic,     ONLY : cpf_t, cvf_t, betaf_t, gammaf_t, b0f_s, vminf_t, b0f_t
USE temperature,    ONLY : temp, ntemp
USE ph_freq_thermodynamics, ONLY : phf_cv
USE constants,      ONLY : ry_kbar

IMPLICIT NONE
INTEGER :: itemp, igeo, igeo1, igeo2

DO itemp=1,ntemp
!
!  find the two volumes closer to omega_t
!
   igeo1=1
   igeo2=2
   DO igeo = 1, ngeo-1
      IF ( vminf_t(itemp) > omega_geo(igeo) ) THEN
         igeo1=igeo
         igeo2=igeo1+1
      END IF
   ENDDO
!
!  and interpolate linearly the specific heat
!
   cvf_t(itemp) = phf_cv(itemp,igeo1) + ( phf_cv(itemp,igeo2) -              &
                  phf_cv(itemp,igeo1) ) * (vminf_t(itemp)-omega_geo(igeo1)) /  &
                                        (omega_geo(igeo2)-omega_geo(igeo1))
ENDDO
!
!  then calculate the other quantities
!
DO itemp=2,ntemp-1

   cpf_t(itemp) = cvf_t(itemp) + temp(itemp) * vminf_t(itemp) * &
                           betaf_t(itemp)**2 * b0f_t(itemp) / ry_kbar

   b0f_s(itemp) = b0f_t(itemp) + temp(itemp) * vminf_t(itemp) * &
                           betaf_t(itemp)**2 * b0f_t(itemp)**2 / cpf_t(itemp) &
                           / ry_kbar

   gammaf_t(itemp) =  betaf_t(itemp) * b0f_t(itemp) * vminf_t(itemp) / ry_kbar /                      cvf_t(itemp)

END DO

RETURN
END SUBROUTINE compute_ph_freq_cp
