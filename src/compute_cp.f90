!
! Copyright (C) 2014-2015 Andrea Dal Corso
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
USE thermo_mod,     ONLY : ngeo, omega_geo, no_ph
USE temperature,    ONLY : temp, ntemp
USE isoentropic,    ONLY : isobaric_heat_capacity, isoentropic_bulk_modulus, &
                           average_gruneisen
USE constants,      ONLY : ry_kbar

IMPLICIT NONE
REAL(DP), INTENT(IN)  :: beta_t(ntemp), vmin_t(ntemp), b0_t(ntemp), &
                         ph_cv(ntemp, ngeo(1))
REAL(DP), INTENT(OUT) :: cv_t(ntemp), cp_t(ntemp), b0_s(ntemp), gamma_t(ntemp)
INTEGER :: itemp, igeo, igeo1, igeo2, jgeo

DO itemp=1,ntemp
!
!  find the two volumes closer to omega_t
!
   igeo1=1
   DO igeo=2,ngeo(1)
      IF (.NOT. no_ph(igeo) ) THEN
         igeo2=igeo
         EXIT
      ENDIF
   ENDDO

   DO igeo = 1, ngeo(1)-1
      IF ( vmin_t(itemp) > omega_geo(igeo) .AND. .NOT. no_ph(igeo)) THEN
         igeo1=igeo
internal:DO jgeo=igeo1+1, ngeo(1)
            IF (.NOT. no_ph(jgeo)) THEN
               igeo2=jgeo
               EXIT internal
            ENDIF
         ENDDO internal
      END IF
   ENDDO
!
!  and interpolate linearly the specific heat
!
   cv_t(itemp) = ph_cv(itemp,igeo1) + ( ph_cv(itemp,igeo2) -              &
                 ph_cv(itemp,igeo1) ) * (vmin_t(itemp)-omega_geo(igeo1)) /  &
                                        (omega_geo(igeo2)-omega_geo(igeo1))
ENDDO
!
CALL isobaric_heat_capacity(vmin_t,b0_t,beta_t,temp,cp_t,ntemp)
!
!  The routine gives only the difference
!
cp_t(1:ntemp) = cp_t(1:ntemp)+cv_t(1:ntemp)

CALL isoentropic_bulk_modulus(vmin_t,b0_t,beta_t,cp_t,temp,b0_s,ntemp)
!
!  The routine gives only the difference
!
b0_s(1:ntemp) = b0_s(1:ntemp)+b0_t(1:ntemp)

CALL average_gruneisen(vmin_t,b0_t,beta_t,cv_t,temp,gamma_t,ntemp)

RETURN
END SUBROUTINE compute_cp
