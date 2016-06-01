!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE interpolate_cv(vmin_t, ph_cv, cv_t)
!
!  This subroutine receives the constant volume heat capacity at
!  several volumes, the equilibrium volume as a function of temperature,
!  and interpolate the isobaric heat capacity at the volume that corresponds
!  to a given temperature.
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : ngeo, omega_geo, no_ph, tot_ngeo
USE temperature,    ONLY : ntemp

IMPLICIT NONE
REAL(DP), INTENT(IN)  :: vmin_t(ntemp), ph_cv(ntemp, tot_ngeo)
REAL(DP), INTENT(OUT) :: cv_t(ntemp)
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
RETURN
END SUBROUTINE interpolate_cv
