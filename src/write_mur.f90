!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_mur(omega0, b0in, b01, emin)
!
!  This routine writes on file the Murnaghan energy versus volume
!  curve, together with the pressure versus volume curve. 
!  It receives the parameters of the Murnaghan equation:
!
! in input emin in Ry, omega0 in (a.u.)**3, b0in in kbar, b01 adimensional
! in output omega in (a.u.)**3, p in kbar, e in Ry
!
USE kinds,            ONLY : DP
USE data_files,       ONLY : flevdat
USE thermo_mod,       ONLY : ngeo, omega_geo
USE control_mur,      ONLY : nvol, vmin_input, vmax_input, deltav
USE control_pressure, ONLY : pressure_kb
USE mp_images,        ONLY : root_image, my_image_id
USE constants,        ONLY : ry_kbar
USE io_global,        ONLY : ionode

IMPLICIT NONE

REAL(DP), INTENT(IN) :: emin, omega0, b0in, b01
CHARACTER(LEN=256)   :: filename
REAL(DP) :: omega, e, p, b0
INTEGER  :: i, iu_mur
INTEGER  :: find_free_unit

IF (my_image_id /= root_image) RETURN

filename="energy_files/"//TRIM(flevdat)//'_mur'
CALL add_pressure(filename)

b0 = b0in / ry_kbar
IF (vmin_input == 0.0_DP) vmin_input=omega_geo(1) * 0.98_DP
IF (vmax_input == 0.0_DP) vmax_input=omega_geo(ngeo(1)) * 1.02_DP
IF (nvol > 1) THEN
   deltav = (vmax_input - vmin_input)/(nvol-1)
ELSE
   IF (deltav > 0.0_DP) THEN
      nvol = NINT ( ( vmax_input - vmin_input ) / deltav + 1.51d0 )
   ELSE
      nvol = 51
      deltav = (vmax_input - vmin_input)/(nvol-1)
   ENDIF
END IF

IF (ionode) THEN
   iu_mur=find_free_unit()
   OPEN(UNIT=iu_mur, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   IF (pressure_kb /= 0.0_DP) THEN
      WRITE(iu_mur,'( "# omega (a.u.)**3      enthalpy (Ry)   pressure (kbar)" )')
   ELSE
      WRITE(iu_mur,'( "# omega (a.u.)**3       energy (Ry)      pressure (kbar)" )')
   END IF
   DO i=1,nvol
      omega= vmin_input + deltav * (i-1)
      e=emin+(omega0 * b0) / b01 *( (1.0_DP / (b01-1.0_DP))*((omega0 / omega)** &
                (b01-1.0_DP)) + omega / omega0 ) - (omega0 *b0 / (b01-1.0_DP))
      p= b0in * ((omega0 / omega)**b01 - 1.0_DP) / b01
      WRITE(iu_mur,'(3f20.10)') omega, e, p + pressure_kb
   ENDDO 
   CLOSE(UNIT=iu_mur, STATUS='KEEP')

END IF

RETURN
END SUBROUTINE write_mur
