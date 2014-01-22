!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE mur(omega0, b0in, b01, emin)
!
! in input emin in Ry, omega0 in (a.u.)**3, b0in in kbar
! in output omega in (a.u.)**3, p in kbar, e in Ry
!
USE kinds,          ONLY : DP
USE control_thermo, ONLY : flevdat
USE thermo_mod,     ONLY : vmin_input, vmax_input, deltav, nvol, energy_geo
USE mp_images,      ONLY : root_image, my_image_id
USE thermodynamics, ONLY : omegav, ngeo
USE constants,      ONLY : ry_kbar
USE io_global,      ONLY : ionode

IMPLICIT NONE

REAL(DP), INTENT(IN) :: emin, omega0, b0in, b01
CHARACTER(LEN=256) :: filename, filename1
REAL(DP) :: omega, e, p, b0
INTEGER :: i, iu_mur

IF (my_image_id /= root_image) RETURN

filename=TRIM(flevdat)//'_mur'
filename1=TRIM(flevdat)//'_mur1'
b0 = b0in / ry_kbar
IF (vmin_input == 0.0_DP) vmin_input=omegav(1) * 0.98_DP
IF (vmax_input == 0.0_DP) vmax_input=omegav(ngeo) * 1.02_DP
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
   iu_mur=2
   OPEN(UNIT=iu_mur, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_mur,'( "# omega (a.u.)**3       energy (Ry)      pressure (kbar)" )')
   DO i=1,nvol
      omega= vmin_input + deltav * (i-1)
      e=emin+(omega0 * b0) / b01 *( (1.0_DP / (b01-1.0_DP))*((omega0 / omega)** &
                (b01-1.0_DP)) + omega / omega0 ) - (omega0 *b0 / (b01-1.0_DP))
      p= b0in * ((omega0 / omega)**b01 - 1.0_DP) / b01
      WRITE(iu_mur,'(3f20.10)') omega, e, p
   ENDDO 
   CLOSE(UNIT=iu_mur, STATUS='KEEP')
   OPEN(UNIT=iu_mur, FILE=TRIM(filename1), STATUS='UNKNOWN', FORM='FORMATTED')
   DO i=1,ngeo
      WRITE(iu_mur,'(2f20.10)') omegav(i), energy_geo(i)
   ENDDO
   CLOSE(UNIT=iu_mur, STATUS='KEEP')
END IF

RETURN
END SUBROUTINE mur
