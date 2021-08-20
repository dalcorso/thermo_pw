!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE write_mur(omega0, b0in, b01, emin, itemp)
!----------------------------------------------------------------------
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
INTEGER, INTENT(IN) :: itemp
CHARACTER(LEN=256)   :: filename
CHARACTER(LEN=6) :: int_to_char
REAL(DP) :: omega, e, p, b0
INTEGER  :: i, iu_mur
INTEGER  :: find_free_unit

IF (my_image_id /= root_image) RETURN

IF (itemp > 0) THEN
   filename="therm_files/"//TRIM(flevdat)//'_mur'//TRIM(int_to_char(itemp))
ELSE
   filename="energy_files/"//TRIM(flevdat)//'_mur'
ENDIF
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
   CALL write_mur_start_line(itemp, iu_mur)
   DO i=1,nvol
      omega= vmin_input + deltav * (i-1)
      e=emin+(omega0 * b0) / b01 *( (1.0_DP / (b01-1.0_DP))*((omega0 / omega)** &
                (b01-1.0_DP)) + omega / omega0 ) - (omega0 *b0 / (b01-1.0_DP))
      p= b0in * ((omega0 / omega)**b01 - 1.0_DP) / b01
      WRITE(iu_mur,'(f18.10,3f20.10)') omega, e, e+p*omega/ry_kbar, &
                                                 p + pressure_kb
   ENDDO 
   CLOSE(UNIT=iu_mur, STATUS='KEEP')
ENDIF

RETURN
END SUBROUTINE write_mur

!-----------------------------------------------------------------------
SUBROUTINE write_mur_start_line(itemp, iu_mur)
!-----------------------------------------------------------------------
USE kinds,            ONLY : DP
USE control_pressure, ONLY : pressure_kb
USE temperature,      ONLY : temp
IMPLICIT NONE
INTEGER :: itemp, iu_mur
REAL(DP) :: temper

IF (itemp==0) THEN
   IF (pressure_kb /= 0.0_DP) THEN
      WRITE(iu_mur,'( "#",2x,"omega (a.u.)**3",6x,"enthalpy (Ry)",6x,&
       & "enthalpy(p) (Ry)",5x,"pressure (kbar)")') 
   ELSE
      WRITE(iu_mur,'( "#",2x,"omega (a.u.)**3",7x,"energy (Ry)",8x,&
       & "enthalpy(p) (Ry)",4x,"pressure (kbar)")') 
   ENDIF
ELSE
   temper=temp(itemp)
   IF (pressure_kb /= 0.0_DP) THEN
      WRITE(iu_mur,'( "#","omega(a.u.)**3",7x,"Gibbs energy (Ry)",  &
       & " Gibbs energy(p,T)(Ry)"," pressure (kbar), T=",f8.2," K")') temper
   ELSE
      WRITE(iu_mur,'( "#","omega(a.u.)**3",3x,"Helm.&
       & free-energy(Ry)"," Gibbs energy(p,T) (Ry)",  &
       & " pressure (kbar), T=",f8.2," K")') temper
   ENDIF
ENDIF

RETURN
END SUBROUTINE write_mur_start_line
