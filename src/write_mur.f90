!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE write_mur(omega0, b0, b01, b02, emin)
!----------------------------------------------------------------------
!
!  This routine writes on file the energy versus volume
!  curve, together with the pressure versus volume curve. 
!  It receives the parameters of the Murnaghan (or of a Birch
!  Murnaghan) equation:
!
! in input emin in Ry, omega0 in (a.u.)**3, b0 in kbar, b01 adimensional
! b02 in 1/kbar
! in output omega in (a.u.)**3, p in kbar, e in Ry
!
USE kinds,            ONLY : DP
USE data_files,       ONLY : flevdat
USE thermo_mod,       ONLY : ngeo, omega_geo
USE control_mur,      ONLY : p0
USE control_vol,      ONLY : nvol, vmin_input, vmax_input, deltav
USE control_ev,       ONLY : ieos
USE eos,              ONLY : eos_energy, eos_press
USE control_pressure, ONLY : pressure_kb, pressure
USE mp_images,        ONLY : root_image, my_image_id
USE constants,        ONLY : ry_kbar
USE io_global,        ONLY : ionode

IMPLICIT NONE

REAL(DP), INTENT(IN) :: emin, omega0, b0, b01, b02
CHARACTER(LEN=256)   :: filename
CHARACTER(LEN=6) :: int_to_char
REAL(DP) :: omega, e, p
INTEGER  :: i, iu_mur
INTEGER  :: find_free_unit

IF (my_image_id /= root_image) RETURN

filename="energy_files/"//TRIM(flevdat)//'_mur'
CALL add_pressure(filename)

IF (ionode) THEN
   iu_mur=find_free_unit()
   OPEN(UNIT=iu_mur, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   CALL write_mur_start_line(0, iu_mur)
   DO i=1,nvol
      omega= vmin_input + deltav * (i-1)
      CALL eos_energy(ieos, omega, e, omega0, b0/ry_kbar, b01, b02*ry_kbar)
      e=e+emin
      CALL eos_press(ieos, omega, p, omega0, b0/ry_kbar, b01, b02*ry_kbar)
      WRITE(iu_mur,'(f18.10,3f20.10)') omega, e, e+p*omega, &
                                            p*ry_kbar + pressure_kb
      p0(i)= (p + pressure) * ry_kbar
   ENDDO 
   CLOSE(UNIT=iu_mur, STATUS='KEEP')
ENDIF

RETURN
END SUBROUTINE write_mur
!
!-----------------------------------------------------------------------
SUBROUTINE write_mur_p()
!-----------------------------------------------------------------------
!
!  This routine writes on file the parameters of the equation of state
!  as a function of pressure
!
USE data_files,       ONLY : flevdat
USE temperature,      ONLY : ntemp_plot
USE control_pressure, ONLY : press, npress, npress_plot
USE control_mur_p,    ONLY : vmin_p, b0_p, b01_p, b02_p, emin_p
USE control_ev,       ONLY : ieos
USE io_global,        ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: iu_mur, ipress
INTEGER :: find_free_unit

IF (ntemp_plot==0.AND.npress_plot==0) RETURN

IF (ionode) THEN
   iu_mur=find_free_unit()
   filename="energy_files/"//TRIM(flevdat)//'_mur_press'
   CALL add_pressure(filename)
   OPEN(UNIT=iu_mur, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   IF (ieos==2) THEN
      WRITE(iu_mur,'("#",2x,"press (kbar)",6x,"v0 (a.u.)^3",6x,&
       & "b0 (kbar) ",5x,"b01 ",8x," b02 (1/kbar)")') 
      DO ipress=1,npress
         WRITE(iu_mur,'(5e16.8)') press(ipress), vmin_p(ipress), &
                           b0_p(ipress), b01_p(ipress), b02_p(ipress)     
      ENDDO
   ELSE
      WRITE(iu_mur,'("#",2x,"press (kbar)",6x,"v0 (a.u.)^3",6x,&
       & "b0 (kbar) ",5x,"b01 ")') 
      DO ipress=1,npress
         WRITE(iu_mur,'(4e16.8)') press(ipress), vmin_p(ipress), &
                                   b0_p(ipress), b01_p(ipress)
      ENDDO
   ENDIF
   CLOSE(iu_mur)
ENDIF

RETURN
END SUBROUTINE write_mur_p
!
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
       & " Gibbs energy(p,T)(Ry)"," pressure (kbar),   p(T=0K), T=",f8.2," K")') temper
   ELSE
      WRITE(iu_mur,'( "#","omega(a.u.)**3",3x,"Helm.&
       & free-energy(Ry)"," Gibbs energy(p,T) (Ry)",  &
       & " pressure (kbar), p(T=0K), T=",f8.2," K")') temper
   ENDIF
ENDIF

RETURN
END SUBROUTINE write_mur_start_line
