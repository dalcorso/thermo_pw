!
! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_quadratic()
!
! This routine writes the energy as a function of the celldm(1) for
! the polynomial fit of cubic systems.
! in input emin in Ry, omega0 in (a.u.)**3, b0in in kbar
! in output omega in (a.u.)**3, p in kbar, e in Ry
!
USE kinds,          ONLY : DP
USE data_files,     ONLY : flevdat
USE thermo_mod,     ONLY : ngeo, energy_geo, celldm_geo
USE cell_base,      ONLY : ibrav
USE control_mur,    ONLY : nvol
USE control_quadratic_energy, ONLY : coeff
USE control_pressure, ONLY : pressure, pressure_kb
USE mp_images,      ONLY : root_image, my_image_id
USE constants,      ONLY : ry_kbar
USE io_global,      ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: filename, filename1
CHARACTER(LEN=8) :: float_to_char
REAL(DP) :: a, e, p, fact, xmax, xmin, deltaa
INTEGER :: i, iu_mur, iwork

IF (my_image_id /= root_image) RETURN

IF (ibrav<1.OR.ibrav>3) RETURN
IF (ibrav==1) fact=1.0_DP
IF (ibrav==2) fact=0.25_DP
IF (ibrav==3) fact=0.5_DP

filename=TRIM(flevdat)//'_quadratic'
IF (pressure /= 0.0_DP) &
   filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))

xmax=0.0_DP
xmin=100000.0_DP
DO iwork = 1, ngeo(1)
   IF (celldm_geo(1,iwork) > xmax) xmax=celldm_geo(1,iwork)
   IF (celldm_geo(1,iwork) < xmin) xmin=celldm_geo(1,iwork)
ENDDO
xmin=xmin*0.99_DP
xmax=xmax*1.01_DP
IF (nvol <= 1) nvol=50
deltaa = (xmax - xmin)/(nvol-1)

IF (ionode) THEN
   iu_mur=2
   OPEN(UNIT=iu_mur, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   IF (pressure /= 0.0_DP) THEN
      WRITE(iu_mur,'( "# a (a.u.)       enthalpy (Ry)    pressure (kbar)" )')
   ELSE
      WRITE(iu_mur,'( "# a (a.u.)       energy (Ry)      pressure (kbar)" )')
   END IF

   DO i=1,nvol
      a= xmin + deltaa * (i-1)
      e= coeff(1) + coeff(2)*a + coeff(3)*a**2
      p= -(coeff(2) + 2.0_DP * coeff(3) * a ) * ry_kbar / (fact*3.0_DP * a**2)
      WRITE(iu_mur,'(3f20.10)') a, e, p + pressure_kb
   ENDDO 

   CLOSE(UNIT=iu_mur, STATUS='KEEP')
END IF

RETURN
END SUBROUTINE write_quadratic
