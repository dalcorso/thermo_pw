!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_e_ke()
USE kinds,          ONLY : DP
USE control_conv,   ONLY : nke, ke, nkeden
USE thermo_mod,     ONLY : energy_geo
USE data_files,     ONLY : flkeconv
USE io_global,      ONLY : ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER :: ike, iden, iu_eke
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=256) :: filename
LOGICAL :: exst, parallelfs

IF (my_image_id /= root_image) RETURN

iu_eke=2
DO iden=1, nkeden
   IF (nkeden > 1) THEN
      filename='energy_files/'//TRIM(flkeconv)//int_to_char(iden)
      CALL check_tempdir ( filename, exst, parallelfs )     
      filename=TRIM(filename)//'/'//TRIM(flkeconv)
   ELSE
      filename='energy_files/'//TRIM(flkeconv)
   END IF
   IF (ionode) THEN
      OPEN(UNIT=iu_eke, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
      WRITE(iu_eke,'("#   E_kin (Ry)       E_tot (Ry) ")' )
      DO ike = 1, nke
         WRITE(iu_eke, '(2e20.10)') ke(ike), energy_geo( ike + (iden-1) * nke ) 
      END DO
      CLOSE(iu_eke)
   END IF
END DO

RETURN
END SUBROUTINE write_e_ke
