!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_e_ke()
USE kinds,          ONLY : DP
USE control_conv,   ONLY : nke, ke, nkeden, deltakeden, deltake
USE initial_param,  ONLY : ecutrho0, ecutwfc0
USE thermo_mod,     ONLY : energy_geo
USE data_files,     ONLY : flkeconv
USE io_global,      ONLY : ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER :: ike, iden, icount, iu_eke
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=256) :: filename
REAL(DP) :: kev, kedenv
LOGICAL :: exst, parallelfs

IF (my_image_id /= root_image) RETURN

iu_eke=2
icount = 0
DO iden=1, nkeden
   kedenv=ecutrho0 + (iden-1) * deltakeden
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
         kev = ecutwfc0 + (ike-1) * deltake
         IF (kedenv/kev > 3.9999_DP) THEN
            icount = icount + 1
            WRITE(iu_eke, '(2e20.10)') ke(icount), energy_geo( icount ) 
         ENDIF
      ENDDO
      CLOSE(UNIT=iu_eke, STATUS='KEEP')
   ENDIF
ENDDO

RETURN
END SUBROUTINE write_e_ke
