!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_e_nk()
USE kinds,          ONLY : DP
USE control_conv,   ONLY : nnk, nk_test, nsigma
USE thermo_mod,     ONLY : energy_geo
USE data_files,     ONLY : flnkconv
USE io_global,      ONLY : ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER :: ink, isigma, iu_enk
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=256) :: filename
LOGICAL :: exst, parallelfs

IF (my_image_id /= root_image) RETURN

iu_enk=2
DO isigma=1, nsigma
   IF (nsigma > 1) THEN
      filename=TRIM(flnkconv)//int_to_char(isigma)
      CALL check_tempdir ( filename, exst, parallelfs )     
      filename=TRIM(filename)//'/'//TRIM(flnkconv)
   ELSE
      filename=TRIM(flnkconv)
   END IF
   IF (ionode) THEN
      OPEN(UNIT=iu_enk, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
      WRITE(iu_enk,'("#   nk       E_tot (Ry) ")' )
      DO ink = 1, nnk
         WRITE(iu_enk, '(i5,e20.10)') nk_test(ink), energy_geo(ink + &
                                     (isigma -1) * nnk ) 
      END DO
      CLOSE(iu_enk)
   END IF
END DO

RETURN
END SUBROUTINE write_e_nk
