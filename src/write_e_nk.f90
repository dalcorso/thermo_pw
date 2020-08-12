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
USE io_files,       ONLY : check_tempdir
USE io_global,      ONLY : ionode
USE mp_world,       ONLY : world_comm
USE mp_images,      ONLY : my_image_id, root_image, nproc_image
USE mp,             ONLY : mp_sum

IMPLICIT NONE
INTEGER :: ink, isigma, iu_enk
INTEGER :: find_free_unit
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=256) :: filename
LOGICAL :: exst, parallelfs
!
!  First collect the total energies
!
CALL mp_sum(energy_geo, world_comm)
energy_geo=energy_geo / nproc_image

IF (my_image_id /= root_image) RETURN

iu_enk=find_free_unit()
DO isigma=1, nsigma
   IF (nsigma > 1) THEN
      filename='energy_files/'//TRIM(flnkconv)//int_to_char(isigma)
      CALL check_tempdir ( filename, exst, parallelfs )     
      filename=TRIM(filename)//'/'//TRIM(flnkconv)
   ELSE
      filename='energy_files/'//TRIM(flnkconv)
   ENDIF
   IF (ionode) THEN
      OPEN(UNIT=iu_enk, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
      WRITE(iu_enk,'("#   nk1   nk2    nk3          E_tot (Ry) ")' )
      DO ink = 1, nnk
         WRITE(iu_enk, '(3i5,e20.10)') nk_test(1, ink), nk_test(2, ink), &
                                       nk_test(3, ink), energy_geo(ink + &
                                     (isigma -1) * nnk ) 
      ENDDO
      CLOSE(UNIT=iu_enk, STATUS='KEEP')
   ENDIF
ENDDO

RETURN
END SUBROUTINE write_e_nk
