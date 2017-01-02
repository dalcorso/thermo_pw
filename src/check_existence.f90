!
! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE check_existence(iwork, part, iaux, run)
!
!   This routine can be called separately by each image
!
  USE kinds,           ONLY : DP
  USE constants,       ONLY : ry_kbar

  USE control_thermo,  ONLY : lpwscf, lstress

  USE thermo_mod,      ONLY : energy_geo
  USE elastic_constants, ONLY : sigma_geo

  USE io_global,       ONLY : ionode, ionode_id, stdout
  USE mp_images,       ONLY : intra_image_comm
  USE mp,              ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iwork, part, iaux
  LOGICAL, INTENT(OUT) :: run

  INTEGER  :: iu_ene, ios, ipol, jpol
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=256) :: filename
  REAL(DP) :: energy, stress(3,3)
  LOGICAL :: exst
  run=.TRUE.
  
  iu_ene=2
  IF (ionode) THEN
     filename='restart/e_work_part.'//TRIM(int_to_char(iwork))//'.'//&
                              TRIM(int_to_char(part+iaux))
     INQUIRE(FILE=TRIM(filename),EXIST=exst)
     IF (exst) THEN
        WRITE(stdout,'(5x,"Data found on file")')
        OPEN(UNIT=iu_ene, FILE=TRIM(filename), STATUS='OLD', FORM='FORMATTED', &
            ERR=20, IOSTAT=ios)

        READ(iu_ene,*) energy
        WRITE(stdout,'(5x,"Total energy = ",f20.8," Ry")') energy
        IF (iwork > 0) THEN
           IF (lstress(iwork)) THEN
              WRITE(stdout,'(5x,"Stress (kbar)")')
              DO ipol=1,3
                 READ(iu_ene,*) (stress(ipol,jpol),jpol=1,3)
                 WRITE(stdout,'(3f15.7)') (stress(ipol,jpol)*ry_kbar,jpol=1,3)
              ENDDO
           END IF
        ELSE
           WRITE(stdout,'(5x,"Stress (kbar)")')
           DO ipol=1,3
              READ(iu_ene,*) (stress(ipol,jpol),jpol=1,3)
              WRITE(stdout,'(3f15.7)') (stress(ipol,jpol)*ry_kbar,jpol=1,3)
           ENDDO
        END IF
        run=.FALSE.
        CLOSE(iu_ene)
     ENDIF
!
!  If the file cannot be opened it is considered unexistent and the program
!  continue
!
20   CONTINUE
  ENDIF
  CALL mp_bcast(run, ionode_id, intra_image_comm)
  IF (.NOT.run .AND. iwork > 0) THEN
     CALL mp_bcast(energy, ionode_id, intra_image_comm)
     energy_geo(iwork)=energy
     IF (lstress(iwork)) THEN
        CALL mp_bcast(stress, ionode_id, intra_image_comm)
        sigma_geo(:,:,iwork)=stress(:,:)
     END IF
  END IF

  RETURN
END SUBROUTINE check_existence

SUBROUTINE save_existence(iwork, part, iaux)

  USE control_thermo,  ONLY : lpwscf, lstress

  USE thermo_mod,      ONLY : energy_geo
  USE elastic_constants, ONLY : sigma_geo

  USE ener,            ONLY : etot
  USE force_mod,       ONLY : sigma

  USE io_global,       ONLY : ionode, ionode_id, meta_ionode_id
  USE mp_images,       ONLY : intra_image_comm
  USE mp,              ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iwork, part, iaux

  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=256) :: filename
  INTEGER :: iu_ene, ios, ipol, jpol

  iu_ene=2
  IF (ionode) THEN
     filename='restart/e_work_part.'//TRIM(int_to_char(iwork))//'.'//&
                              TRIM(int_to_char(part+iaux))
     OPEN(UNIT=iu_ene, FILE=TRIM(filename), STATUS='UNKNOWN', &
                       FORM='FORMATTED', ERR=20, IOSTAT=ios)

     IF (iwork > 0) THEN
        WRITE(iu_ene,*) energy_geo(iwork)
        IF (lstress(iwork)) THEN
           DO ipol=1,3
              WRITE(iu_ene,*) (sigma_geo(ipol,jpol,iwork),jpol=1,3)
           ENDDO
        ENDIF
     ELSE
        WRITE(iu_ene,*) etot
        DO ipol=1,3
           WRITE(iu_ene,*) (sigma(ipol,jpol),jpol=1,3)
        ENDDO
     ENDIF
     CLOSE(iu_ene)
!
!  If the file cannot be written we do nothing
!
20   CONTINUE
  END IF
  
  RETURN
END SUBROUTINE save_existence
