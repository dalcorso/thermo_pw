!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------
SUBROUTINE read_state_densities()
!----------------------------------------------------

USE kinds,             ONLY : DP
USE ions_base,         ONLY : nat
USE control_2d_bands,  ONLY : averag, vacuum, nlayers, identify_sur, &
                              surface1, surface2
USE data_files,        ONLY : flprojlayer
USE io_global,         ONLY : ionode, ionode_id
USE mp_images,         ONLY : intra_image_comm
USE mp,                ONLY : mp_bcast

IMPLICIT NONE

LOGICAL :: exst
INTEGER :: iun, ios, idum, ilayer, ik, ibnd, nspin, nat_, nbnd_, nkstot_
INTEGER :: ispin
INTEGER :: find_free_unit

IF (identify_sur) THEN
   IF (ionode) &
        INQUIRE( FILE = TRIM(flprojlayer), EXIST = exst ) 
     CALL mp_bcast(exst,ionode_id,intra_image_comm)
     
     IF (exst) THEN
        iun=find_free_unit()
        IF (ionode) THEN
           OPEN(UNIT=iun,FILE=TRIM(flprojlayer),STATUS='old',ERR=300,&
                                                              IOSTAT=ios)
           READ(iun, '(5i8)') nat, nlayers, nbnd_, nkstot_, nspin     
        ENDIF
300     CALL mp_bcast(ios,ionode_id,intra_image_comm)
        IF (ios /= 0) CALL errore('read_state_densities','problems with flprojlayer',ABS(ios))
        CALL mp_bcast(nat, ionode_id, intra_image_comm)
        CALL mp_bcast(nlayers, ionode_id, intra_image_comm)
        CALL mp_bcast(nbnd_, ionode_id, intra_image_comm)
        CALL mp_bcast(nkstot_, ionode_id, intra_image_comm)
        CALL mp_bcast(nspin, ionode_id, intra_image_comm)

        ALLOCATE(averag(nat, nspin, nbnd_, nkstot_))
        ALLOCATE(vacuum(nspin, nbnd_, nkstot_))
        IF (ionode) THEN
           READ(iun,'(2i8)') surface1, surface2
           DO ik = 1, nkstot_
              DO ibnd = 1, nbnd_
                 READ(iun,'(2i8)') idum, idum
                 DO ilayer=1,nlayers
                    READ(iun,'(i8,4f17.12)') idum, &
                        (averag(ilayer, ispin, ibnd, ik), ispin=1,nspin)
                 ENDDO
                 READ(iun,'(4f20.12)') (vacuum(ispin, ibnd, ik), ispin=1,nspin)
              ENDDO
           ENDDO
           CLOSE(iun)
        ENDIF
        CALL mp_bcast(surface1, ionode_id, intra_image_comm)
        CALL mp_bcast(surface2, ionode_id, intra_image_comm)
        CALL mp_bcast(averag, ionode_id, intra_image_comm)
        CALL mp_bcast(vacuum, ionode_id, intra_image_comm)
     ELSE
100     CALL errore('read_state_densities','problem finding projection file',1)
     ENDIF   
  END IF

  RETURN
  END SUBROUTINE read_state_densities
