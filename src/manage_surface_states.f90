!
! Copyright (C) 2018 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE manage_surface_states()
  !-----------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE ions_base,          ONLY : nat
  USE lsda_mod,           ONLY : nspin
  USE wvfct,              ONLY : nbnd
  USE klist,              ONLY : nkstot
  USE control_2d_bands,   ONLY : averag, vacuum, nlayers, surface1, surface2
  USE data_files,         ONLY : flprojlayer
  USE control_thermo,     ONLY : spin_component
  USE io_global,          ONLY : ionode, ionode_id, stdout
  USE mp,                 ONLY : mp_bcast
  USE mp_images,          ONLY : intra_image_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ibnd, ilayer, ios, iun
  INTEGER :: nks1, nks2, nks1tot, nks2tot
  LOGICAL :: exst
  INTEGER :: find_free_unit
!
!   Find which k points must be done by this pool
!
  CALL find_nks1nks2(1,nkstot,nks1tot,nks1,nks2tot,nks2,spin_component)

  IF (ionode) INQUIRE( FILE = TRIM(flprojlayer), EXIST = exst )
  CALL mp_bcast(exst, ionode_id, intra_image_comm)
!
!   the file with the projections is created here if it does not exist,
!   otherwise we assume that it has been already calculated in a previous run
!   and exit
!
  IF (exst) RETURN

  ALLOCATE(averag(nat, nspin, nbnd, nkstot))
  ALLOCATE(vacuum(nspin, nbnd, nkstot))
  CALL plan_avg_sub(averag, vacuum, nat, nbnd, nkstot, nlayers, surface1, &
                                                                surface2)
  IF (ionode) THEN
     iun=find_free_unit()
     OPEN(UNIT=iun, FILE=TRIM(flprojlayer), STATUS='unknown', ERR=400, &
                                                            IOSTAT=ios)
     WRITE(iun, '(5i8)') nat, nlayers, nbnd, nkstot, nspin     
     WRITE(iun, '(4i8)') surface1, surface2    
     DO ik=nks1tot, nks2tot
        DO ibnd=1, nbnd
           WRITE(iun,'(2i8)') ik, ibnd
           DO ilayer=1,nlayers
              WRITE(iun,'(i8,4f17.12)') ilayer, averag(ilayer, 1:nspin, &
                                                                  ibnd, ik)
           ENDDO
           WRITE(iun,'(4f20.12)') vacuum(1:nspin,ibnd, ik)
        ENDDO
     ENDDO
     CLOSE(iun)
  ENDIF
400  CALL mp_bcast(ios,ionode_id,intra_image_comm)
  IF (ios /= 0) CALL errore('sym_band_sub','problems with flprojlayer',1)

  DEALLOCATE (vacuum)
  DEALLOCATE (averag)
  !
  RETURN
END SUBROUTINE manage_surface_states
