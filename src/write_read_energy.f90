!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE read_energy(nwork, filename)
  !-----------------------------------------------------------------------
  !
  !  This routine reads the energy for each run if this is
  !  already on file. It must be called by all processors, only the
  !  meta_ionode reads the data and sends them to all the others
  !
  !
  USE kinds,      ONLY : DP
  USE thermo_mod, ONLY : energy_geo
  USE io_global,  ONLY : meta_ionode, meta_ionode_id
  USE mp_world,   ONLY : world_comm
  USE mp,         ONLY : mp_bcast

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nwork
  CHARACTER(LEN=256), INTENT(IN) :: filename
  INTEGER :: iu_ev, iwork, ios
  !
  IF (meta_ionode) THEN
     iu_ev=2
     OPEN(UNIT=iu_ev, FILE=TRIM(filename), STATUS='OLD', FORM='FORMATTED', &
         ERR=20, IOSTAT=ios)
     DO iwork=1,nwork
        READ(iu_ev,'(e30.15)', ERR=20, IOSTAT=ios) energy_geo(iwork)
     ENDDO
     !
     CLOSE(iu_ev)
     !
  END IF

20 CALL mp_bcast(ios, meta_ionode_id, world_comm)
  IF (ios /= 0 ) CALL errore('read_ev_dat','opening or reading input file',1)
  CALL mp_bcast(energy_geo, meta_ionode_id, world_comm)
  !
  RETURN
  !
END SUBROUTINE read_energy
!
!-------------------------------------------------------------------------
SUBROUTINE write_energy(nwork, filename)
  !-----------------------------------------------------------------------
  !
  !  This routine reads the energy at each lattice constant if this is
  !  already on file. It must be called by all processors, only the
  !  meta_ionode reads the data and sends them to all the others
  !
  !
  USE kinds,      ONLY : DP
  USE thermo_mod, ONLY : energy_geo
  USE io_global,  ONLY : meta_ionode, meta_ionode_id
  USE mp_world,   ONLY : world_comm
  USE mp,         ONLY : mp_bcast

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nwork
  CHARACTER(LEN=256), INTENT(IN) :: filename
  INTEGER :: iu_ev, iwork, ios
  !
  IF (meta_ionode) THEN
     iu_ev=2
     OPEN(UNIT=iu_ev, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED', &
         ERR=20, IOSTAT=ios)
     DO iwork=1,nwork
        WRITE(iu_ev,'(e30.15)', ERR=20, IOSTAT=ios) energy_geo(iwork)
     ENDDO
     !
     CLOSE(iu_ev)
     !
  END IF
  !
20 CALL mp_bcast(ios, meta_ionode_id, world_comm)
  IF (ios /= 0 ) CALL errore('write_ev_dat','opening or writing output file',1)
  !
  RETURN
  !
END SUBROUTINE write_energy
