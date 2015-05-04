!
! Copyright (C) 2014-2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE write_gnuplot_energy(nwork)
  !-----------------------------------------------------------------------
  !
  !  This routine writes the energy at each geometry in a form that
  !  can be read by gnuplot to make two dimensional plots of the 
  !  energy.
  !  It must be called by all processors, only the meta_ionode writes 
  !  the data.
  !
  !
  USE kinds,      ONLY : DP
  USE thermo_mod, ONLY : energy_geo, ngeo, celldm_geo, omega_geo
  USE cell_base,  ONLY : ibrav
  USE control_mur, ONLY : lmurn
  USE data_files, ONLY : flenergy
  USE io_global,  ONLY : stdout, meta_ionode, meta_ionode_id
  USE mp_world,   ONLY : world_comm
  USE mp,         ONLY : mp_bcast

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nwork
  CHARACTER(LEN=256) :: fileout
  INTEGER :: iu_ev, iwork, ifiles, nfiles, ios
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=8) :: float_to_char
  !
  SELECT CASE (ibrav)
     CASE(1,2,3)
        nfiles=1
     CASE(4,5,6,7)
        nfiles=1
     CASE(8,9,91,10,11)
        nfiles=ngeo(3)
     CASE DEFAULT
     RETURN
  END SELECT

  IF (meta_ionode) THEN
     iu_ev=2
     DO ifiles = 1, nfiles
        fileout=TRIM(flenergy)//int_to_char(ifiles)
        OPEN(UNIT=iu_ev, FILE=TRIM(fileout), STATUS='UNKNOWN', &
                         FORM='FORMATTED', ERR=20, IOSTAT=ios)

        IF (lmurn) THEN
           DO iwork=1,ngeo(1)
              WRITE(iu_ev,'(2e30.15)', ERR=20, IOSTAT=ios) &
                              celldm_geo(1,iwork), energy_geo(iwork) 
           ENDDO
        ELSE
           SELECT CASE (ibrav)
              CASE(1,2,3)
                 DO iwork=1,ngeo(1)
                    WRITE(iu_ev,'(2e30.15)', ERR=20, IOSTAT=ios) &
                              celldm_geo(1,iwork), energy_geo(iwork) 
                 ENDDO
              CASE(4,6,7)
                 DO iwork=1,ngeo(1)*ngeo(3)
                     WRITE (iu_ev,'(3e25.12)', ERR=20, IOSTAT=ios) &
                           celldm_geo(1,iwork), &
                           celldm_geo(3,iwork), & 
                           energy_geo(iwork) 
                 ENDDO
              CASE(5)
                 DO iwork=1,ngeo(1)*ngeo(4)
                    WRITE(iu_ev,'(3e25.12)', ERR=20, IOSTAT=ios) &
                         celldm_geo(1,iwork), &
                         celldm_geo(4,iwork), &
                         energy_geo(iwork) 
                 ENDDO
              CASE(8,9,91,10,11)
                 DO iwork=1+(ifiles-1)*ngeo(1)*ngeo(2), ifiles*ngeo(1)*ngeo(2)
                    WRITE(iu_ev,'(3e25.12)', ERR=20, IOSTAT=ios) &
                        celldm_geo(1,iwork), &
                        celldm_geo(2,iwork), energy_geo(iwork)  
                                          
                 ENDDO
           END SELECT
        ENDIF
        CLOSE(iu_ev)
     ENDDO
     !
  END IF
  !
20 CALL mp_bcast(ios, meta_ionode_id, world_comm)
  IF (ios /= 0 ) CALL errore('write_ev_dat','opening or writing output file',1)
  !
  RETURN
  !
END SUBROUTINE write_gnuplot_energy

