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
  !  can be read by gnuplot to make two dimensional countour plots of the 
  !  energy.
  !  It must be called by all processors, only the meta_ionode writes 
  !  the data.
  !
  !
  USE kinds,      ONLY : DP
  USE thermo_mod, ONLY : energy_geo_eos, ngeo, celldm_geo_eos, omega_geo_eos
  USE cell_base,  ONLY : ibrav
  USE data_files, ONLY : flenergy
  USE control_pressure, ONLY : pressure
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE mp_images,  ONLY : my_image_id, root_image, intra_image_comm
  USE mp,         ONLY : mp_bcast

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nwork
  CHARACTER(LEN=256)  :: fileout
  INTEGER :: iu_ev, iwork, ifiles, nfiles, ios
  INTEGER :: find_free_unit
  CHARACTER(LEN=6) :: int_to_char
  !
  IF (my_image_id /= root_image) RETURN

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

  IF (ionode) THEN
     iu_ev=find_free_unit()
     DO ifiles = 1, nfiles
        fileout='energy_files/'//TRIM(flenergy)//int_to_char(ifiles)
        CALL add_pressure(fileout)
        OPEN(UNIT=iu_ev, FILE=TRIM(fileout), STATUS='UNKNOWN', &
                         FORM='FORMATTED', ERR=20, IOSTAT=ios)

        SELECT CASE (ibrav)
           CASE(1,2,3)
              DO iwork=1,nwork
                 WRITE(iu_ev,'(2e30.15)', ERR=20, IOSTAT=ios) &
                           celldm_geo_eos(1,iwork), energy_geo_eos(iwork) &
                                              + pressure * omega_geo_eos(iwork)
              ENDDO
           CASE(4,6,7)
              DO iwork=1,nwork
                  WRITE (iu_ev,'(3e25.12)', ERR=20, IOSTAT=ios) &
                        celldm_geo_eos(1,iwork), &
                        celldm_geo_eos(3,iwork), & 
                        energy_geo_eos(iwork) + pressure * omega_geo_eos(iwork)
              ENDDO
           CASE(5)
              DO iwork=1,nwork
                 WRITE(iu_ev,'(3e25.12)', ERR=20, IOSTAT=ios) &
                      celldm_geo_eos(1,iwork), &
                      celldm_geo_eos(4,iwork), &
                      energy_geo_eos(iwork) + pressure * omega_geo_eos(iwork)
              ENDDO
           CASE(8,9,91,10,11)
              DO iwork=1+(ifiles-1)*ngeo(1)*ngeo(2), ifiles*ngeo(1)*ngeo(2)
                 WRITE(iu_ev,'(3e25.12)', ERR=20, IOSTAT=ios) &
                     celldm_geo_eos(1,iwork), &
                     celldm_geo_eos(2,iwork), energy_geo_eos(iwork) &  
                              + pressure * omega_geo_eos(iwork)
              ENDDO
           CASE DEFAULT
              DO iwork = 1, nwork
                 WRITE(iu_ev,'(7e20.12)', ERR=20, IOSTAT=ios)  &
                     celldm_geo_eos(1,iwork), celldm_geo_eos(2,iwork), &
                     celldm_geo_eos(3,iwork), celldm_geo_eos(4,iwork), &
                     celldm_geo_eos(4,iwork), celldm_geo_eos(5,iwork), &
                     celldm_geo_eos(6,iwork), energy_geo_eos(iwork)    &  
                              + pressure * omega_geo_eos(iwork)
              ENDDO
        END SELECT
        CLOSE(iu_ev)
     ENDDO
     !
  END IF
  !
20 CALL mp_bcast(ios, ionode_id, intra_image_comm)
  IF (ios /= 0 ) CALL errore('write_gnuplot_energy',&
                                  'opening or writing output file',1)
  !
  RETURN
  !
END SUBROUTINE write_gnuplot_energy
