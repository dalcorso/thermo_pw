!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

SUBROUTINE read_state_densities()

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

  SUBROUTINE read_minimal_info(read_not_write,ierr)
  !
  !  This routine writes and reads a minimal amount of information necessary
  !  to make a plot of the bands without redoing the self consistent calculation
  !  and the band structure calculation.
  !  When it writes it assumes that a scf calculation has been done so it
  !  has all the information on the crystal and on the scf calculation.
  !  When it reads, it assumes that the code has the information read 
  !  on the input of pw and on the input of thermo_pw, it has all the parallelization
  !  structure already in place but it has not the k points
  !  paths and the Fermi energy, so it read it here. The name of the file
  !  used here is internal and not available to the user
  !
  USE kinds,             ONLY : DP
  USE ener,              ONLY : ef
  USE klist,             ONLY : nelec
  USE control_paths,     ONLY : nqaux, disp_nqs, label_disp_q, nrap_plot, &
                                rap_plot, q2d
  USE klist,             ONLY : nkstot
  USE control_2d_bands,  ONLY : nkz, aux_ind_sur, sym_divide
  USE io_global,         ONLY : ionode, ionode_id
  USE mp_images,         ONLY : intra_image_comm
  USE mp,                ONLY : mp_bcast

  USE kinds, ONLY : DP
  IMPLICIT NONE
  LOGICAL :: read_not_write
  INTEGER :: ierr
  LOGICAL :: exst
  INTEGER :: iun, ios, n, i, ik, ikz, nqaux_, nks_
  INTEGER :: find_free_unit

  ierr=0
  IF (read_not_write) THEN
     IF (ionode) &
        INQUIRE( FILE = 'band_files/info_data', EXIST = exst )
     CALL mp_bcast(exst, ionode_id, intra_image_comm)
     IF (.NOT. exst) THEN
         CALL errore('read_minimal_info','info_data not present',1)
         ef=0.0_DP
         ierr=-1
         RETURN
     ENDIF
     IF (ionode) THEN
        iun=find_free_unit()
        OPEN(UNIT=iun, FILE='band_files/info_data', STATUS='old', ERR=300,  &
                                                                IOSTAT=ios)
        READ(iun, *) ef, nelec
        READ(iun, *) nqaux_
        READ(iun, *) q2d
        IF (.NOT.q2d) READ(iun, '(9i8)') (label_disp_q(n), n=1,nqaux)
        READ(iun, *) disp_nqs
        IF (disp_nqs > 0.AND..NOT.q2d) THEN
           ALLOCATE(nrap_plot(disp_nqs))
           ALLOCATE(rap_plot(12,disp_nqs))
           nrap_plot=0
           rap_plot=0
        ENDIF
        DO n=1,disp_nqs
           READ(iun,*) nrap_plot(n)
           IF (nrap_plot(n) > 0)  &
              READ(iun,'(9i8)') (rap_plot(i,n), i=1,nrap_plot(n))
        ENDDO
        IF (nkz>1 .AND. sym_divide) THEN
           READ(iun, '(2i8)') nkz, nks_
           ALLOCATE(aux_ind_sur(nks_, nkz))
           DO ik=1, nks_
              READ(iun, '(15i4)') (aux_ind_sur(ik,ikz), ikz=1,nkz)
           ENDDO
        ENDIF

        CLOSE(iun)
     ENDIF
300  CALL mp_bcast(ios,ionode_id,intra_image_comm)
     IF (ios /= 0) CALL errore('read_minimal_info','problems with file',ABS(ios))
     CALL mp_bcast(nelec, ionode_id, intra_image_comm)
     CALL mp_bcast(ef, ionode_id, intra_image_comm)
     IF (nqaux>0) CALL mp_bcast(label_disp_q, ionode_id, intra_image_comm)
     CALL mp_bcast(disp_nqs, ionode_id, intra_image_comm)
     IF (nkz > 1 .AND. sym_divide) &
        CALL mp_bcast(nks_, ionode_id, intra_image_comm)
     IF (.NOT.ionode) THEN
        IF (disp_nqs>0) THEN
           ALLOCATE(nrap_plot(disp_nqs))
           ALLOCATE(rap_plot(12,disp_nqs))
        ENDIF
        IF (nkz>1 .AND. sym_divide) ALLOCATE(aux_ind_sur(nks_,nkz))
     ENDIF
     IF (disp_nqs>0) THEN
        CALL mp_bcast(nrap_plot, ionode_id, intra_image_comm)
        CALL mp_bcast(rap_plot, ionode_id, intra_image_comm)
     ENDIF
     IF (nkz > 1 .AND. sym_divide) &
        CALL mp_bcast(aux_ind_sur, ionode_id, intra_image_comm)
  ELSE
     IF (ionode) THEN
        iun=find_free_unit()
        OPEN(UNIT=iun,FILE='band_files/info_data',STATUS='unknown',&
                                                          ERR=400,IOSTAT=ios)
        WRITE(iun, '(2f16.11)') ef, nelec
        WRITE(iun, '(i8)') nqaux
        WRITE(iun, '(l8)') q2d
        IF (.NOT.q2d) WRITE(iun, '(9i8)') (label_disp_q(n), n=1,nqaux)
        WRITE(iun, '(i8)') disp_nqs
        IF (.NOT.q2d) THEN
           DO n=1,disp_nqs
              WRITE(iun,'(i8)') nrap_plot(n)
              IF (nrap_plot(n) > 0)  &
                 WRITE(iun,'(9i8)') (rap_plot(i,n), i=1,nrap_plot(n))
           ENDDO
        ENDIF
        IF (nkz>1.AND.sym_divide) THEN
           WRITE(iun, '(2i8)') nkz, nkstot / nkz
           DO ik=1, nkstot / nkz
              WRITE(iun, '(15i4)') (aux_ind_sur(ik,ikz), ikz=1,nkz)
           ENDDO
        ENDIF 
        CLOSE(iun)
     ENDIF
400  CALL mp_bcast(ios,ionode_id,intra_image_comm)
     IF (ios /= 0) CALL errore('read_minimal_info','problems writing file',ABS(ios))
  ENDIF

  RETURN
  END SUBROUTINE read_minimal_info

