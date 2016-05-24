! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
SUBROUTINE write_phdos(igeom)
  !-----------------------------------------------------------------------
  !
  !  This routine computes the phonon density of states. 
  !  It assumes that the frequencies are already known and are in the
  !  variable freq_save. It needs also the weight of each q point.
  !  After the call to this routine the variable phdos_save(igeom) 
  !  will contain the phonon dos.
  !  The dos are also written on file fldos. If the file is found in
  !  the phdisp_files directory the routine will read the dos, set the 
  !  phdos_save(igeom) variables and exits.
  !
  USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_sum
  USE mp_images,  ONLY : my_image_id, root_image, intra_image_comm
  USE io_global,  ONLY : ionode, stdout
  USE ions_base,  ONLY : nat
  USE control_dosq, ONLY : phdos_sigma, deltafreq, freqmin, freqmax, &
                         ndos_input, &
                         freqmin_input, freqmax_input, nq1_d, nq2_d, nq3_d
  USE phonon_save,    ONLY : freq_save
  USE thermo_mod,     ONLY : tot_ngeo
  USE thermodynamics, ONLY : phdos_save
  USE control_dosq,   ONLY : dos_wq, dos_nqs
  USE data_files,     ONLY : fldos
  USE phdos_module,   ONLY : set_phdos, read_phdos_data, find_minimum_maximum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: igeom

  CHARACTER(LEN=256) :: filedos
  REAL(DP) :: e, emin, emax, dosofe(2)
  INTEGER :: n, i, ndos, nq, nstart, nlast
  LOGICAL :: check_file_exists
!
!  If phdos is on file it is read
!
  filedos="phdisp_files/"//TRIM(fldos)
  IF (.NOT.ALLOCATED(phdos_save) ) ALLOCATE(phdos_save(tot_ngeo))
  IF ( check_file_exists(filedos) ) THEN
     IF ( my_image_id == root_image ) THEN
        WRITE(stdout,'(/,2x,76("-"))')
        WRITE(stdout,'(5x,"Readin phdos from file ")') 
        WRITE(stdout,'(5x,a)') TRIM(filedos)
        WRITE(stdout,'(2x,76("-"),/)')
        CALL read_phdos_data(phdos_save(igeom),filedos)
        CALL find_minimum_maximum(phdos_save(igeom), freqmin, freqmax)
        RETURN
     END IF
  END IF

  IF ( my_image_id /= root_image ) RETURN

  WRITE(stdout,'(/,2x,76("+"))')
  WRITE(stdout,'(5x,"Writing phdos on file ",a)') 
  WRITE(stdout,'(5x,a)') TRIM(filedos)
  WRITE(stdout,'(2x,76("+"),/)')
  nq=dos_nqs
!
! compute the dos
!
  emin = 0.0d0
  emax = 0.0d0
  DO n=1,nq
     DO i=1, 3*nat
        emin = MIN (emin, freq_save(i,n))
        emax = MAX (emax, freq_save(i,n))
     END DO
  END DO
  emax=emax*1.02_DP
  !
  IF (freqmin_input > 0.0_DP) emin=freqmin_input
  freqmin=emin
  !
  IF (freqmax_input > 0.0_DP) emax=freqmax_input
  freqmax=NINT(emax)
  !
  IF (ndos_input > 1) THEN
     deltafreq = (emax - emin)/(ndos_input-1)
     ndos = ndos_input
  ELSE
     ndos = NINT ( (emax - emin) / deltafreq + 1.51d0 )
     ndos_input = ndos
  END IF
!
! initialize the phdos_save space 
!
  CALL set_phdos(phdos_save(igeom),ndos,deltafreq)
!
!   Divide the calculation of the density of states among processors
!   of one image
!
  CALL divide (intra_image_comm, ndos, nstart, nlast)
  phdos_save(igeom)%nu=0.0_DP
  phdos_save(igeom)%phdos=0.0_DP
  DO n= nstart, nlast
     e = emin + (n - 1) * deltafreq
     !
     CALL dos_g(freq_save, 1, 3*nat, nq, dos_wq, phdos_sigma, 0, e, dosofe)
     !
     phdos_save(igeom)%nu(n) = e
     phdos_save(igeom)%phdos(n) = dosofe(1)
  END DO
!
!  and collect the results
!
  CALL mp_sum(phdos_save(igeom)%nu, intra_image_comm)
  CALL mp_sum(phdos_save(igeom)%phdos, intra_image_comm)

  IF (ionode) THEN
     OPEN (UNIT=2,FILE=filedos,STATUS='unknown',FORM='formatted')
     DO n=1, ndos 
        WRITE (2, '(ES30.15, ES30.15)') phdos_save(igeom)%nu(n),  &
                                             phdos_save(igeom)%phdos(n)
     ENDDO
     CLOSE(unit=2)
  END IF
  !
  RETURN
END SUBROUTINE write_phdos
!
