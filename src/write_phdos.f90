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
  USE kinds,          ONLY : DP
  USE mp,             ONLY : mp_sum, mp_bcast, mp_max, mp_min
  USE mp_world,       ONLY : world_comm
  USE mp_images,      ONLY : my_image_id, root_image 
  USE constants,      ONLY : amu_ry
  USE io_global,      ONLY : meta_ionode, meta_ionode_id, stdout
  USE ions_base,      ONLY : nat, amass
  USE control_dosq,   ONLY : phdos_sigma, deltafreq, freqmin, freqmax, &
                             ndos_input, freqmin_input, freqmax_input
  USE phonon_save,    ONLY : freq_save, z_save
  USE thermo_mod,     ONLY : tot_ngeo
  USE control_thermo, ONLY : with_eigen
  USE thermodynamics, ONLY : phdos_save, gen_phdos_save
  USE ph_freq_thermodynamics, ONLY : ph_freq_save
  USE control_dosq,   ONLY : dos_wq, dos_nqs
  USE data_files,     ONLY : fldos
  USE phdos_module,   ONLY : set_phdos, read_phdos_data, find_minimum_maximum,&
                             set_gen_phdos
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: igeom

  CHARACTER(LEN=256) :: filedos, filename
  REAL(DP) :: e, emin, emax, dosofe(2)
  REAL(DP), ALLOCATABLE :: gen_dos(:,:)
  INTEGER :: n, i, ndos, nq, startq, lastq, nq_eff, na, ijpol, iundos
  INTEGER :: find_free_unit
  LOGICAL :: check_file_exists
  CHARACTER(LEN=6) :: int_to_char
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
     ENDIF
     CALL mp_bcast(phdos_save(igeom)%number_of_points, meta_ionode_id, &
                                                       world_comm)
     CALL mp_bcast(phdos_save(igeom)%de, meta_ionode_id, world_comm)
     IF ( my_image_id /= root_image ) THEN
        IF (.NOT.ALLOCATED(phdos_save(igeom)%nu)) &
          ALLOCATE(phdos_save(igeom)%nu(phdos_save(igeom)%number_of_points)) 
        IF (.NOT.ALLOCATED(phdos_save(igeom)%phdos)) &
          ALLOCATE(phdos_save(igeom)%phdos(phdos_save(igeom)%number_of_points)) 
     END IF
     CALL mp_bcast(phdos_save(igeom)%nu, meta_ionode_id, world_comm)
     CALL mp_bcast(phdos_save(igeom)%phdos, meta_ionode_id, world_comm)

     CALL find_minimum_maximum(phdos_save(igeom), freqmin, freqmax)
     RETURN
  END IF

  WRITE(stdout,'(/,2x,76("+"))')
  WRITE(stdout,'(5x,"Writing phdos on file ",a)') 
  WRITE(stdout,'(5x,a)') TRIM(filedos)
  WRITE(stdout,'(2x,76("+"),/)')
  nq=dos_nqs
!
! compute the dos
!
  CALL divide(world_comm, nq, startq, lastq)
  emin = 0.0d0
  emax = 0.0d0
  DO n=startq, lastq
     DO i=1, 3*nat
        emin = MIN (emin, freq_save(i,n))
        emax = MAX (emax, freq_save(i,n))
     END DO
  END DO
  emax=emax*1.02_DP
  CALL mp_max(emax, world_comm)
  CALL mp_min(emin, world_comm)
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
  phdos_save(igeom)%nu=0.0_DP
  phdos_save(igeom)%phdos=0.0_DP
  nq_eff=ph_freq_save(igeom)%nq_eff
  startq=ph_freq_save(igeom)%startq
  DO n= 1, ndos
     e = emin + (n - 1) * deltafreq
     !
     CALL dos_g(freq_save(1,startq), 1, 3*nat, nq_eff, dos_wq(startq), &
                                         phdos_sigma, 0, e, dosofe)
     !
     phdos_save(igeom)%nu(n) = e
     phdos_save(igeom)%phdos(n) = dosofe(1)
  END DO

  IF (with_eigen) THEN
     ALLOCATE(gen_dos(6,nat))
     CALL set_gen_phdos(gen_phdos_save,ndos,nat,deltafreq)
     DO n= 1, ndos
        e = emin + (n - 1) * deltafreq
        CALL generalized_phdos(freq_save(1,startq), 3*nat, nq_eff, &
                              dos_wq(startq), phdos_sigma, 0, e, &
                              z_save(1,1,startq), nat, gen_dos)
        gen_phdos_save%nu(n) = e
        gen_phdos_save%phdos(:,:,n) = gen_dos(:,:)
     ENDDO
     CALL mp_sum(gen_phdos_save%phdos, world_comm)
     DO n= 1, ndos
        CALL symmetrize_gen_phdos(gen_phdos_save%phdos(1,1,n),nat)
     ENDDO
     DEALLOCATE(gen_dos)
  ENDIF
!
!  and collect the results
!
  CALL mp_sum(phdos_save(igeom)%phdos, world_comm)

  IF (meta_ionode) THEN
     iundos=find_free_unit()
     OPEN (UNIT=iundos,FILE=filedos,STATUS='unknown',FORM='formatted')
     DO n=1, ndos 
        WRITE (iundos, '(ES30.15, ES30.15)') phdos_save(igeom)%nu(n),  &
                                             phdos_save(igeom)%phdos(n)
     ENDDO
     CLOSE(UNIT=iundos,STATUS='keep')

     IF (with_eigen) THEN
        DO na=1,nat
           filename=TRIM(filedos)//'.'//TRIM(int_to_char(na))
           OPEN (UNIT=iundos,FILE=filename,STATUS='unknown',FORM='formatted')
           DO n=1, ndos 
              WRITE (iundos, '(ES20.11, 6ES20.11)') gen_phdos_save%nu(n),  &
                                (gen_phdos_save%phdos(ijpol,na,n),ijpol=1,6)
           ENDDO
           CLOSE(UNIT=iundos,STATUS='keep')
        ENDDO
     ENDIF
  ENDIF
  !
  RETURN
END SUBROUTINE write_phdos
!
SUBROUTINE symmetrize_gen_phdos(phdos, nat)
USE kinds, ONLY : DP
USE symme, ONLY : symtensor

IMPLICIT NONE
INTEGER :: nat
REAL(DP) :: phdos(6,nat)

REAL(DP) :: tot_phdos(3,3,nat)

INTEGER :: ipol, jpol, ijpol, na

DO na=1,nat
   ijpol=0
   DO ipol=1,3
      DO jpol=ipol,3
         ijpol=ijpol+1
         tot_phdos(ipol,jpol, na)= phdos(ijpol,na)
         IF (ipol/=jpol) tot_phdos(jpol, ipol, na)=tot_phdos(ipol, jpol, na)
      ENDDO
   ENDDO
ENDDO
CALL symtensor(nat, tot_phdos)

DO na=1,nat
   ijpol=0
   DO ipol=1,3
      DO jpol=ipol,3
         ijpol=ijpol+1
         phdos(ijpol, na)= tot_phdos(ipol,jpol,na)
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE symmetrize_gen_phdos
