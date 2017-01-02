!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_gruneisen_band(file_disp, file_vec)
  !
  ! read data files produced by "bands_sub" for ngeo geometries, writes
  ! a file with the gruneisen parameters: the derivatives of the phonon 
  ! dispersions with respect to the volume (gruneisen parameters)
  ! This is calculated at the volume given in input or at the volume
  ! that corresponds to the temperature given in input.
  ! 
  USE kinds,          ONLY : DP
  USE ions_base,      ONLY : nat, ntyp => nsp
  USE control_thermo, ONLY : with_eigen
  USE data_files,     ONLY : flgrun
  USE thermo_mod,     ONLY : ngeo, omega_geo, no_ph
  USE anharmonic,     ONLY : vmin_t
  USE ph_freq_anharmonic, ONLY : vminf_t
  USE control_grun,   ONLY : temp_ph, volume_ph
  USE initial_conf,   ONLY : amass_save, ityp_save
  USE control_mur,    ONLY : vmin
  USE control_thermo, ONLY : ltherm_dos, ltherm_freq
  USE point_group,    ONLY : nsym_group
  USE temperature,    ONLY : temp, ntemp
  USE io_bands,       ONLY : read_bands, read_parameters, &
                             read_representations, write_bands
  USE mp,             ONLY : mp_bcast
  USE io_global,      ONLY : stdout, ionode, ionode_id
  USE mp_images,      ONLY : intra_image_comm, root_image, my_image_id

  IMPLICIT NONE

  CHARACTER(LEN=256), INTENT(IN) :: file_disp

  REAL(DP) :: eps=1.d-4
  REAL(DP), ALLOCATABLE :: freq_geo(:,:,:), k(:,:), k_rap(:,:), gaugek(:,:)
  REAL(DP), ALLOCATABLE :: frequency_geo(:,:)
  COMPLEX(DP), ALLOCATABLE :: displa_geo(:,:,:,:)
  INTEGER, ALLOCATABLE :: rap_geo(:,:,:), repres_geo(:,:), gcodek(:), &
                          aux_ind(:), gcodek_ext(:), ptypek(:,:), lprojk(:)
  LOGICAL, ALLOCATABLE :: same_next(:)
  INTEGER :: nks, nbnd, nks_rap, nbnd_rap 
  INTEGER :: ibnd, jbnd, irap, ios, i, n, ierr, igeo
  INTEGER :: iu_grun, iumode
  INTEGER :: poly_order
  REAL(DP), ALLOCATABLE :: poly_grun(:,:), frequency(:,:), gruneisen(:,:)
  REAL(DP) :: vm
  LOGICAL, ALLOCATABLE :: high_symmetry(:), is_gamma(:)
  LOGICAL :: copy_before, exist_rap, allocated_variables
  CHARACTER(LEN=256) :: filename, filedata, file_vec, filegrun, filefreq
  CHARACTER(LEN=6), EXTERNAL :: int_to_char


  IF ( my_image_id /= root_image ) RETURN

  IF (flgrun == ' ') RETURN

  iumode=23
  WRITE(stdout,*)
  exist_rap=.TRUE.
  allocated_variables=.FALSE.
  DO igeo = 1, ngeo(1)

     IF (no_ph(igeo)) CYCLE
     filedata = "phdisp_files/"//TRIM(file_disp)//'.g'//TRIM(int_to_char(igeo))
     CALL read_parameters(nks, nbnd, filedata)

     IF (nks <= 0 .or. nbnd <= 0) THEN
        CALL errore('write_gruneisen_band','reading plot namelist',ABS(ios))
     ELSE
        WRITE(stdout, '(5x,"Reading ",i4," dispersions at ",i6," k-points for&
                       & geometry",i4)') nbnd, nks, igeo
     ENDIF
     IF (.NOT.allocated_variables) THEN
        ALLOCATE (freq_geo(nbnd,nks,ngeo(1)))
        IF (with_eigen) ALLOCATE (displa_geo(nbnd,nbnd,ngeo(1),nks))
        ALLOCATE (rap_geo(nbnd,nks,ngeo(1)))
        ALLOCATE (k(3,nks)) 
        ALLOCATE (k_rap(3,nks))
        ALLOCATE (high_symmetry(nks))
        ALLOCATE (gcodek(nks))
        ALLOCATE (gcodek_ext(nks))
        ALLOCATE (aux_ind(nks))
        ALLOCATE (ptypek(3,nks))
        ALLOCATE (lprojk(nks))
        ALLOCATE (same_next(nks))
        ALLOCATE (gaugek(48,nks))
        ALLOCATE (is_gamma(nks))
        allocated_variables=.TRUE.
     ENDIF
     CALL read_bands(nks, nbnd, k, freq_geo(1,1,igeo), filedata)

     filename=TRIM(filedata)//".rap"
     k_rap(:,:)=k(:,:)
     rap_geo(:,:,igeo)=-1
     CALL read_representations(nks, nbnd, k_rap, rap_geo(1,1,igeo),     &
                          high_symmetry, gcodek, aux_ind, gcodek_ext,   &
                          ptypek, lprojk, same_next, gaugek, exist_rap, &
                          filename)
     nks_rap=nks
     nbnd_rap=nbnd
     DO n=1,nks
        is_gamma(n) = (( k(1,n)**2 + k(2,n)**2 + k(3,n)**2) < 1.d-12)
     ENDDO

     IF (with_eigen) THEN
        filename="phdisp_files/"//TRIM(file_vec)//".g"//TRIM(int_to_char(igeo))
        IF (ionode) OPEN(UNIT=iumode, FILE=TRIM(filename), FORM='formatted', &
                   STATUS='old', ERR=210, IOSTAT=ios)
210     CALL mp_bcast(ios, ionode_id, intra_image_comm)
        CALL errore('write_gruneisen_band','modes are needed',ABS(ios))
        IF (ionode) THEN
           CALL readmodes(nat,nks,k,displa_geo,ngeo(1),igeo,ntyp,ityp_save,  &
                                                            amass_save,iumode)

           CLOSE(UNIT=iumode, STATUS='KEEP')
        ENDIF
     ENDIF
  ENDDO
  CALL mp_bcast(displa_geo, ionode_id, intra_image_comm)
!
!  Part two: Compute the Gruneisen parameters
!
  copy_before=.FALSE.
  poly_order=5
  ALLOCATE(frequency_geo(nbnd,ngeo(1)))
  ALLOCATE(repres_geo(nbnd,ngeo(1)))
  ALLOCATE(poly_grun(poly_order,nbnd))
  ALLOCATE(frequency(nbnd,nks))
  ALLOCATE(gruneisen(nbnd,nks))
  frequency(:,:)= 0.0_DP
  gruneisen(:,:)= 0.0_DP
  IF (volume_ph==0.0_DP) THEN
     IF (ltherm_freq) THEN
        CALL evaluate_vm(temp_ph, vm, ntemp, temp, vminf_t)
     ELSEIF (ltherm_dos) THEN
        CALL evaluate_vm(temp_ph, vm, ntemp, temp, vmin_t)
     ELSE
        vm=vmin
     ENDIF
  ELSE
     vm=volume_ph
  ENDIF

  WRITE(stdout,'(/,5x,"Plotting Gruneisen parameters at volume",f17.8,&
                                                    &" (a.u.)^3")') vm
  IF (volume_ph==0.0_DP.AND.(ltherm_freq.OR.ltherm_dos)) &
            WRITE(stdout,'(5x,"Corresponding to T=",f17.8)') temp_ph

  DO n = 1,nks
     IF (is_gamma(n)) THEN
!
!    In the gamma point the Gruneisen parameters are not defined.
!    In order to have a continuous plot we take the same parameters
!    of the previous point, if this point exists and is not gamma.
!    Otherwise at the next point we copy the parameters in the present one
!
        copy_before=.FALSE.
        IF (n==1) THEN
           copy_before=.TRUE.
        ELSEIF (is_gamma(n-1)) THEN
           copy_before=.TRUE.
        ELSE
           DO ibnd=1,nbnd
              gruneisen(ibnd,n)=gruneisen(ibnd,n-1)
              frequency(ibnd,n)=frequency(ibnd,n-1)
           ENDDO
!
!  At the gamma point the first three frequencies vanishes
!
           DO ibnd=1,3
              frequency(ibnd,n)=0.0_DP
           ENDDO
        ENDIF
     ELSE
        frequency_geo(:,:)=freq_geo(:,n,:)
        IF (with_eigen) THEN
           CALL compute_freq_derivative_eigen(ngeo(1),frequency_geo,   &
                        omega_geo, displa_geo(1,1,1,n),no_ph,poly_order, &
                        poly_grun)
        ELSE 
           repres_geo(:,:)=rap_geo(:,n,:)
           CALL compute_freq_derivative(ngeo,frequency_geo,repres_geo, &
                                   omega_geo,no_ph,poly_order,poly_grun)
        ENDIF
!
!  frequencies and gruneisen parameters are calculated at the chosen
!  volume
!
        DO ibnd=1,nbnd
           DO i=1,poly_order
              frequency(ibnd,n) = frequency(ibnd,n) + &
                    poly_grun(i,ibnd) * vm**(i-1)
              gruneisen(ibnd,n) = gruneisen(ibnd,n) - &
                    poly_grun(i,ibnd) * vm**(i-1) * (i-1.0_DP)
           END DO
           IF (frequency(ibnd,n) > 0.0_DP ) THEN
              gruneisen(ibnd,n) = gruneisen(ibnd,n) / frequency(ibnd,n)
           ELSE
              gruneisen(ibnd,n) = 0.0_DP
           ENDIF
           IF (copy_before) THEN
              gruneisen(ibnd,n-1) = gruneisen(ibnd,n)
              frequency(ibnd,n-1) = frequency(ibnd,n)
           ENDIF 
        ENDDO
        copy_before=.FALSE.
!
!  At the gamma point the first three frequencies vanishes
!
        IF (n>1.AND.is_gamma(n-1)) THEN
           DO ibnd=1,3
              frequency(ibnd,n-1)=0.0_DP
           ENDDO
        ENDIF
     ENDIF
  ENDDO
!
!  Third part: writes Gruneisen parameters on file
!
   filegrun="anhar_files/"//TRIM(flgrun)
   CALL write_bands(nks, nbnd, k, gruneisen, 1.0_DP, filegrun)
!
!  writes frequencies at the chosen volume on file
!
   filefreq=TRIM(filegrun)//'_freq'
   CALL write_bands(nks, nbnd, k, frequency, 1.0_DP, filefreq)

   DEALLOCATE ( is_gamma )
   DEALLOCATE ( high_symmetry )
   DEALLOCATE ( k_rap )
   DEALLOCATE ( k ) 
   DEALLOCATE ( rap_geo )
   DEALLOCATE ( freq_geo )
   DEALLOCATE ( gcodek )
   DEALLOCATE ( gcodek_ext )
   DEALLOCATE ( aux_ind )
   DEALLOCATE ( ptypek )
   DEALLOCATE ( lprojk )
   DEALLOCATE ( same_next )
   DEALLOCATE ( gaugek )
   DEALLOCATE ( poly_grun )
   DEALLOCATE ( frequency_geo )
   DEALLOCATE ( repres_geo )
   DEALLOCATE ( frequency )
   DEALLOCATE ( gruneisen )

   IF (with_eigen) DEALLOCATE ( displa_geo )

   RETURN
END SUBROUTINE write_gruneisen_band


SUBROUTINE evaluate_vm(temp_ph, vm, ntemp, temp, vminf_t)

USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), vminf_t(ntemp)
REAL(DP), INTENT(INOUT) :: temp_ph, vm

INTEGER :: itemp0, itemp1, itemp

itemp0=1
DO itemp=1,ntemp
   IF (temp(itemp) < temp_ph) itemp0=itemp
ENDDO

IF (itemp0 == ntemp) THEN
   WRITE(stdout,'(5x,"temp_ph too large setting to",f15.8 )') temp(ntemp-1)
   temp_ph=temp(ntemp-1)
   vm=vminf_t(ntemp-1)
   RETURN
ENDIF

IF (itemp0 == 1) THEN
   WRITE(stdout,'(5x,"temp_ph too small setting to",f15.8 )') temp(2)
   temp_ph=temp(2)
   vm=vminf_t(2)
   RETURN
ENDIF

itemp1=itemp0+1

vm = vminf_t(itemp0) + (temp_ph - temp(itemp0)) *          &
                       (vminf_t(itemp1)-vminf_t(itemp0)) / &
                       (temp(itemp1)-temp(itemp0))

RETURN
END SUBROUTINE evaluate_vm
