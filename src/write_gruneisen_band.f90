!
! Copyright (C) 2015-2018 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE write_gruneisen_band(file_disp, file_vec)
!-------------------------------------------------------------------------
  !
  ! reads data files produced by write_ph_dispersions for ngeo geometries, 
  ! interpolates them with a polynomial and computes and writes a file with 
  ! the mode gruneisen parameters: minus the derivatives of the logarithm of 
  ! the phonon frequencies with respect to the logarithm of the volume.
  ! This is calculated at the volume given in input or at the volume
  ! that corresponds to the temperature given in input. This routine is 
  ! used when lmurn=.TRUE..
  ! 
  USE kinds,          ONLY : DP
  USE ions_base,      ONLY : nat, ntyp => nsp, amass
  USE data_files,     ONLY : flgrun
  USE thermo_mod,     ONLY : ngeo, omega_geo, no_ph
  USE anharmonic,     ONLY : vmin_t
  USE ph_freq_anharmonic, ONLY : vminf_t
  USE grun_anharmonic, ONLY : poly_degree_grun
  USE control_grun,   ONLY : temp_ph, volume_ph
  USE initial_conf,   ONLY : ityp_save, ibrav_save
  USE control_mur,    ONLY : vmin
  USE control_thermo, ONLY : ltherm_dos, ltherm_freq
  USE temperature,    ONLY : temp, ntemp
  USE freq_interpolate, ONLY : interp_freq_eigen 
  USE polyfit_mod,    ONLY : compute_poly, compute_poly_deriv
  USE io_bands,       ONLY : read_bands, read_parameters, write_bands
  USE mp,             ONLY : mp_bcast
  USE io_global,      ONLY : stdout, ionode, ionode_id
  USE mp_images,      ONLY : intra_image_comm, root_image, my_image_id

  IMPLICIT NONE

  CHARACTER(LEN=256), INTENT(IN) :: file_disp, file_vec

  REAL(DP), ALLOCATABLE :: freq_geo(:,:,:), k(:,:), frequency_geo(:,:), &
                           omega_data(:), poly_grun(:,:), frequency(:,:), &
                           gruneisen(:,:)
  COMPLEX(DP), ALLOCATABLE :: displa_geo(:,:,:,:), displa(:,:,:)
  INTEGER :: nks, nmodes, cgeo_eff, central_geo, imode, ios, n, igeo, ndata, &
             iumode
  INTEGER :: find_free_unit
  REAL(DP) :: vm, f, g
  LOGICAL, ALLOCATABLE :: is_gamma(:)
  LOGICAL :: copy_before, allocated_variables
  CHARACTER(LEN=256) :: filename, filedata, filegrun, filefreq
  CHARACTER(LEN=6), EXTERNAL :: int_to_char

  IF ( my_image_id /= root_image ) RETURN
  IF (flgrun == ' ') RETURN
!
!  Part one: read the frequencies and the mode eigenvectors.
!
  IF (ionode) iumode=find_free_unit()
  allocated_variables=.FALSE.
  DO igeo = 1, ngeo(1)
     IF (no_ph(igeo)) CYCLE
     filedata = "phdisp_files/"//TRIM(file_disp)//'.g'//TRIM(int_to_char(igeo))
     CALL read_parameters(nks, nmodes, filedata)
     IF (nks <= 0 .OR. nmodes <= 0) THEN
        CALL errore('write_gruneisen_band','reading plot namelist',ABS(ios))
     ELSE
        WRITE(stdout, '(5x,"Reading ",i4," dispersions at ",i6," k-points for&
                       & geometry",i4)') nmodes, nks, igeo
     ENDIF
     IF (.NOT.allocated_variables) THEN
        ALLOCATE (freq_geo(nmodes,nks,ngeo(1)))
        ALLOCATE (displa_geo(nmodes,nmodes,ngeo(1),nks))
        ALLOCATE (k(3,nks)) 
        allocated_variables=.TRUE.
     ENDIF
     CALL read_bands(nks, nmodes, k, freq_geo(1,1,igeo), filedata)

     filename="phdisp_files/"//TRIM(file_vec)//".g"//TRIM(int_to_char(igeo))
     IF (ionode) OPEN(UNIT=iumode, FILE=TRIM(filename), FORM='formatted', &
                STATUS='old', ERR=210, IOSTAT=ios)
210  CALL mp_bcast(ios, ionode_id, intra_image_comm)
     CALL errore('write_gruneisen_band','modes are needed',ABS(ios))
     IF (ionode) THEN
        CALL readmodes(nat,nks,k,displa_geo,ngeo(1),igeo,ntyp,ityp_save,  &
                                                         amass,-1,iumode)
        CLOSE(UNIT=iumode, STATUS='KEEP')
     ENDIF
  ENDDO
  CALL mp_bcast(displa_geo, ionode_id, intra_image_comm)
  ALLOCATE (is_gamma(nks))
  DO n=1,nks
     is_gamma(n) = (( k(1,n)**2 + k(2,n)**2 + k(3,n)**2) < 1.d-12)
  ENDDO
!
!  Part two: Compute the Gruneisen parameters
!
!  find how many data have been really calculated 
!
  ndata=0
  DO igeo=1, ngeo(1)
     IF (.NOT.no_ph(igeo)) ndata=ndata+1
  ENDDO
!
!   find the central geometry of this set of data and save the volume of
!   each computed geometry
!
  ALLOCATE(omega_data(ndata))
  CALL find_central_geo(ngeo,no_ph,central_geo)
  ndata=0
  DO igeo=1, ngeo(1)
     IF (no_ph(igeo)) CYCLE
     ndata=ndata+1
     IF (central_geo==igeo) cgeo_eff=ndata
     omega_data(ndata)=omega_geo(igeo)
  ENDDO

!
!  Compute the volume at which the Gruneisen parameters and the frequencies
!  are interpolated
!
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
            WRITE(stdout,'(5x,"Corresponding to T=",f17.4," K")') temp_ph
!
!  Allocate space for interpolating the frequencies
!
  ALLOCATE(frequency_geo(nmodes,ndata))
  ALLOCATE(displa(nmodes,nmodes,ndata))
  ALLOCATE(poly_grun(poly_degree_grun+1,nmodes))
  ALLOCATE(frequency(nmodes,nks))
  ALLOCATE(gruneisen(nmodes,nks))

  copy_before=.FALSE.
  frequency(:,:)= 0.0_DP
  gruneisen(:,:)= 0.0_DP
  DO n = 1,nks
     IF (is_gamma(n)) THEN
!
!    At the gamma point the Gruneisen parameters are not defined.
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
           gruneisen(1:nmodes,n)=gruneisen(1:nmodes,n-1)
           frequency(1:nmodes,n)=frequency(1:nmodes,n-1)
!
!  At the gamma point the first three frequencies vanishes
!
           DO imode=1,3
              frequency(imode,n)=0.0_DP
           ENDDO
        ENDIF
     ELSE
!
!  interpolates the frequencies with a polynomial and gives the coefficients
!  poly_grun
!
        ndata=0
        DO igeo=1, ngeo(1)
           IF (no_ph(igeo)) CYCLE
           ndata=ndata+1
           frequency_geo(1:nmodes,ndata)=freq_geo(1:nmodes,n,igeo)
           displa(1:nmodes,1:nmodes,ndata) = displa_geo(1:nmodes,1:nmodes,igeo,n)
        ENDDO
        CALL interp_freq_eigen(ndata, frequency_geo, omega_data, &
                          cgeo_eff, displa, poly_degree_grun, poly_grun)
!
!  frequencies and gruneisen parameters are calculated at the chosen
!  volume using the intepolating polynomial
!
        DO imode=1,nmodes
           CALL compute_poly(vm, poly_degree_grun, poly_grun(:,imode),f)
           CALL compute_poly_deriv(vm, poly_degree_grun, poly_grun(:,imode),g)
           frequency(imode,n)=f
!
!     g here is V d w / d V. We change sign and divide by the frequency w 
!     to get the gruneisen parameter.
!
           IF (f > 0.0_DP ) THEN
              gruneisen(imode,n) = - vm * g / f
           ELSE
              gruneisen(imode,n) = 0.0_DP
           ENDIF
           IF (copy_before) THEN
              gruneisen(imode,n-1) = gruneisen(imode,n)
              frequency(imode,n-1) = frequency(imode,n)
           ENDIF 
        ENDDO
        copy_before=.FALSE.
!
!  At the gamma point the first three frequencies vanishes
!
        IF (n>1.AND.is_gamma(n-1)) THEN
           DO imode=1,3
              frequency(imode,n-1)=0.0_DP
           ENDDO
        ENDIF
     ENDIF
  ENDDO
!
!  Third part: writes Gruneisen parameters on file
!
   filegrun="anhar_files/"//TRIM(flgrun)
   CALL write_bands(nks, nmodes, k, gruneisen, 1.0_DP, filegrun)
!
!  writes frequencies at the chosen volume on file
!
   filefreq=TRIM(filegrun)//'_freq'
   CALL write_bands(nks, nmodes, k, frequency, 1.0_DP, filefreq)

   DEALLOCATE( gruneisen )
   DEALLOCATE( frequency )
   DEALLOCATE( poly_grun )
   DEALLOCATE( displa )
   DEALLOCATE( frequency_geo )
   DEALLOCATE( omega_data )
   DEALLOCATE( is_gamma )
   DEALLOCATE( k ) 
   DEALLOCATE( displa_geo )
   DEALLOCATE( freq_geo )

   RETURN
END SUBROUTINE write_gruneisen_band

