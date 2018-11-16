!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_gruneisen_band_anis(file_disp, file_vec)
  !
  ! reads data files produced by write_ph_dispersions for ngeo geometries, 
  ! interpolates them with a quadratic polynomial, computes an writes 
  ! several files with the Gruneisen parameters: the derivatives of the 
  ! logarithm of the phonon frequencies with respect to the logarithm of
  ! the relevant celldm parameters.
  ! This is calculated at the celldm given in input or at the celldm
  ! that corresponds to the temperature given in input.
  ! 
  USE kinds,          ONLY : DP
  USE ions_base,      ONLY : nat, ntyp => nsp
  USE cell_base,      ONLY : celldm
  USE data_files,     ONLY : flgrun
  USE control_paths,  ONLY : nqaux, disp_q, disp_nqs
  USE equilibrium_conf, ONLY : celldm0
  USE thermo_mod,     ONLY : celldm_geo, no_ph
  USE anharmonic,     ONLY : celldm_t
  USE ph_freq_anharmonic,  ONLY : celldmf_t
  USE control_grun,   ONLY : temp_ph, celldm_ph
  USE initial_conf,   ONLY : ibrav_save, amass_save, ityp_save
  USE control_thermo, ONLY : ltherm_dos, ltherm_freq, set_internal_path
  USE freq_interpolate, ONLY : interp_freq_anis_eigen
  USE temperature,    ONLY : temp, ntemp
  USE quadratic_surfaces, ONLY : evaluate_fit_quadratic, &
                                 evaluate_fit_grad_quadratic, quadratic_var
  USE lattices,       ONLY : compress_celldm, crystal_parameters
  USE io_bands,       ONLY : read_bands, read_parameters, write_bands
  USE mp,             ONLY : mp_bcast
  USE io_global,      ONLY : stdout, ionode, ionode_id
  USE mp_images,      ONLY : intra_image_comm, root_image, my_image_id

  IMPLICIT NONE

  CHARACTER(LEN=256), INTENT(IN) :: file_disp, file_vec

  REAL(DP), ALLOCATABLE :: freq_geo(:,:,:), k(:,:), frequency_geo(:,:), &
                           x_data(:,:), poly_grun(:,:), frequency(:,:), &
                           gruneisen(:,:,:), grad(:), x(:)
  COMPLEX(DP), ALLOCATABLE :: displa_geo(:,:,:,:), displa(:,:,:)
  INTEGER :: nks, nbnd, cgeo_eff, central_geo, ibnd, ios, i, n, igeo, &
             nwork, idata, ndata, iumode, nvar, degree, icrys
  INTEGER :: find_free_unit, compute_nwork
  REAL(DP) :: cm(6), f
  LOGICAL, ALLOCATABLE :: is_gamma(:)
  LOGICAL :: copy_before, allocated_variables
  CHARACTER(LEN=256) :: filename, filedata, filegrun, filefreq
  CHARACTER(LEN=6), EXTERNAL :: int_to_char

  IF ( my_image_id /= root_image ) RETURN
  IF (flgrun == ' ') RETURN

  IF (ionode) iumode=find_free_unit()
  nwork=compute_nwork()
!
!  Part one: Reads the frequencies and the displacements from file.
!
  allocated_variables=.FALSE.
  DO igeo = 1, nwork
     
     IF (no_ph(igeo)) CYCLE
     filedata = 'phdisp_files/'//TRIM(file_disp)//'.g'//TRIM(int_to_char(igeo))
     CALL read_parameters(nks, nbnd, filedata)
     !
     IF (nks <= 0 .or. nbnd <= 0) THEN
        CALL errore('write_gruneisen_band_anis','reading plot namelist',&
                                                                 ABS(ios))
     ELSE
        WRITE(stdout, '(5x, "Reading ",i4," dispersions at ",i6," k-points for&
                       & geometry",i4)') nbnd, nks, igeo
     ENDIF

     IF (.NOT.allocated_variables) THEN
        ALLOCATE (freq_geo(nbnd,nks,nwork))
        ALLOCATE (displa_geo(nbnd,nbnd,nwork,nks))
        ALLOCATE (k(3,nks))
        allocated_variables=.TRUE.
     ENDIF
     CALL read_bands(nks, nbnd, k, freq_geo(1,1,igeo), filedata)

     filename='phdisp_files/'//TRIM(file_vec)//".g"//TRIM(int_to_char(igeo))
     IF (ionode) OPEN(UNIT=iumode, FILE=TRIM(filename), FORM='formatted', &
                   STATUS='old', ERR=200, IOSTAT=ios)
200  CALL mp_bcast(ios, ionode_id, intra_image_comm)
     CALL errore('write_gruneisen_band_anis','modes are needed',ABS(ios))
     !
     IF (ionode) THEN
!
!  readmodes reads a file with the displacements, but writes in displa_geo
!  the normalized eigenvectors of the dynamical matrix
!
        CALL readmodes(nat,nks,k,displa_geo,nwork,igeo,ntyp,ityp_save, &
                                       amass_save, iumode)
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
  degree=crystal_parameters(ibrav_save)
  nvar=quadratic_var(degree)
!
!  find how many geometries have been really calculated
! 
  ndata=0
  DO idata=1, nwork
     IF (.NOT.no_ph(idata)) ndata=ndata+1
  ENDDO
!
!  And find the central geometry in this list of phonon. Prepare also the
!  x data for the computed geometries
!
  CALL find_central_geo(nwork,no_ph,central_geo)
  ALLOCATE(x_data(degree,ndata))
  ndata=0
  DO idata=1, nwork
     IF (.NOT.no_ph(idata)) THEN
        ndata=ndata+1
        IF (central_geo==idata) cgeo_eff=ndata
        CALL compress_celldm(celldm_geo(1,idata),x_data(1,ndata),degree, &
                                                                   ibrav_save)
     ENDIF
  ENDDO 
!
!  find the celldm at which we compute the Gruneisen parameters and compress
!  it
  ALLOCATE(x(degree))
  IF (celldm_ph(1)==0.0_DP) THEN
     IF (ltherm_freq) THEN
        CALL evaluate_celldm(temp_ph, cm, ntemp, temp, celldmf_t)
     ELSEIF(ltherm_dos) THEN
        CALL evaluate_celldm(temp_ph, cm, ntemp, temp, celldm_t)
     ELSE
        cm(:)=celldm0(:)
     ENDIF
  ELSE
     cm=celldm_ph(:)
  ENDIF
  CALL compress_celldm(cm,x,degree,ibrav_save)
!
!  calculate the BZ path that corresponds to the cm parameters
!
  celldm(:)=cm(:)
  IF (set_internal_path) CALL set_bz_path()
  IF (nqaux > 0) CALL set_paths_disp()
  IF (disp_nqs /= nks) &
     CALL errore('write_gruneisen_band_anis','Problem with the path',1)

  WRITE(stdout,'(/,5x,"Plotting Gruneisen parameters at celldm:")') 
  WRITE(stdout,'(5x,6f12.7)') cm 
  IF (celldm_ph(1)==0.0_DP.AND.(ltherm_freq.OR.ltherm_dos)) &
            WRITE(stdout,'(5x,"Corresponding to T=",f17.8)') temp_ph
!
! Now allocate space to save the frequencies to interpolate and the
! coeffieints of the polynomial. Allocate also space to save the
! interpolated frequencies and the gruneisen parameters
!
  ALLOCATE(frequency_geo(nbnd,ndata))
  ALLOCATE(displa(nbnd,nbnd,ndata))
  ALLOCATE(poly_grun(nvar,nbnd))
  ALLOCATE(frequency(nbnd,nks))
  ALLOCATE(gruneisen(nbnd,nks,degree))
  ALLOCATE(grad(degree))
!
!  then intepolates and computes the derivatives of the polynomial.
!  At the gamma point the gruneisen parameters are not defined, so
!  we copy those of the previous or the next point
!
  copy_before=.FALSE.
  frequency(:,:)= 0.0_DP
  gruneisen(:,:,:)= 0.0_DP
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
        ELSEIF(is_gamma(n-1)) THEN
           copy_before=.TRUE.
        ELSE
           DO ibnd=1,nbnd
              gruneisen(ibnd,n,:)=gruneisen(ibnd,n-1,:)
              frequency(ibnd,n)=frequency(ibnd,n-1)
           END DO
!
!  At the gamma point the first three frequencies vanish
!
           DO ibnd=1,3
              frequency(ibnd,n)=0.0_DP
           END DO
        ENDIF
     ELSE
        ndata=0
        DO idata=1,nwork
           IF (.NOT.no_ph(idata)) THEN
              ndata=ndata+1
              frequency_geo(1:nbnd,ndata)=freq_geo(1:nbnd,n,idata)
              displa(1:nbnd,1:nbnd,ndata)=displa_geo(1:nbnd,1:nbnd,idata,n)
           ENDIF
        ENDDO
        CALL interp_freq_anis_eigen(nwork, frequency_geo, x_data, cgeo_eff, &
                                    displa, degree, nvar, poly_grun)
!
!  frequencies and gruneisen parameters are calculated at the chosen volume
!
        DO ibnd=1,nbnd
           CALL evaluate_fit_quadratic(degree,nvar,x,f,poly_grun(1,ibnd))
           CALL evaluate_fit_grad_quadratic(degree,nvar,x,grad, &
                                                        poly_grun(1,ibnd))
           frequency(ibnd,n) = f 
           gruneisen(ibnd,n,:) = -grad(:)
           IF (frequency(ibnd,n) > 0.0_DP ) THEN
              DO i=1,degree
                 gruneisen(ibnd,n,i) = x(i) * gruneisen(ibnd,n,i) /  &
                                              frequency(ibnd,n)
              END DO
           ELSE
              gruneisen(ibnd,n,:) = 0.0_DP
           ENDIF
           IF (copy_before) THEN
              gruneisen(ibnd,n-1,:) = gruneisen(ibnd,n,:)
              frequency(ibnd,n-1) = frequency(ibnd,n)
           ENDIF 
        ENDDO
        copy_before=.FALSE.
!
!  At the gamma point the first three frequencies vanish
!
        IF (n>1.AND.is_gamma(n-1)) THEN
           DO ibnd=1,3
              frequency(ibnd,n-1)=0.0_DP
           ENDDO
        ENDIF
     ENDIF
  ENDDO
!
!  Writes Gruneisen parameters on file
!
  DO icrys=1, degree
     filegrun='anhar_files/'//TRIM(flgrun)//'_'//TRIM(INT_TO_CHAR(icrys))
     CALL write_bands(nks, nbnd, disp_q, gruneisen(1,1,icrys), 1.0_DP, filegrun)
  END DO
!
!  writes frequencies at the chosen geometry on file
!
  filefreq='anhar_files/'//TRIM(flgrun)//'_freq'
  CALL write_bands(nks, nbnd, disp_q, frequency, 1.0_DP, filefreq)

  DEALLOCATE( grad )
  DEALLOCATE( gruneisen )
  DEALLOCATE( frequency )
  DEALLOCATE( poly_grun )
  DEALLOCATE( displa )
  DEALLOCATE( frequency_geo )
  DEALLOCATE( x )
  DEALLOCATE( x_data )
  DEALLOCATE( is_gamma )
  DEALLOCATE( k ) 
  DEALLOCATE( displa_geo )
  DEALLOCATE( freq_geo )

RETURN
END SUBROUTINE write_gruneisen_band_anis

SUBROUTINE evaluate_celldm(temp_ph, cm, ntemp, temp, celldmf_t)

USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), celldmf_t(6,ntemp)
REAL(DP), INTENT(INOUT) :: temp_ph, cm(6)

INTEGER :: itemp0, itemp1, itemp

itemp0=1
DO itemp=1,ntemp
   IF (temp(itemp) < temp_ph) itemp0=itemp
ENDDO

IF (itemp0 == ntemp) THEN
   WRITE(stdout,'(5x,"temp_ph too large setting to",f15.8 )') temp(ntemp-1)
   temp_ph=temp(ntemp-1)
   cm(:)=celldmf_t(:,ntemp-1)
   RETURN
ENDIF

IF (itemp0 == 1) THEN
   WRITE(stdout,'(5x,"temp_ph too small setting to",f15.8 )') temp(2)
   temp_ph=temp(2)
   cm(:)=celldmf_t(:,2)
   RETURN
ENDIF

itemp1=itemp0+1

cm(:) = celldmf_t(:,itemp0) + (temp_ph - temp(itemp0)) *           &
                       (celldmf_t(:,itemp1)-celldmf_t(:,itemp0)) / &
                       (temp(itemp1)-temp(itemp0))

RETURN
END SUBROUTINE evaluate_celldm

