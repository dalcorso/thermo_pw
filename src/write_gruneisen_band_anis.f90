!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE write_gruneisen_band_anis(file_disp, file_vec)
!------------------------------------------------------------------------
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
  USE ions_base,      ONLY : nat, ntyp => nsp, amass
  USE cell_base,      ONLY : celldm
  USE data_files,     ONLY : flgrun
  USE control_paths,  ONLY : nqaux, disp_q, disp_nqs
  USE equilibrium_conf, ONLY : celldm0
  USE thermo_mod,     ONLY : celldm_geo, no_ph, in_degree, reduced_grid, red_central_geo
  USE anharmonic,     ONLY : celldm_t
  USE ph_freq_anharmonic,  ONLY : celldmf_t
  USE grun_anharmonic, ONLY : poly_degree_grun
  USE control_grun,   ONLY : temp_ph, celldm_ph
  USE initial_conf,   ONLY : ibrav_save, ityp_save
  USE control_thermo, ONLY : ltherm_dos, ltherm_freq, set_internal_path
  USE control_quartic_energy, ONLY : lsolve
  USE freq_interpolate, ONLY : interp_freq_anis_eigen, interp_freq_eigen
  USE polyfit_mod,    ONLY : compute_poly, compute_poly_deriv
  USE temperature,    ONLY : temp, ntemp
  USE quadratic_surfaces, ONLY : evaluate_fit_quadratic, &
                                 evaluate_quadratic_grad
  USE lattices,       ONLY : compress_celldm, crystal_parameters
  USE polynomial,     ONLY : poly2, init_poly, clean_poly
  USE io_bands,       ONLY : read_bands, read_parameters, write_bands
  USE mp,             ONLY : mp_bcast
  USE io_global,      ONLY : stdout, ionode, ionode_id
  USE mp_images,      ONLY : intra_image_comm, root_image, my_image_id

  IMPLICIT NONE

  CHARACTER(LEN=256), INTENT(IN) :: file_disp, file_vec

  REAL(DP), ALLOCATABLE :: freq_geo(:,:,:), k(:,:), frequency_geo(:,:), &
                           x_data(:,:), poly_grun(:,:), frequency(:,:), &
                           gruneisen(:,:,:), grad(:), x(:), xd(:)
  TYPE(poly2), ALLOCATABLE :: p_grun_p2(:)
  COMPLEX(DP), ALLOCATABLE :: displa_geo(:,:,:,:), displa(:,:,:)
  INTEGER :: nks, nmodes, cgeo_eff, central_geo, imode, ios, n, igeo, &
             nwork, idata, ndata, iumode, nvar, icrys
  INTEGER :: find_free_unit, compute_nwork
  REAL(DP) :: cm(6), f, g
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
     CALL read_parameters(nks, nmodes, filedata)
     !
     IF (nks <= 0 .or. nmodes <= 0) THEN
        CALL errore('write_gruneisen_band_anis','reading plot namelist',&
                                                                 ABS(ios))
     ELSE
        WRITE(stdout, '(5x, "Reading ",i4," dispersions at ",i6," k-points for&
                       & geometry",i4)') nmodes, nks, igeo
     ENDIF

     IF (.NOT.allocated_variables) THEN
        ALLOCATE (freq_geo(nmodes,nks,nwork))
        ALLOCATE (displa_geo(nmodes,nmodes,nwork,nks))
        ALLOCATE (k(3,nks))
        allocated_variables=.TRUE.
     ENDIF
     CALL read_bands(nks, nmodes, k, freq_geo(1,1,igeo), filedata)

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
                                       amass, -1, iumode)
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
  nvar=crystal_parameters(ibrav_save)
!
!  find the central geometry in this list of phonons. Prepare also the
!  x data for the computed geometries
!
  CALL find_central_geo(nwork,no_ph,central_geo)
  ALLOCATE(x_data(nvar,nwork))
  ndata=0
  DO idata=1, nwork
     IF (no_ph(idata)) CYCLE
     ndata=ndata+1
     IF (central_geo==idata) cgeo_eff=ndata
     IF (reduced_grid) THEN
        CALL compress_celldm(celldm_geo(1,idata),x_data(1,idata),nvar, &
                                                                ibrav_save)
     ELSE
        CALL compress_celldm(celldm_geo(1,idata),x_data(1,ndata),nvar, &
                                                                ibrav_save)
     ENDIF
  ENDDO 
!
!  find the celldm at which we compute the Gruneisen parameters and compress it
!
  ALLOCATE(x(nvar))
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
  CALL compress_celldm(cm,x,nvar,ibrav_save)
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
            WRITE(stdout,'(5x,"Corresponding to T=",f17.4," K")') temp_ph
!
! Now allocate space to save the frequencies to interpolate and the
! coeffieints of the polynomial. Allocate also space to save the
! interpolated frequencies and the gruneisen parameters
!
  ALLOCATE(frequency_geo(nmodes,ndata))
  ALLOCATE(displa(nmodes,nmodes,ndata))
  ALLOCATE(p_grun_p2(nmodes))
  ALLOCATE(frequency(nmodes,nks))
  ALLOCATE(gruneisen(nmodes,nks,nvar))
  ALLOCATE(xd(ndata))
  ALLOCATE(grad(nvar))
  DO imode=1,3*nat
     CALL init_poly(nvar,p_grun_p2(imode))
  ENDDO
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
           DO imode=1,nmodes
              gruneisen(imode,n,:)=gruneisen(imode,n-1,:)
              frequency(imode,n)=frequency(imode,n-1)
           END DO
!
!  At the gamma point the first three frequencies vanish
!
           DO imode=1,3
              frequency(imode,n)=0.0_DP
           END DO
        ENDIF
     ELSE
        IF (reduced_grid) THEN
!
!  with reduced grid each degree of freedom is calculated separately
!
           DO icrys=1, nvar
              ndata=0
              cgeo_eff=0
              DO idata=1,nwork
                 IF (in_degree(idata)==icrys.OR.idata==red_central_geo) THEN
                    ndata=ndata+1
                    xd(ndata)=x_data(icrys,idata)
                    frequency_geo(1:nmodes,ndata)=freq_geo(1:nmodes,n,idata)
                    displa(1:nmodes,1:nmodes,ndata)=&
                                   displa_geo(1:nmodes,1:nmodes,idata,n)
                    IF (idata==red_central_geo) cgeo_eff=ndata
                 ENDIF
              ENDDO
              IF (cgeo_eff==0) cgeo_eff=ndata/2
              CALL interp_freq_eigen(ndata, frequency_geo, xd, cgeo_eff, &
                                          displa, poly_degree_grun, poly_grun)
              DO imode=1,nmodes
                 CALL compute_poly(x(icrys), poly_degree_grun, &
                                                         poly_grun(:,imode),f)
!
!  this function gives the derivative with respect to x(i) multiplied by x(i)
!
                 CALL compute_poly_deriv(x(icrys), poly_degree_grun, &
                                                   poly_grun(:,imode),g)
                 frequency(imode,n) = f
                 gruneisen(imode,n,icrys) = - x(icrys) * g
                 IF (f > 0.0_DP ) THEN
                    gruneisen(imode,n,icrys) = gruneisen(imode,n,icrys)/f
                 ELSE
                    gruneisen(imode,n,icrys) = 0.0_DP
                 ENDIF
              ENDDO
           ENDDO
        ELSE
!
!  In this case the frequencies are fitted with a multidimensional polynomial
!
           ndata=0
           DO idata=1,nwork
              IF (no_ph(idata)) CYCLE
              ndata=ndata+1
              frequency_geo(1:nmodes,ndata)=freq_geo(1:nmodes,n,idata)
              displa(1:nmodes,1:nmodes,ndata)=displa_geo(1:nmodes,1:nmodes,idata,n)
           ENDDO
           CALL interp_freq_anis_eigen(ndata, frequency_geo, lsolve, &
                     x_data, cgeo_eff, displa, nvar, p_grun_p2 )
!
!  frequencies and gruneisen parameters are calculated at the chosen volume
!
           DO imode=1, nmodes
              CALL evaluate_fit_quadratic(nvar,x,f,p_grun_p2(imode))
              CALL evaluate_quadratic_grad(nvar,x,grad, p_grun_p2(imode))
              frequency(imode,n) = f 
              gruneisen(imode,n,:) = -grad(:)
              IF (frequency(imode,n) > 0.0_DP ) THEN
                 DO icrys=1,nvar
                    gruneisen(imode,n,icrys) = x(icrys) * gruneisen(imode,n,icrys) /  &
                                                 frequency(imode,n)
                 END DO
              ELSE
                 gruneisen(imode,n,:) = 0.0_DP
              ENDIF
           ENDDO
        ENDIF
        IF (copy_before) THEN
           gruneisen(1:nmodes,n-1,1:nvar) = gruneisen(1:nmodes,n,1:nvar)
           frequency(1:nmodes,n-1) = frequency(1:nmodes,n)
        ENDIF
        copy_before=.FALSE.
!
!  At the gamma point the first three frequencies vanish
!
        IF (n>1.AND.is_gamma(n-1)) THEN
           DO imode=1,3
              frequency(imode,n-1)=0.0_DP
           ENDDO
        ENDIF
     ENDIF
  ENDDO
!
!  Writes Gruneisen parameters on file
!
  DO icrys=1, nvar
     filegrun='anhar_files/'//TRIM(flgrun)//'_'//TRIM(INT_TO_CHAR(icrys))
     CALL write_bands(nks, nmodes, disp_q, gruneisen(1,1,icrys), 1.0_DP, filegrun)
  END DO
!
!  writes frequencies at the chosen geometry on file
!
  filefreq='anhar_files/'//TRIM(flgrun)//'_freq'
  CALL write_bands(nks, nmodes, disp_q, frequency, 1.0_DP, filefreq)


  DO imode=1,3*nat
     CALL clean_poly(p_grun_p2(imode))
  ENDDO

  DEALLOCATE( grad )
  DEALLOCATE( gruneisen )
  DEALLOCATE( frequency )
  DEALLOCATE( p_grun_p2 )
  DEALLOCATE( displa )
  DEALLOCATE( frequency_geo )
  DEALLOCATE( xd )
  DEALLOCATE( x )
  DEALLOCATE( x_data )
  DEALLOCATE( is_gamma )
  DEALLOCATE( k ) 
  DEALLOCATE( displa_geo )
  DEALLOCATE( freq_geo )

RETURN
END SUBROUTINE write_gruneisen_band_anis
