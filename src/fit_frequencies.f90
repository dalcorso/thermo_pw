!
! Copyright (C) 2014-2018 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE fit_frequencies()
!--------------------------------------------------------------------------
  !
  ! Fits with a polynomial the phonon frequencies calculated in a uniform 
  ! mesh of ngeo geometries.
  ! 
  USE kinds,                  ONLY : DP
  USE thermo_mod,             ONLY : ngeo, omega_geo, no_ph
  USE ions_base,              ONLY : nat
  USE ph_freq_thermodynamics, ONLY : ph_freq_save
  USE grun_anharmonic,        ONLY : poly_grun, poly_order
  USE control_thermo,         ONLY : with_eigen
  USE freq_interpolate,       ONLY : interp_freq_eigen, interp_freq

  IMPLICIT NONE

  REAL(DP),    ALLOCATABLE :: freq_geo(:,:), omega_data(:)
  COMPLEX(DP), ALLOCATABLE :: displa_geo(:,:,:)
  INTEGER :: iq_eff, n, igeo, startq, lastq, cgeo_eff, ndata, central_geo
!
!  Finds how many real phonon data are available 
!
  ndata=0
  DO igeo=1,ngeo(1)
     IF (.NOT. no_ph(igeo)) ndata=ndata+1
  ENDDO
  ALLOCATE(omega_data(ndata))
!
!  finds the central geometry between those calculated
!  and collects the volumes omega for these geometries
!
  CALL find_central_geo(ngeo, no_ph, central_geo)
  ndata=0
  DO igeo=1,ngeo(1)
     IF (no_ph(igeo)) CYCLE
     ndata=ndata+1
     IF (igeo==central_geo) cgeo_eff=ndata
     omega_data(ndata)=omega_geo(igeo)
  ENDDO
!
!  divides the q vectors among all processors. Each processor computes
!  the polynomials for the wave vectors q that belong to it
!
  n=ph_freq_save(central_geo)%nq
  startq=ph_freq_save(central_geo)%startq
  lastq=ph_freq_save(central_geo)%lastq
!
!  allocates space for the fit of the frequencies with respect to the
!  volume. Note that poly_grun is not deallocated and is the output of 
!  this routine, used in the following ones.
!
  ALLOCATE(poly_grun(poly_order,3*nat,startq:lastq))
  ALLOCATE(freq_geo(3*nat,ndata))
  IF (with_eigen) ALLOCATE(displa_geo(3*nat,3*nat,ndata))
 
  poly_grun=0.0_DP
  iq_eff=0
  DO n = startq, lastq
     iq_eff=iq_eff+1
!
!  for each q point selects only the frequencies and eigenvectors
!  of the geometries that have been calculated
!
     ndata=0
     DO igeo=1,ngeo(1)
        IF (no_ph(igeo)) CYCLE
        ndata=ndata+1
        freq_geo(1:3*nat,ndata)=ph_freq_save(igeo)%nu(1:3*nat,iq_eff)
        IF (with_eigen) displa_geo(1:3*nat, 1:3*nat, ndata)= &
                           ph_freq_save(igeo)%displa(1:3*nat,1:3*nat,iq_eff)
     ENDDO
!
!  and interpolates the data
!
     IF (with_eigen) THEN
        CALL interp_freq_eigen(ndata,freq_geo,omega_data,cgeo_eff, &
                               displa_geo,poly_order,poly_grun(1,1,n))
     ELSE
        CALL interp_freq(ndata,freq_geo,omega_data, &
                                             poly_order,poly_grun(1,1,n))
     END IF
  ENDDO
!
!  deallocates space
!
  IF (with_eigen) DEALLOCATE ( displa_geo )
  DEALLOCATE ( freq_geo )
  DEALLOCATE ( omega_data )

  RETURN
END SUBROUTINE fit_frequencies
!
!--------------------------------------------------------------------------
SUBROUTINE fit_frequencies_anis()
!--------------------------------------------------------------------------
  !
  ! Computes the coefficients of a second order polynomial of
  ! the crystal parameters and fits the phonon frequencies calculated 
  ! for nwork geometries. 
  ! 
  USE kinds,                  ONLY : DP
  USE thermo_mod,             ONLY : celldm_geo, no_ph
  USE ions_base,              ONLY : nat
  USE cell_base,              ONLY : ibrav
  USE ph_freq_thermodynamics, ONLY : ph_freq_save
  USE grun_anharmonic,        ONLY : poly_grun
  USE quadratic_surfaces,     ONLY : quadratic_ncoeff
  USE control_thermo,         ONLY : with_eigen
  USE freq_interpolate,       ONLY : interp_freq_anis, interp_freq_anis_eigen
  USE lattices,               ONLY : compress_celldm, crystal_parameters

  IMPLICIT NONE

  REAL(DP),    ALLOCATABLE :: freq_geo(:,:), x(:,:)
  COMPLEX(DP), ALLOCATABLE :: displa_geo(:,:,:)
  INTEGER :: n, igeo, nvar, ncoeff, nwork, startq, lastq, iq_eff, ndata, &
             cgeo_eff, central_geo
  INTEGER :: compute_nwork
!
!  Finds how many data have been calculated
!
  nwork=compute_nwork()
  nvar=crystal_parameters(ibrav)
  ncoeff=quadratic_ncoeff(nvar)
  ndata=0
  DO igeo=1,nwork
     IF (.NOT.no_ph(igeo)) ndata=ndata+1
  ENDDO
  ALLOCATE(x(nvar,ndata))
!
!  finds the central geometry among those calculated and extracts 
!  from celldm_geo the relevant crystal parameters
!
  CALL find_central_geo(nwork,no_ph,central_geo)
  ndata=0
  DO igeo=1,nwork
     IF (no_ph(igeo)) CYCLE
     ndata=ndata+1
     IF (central_geo==igeo) cgeo_eff=ndata
     CALL compress_celldm(celldm_geo(1,igeo),x(1,ndata),nvar,ibrav)
  ENDDO
!
!  divides the q vectors among all processors. Each processor computes
!  the polynomials for the wave vectors q that belong to it
!
  n=ph_freq_save(central_geo)%nq
  startq=ph_freq_save(central_geo)%startq
  lastq=ph_freq_save(central_geo)%lastq
!
!  allocates space for the fit of the frequencies with respect to the
!  crystal parameters. Note that poly_grun is not deallocated and 
!  is the output of this routine, used in the following ones.
!
  ALLOCATE(poly_grun(ncoeff,3*nat,startq:lastq))
  ALLOCATE(freq_geo(3*nat,ndata))
  IF (with_eigen) ALLOCATE(displa_geo(3*nat,3*nat,ndata))
!
  iq_eff=0
  DO n = startq, lastq
     iq_eff=iq_eff+1
!
!  for each q point selects only the frequencies and eigenvectors
!  of the geometries that have been calculated
!
     ndata=0
     DO igeo=1,nwork
        IF (no_ph(igeo)) CYCLE
        ndata=ndata+1
        freq_geo(1:3*nat,ndata)=ph_freq_save(igeo)%nu(1:3*nat,iq_eff)
        IF (with_eigen) displa_geo(1:3*nat, 1:3*nat, ndata)= &
                           ph_freq_save(igeo)%displa(1:3*nat,1:3*nat,iq_eff)
     ENDDO
!
!   and interpolates the data
!
     IF (with_eigen) THEN
        CALL interp_freq_anis_eigen(ndata,freq_geo,x,cgeo_eff,displa_geo, &
                           nvar,ncoeff,poly_grun(1,1,n))
     ELSE
        CALL interp_freq_anis(ndata,freq_geo,x,nvar,ncoeff,poly_grun(1,1,n))
     ENDIF
  ENDDO
!
!  deallocates space
!
  IF (with_eigen) DEALLOCATE (displa_geo)
  DEALLOCATE(freq_geo)
  DEALLOCATE(x)

  RETURN
END SUBROUTINE fit_frequencies_anis
!
!--------------------------------------------------------------------------
SUBROUTINE fit_frequencies_anis_reduced()
!--------------------------------------------------------------------------
  !
  ! Computes the coefficients of a polynomial of the crystal parameters 
  ! (one at a time) and fits the phonon frequencies calculated for 
  ! nwork geometries. 
  ! 
  USE kinds,                  ONLY : DP
  USE thermo_mod,             ONLY : celldm_geo, no_ph, in_degree, &
                                     red_central_geo
  USE ions_base,              ONLY : nat
  USE cell_base,              ONLY : ibrav
  USE ph_freq_thermodynamics, ONLY : ph_freq_save
  USE grun_anharmonic,        ONLY : poly_grun_red, poly_order
  USE control_thermo,         ONLY : with_eigen
  USE freq_interpolate,       ONLY : interp_freq, interp_freq_eigen
  USE lattices,               ONLY : compress_celldm, crystal_parameters

  IMPLICIT NONE

  REAL(DP),    ALLOCATABLE :: freq_geo(:,:), x(:,:), xd(:)
  COMPLEX(DP), ALLOCATABLE :: displa_geo(:,:,:)
  INTEGER :: n, igeo, nvar, nwork, startq, lastq, iq_eff, ndata, &
             cgeo_eff, i
  INTEGER :: compute_nwork
!
!  Finds how many data have been calculated
!
  nwork=compute_nwork()
  nvar=crystal_parameters(ibrav)
  ndata=0
  DO igeo=1,nwork
     IF (.NOT.no_ph(igeo)) ndata=ndata+1
  ENDDO
  ALLOCATE(x(nvar,nwork))
!
!  finds the central geometry among those calculated and extracts 
!  from celldm_geo the relevant crystal parameters
!
  DO igeo=1,nwork
     IF (no_ph(igeo)) CYCLE
     CALL compress_celldm(celldm_geo(1,igeo),x(1,igeo),nvar,ibrav)
  ENDDO
!
!  divides the q vectors among all processors. Each processor computes
!  the polynomials for the wave vectors q that belong to it
!
  n=ph_freq_save(red_central_geo)%nq
  startq=ph_freq_save(red_central_geo)%startq
  lastq=ph_freq_save(red_central_geo)%lastq
!
!  allocates space for the fit of the frequencies with respect to the
!  crystal parameters. Note that poly_grun is not deallocated and 
!  is the output of this routine, used in the following ones.
!
  ALLOCATE(poly_grun_red(poly_order,3*nat,nvar,startq:lastq))
  ALLOCATE(xd(ndata))
  ALLOCATE(freq_geo(3*nat,ndata))
  IF (with_eigen) ALLOCATE(displa_geo(3*nat,3*nat,ndata))
!
  iq_eff=0
  DO n = startq, lastq
     iq_eff=iq_eff+1
!
!  for each q point selects only the frequencies and eigenvectors
!  of the geometries that have been calculated
!
     DO i=1, nvar
        ndata=0
        cgeo_eff=0
        DO igeo=1, nwork
           IF (in_degree(igeo)==i.OR.igeo==red_central_geo) THEN
              ndata=ndata+1
              freq_geo(1:3*nat,ndata)=ph_freq_save(igeo)%nu(1:3*nat,iq_eff)
              IF (with_eigen) displa_geo(1:3*nat, 1:3*nat, ndata)= &
                           ph_freq_save(igeo)%displa(1:3*nat,1:3*nat,iq_eff)

              xd(ndata)=x(i,igeo)
              IF (igeo==red_central_geo) cgeo_eff=ndata
           ENDIF
        ENDDO
        IF (cgeo_eff==0) cgeo_eff=ndata/2
!
!   and interpolates the data
!
        IF (with_eigen) THEN
           CALL interp_freq_eigen(ndata,freq_geo,xd,cgeo_eff,displa_geo, &
                           poly_order,poly_grun_red(1,1,i,n))
        ELSE
           CALL interp_freq(ndata,freq_geo,xd,poly_order,&
                                      poly_grun_red(1,1,i,n))
        ENDIF
     ENDDO
  ENDDO
!
!  deallocates space
!
  IF (with_eigen) DEALLOCATE (displa_geo)
  DEALLOCATE(freq_geo)
  DEALLOCATE(x)
  DEALLOCATE(xd)

  RETURN
END SUBROUTINE fit_frequencies_anis_reduced
