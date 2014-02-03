!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE compute_gruneisen()
  !
  ! use the phonon frequencies calculated for a uniform mesh at ngeo 
  ! geometries and computes the derivatives of the phonon dispersions 
  ! with respect to the volume (gruneisen parameters) on the same mesh.
  ! Only three geometries are used if ngeo is odd or two if ngeo is
  ! even.
  ! NB: we do not multiply by the volume, so these are actualy the
  ! gruneisen parameters divided by the cell volume (in (a.u.)^3). 
  !
  ! 
  USE kinds,                  ONLY : DP
  USE thermo_mod,             ONLY : ngeo, omega_geo
  USE ions_base,              ONLY : nat
  USE ph_freq_thermodynamics, ONLY : ph_freq_save
  USE ph_freq_anharmonic,     ONLY : vminf_t
  USE grun_anharmonic,        ONLY : ph_grun
  USE ph_freq_module,         ONLY : init_ph_freq
  USE ifc,                    ONLY : nq1_d, nq2_d, nq3_d
  USE mp_images,              ONLY : root_image, my_image_id

  IMPLICIT NONE


  INTEGER, PARAMETER :: m1 = 3   ! number of polynomial coefficients
  REAL(DP) :: alpha(m1)          ! the polynomial coefficients
  REAL(DP), ALLOCATABLE :: frequences(:)

  REAL(DP) :: freq
  INTEGER, ALLOCATABLE :: level(:,:)
  INTEGER :: ibnd, jbnd, n, ipos, jpos, irap, igeo
  INTEGER :: geo, use_geo, nq

  IF ( my_image_id /= root_image ) RETURN

  IF ( .NOT. ALLOCATED( level ) ) ALLOCATE (level(12,ngeo))
  IF ( .NOT. ALLOCATED( frequences ) ) ALLOCATE (frequences(ngeo))
!
!  allocate space for the gruneisen parameters
!
!
!  Compute the Gruneisen parameters
!
  geo=ngeo/2-1
  IF (mod(ngeo,2)==0) THEN
     use_geo=2
  ELSE
     use_geo=3
  ENDIF
  nq=ph_freq_save(geo+2)%nq
  CALL init_ph_freq(ph_grun, nat, nq1_d, nq2_d, nq3_d, nq)
  DO n = 1, ph_grun%nq
     level=1
     DO ibnd=1, 3*nat
!
!   there are several representation files, we choose to order the
!   Gruneisen parameters on file as those of the central geometry
!
        ipos = ibnd + (n-1) * 3 * nat
        irap = ph_freq_save(geo+2)%rap(ibnd,n)
          
        IF (irap == -1) THEN
           DO igeo=1,use_geo
              frequences(igeo)=ph_freq_save(geo+igeo)%nu(ibnd,n)
           END DO
        ELSE
           DO igeo=1,use_geo
              DO jbnd=level(irap,igeo), 3*nat
                 jpos = jbnd + (n-1) * 3 * nat
                 IF (ph_freq_save(geo+igeo)%rap(jbnd,n)==irap) THEN
                    frequences(igeo)=ph_freq_save(geo+igeo)%nu(jbnd,n)
                    level(irap,igeo)=jbnd+1
                    GOTO 20
                 ENDIF
              ENDDO
              CALL errore('compute_gruneisen','representation not found',1)
20            CONTINUE
           ENDDO
        ENDIF
!
!    Use the frequencies of the central geometry if there is 
!    an odd number of geometries, or the average of the two central
!    geometries if there is an even number
!
        IF (mod(ngeo,2)==0) THEN
           freq  = 0.5_DP * ( ph_freq_save(ngeo/2)%nu(ibnd,n) + &
                              ph_freq_save(ngeo/2+1)%nu(ibnd,n) )
        ELSE
           freq = ph_freq_save(ngeo/2+1)%nu(ibnd,n)
        ENDIF
        CALL polifit( omega_geo(geo+1), frequences, use_geo, alpha, m1 )
        IF ( freq > 1.d-4 ) THEN
           freq = alpha(1) + alpha(2) * vminf_t(1) + alpha(3) * vminf_t(1)**2
           ph_grun%nu(ibnd,n) = - (alpha(2) + 2.0_DP * alpha(3) * &
                                              vminf_t(1)) / freq
           ph_freq_save(ngeo/2+1)%nu(ibnd,n) = freq
        ELSE
           ph_grun%nu(ibnd,n) = 0.0_DP
        END IF
     ENDDO
  ENDDO
  DEALLOCATE ( level )
  DEALLOCATE ( frequences )

   RETURN
END SUBROUTINE compute_gruneisen
