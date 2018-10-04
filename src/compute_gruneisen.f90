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
  ! with respect to the volume (Gruneisen parameters) on the same mesh.
  !
  ! NB: we do not multiply by the volume, so these are actualy the
  ! Gruneisen parameters divided by the cell volume (in (a.u.)^3). 
  !
  ! 
  USE kinds,                  ONLY : DP
  USE thermo_mod,             ONLY : ngeo, omega_geo, no_ph
  USE ions_base,              ONLY : nat
  USE ph_freq_thermodynamics, ONLY : ph_freq_save
  USE grun_anharmonic,        ONLY : poly_grun, poly_order
  USE control_thermo,         ONLY : with_eigen
  USE mp_world,               ONLY : world_comm
  USE mp,                     ONLY : mp_sum

  IMPLICIT NONE

  REAL(DP), ALLOCATABLE :: freq_geo(:,:)
  COMPLEX(DP), ALLOCATABLE ::  displa_geo(:,:,:)
  INTEGER, ALLOCATABLE :: rap_geo(:,:)
  INTEGER :: iq_eff, n, igeo, nq, startq, lastq
!
!  allocate space for the fit of the frequencies with respect to the
!  volume
!
  nq=ph_freq_save(1)%nq
  poly_order=5

  CALL divide (world_comm, nq, startq, lastq)
  IF ( .NOT. ALLOCATED( poly_grun ) ) ALLOCATE(poly_grun(poly_order,3*nat,&
                                                      startq:lastq))

  ALLOCATE(freq_geo(3*nat,ngeo(1)))
  ALLOCATE(rap_geo(3*nat,ngeo(1)))
  IF (with_eigen) ALLOCATE(displa_geo(3*nat,3*nat,ngeo(1)))
!
!  representations are not used here
!
  rap_geo=-1
  poly_grun=0.0_DP
  iq_eff=0
  DO n = startq, lastq
     iq_eff=iq_eff+1
     DO igeo=1,ngeo(1)
        IF (.NOT. no_ph(igeo)) THEN
           freq_geo(1:3*nat,igeo)=ph_freq_save(igeo)%nu(1:3*nat,iq_eff)
           IF (with_eigen) displa_geo(1:3*nat, 1:3*nat, igeo)= &
                           ph_freq_save(igeo)%displa(1:3*nat,1:3*nat,iq_eff)
        ENDIF
     ENDDO
     IF (with_eigen) THEN
        CALL compute_freq_derivative_eigen(ngeo(1),freq_geo,omega_geo,     &
                               displa_geo,no_ph, poly_order,poly_grun(1,1,n))
     ELSE
        CALL compute_freq_derivative(ngeo(1),freq_geo,rap_geo,omega_geo,   &
                                             no_ph, poly_order,poly_grun(1,1,n))
     END IF
  ENDDO
!
!  the poly_grun are still distributed
!
  IF (with_eigen) DEALLOCATE ( displa_geo )
  DEALLOCATE ( freq_geo )
  DEALLOCATE ( rap_geo )

  RETURN
END SUBROUTINE compute_gruneisen

SUBROUTINE compute_freq_derivative(ngeo,freq_geo,rap_geo,omega_geo,no_ph,&
                                                        poly_order,poly_grun)
USE kinds, ONLY : DP
USE ions_base, ONLY : nat
USE quadratic_surfaces, ONLY : polifit
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ngeo, poly_order
INTEGER,  INTENT(IN) :: rap_geo(3*nat,ngeo)
REAL(DP), INTENT(IN) :: freq_geo(3*nat,ngeo), omega_geo(ngeo)
LOGICAL, INTENT(IN)  :: no_ph(ngeo)
REAL(DP), INTENT(INOUT) :: poly_grun(poly_order,3*nat)

REAL(DP), ALLOCATABLE :: frequences(:), omega_data(:)
INTEGER,  ALLOCATABLE :: level(:,:)
INTEGER :: ibnd, jbnd, igeo, irap, central_geo, ndata

ALLOCATE (level(12,ngeo))
ALLOCATE (frequences(ngeo))
ALLOCATE (omega_data(ngeo))

CALL find_central_geo(ngeo, no_ph, central_geo)

ndata=0
DO igeo=1,ngeo
   IF (.NOT. no_ph(igeo)) THEN
      ndata=ndata+1
      omega_data(ndata)=omega_geo(igeo)
   END IF
END DO
level=1
DO ibnd=1, 3*nat
!
!   there are several representation files, we choose to order the
!   Gruneisen parameters on file as those of the central geometry
!
   irap = rap_geo(ibnd,central_geo)

   IF (irap == -1) THEN
      ndata=0
      DO igeo=1,ngeo
         IF (.NOT. no_ph(igeo)) THEN
            ndata=ndata+1
            frequences(ndata)=freq_geo(ibnd,igeo)
         ENDIF
      END DO
   ELSE
      ndata=0
      DO igeo=1,ngeo
         IF (.NOT. no_ph(igeo)) THEN
            ndata=ndata+1
            DO jbnd=level(irap,igeo), 3*nat
               IF (rap_geo(jbnd,igeo)==irap) THEN
                  frequences(ndata)=freq_geo(jbnd,igeo)
                  level(irap,igeo)=jbnd+1
                  GOTO 20
               ENDIF
            ENDDO
            CALL errore('compute_freq_derivative','representation not found',1)
20          CONTINUE
         ENDIF
      ENDDO
   ENDIF
!
!    Fit the frequencies as a function of the volume
!
   CALL polifit( omega_data, frequences, ndata,  &
                                     poly_grun(:,ibnd), poly_order )
ENDDO

DEALLOCATE ( level )
DEALLOCATE ( frequences )
DEALLOCATE (omega_data)
RETURN
END SUBROUTINE compute_freq_derivative

SUBROUTINE compute_freq_derivative_eigen(ngeo,freq_geo,omega_geo,&
                                         displa_geo, no_ph, poly_order,poly_grun)

USE kinds, ONLY : DP
USE ions_base, ONLY : nat
USE quadratic_surfaces, ONLY : polifit
IMPLICIT NONE
INTEGER, INTENT(IN) :: ngeo, poly_order
REAL(DP), INTENT(IN) :: freq_geo(3*nat,ngeo), omega_geo(ngeo)
COMPLEX(DP), INTENT(IN) :: displa_geo(3*nat,3*nat,ngeo)
LOGICAL, INTENT(IN) :: no_ph(ngeo)
REAL(DP), INTENT(INOUT) :: poly_grun(poly_order,3*nat)

REAL(DP), ALLOCATABLE :: frequences(:), omega_data(:)
REAL(DP) :: overlap
COMPLEX(DP) :: ZDOTC
INTEGER :: ibnd, igeo, central_geo, jmode, ndata

ALLOCATE (frequences(ngeo))
ALLOCATE(omega_data(ngeo))

CALL find_central_geo(ngeo, no_ph, central_geo)

ndata=0
DO igeo=1,ngeo
   IF (.NOT. no_ph(igeo)) THEN
      ndata=ndata+1
      omega_data(ndata)=omega_geo(igeo)
   ENDIF
ENDDO


DO ibnd=1, 3*nat
   frequences=0.0_DP
   ndata=0
   DO igeo=1, ngeo
      IF (.NOT. no_ph(igeo)) THEN
         ndata=ndata+1
         IF (igeo /= central_geo) THEN
!
!   It can be shown that the frequency function defined below has the same 
!   derivatives of the real frequency at the central geometry.
!   For other geometries these frequencies often differ from the real 
!   frequencies for a term quadratic in the difference between the 
!   central geometry and the actual geometry, so be careful when you 
!   plot the Gruneisen parameters for volumes too far from the central 
!   geometry. In some particular points the difference can be larger.
!   Gruneisen parameters calculated using this formula should be the
!   same as those calculated as the expectation value of the derivative 
!   of the dynamical matrix on the central geometry eigenvectors.
!
            DO jmode=1,3*nat
               overlap=ABS(ZDOTC(3*nat,displa_geo(1,jmode,igeo),1,&
                                        displa_geo(1,ibnd,central_geo),1))**2
               frequences(ndata)=frequences(ndata) + freq_geo(jmode,igeo)* overlap
            ENDDO
         ELSE
            frequences(ndata)=freq_geo(ibnd,igeo)
         ENDIF
      END IF
   END DO

   CALL polifit( omega_data, frequences, ndata,  &
                                     poly_grun(:,ibnd), poly_order )
ENDDO

DEALLOCATE ( frequences )
DEALLOCATE ( omega_data )

RETURN
END SUBROUTINE compute_freq_derivative_eigen
