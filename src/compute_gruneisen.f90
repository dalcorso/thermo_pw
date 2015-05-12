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
  USE grun_anharmonic,        ONLY : poly_grun, poly_order
  USE mp_images,              ONLY : root_image, my_image_id

  IMPLICIT NONE

  REAL(DP), ALLOCATABLE :: freq_geo(:,:)
  INTEGER, ALLOCATABLE :: rap_geo(:,:)
  INTEGER :: n, igeo, nq

  IF ( my_image_id /= root_image ) RETURN
!
!  allocate space for the fit of the frequencies with respect to the
!  volume
!
  nq=ph_freq_save(1)%nq
  poly_order=3

  IF ( .NOT. ALLOCATED( poly_grun ) ) ALLOCATE(poly_grun(poly_order,3*nat,nq))

  ALLOCATE(freq_geo(3*nat,ngeo(1)))
  ALLOCATE(rap_geo(3*nat,ngeo(1)))

  DO n = 1, nq
     DO igeo=1,ngeo(1)
        freq_geo(1:3*nat,igeo)=ph_freq_save(igeo)%nu(1:3*nat,n)
        rap_geo(1:3*nat,igeo)=-1
     ENDDO
     CALL compute_freq_derivative(ngeo(1),freq_geo,rap_geo,omega_geo,&
                                         poly_order,poly_grun(1,1,n))
  ENDDO
  DEALLOCATE ( freq_geo )
  DEALLOCATE ( rap_geo )

  RETURN
END SUBROUTINE compute_gruneisen

SUBROUTINE compute_freq_derivative(ngeo,freq_geo,rap_geo,omega_geo,&
                                                        poly_order,poly_grun)
USE kinds, ONLY : DP
USE ions_base, ONLY : nat
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ngeo, poly_order
INTEGER,  INTENT(IN) :: rap_geo(3*nat,ngeo)
REAL(DP), INTENT(IN) :: freq_geo(3*nat,ngeo), omega_geo(ngeo)
REAL(DP), INTENT(INOUT) :: poly_grun(poly_order,3*nat)

REAL(DP), ALLOCATABLE :: frequences(:)
INTEGER,  ALLOCATABLE :: level(:,:)
INTEGER :: ibnd, jbnd, igeo, irap, central_geo

ALLOCATE (level(12,ngeo))
ALLOCATE (frequences(ngeo))

central_geo=ngeo/2
IF (MOD(ngeo,2)==1) central_geo=central_geo+1
level=1
DO ibnd=1, 3*nat
!
!   there are several representation files, we choose to order the
!   Gruneisen parameters on file as those of the central geometry
!
   irap = rap_geo(ibnd,central_geo)

   IF (irap == -1) THEN
      DO igeo=1,ngeo
         frequences(igeo)=freq_geo(ibnd,igeo)
      END DO
   ELSE
      DO igeo=1,ngeo
         DO jbnd=level(irap,igeo), 3*nat
            IF (rap_geo(jbnd,igeo)==irap) THEN
               frequences(igeo)=freq_geo(jbnd,igeo)
               level(irap,igeo)=jbnd+1
               GOTO 20
            ENDIF
         ENDDO
         CALL errore('compute_freq_derivative','representation not found',1)
20       CONTINUE
      ENDDO
   ENDIF
!
!    Fit the frequencies as a function of the volume
!
   CALL polifit( omega_geo, frequences, ngeo,  &
                                     poly_grun(:,ibnd), poly_order )
ENDDO

DEALLOCATE ( level )
DEALLOCATE ( frequences )
RETURN
END SUBROUTINE compute_freq_derivative
