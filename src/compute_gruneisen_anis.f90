!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE fit_frequencies_anis()
  !
  ! use the phonon frequencies calculated for a uniform mesh at ngeo 
  ! geometries and computes the coefficients of a second order polynomial 
  ! that fit the frequencies as a function of the crystal parameters
  ! 
  USE kinds,                  ONLY : DP
  USE thermo_mod,             ONLY : ngeo, celldm_geo, no_ph
  USE ions_base,              ONLY : nat
  USE cell_base,              ONLY : ibrav
  USE ph_freq_thermodynamics, ONLY : ph_freq_save
  USE grun_anharmonic,        ONLY : poly_grun
  USE control_thermo,         ONLY : with_eigen

  IMPLICIT NONE

  REAL(DP), ALLOCATABLE :: freq_geo(:,:)
  COMPLEX(DP), ALLOCATABLE ::  displa_geo(:,:,:)
  INTEGER :: n, igeo, nq, degree, nvar, nwork, startq, lastq, iq_eff
  INTEGER :: compute_nwork
!
!  allocate space for the fit of the frequencies with respect to the
!  celldm parameters
!
  nq=ph_freq_save(1)%nq
  startq=ph_freq_save(1)%startq
  lastq=ph_freq_save(1)%lastq
  CALL compute_degree(ibrav,degree,nvar)
  nwork=compute_nwork()

  IF ( .NOT. ALLOCATED( poly_grun ) ) ALLOCATE(poly_grun(nvar,3*nat,&
                                                              startq:lastq))

  ALLOCATE(freq_geo(3*nat,nwork))
  IF (with_eigen) ALLOCATE(displa_geo(3*nat,3*nat,nwork))
!
  iq_eff=0
  DO n = startq, lastq
     iq_eff=iq_eff+1
     DO igeo=1,nwork
        IF (.NOT. no_ph(igeo)) THEN
           freq_geo(1:3*nat,igeo)=ph_freq_save(igeo)%nu(1:3*nat,iq_eff)
           IF (with_eigen) displa_geo(1:3*nat, 1:3*nat, igeo)= &
                           ph_freq_save(igeo)%displa(1:3*nat,1:3*nat,iq_eff)
        ENDIF
     ENDDO
     IF (with_eigen) THEN
        CALL compute_freq_derivative_anis_eigen(nwork,freq_geo, &
                               celldm_geo,displa_geo,degree,nvar,ibrav, &
                               no_ph,poly_grun(1,1,n))
     ELSE
        CALL compute_freq_derivative_anis(nwork,freq_geo, &
                    celldm_geo,degree,nvar,ibrav,no_ph,poly_grun(1,1,n))
     END IF
  ENDDO

  IF (with_eigen) DEALLOCATE ( displa_geo )
  DEALLOCATE ( freq_geo )

  RETURN
END SUBROUTINE fit_frequencies_anis

SUBROUTINE compute_freq_derivative_anis_eigen(ngeo,freq_geo,celldm_geo,&
                     displa_geo,degree,nvar,ibrav,no_ph,poly_grun)

USE kinds, ONLY : DP
USE ions_base, ONLY : nat
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE lattices,  ONLY : compress_celldm
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ngeo, degree, nvar, ibrav
REAL(DP), INTENT(IN) :: freq_geo(3*nat,ngeo), celldm_geo(6,ngeo)
COMPLEX(DP), INTENT(IN) :: displa_geo(3*nat,3*nat,ngeo)
LOGICAL,  INTENT(IN) :: no_ph(ngeo)
REAL(DP), INTENT(INOUT) :: poly_grun(nvar,3*nat)

REAL(DP), ALLOCATABLE :: x(:,:), f(:)
REAL(DP) :: overlap
INTEGER :: imode, jmode, igeo, idata, central_geo, ndata
COMPLEX(DP), EXTERNAL :: ZDOTC

ALLOCATE (x(degree,ngeo))
ALLOCATE (f(ngeo))

central_geo=ngeo/2+1
IF (no_ph(central_geo)) THEN
   DO igeo=1,ngeo/2
      central_geo=central_geo-igeo
      IF (.NOT. no_ph(central_geo)) EXIT
      central_geo=central_geo+2*igeo
      IF (.NOT. no_ph(central_geo)) EXIT
      central_geo=central_geo-igeo
   ENDDO
ENDIF

ndata=0
DO idata=1,ngeo
   IF (.NOT. no_ph(idata)) THEN
      ndata=ndata+1
      CALL compress_celldm(celldm_geo(1,idata),x(1,idata),degree,ibrav)
   ENDIF
ENDDO

DO imode=1, 3*nat
   f=0.0_DP
   ndata=0
   DO igeo=1, ngeo
      IF (.NOT. no_ph(igeo)) THEN
         ndata=ndata+1
         IF (igeo /= central_geo) THEN
            DO jmode=1,3*nat
               overlap=ABS(ZDOTC(3*nat,displa_geo(1,jmode,igeo),1,&
                                         displa_geo(1,imode,central_geo),1))**2
               f(ndata)=f(ndata) + freq_geo(jmode,igeo)* overlap
            ENDDO
         ELSE
            f(ndata)=freq_geo(imode,igeo)
         ENDIF
      ENDIF
   ENDDO
   CALL fit_multi_quadratic(ndata,degree,nvar,x,f,poly_grun(:,imode))
ENDDO

DEALLOCATE ( x )
DEALLOCATE ( f )

RETURN
END SUBROUTINE compute_freq_derivative_anis_eigen

SUBROUTINE compute_freq_derivative_anis(ngeo,freq_geo,celldm_geo,degree,nvar,&
                                             ibrav,no_ph,poly_grun)
USE kinds,     ONLY : DP
USE ions_base, ONLY : nat
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE lattices,  ONLY : compress_celldm

IMPLICIT NONE
INTEGER,  INTENT(IN) :: ngeo, degree, nvar, ibrav
REAL(DP), INTENT(IN) :: freq_geo(3*nat,ngeo), celldm_geo(6,ngeo)
LOGICAL,  INTENT(IN) :: no_ph(ngeo)
REAL(DP), INTENT(INOUT) :: poly_grun(nvar,3*nat)

REAL(DP), ALLOCATABLE :: x(:,:), f(:)
INTEGER :: imode, igeo, idata, ndata

ALLOCATE (x(degree,ngeo))
ALLOCATE (f(ngeo))

ndata=0
DO idata=1,ngeo
   IF (.NOT. no_ph(idata)) THEN
      ndata=ndata+1
      CALL compress_celldm(celldm_geo(1,idata),x(1,idata),degree,ibrav)
   ENDIF
ENDDO

DO imode=1, 3*nat
   f=0.0_DP
   ndata=0
   DO igeo=1, ngeo
      IF (.NOT. no_ph(igeo)) THEN
         ndata=ndata+1
         f(ndata)=freq_geo(imode,igeo)
      ENDIF
   ENDDO
   CALL fit_multi_quadratic(ndata,degree,nvar,x,f,poly_grun(:,imode))
ENDDO

DEALLOCATE ( x )
DEALLOCATE ( f )

RETURN
END SUBROUTINE compute_freq_derivative_anis
