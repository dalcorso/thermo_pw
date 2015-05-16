!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE compute_freq_derivative_anis_eigen(ngeo,freq_geo,celldm_geo,&
                     displa_geo, degree,nvar,ibrav,poly_grun)

USE kinds, ONLY : DP
USE ions_base, ONLY : nat
USE quadratic_surfaces, ONLY : fit_multi_quadratic
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ngeo, degree, nvar, ibrav
REAL(DP), INTENT(IN) :: freq_geo(3*nat,ngeo), celldm_geo(6,ngeo)
COMPLEX(DP), INTENT(IN) :: displa_geo(3*nat,3*nat,ngeo)
REAL(DP), INTENT(INOUT) :: poly_grun(nvar,3*nat)

REAL(DP), ALLOCATABLE :: frequences(:), x(:,:), f(:)
REAL(DP) :: overlap
INTEGER :: imode, jmode, igeo, idata, central_geo
COMPLEX(DP) :: ZDOTC

ALLOCATE (frequences(ngeo))
ALLOCATE (x(degree,ngeo))
ALLOCATE (f(ngeo))

central_geo=ngeo/2+1

DO imode=1, 3*nat
   frequences=0.0_DP
   DO igeo=1, ngeo
      IF (igeo /= central_geo) THEN
         DO jmode=1,3*nat
            overlap=ABS(ZDOTC(3*nat,displa_geo(1,jmode,igeo),1,&
                                      displa_geo(1,imode,central_geo),1))**2
            frequences(igeo)=frequences(igeo) + freq_geo(jmode,igeo)* overlap
         ENDDO
      ELSE
         frequences(igeo)=freq_geo(imode,igeo)
      ENDIF
   END DO
!
   SELECT CASE (ibrav)
      CASE(1,2,3)
         DO idata=1,ngeo
            x(1,idata)=celldm_geo(1,idata)
            f(idata)=frequences(idata)
         ENDDO
      CASE(4,5,6,7)
         DO idata=1,ngeo
            x(1,idata)=celldm_geo(1,idata)
            IF (ibrav==5) THEN
               x(2,idata)=ACOS(celldm_geo(4,idata))
            ELSE
               x(2,idata)=celldm_geo(3,idata)
            ENDIF
            f(idata)=frequences(idata)
         END DO
      CASE(8,9,91,10,11)
         DO idata=1,ngeo
            x(1,idata)=celldm_geo(1,idata)
            x(2,idata)=celldm_geo(2,idata)
            x(3,idata)=celldm_geo(3,idata)
            f(idata)=frequences(idata)
         ENDDO
      CASE(12,-12,13,-13)
         DO idata=1,ngeo
            x(1,idata)=celldm_geo(1,idata)
            x(2,idata)=celldm_geo(2,idata)
            x(3,idata)=celldm_geo(3,idata)
            IF (ibrav>0) THEN
!
!   c unique
!
               x(4,idata)=ACOS(celldm_geo(4,idata))
            ELSE
!
!   b unique
!
               x(4,idata)=ACOS(celldm_geo(5,idata))
            ENDIF
            f(idata)=frequences(idata)
         END DO
      CASE DEFAULT
         DO idata=1,ngeo
            x(1,idata)=celldm_geo(1,idata)
            x(2,idata)=celldm_geo(2,idata)
            x(3,idata)=celldm_geo(3,idata)
            x(4,idata)=ACOS(celldm_geo(4,idata))
            x(5,idata)=ACOS(celldm_geo(5,idata))
            x(6,idata)=ACOS(celldm_geo(6,idata))
            f(idata)=frequences(idata) 
         ENDDO
   END SELECT
   CALL fit_multi_quadratic(ngeo,degree,nvar,x,f,poly_grun(:,imode))
ENDDO

DEALLOCATE ( x )
DEALLOCATE ( f )
DEALLOCATE ( frequences )

RETURN
END SUBROUTINE compute_freq_derivative_anis_eigen
