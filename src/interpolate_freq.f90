!
! Copyright (C) 2018 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE freq_interpolate
!
!  This module offers five routines to interpolate the phonon frequencies 
!  as a function of external parameters. It assumes that there 
!  are 3*nat frequencies given as a function of the volume or of
!  the crystal parameters and gives as output, for each frequency, the 
!  coefficients of an interpolating polynomial. In input there are ngeo 
!  frequencies for each mode and ngeo values of the parameter(s). For 
!  one-dimensional fits the order of the interpolating polinomial is given 
!  as input, for the many parameters case it is a quadratic polynomial.
!  The routines are:
!
!  interp_freq frequencies are given as a function of the volume and
!              the routine interpolates the frequencies as they are.
!
!  interp_freq_rap frequencies are given as a function of the volume.
!              The routine receives the symmetry representation 
!              of each mode and uses them to group frequencies with the 
!              same symmetry.
!
!  interp_freq_eigen frequencies are given as a function of the volume.
!              The routine receives the mode eigenvectors and
!              uses them to build a quantity whose derivatives 
!              are the expectation values of the derivative of the
!              dynamical matrix on the eigenvectors of the central geometry.
!
! interp_freq_anis: the frequencies are given as a function of crystal 
!              parameters and the routines interpolates them as they are. 
!
! interp_freq_anis_eigen it is similar to interp_freq_eigen but
!              frequencies are given as a function of the crystal parameters. 
!
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : nat
  !
  SAVE

  PRIVATE

  PUBLIC interp_freq, interp_freq_eigen, interp_freq_rap, interp_freq_anis, &
         interp_freq_anis_eigen

CONTAINS

!-------------------------------------------------------------------------
SUBROUTINE interp_freq(ngeo,freq,omega,poly_degree_grun,poly_grun)
!-------------------------------------------------------------------------

USE polyfit_mod, ONLY : polyfit

IMPLICIT NONE
INTEGER,  INTENT(IN) :: ngeo, poly_degree_grun
REAL(DP), INTENT(IN) :: freq(3*nat,ngeo), omega(ngeo)
REAL(DP), INTENT(INOUT) :: poly_grun(poly_degree_grun+1,3*nat)

INTEGER :: imode

DO imode=1, 3*nat
!
!    Fits the frequencies as a function of the volume
!
   CALL polyfit( omega, freq(imode,:), ngeo, poly_grun(:,imode), &
                                                        poly_degree_grun )
ENDDO

RETURN
END SUBROUTINE interp_freq
!
!-------------------------------------------------------------------------
SUBROUTINE interp_freq_eigen(ngeo, freq, omega, central_geo, &
                                      displa, poly_degree_grun, poly_grun)
!-------------------------------------------------------------------------
USE polyfit_mod, ONLY : polyfit
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ngeo, central_geo, poly_degree_grun
REAL(DP), INTENT(IN) :: freq(3*nat,ngeo), omega(ngeo)
COMPLEX(DP), INTENT(IN) :: displa(3*nat,3*nat,ngeo)
REAL(DP), INTENT(INOUT) :: poly_grun(poly_degree_grun+1,3*nat)

REAL(DP), ALLOCATABLE :: frequences(:)
REAL(DP) :: overlap
COMPLEX(DP) :: ZDOTC
INTEGER :: imode, igeo, jmode

ALLOCATE(frequences(ngeo))

DO imode=1, 3*nat
   frequences=0.0_DP
   DO igeo=1, ngeo
      IF (igeo /= central_geo) THEN
!
!   The frequency function defined below has the same 
!   derivatives of the real frequency at the central geometry.
!   For other geometries these frequencies often differ from the real 
!   frequencies for a term quadratic in the difference between the 
!   central geometry and the actual geometry, so be careful when you 
!   plot the Gruneisen parameters for volumes too far from the central 
!   geometry. In some particular points the difference can be large.
!   Gruneisen parameters calculated using this formula should be the
!   same as those calculated as the expectation value of the derivative 
!   of the square root of the dynamical matrix on the central geometry 
!   eigenvectors.
!
         DO jmode=1,3*nat
            overlap=ABS(ZDOTC(3*nat,displa(1,jmode,igeo),1,&
                                    displa(1,imode,central_geo),1))**2
            frequences(igeo)=frequences(igeo) + freq(jmode,igeo)*overlap
         ENDDO
      ELSE
         frequences(igeo)=freq(imode,igeo)
      ENDIF
   END DO
!
!    Fits the frequencies as a function of the volume
!
   CALL polyfit( omega, frequences, ngeo, poly_grun(:,imode), poly_degree_grun)
ENDDO

DEALLOCATE(frequences)

RETURN
END SUBROUTINE interp_freq_eigen

!-------------------------------------------------------------------------
SUBROUTINE interp_freq_rap(ngeo, freq, omega, central_geo, &
                                         rap, poly_degree_grun, poly_grun)
!-------------------------------------------------------------------------
USE polyfit_mod, ONLY : polyfit

IMPLICIT NONE
INTEGER,  INTENT(IN) :: ngeo, poly_degree_grun, central_geo
INTEGER,  INTENT(IN) :: rap(3*nat,ngeo)
REAL(DP), INTENT(IN) :: freq(3*nat,ngeo), omega(ngeo)
REAL(DP), INTENT(INOUT) :: poly_grun(poly_degree_grun+1,3*nat)

REAL(DP), ALLOCATABLE :: frequences(:)
INTEGER,  ALLOCATABLE :: level(:,:)
INTEGER :: imode, jbnd, igeo, irap

ALLOCATE (level(12,ngeo))
ALLOCATE (frequences(ngeo))

level=1
DO imode=1, 3*nat
!
!   there are several representation files, we choose to order the
!   Gruneisen parameters on file as those of the central geometry
!
   irap = rap(imode,central_geo)
   IF (irap > 12) CALL errore('interp_freq_rap',' Too many representations',1)

   IF (irap == -1) THEN
      frequences(:)=freq(imode,:)
   ELSE
      DO igeo=1,ngeo
         DO jbnd=level(irap,igeo), 3*nat
            IF (rap(jbnd,igeo)==irap) THEN
               frequences(igeo)=freq(jbnd,igeo)
               level(irap,igeo)=jbnd+1
               GOTO 20
            ENDIF
         ENDDO
         CALL errore('interp_freq_rap','representation not found',1)
20       CONTINUE
      ENDDO
   ENDIF
!
!    Fits the frequencies as a function of the volume
!
   CALL polyfit(omega, frequences, ngeo, poly_grun(:,imode), poly_degree_grun)
ENDDO

DEALLOCATE( level )
DEALLOCATE( frequences )
RETURN
END SUBROUTINE interp_freq_rap

!-------------------------------------------------------------------------
SUBROUTINE interp_freq_anis(ngeo,freq,lsolve,x,nvar,p_grun_p2)
!-------------------------------------------------------------------------

USE polynomial, ONLY : poly2
USE quadratic_surfaces, ONLY : fit_multi_quadratic

IMPLICIT NONE
INTEGER,  INTENT(IN) :: ngeo, nvar, lsolve
REAL(DP), INTENT(IN) :: freq(3*nat,ngeo), x(nvar,ngeo)
TYPE(poly2), INTENT(INOUT) :: p_grun_p2(3*nat)

INTEGER :: imode

DO imode=1,3*nat
   CALL fit_multi_quadratic(ngeo, nvar, lsolve, x, freq(imode,:), p_grun_p2(imode))
ENDDO

RETURN
END SUBROUTINE interp_freq_anis

!--------------------------------------------------------------------------
SUBROUTINE interp_freq_anis_eigen(ngeo,freq,lsolve,x,central_geo,displa,nvar, &
                                                               p_grun_p2)
!--------------------------------------------------------------------------
!
USE polynomial, ONLY : poly2
USE quadratic_surfaces, ONLY : fit_multi_quadratic

IMPLICIT NONE
INTEGER,  INTENT(IN) :: ngeo, nvar, central_geo, lsolve
REAL(DP), INTENT(IN) :: freq(3*nat,ngeo), x(nvar,ngeo)
COMPLEX(DP), INTENT(IN) :: displa(3*nat,3*nat,ngeo)
TYPE(poly2), INTENT(INOUT) :: p_grun_p2(3*nat) 

REAL(DP), ALLOCATABLE :: f(:)
REAL(DP) :: overlap
INTEGER  :: imode, jmode, igeo
COMPLEX(DP), EXTERNAL :: ZDOTC

ALLOCATE(f(ngeo))

DO imode=1, 3*nat
   f=0.0_DP
   DO igeo=1, ngeo
      IF (igeo /= central_geo) THEN
         DO jmode=1,3*nat
            overlap=ABS(ZDOTC(3*nat, displa(1,jmode,igeo), 1, &
                                     displa(1,imode,central_geo), 1))**2
            f(igeo)=f(igeo)+freq(jmode,igeo)*overlap
         ENDDO
      ELSE
         f(igeo)=freq(imode,igeo)
      ENDIF
   ENDDO
   CALL fit_multi_quadratic(ngeo, nvar, lsolve, x, f, p_grun_p2(imode)) 
ENDDO

DEALLOCATE(f)

RETURN
END SUBROUTINE interp_freq_anis_eigen

END MODULE freq_interpolate
