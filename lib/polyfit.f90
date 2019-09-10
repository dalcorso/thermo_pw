!
! Copyright (C) 2019 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE polyfit_mod
!
!  This module contains routines for polynomial interpolation.
!  The degree of the polynomial is arbitrary but there is only one 
!  variable. For multivariable polynomial see linear_surfaces, 
!  quadratic_surfaces, cubic_surfaces and quartic_surfaces.
! 
!  Presently the modulus contains the routines
!
!  polyfit : interpolate a set of data with a polynomial of 
!            arbitrary degree
!
!  compute_poly : computes the value of the polynomial in one point  
!
!  compute_poly_deriv : computes the value of the first derivative of the
!                 polynomial in one point
!
!  write_poly : writes on output the coefficients of the polynomial.
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: polyfit, compute_poly, compute_poly_deriv, write_poly

CONTAINS

!--------------------------------------------------------------------------
SUBROUTINE polyfit(x,y,ndati,a,degree)
!--------------------------------------------------------------------------
!
!  This routine fits a set of data with a polynomial of one variable and
!  arbitrary degree. 
!
USE linear_solvers, ONLY : linsolvx
IMPLICIT NONE
INTEGER, INTENT(IN)   ::   ndati, &  ! number of data points
                           degree    ! degree of the polynomial

REAL(DP), INTENT(IN)  :: x(ndati), y(ndati)
REAL(DP), INTENT(OUT) :: a(degree+1)

REAL(DP) :: amat(degree+1,degree+1), bvec(degree+1)
INTEGER  ::  i, j, k, ncoeff 

ncoeff=degree+1

DO k=1,ncoeff
   DO j=1,ncoeff
      amat(k,j)=0.d0
      DO i=1,ndati
         amat(k,j)=amat(k,j)+x(i)**(j-1)*x(i)**(k-1)
      ENDDO
   ENDDO
   bvec(k)=0.d0
   DO i=1,ndati
      bvec(k)=bvec(k)+y(i)*x(i)**(k-1)
   ENDDO
ENDDO

CALL linsolvx(amat,ncoeff,bvec,a)

RETURN
END SUBROUTINE polyfit
!
!--------------------------------------------------------------------------
SUBROUTINE compute_poly(x, degree, poly, f)
!--------------------------------------------------------------------------
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER :: degree
REAL(DP), INTENT(IN) :: poly(degree+1), x
REAL(DP), INTENT(OUT) :: f

INTEGER :: i, ncoeff

ncoeff=degree+1
f=poly(ncoeff)
DO i=ncoeff-1,1,-1
   f = poly(i) + f*x
ENDDO

RETURN
END SUBROUTINE compute_poly

!--------------------------------------------------------------------------
SUBROUTINE compute_poly_deriv(x, degree, poly, g)
!--------------------------------------------------------------------------
!
! gives as output the derivative of the polynomial with respect to x
! 
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree
REAL(DP), INTENT(IN) :: poly(degree+1), x

REAL(DP), INTENT(OUT) :: g

INTEGER :: i, ncoeff

ncoeff=degree+1

g=poly(ncoeff)*(ncoeff-1.0_DP)
DO i=ncoeff-1,2,-1
   g = poly(i)*(i-1.0_DP) + g*x
ENDDO

RETURN
END SUBROUTINE compute_poly_deriv

!--------------------------------------------------------------------------
SUBROUTINE write_poly(a,degree)
!--------------------------------------------------------------------------
!
!   This routine writes on output the coefficients of a polynomial of
!   one variable
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree
REAL(DP), INTENT(IN) :: a(degree+1)

INTEGER :: i, ncoeff
CHARACTER(LEN=6) :: int_to_char

ncoeff=1+degree
WRITE(stdout,'(/,5x,"Polynomial coefficients")')
DO i=1,ncoeff
   WRITE(stdout,'(5x,a,e20.12)') "a"//TRIM(INT_TO_CHAR(i))//"=", a(i)
END DO

RETURN
END SUBROUTINE write_poly

END MODULE polyfit_mod
