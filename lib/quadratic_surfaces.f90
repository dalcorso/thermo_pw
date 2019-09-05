!
! Copyright (C) 2014-2019 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE quadratic_surfaces
!
!   This module contains the support routines for dealing with quadratic
!   surfaces interpolation. 
!
!   It provides the following routines:
!
!   fit_multi_quadratic receives as input the coordinates of some points 
!   and the value of a function in these points and finds the coefficients 
!   of the quadratic polynomial that better interpolates all the points.
!
!   evaluate_fit_quadratic evaluates the quadratic polynomial at a given 
!   input point.
!
!   evaluate_quadratic_grad evaluates the gradient of the quadratic 
!   polynomial at a given input point.
!
!   write_quadratic_hessian gives the Hessian of the quadratic polynomial 
!   and diagonalizes it finding its eigenvalues and eigenvectors.
!
!   find_quadratic_extremum receives as input the coefficients of the 
!   quadratic polynomial and finds the position of the extremum.
!
!   find_quadratic_linear_extremum receives as input the coefficients of 
!   a quadratic and a linear polynomials and finds the position of the 
!   extremum of their sum.
!
!   find_two_quadratic_extremum receives as input the coefficients of 
!   two quadratic polynomials and finds the position of the extremum of 
!   their sum.
!
!   print_quadratic_polynomial writes on output the coefficients of the
!   quadratic polynomial.
!
!   introduce_quadratic_fit writes a message with a few information on the
!   quadratic polynomial and the number of data used to fit it.
!
!   print_chisq_quadratic prints the chi square of the difference between
!   a quadratic interpolating polynomial and the values of the function 
!   in the set of points. 
!
!   compare_quadratic_fit receives a set of data and the coefficients
!   of a fitting quadratic polynomial and writes on output the data, the
!   values of the quadratic polynomial and their difference.
!
!   summarize_fitting_data writes on output the data to be fitted.
!
!
  USE kinds, ONLY : DP
  USE polynomial, ONLY : poly2
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: fit_multi_quadratic, evaluate_fit_quadratic,              &
            evaluate_quadratic_grad, write_quadratic_hessian,         &
            find_quadratic_extremum, find_quadratic_linear_extremum,  &
            find_two_quadratic_extremum, print_quadratic_polynomial,  &
            introduce_quadratic_fit, print_chisq_quadratic,           &
            compare_quadratic_fit, summarize_fitting_data

CONTAINS

SUBROUTINE fit_multi_quadratic(ndata,nvar,lsolve,x,f,p2)
!
!  This routine receives as input a set of vectors x(nvar,ndata) and
!  function values f(ndata) and gives as output the coefficients of
!  a quadratic interpolating polynomial p2. In input ndata is the number 
!  of data points, nvar is the number of variables.
! 
!  The coefficients are organized as follows:
!  a0 + \sum_i v1(i) * x(i,idata) 
!       \sum_i,j (j>=i) phi2(n(i,j)) * x(i,idata) * x(j,idata)  
!
!  n(i,j) is a linear index on phi2 that is computed as:
!  n=0
!  DO i=1,nvar
!     DO j=i,nvar
!        n=n+1
!     ENDDO
!  ENDDO
!
!  The number of coefficients of the polynomial is:
! 
!                        degree
!   number of variables  0   1   2   total number of coefficients
!      1                 1   1   1              3
!      2                 1   2   3              6
!      3                 1   3   6             10 
!      4                 1   4  10             15
!      5                 1   5  15             21
!      6                 1   6  21             28
!      7                 1   7  28             36
!      8                 1   8  36             45
!      ...                    ...
!
!   To interpolate the function it is better to give to fit_multi_quadratic 
!   at least as many data points as the number of coefficients or more. The 
!   routine makes a least square fit of the data.
!
!   lsolve can be 1, 2 or 3. See the routine min_sqr_solve for an
!   explanation of its meaning. 
!
USE linear_solvers,     ONLY : min_sqr_solve
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ndata, lsolve
REAL(DP), INTENT(IN) :: x(nvar,ndata), f(ndata)
TYPE(poly2), INTENT(INOUT) :: p2

REAL(DP), ALLOCATABLE :: amat(:,:), coeff(:)

INTEGER :: i, j, n, idata, ncoeff

ncoeff=1+nvar+p2%ncoeff2

ALLOCATE(amat(ndata,ncoeff))
ALLOCATE(coeff(ncoeff)) 

IF (nvar<1) CALL errore('fit_multi_quadratic','nvar must be larger than 1',1)
IF (ndata < ncoeff) &
   CALL errore('fit_multi_quadratic','Too few sampling data',1)
!
!  prepare the auxiliary matrix
!
amat=0.0_DP

DO idata=1,ndata
   amat(idata,1) = 1.0_DP
   n=0
   DO i=1, nvar
      amat(idata,i+1)=x(i,idata)
      DO j=i,nvar
         n=n+1
         amat(idata,1+nvar+n)=x(i,idata)*x(j,idata)
      ENDDO
   ENDDO
ENDDO
!
CALL min_sqr_solve(ndata, ncoeff, amat, f, coeff, lsolve)
!
!   assign the coefficients to the polynomial
!
p2%a0=coeff(1)
n=0
DO i=1,nvar
   p2%phi1(i)=coeff(1+i)
   DO j=i,nvar
      n=n+1
      p2%phi2(n)=coeff(1+nvar+n)
   ENDDO
ENDDO

DEALLOCATE(amat)
DEALLOCATE(coeff) 

RETURN
END SUBROUTINE fit_multi_quadratic
!
SUBROUTINE evaluate_fit_quadratic(nvar,x,f,p2)
!
!  This routine evaluates the quadratic polynomial at the point x
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
REAL(DP), INTENT(IN) :: x(nvar)
TYPE(poly2) :: p2
REAL(DP), INTENT(OUT) :: f

INTEGER :: n, i, j

n=0
f=p2%a0
DO i=1,nvar
   f=f+p2%phi1(i)*x(i)
   DO j=i, nvar
      n=n+1
      f=f+p2%phi2(n)*x(i)*x(j)
   ENDDO
ENDDO
   
RETURN
END SUBROUTINE evaluate_fit_quadratic

SUBROUTINE evaluate_quadratic_grad(nvar,x,f,p2)
!
!   This routine evaluates the gradient of the quadratic polynomial
!   at the point x.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
REAL(DP), INTENT(IN) :: x(nvar)
REAL(DP), INTENT(OUT) :: f(nvar)
TYPE(poly2), INTENT(IN) :: p2

INTEGER  :: i, j, n

n=0
f=p2%phi1
DO i=1, nvar
   DO j=i, nvar
      n=n+1
      f(i) = f(i) + p2%phi2(n) * x(j)
      f(j) = f(j) + p2%phi2(n) * x(i)
   ENDDO
ENDDO

RETURN
END SUBROUTINE evaluate_quadratic_grad

SUBROUTINE write_quadratic_hessian(nvar, p2, v, e)
!
!  This routine writes the Hessian of the quadratic polynomial and
!  its eigenvalues and eigenvectors.
!
USE diagonalize, ONLY : diagonalize_r
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
TYPE(poly2), INTENT(IN) :: p2
REAL(DP), INTENT(INOUT) :: e(nvar), v(nvar,nvar)
INTEGER ::  ideg, i, j, n

REAL(DP) :: amat(nvar, nvar)

n=0
DO i=1, nvar
   DO j=i, nvar
      n=n+1
      amat(i,j) = p2%phi2(n)
      amat(j,i) = p2%phi2(n)
   ENDDO
ENDDO

WRITE(stdout,'(/,5x,"Hessian:")')
DO ideg=1,nvar
   WRITE(stdout,'(5x,6f18.8)') amat(ideg,:)
ENDDO

IF (nvar > 1) THEN
   CALL diagonalize_r(nvar,nvar,amat,e,v)
   WRITE(stdout,'(/,5x,"Hessian eigenvalues:")')
   WRITE(stdout,'(5x,6f18.8)') e(1:nvar)

   WRITE(stdout,'(/,5x,"Hessian eigenvectors (columns):")')
   DO ideg=1,nvar
      WRITE(stdout,'(5x,6f18.8)') v(ideg,1:nvar)
   ENDDO
ENDIF
WRITE(stdout,*)

RETURN
END SUBROUTINE write_quadratic_hessian
!
SUBROUTINE find_quadratic_extremum(nvar,x,f,p2)
!
!   This routine finds the extremum x of a quadratic polynomial. 
!   On output f contains the value of the polynomial at the extremum. 
!
USE linear_solvers, ONLY : linsolvx
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
TYPE(poly2), INTENT(IN) :: p2
REAL(DP), INTENT(INOUT) :: x(nvar), f

REAL(DP) :: amat(nvar, nvar), v(nvar)
INTEGER :: i, j, n

IF (nvar < 1 ) CALL errore('find_quadratic_extremum','incorrect nvar',1)

amat=0.0_DP

n=0
v=-p2%phi1
DO i=1, nvar
   DO j=i, nvar
      n=n+1
      amat(i,j) = amat(i,j) + p2%phi2(n)
      amat(j,i) = amat(j,i) + p2%phi2(n)
   ENDDO
ENDDO

CALL linsolvx(amat,nvar,v,x)

CALL evaluate_fit_quadratic(nvar,x,f,p2)

RETURN
END SUBROUTINE find_quadratic_extremum

SUBROUTINE find_quadratic_linear_extremum(nvar,x,f,p2,p1)
!
! This routine finds the extremum of the sum of a quadratic and a linear 
! polynomials.
!
USE polynomial, ONLY : poly1, init_poly, clean_poly
IMPLICIT NONE
INTEGER,  INTENT(IN) :: nvar
TYPE(poly2), INTENT(IN) :: p2
TYPE(poly1), INTENT(IN) :: p1
REAL(DP), INTENT(OUT) :: x(nvar), f

TYPE(poly2) :: ps

CALL init_poly(nvar,ps)

ps%a0 = p2%a0 + p1%a0
ps%phi1 = p2%phi1 + p1%phi1
ps%phi2 = p2%phi2 

CALL find_quadratic_extremum(nvar,x,f,ps)

CALL clean_poly(ps)

RETURN
END SUBROUTINE find_quadratic_linear_extremum
!
SUBROUTINE find_two_quadratic_extremum(nvar,x,f,p2,p21)
!
! This routine finds the extremum of the sum of two quadratic polynomials
! p2 and p21. 
!
USE polynomial, ONLY : init_poly, clean_poly
IMPLICIT NONE
INTEGER,  INTENT(IN) :: nvar
TYPE(poly2), INTENT(IN) :: p2, p21
REAL(DP), INTENT(INOUT) :: x(nvar), f

TYPE(poly2) :: ps

CALL init_poly(nvar,ps)

ps%a0 = p2%a0 + p21%a0
ps%phi1 = p2%phi1 + p21%phi1
ps%phi2 = p2%phi2 + p21%phi2

CALL find_quadratic_extremum(nvar,x,f,ps)

CALL clean_poly(ps)

RETURN
END SUBROUTINE find_two_quadratic_extremum
!
SUBROUTINE print_quadratic_polynomial(nvar, p2)
!
!  This routine writes on output the coefficients of the quadratic polynomial. 
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nvar
TYPE(poly2), INTENT(IN) :: p2

CHARACTER(LEN=6) :: int_to_char
INTEGER :: i, j, n

WRITE(stdout,'(/,5x,"Quadratic polynomial:")') 
!
!   term of 0 degree
!
WRITE(stdout,'(5x,e15.7)') p2%a0
!
!   terms of 1 degree
!
DO i=1,nvar
   WRITE(stdout,'(5x," +",e15.7,a)') p2%phi1(i), '  x'//TRIM(int_to_char(i))
ENDDO
!
!   terms of 2 degree
!
n=0
DO i=1, nvar
   DO j=i,nvar
      n=n+1
      WRITE(stdout,'(5x," +",e15.7,a)') p2%phi2(n), &
                     '  x'//TRIM(int_to_char(i))//' x'//TRIM(int_to_char(j)) 
   ENDDO
ENDDO   

RETURN
END SUBROUTINE print_quadratic_polynomial

SUBROUTINE introduce_quadratic_fit(nvar, ncoeff, ndata)
!
!   This routine writes a small message with the information on 
!   the quadratic polynomial
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ncoeff, ndata

  WRITE(stdout,'(/,5x,"Fitting the data with a quadratic polynomial:")')

  WRITE(stdout,'(/,5x,"Number of variables:",10x,i5)')  nvar
  WRITE(stdout,'(5x,"Coefficients of the quadratic polynomial:",i5)')  ncoeff
  WRITE(stdout,'(5x,"Number of fitting data:",7x,i5,/)')  ndata

RETURN
END SUBROUTINE introduce_quadratic_fit

SUBROUTINE print_chisq_quadratic(ndata, nvar, x, f, p2)
!
!   This routine receives as input the values of a function f for ndata
!   values of the independent variables x, a set of ncoeff coefficients
!   of a quadratic interpolating polynomial and writes as output
!   the sum of the squares of the differences between the values of
!   the function and of the interpolating polynomial. 
!
IMPLICIT NONE

INTEGER :: ndata, nvar
REAL(DP) :: x(nvar, ndata), f(ndata)
TYPE(poly2) :: p2

REAL(DP) :: chisq, perc, aux
INTEGER :: idata

chisq=0.0_DP
perc=0.0_DP
DO idata=1,ndata
   CALL evaluate_fit_quadratic(nvar,x(1,idata),aux,p2)
!   WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
   IF (ABS(f(idata))>1.D-12) perc= perc + ABS((f(idata)-aux) / f(idata))
ENDDO
WRITE(stdout,'(5x,"chi square quadratic=",e18.5," relative error",e18.5,&
                                  &" %",/)') chisq/ndata, perc * 100 / ndata
RETURN
END SUBROUTINE print_chisq_quadratic

SUBROUTINE compare_quadratic_fit(ndata, nvar, x, f, p2)
!
!   This routine receives as input the values of a function f for ndata
!   values of the independent variables x, a set of ncoeff coefficients
!   of a quadratic interpolating polynomial and writes as output the 
!   input data, the interpolated data, and their differences. 
!
IMPLICIT NONE

INTEGER :: ndata, nvar
REAL(DP) :: x(nvar, ndata), f(ndata)
TYPE(poly2) :: p2

REAL(DP) :: aux
INTEGER :: idata

DO idata=1,ndata
   CALL evaluate_fit_quadratic(nvar,x(1,idata),aux,p2)
   WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
ENDDO

RETURN
END SUBROUTINE compare_quadratic_fit

SUBROUTINE summarize_fitting_data(nvar, ndata, x, f)
!
!   This routine writes in output the values of the function f in
!   the set of points x. nvar is the number of variables of the vector x.
!
IMPLICIT NONE
INTEGER :: nvar, ndata
REAL(DP) :: x(nvar,ndata), f(ndata)

INTEGER :: idata

WRITE(stdout,'(5x," Fitting the following data")')

DO idata=1,ndata
   WRITE(stdout,'(8f15.8)') x(1:nvar,idata), f(idata)
ENDDO

RETURN
END SUBROUTINE summarize_fitting_data

END MODULE quadratic_surfaces
