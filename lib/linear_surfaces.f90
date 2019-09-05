!
! Copyright (C) 2019 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE linear_surfaces
!
!   This module contains the support routines for dealing with linear
!   surfaces interpolation. 
!
!   It provides the following routines:
!
!   fit_multi_linear receives as input the coordinates of some
!   points and the value of the function in these points and
!   finds the coefficients of the linear polynomial that better
!   interpolates all the points.
!
!   evaluate_fit_linear evaluates the linear polynomial at a given 
!   input point.
!
!   set_linear_gradient gives the gradient of the linear polynomial.
!   (which is constant independent on the input point).
!
!   print_linear_polynomial writes on output the coefficients of the linear
!   polynomial.
!
!   introduce_linear_fit writes a message with a few information on the
!   linear polynomial and the number of data used to fit it.
!
!   print_chisq_linear writes the sum of the squares of the difference 
!   between the input data and the interpolating linear polynomial divided
!   by the number of data points.
!
  USE kinds, ONLY : DP
  USE polynomial, ONLY : poly1
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: fit_multi_linear, evaluate_fit_linear, set_linear_grad,  &
            print_linear_polynomial, introduce_linear_fit, &
            print_chisq_linear                                        

CONTAINS

SUBROUTINE fit_multi_linear(ndata,nvar,lsolve,x,f,p_in)
!
!  This routine receives as input a set of vectors x(nvar,ndata) and
!  function values f(ndata) and gives as output the coefficients of
!  a linear interpolating polynomial a0 and v1(nvar). In input
!  ndata is the number of data points, nvar is the number of variables.
!
!  The coefficients are organized as follows:
!  p_in%a0 + \sum_j p_in%phi1(j)  x(j,i) 
!
!  The number of coefficients necessary for the fit is nvar+1, where
!  nvar is the number of variables of the polynomial.
!  To interpolate the function it is better to give to fit_multi_linear
!  at least as many data points as the number of coefficients or more. The 
!  routine makes a least square fit of the data.
!
!  lsolve can be 1, 2 or 3. See the routine min_sqr_solve for an
!  explanation of its meaning. 
!
USE linear_solvers,     ONLY : min_sqr_solve
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ndata, lsolve
REAL(DP), INTENT(IN) :: x(nvar,ndata), f(ndata)
TYPE(poly1), INTENT(INOUT) :: p_in

REAL(DP), ALLOCATABLE :: amat(:,:), coeff(:) 

INTEGER :: i, idata, ncoeff

IF (nvar<1) CALL errore('fit_multi_linear','nvar must be larger than 1',1)
IF (ndata < nvar+1) &
   CALL errore('fit_multi_linear','Be careful: there are too few &
                                                        &sampling data',1)
ncoeff=nvar+1
ALLOCATE(amat(ndata,ncoeff))
ALLOCATE(coeff(ncoeff))
!
!  prepare the auxiliary matrix
!
amat=0.0_DP

DO idata=1,ndata
   amat(idata,1) = 1.0_DP
   DO i=1,nvar
      amat(idata,1+i) = x(i,idata)
   ENDDO
ENDDO
!
CALL min_sqr_solve(ndata, ncoeff, amat, f, coeff, lsolve)
!
!   assign the coefficients to the polynomial 
!
p_in%a0=coeff(1)
DO i=1,nvar
   p_in%phi1(i)=coeff(1+i)
ENDDO

DEALLOCATE(amat)
DEALLOCATE(coeff)

RETURN
END SUBROUTINE fit_multi_linear
!
SUBROUTINE evaluate_fit_linear(nvar,x,f,p_in)
!
!  This routine evaluates the linear polynomial at the point x
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
REAL(DP), INTENT(IN) :: x(nvar)
TYPE(poly1), INTENT(IN) :: p_in
REAL(DP), INTENT(OUT) :: f

INTEGER  :: i

f=p_in%a0
DO i=1, nvar
   f=f+p_in%phi1(i)*x(i)
ENDDO

RETURN
END SUBROUTINE evaluate_fit_linear

SUBROUTINE set_linear_grad(nvar,f,p_in)
!
!   This routine gives the gradient of the linear polynomial.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
TYPE(poly1), INTENT(IN) :: p_in
REAL(DP), INTENT(INOUT) :: f(nvar)

f=p_in%phi1

RETURN
END SUBROUTINE set_linear_grad
!
SUBROUTINE print_linear_polynomial(nvar, p_in)
!
!  This routine writes on output the coefficients of the linear polynomial. 
!
IMPLICIT NONE

INTEGER, INTENT(IN)     :: nvar
TYPE(poly1), INTENT(IN) :: p_in

CHARACTER(LEN=6) :: int_to_char
INTEGER :: i

WRITE(stdout,'(/,5x,"Linear polynomial:")') 
!
!  term of 0 degree
!
WRITE(stdout,'(5x,e15.7)') p_in%a0
!
!  terms of 1 degree
!
DO i=1,nvar
   WRITE(stdout,'(5x," +",e15.7,a)') p_in%phi1(i), '  x'//TRIM(int_to_char(i))
ENDDO

RETURN
END SUBROUTINE print_linear_polynomial

SUBROUTINE introduce_linear_fit(nvar, ndata)
!
!   This routine writes a small message with the information on 
!   the fit with a linear polynomial.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ndata

  WRITE(stdout,'(/,5x,"Fitting the data with a linear polynomial:")')

  WRITE(stdout,'(/,5x,"Number of variables:",10x,i5)')  nvar
  WRITE(stdout,'(5x,"Coefficients of the linear polynomial:",i5)')  nvar+1
  WRITE(stdout,'(5x,"Number of fitting data:",7x,i5,/)')  ndata

RETURN
END SUBROUTINE introduce_linear_fit

SUBROUTINE print_chisq_linear(ndata, nvar, x, f, p_in)
!
!   This routine receives as input the values of a function f for ndata
!   values of the independent variables x, a set of coefficients
!   of a linear interpolating polynomial p_in and writes as output
!   the sum of the squares of the differences between the values of
!   the function and of the interpolating polynomial divided by the
!   number of data.
!
IMPLICIT NONE

INTEGER     :: ndata, nvar
REAL(DP)    :: x(nvar, ndata), f(ndata)
TYPE(poly1) :: p_in

REAL(DP) :: chisq, perc, aux
INTEGER  :: idata

chisq=0.0_DP
perc=0.0_DP
DO idata=1,ndata
   CALL evaluate_fit_linear(nvar,x(1,idata),aux,p_in)
!  WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
   IF (ABS(f(idata))>1.D-12) perc= perc + ABS((f(idata)-aux) / f(idata))
ENDDO
WRITE(stdout,'(5x,"chi square linear=",e18.5," relative error",e18.5,&
                                     &" %",/)') chisq/ndata, perc*100 / ndata
RETURN
END SUBROUTINE print_chisq_linear

END MODULE linear_surfaces
