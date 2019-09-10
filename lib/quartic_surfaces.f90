!
! Copyright (C) 2019 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE quartic_surfaces
!
!   This module contains the support routines for dealing with quartic
!   surfaces interpolation. 
!
!   It provides the following routines:
!
!   fit_multi_quartic receives as input the coordinates of some
!   points and the value of a function in these points and
!   finds the coefficients of the quartic polynomial that better fits
!   the data.
!
!   evaluate_fit_quartic evaluates the quartic polynomial on an given input
!   point.
!
!   evaluate_quartic_grad evaluates the gradient of the quartic polynomial
!   on a given input point.
!
!   evaluate_quartic_hessian evaluates the Hessian of the quartic polynomial
!   on a given input point.
!
!   find_quartic_extremum finds the extremum of the quartic polynomial 
!   closest to a given input point.
!
!   print_quartic_polynomial writes on output the coefficients of the
!   cubic polynomial.
!
!   introduce_quartic_fit writes on output a few information
!   on the quartic polynomial and the number of data used for the fit. 
!
!   print_chisq_quartic writes on output the chi square of a given 
!   quartic polynomial interpolation.
!   
!   evaluate_two_quartic evaluates the sum of two quartic polynomials
!   on a given input point.
!  
!   find_two_quartic_extremum finds the extremum of the sum of two quartic
!   polynomial closest to a given input point.
!
!   print_chisq_two_quartic writes on output the chi square of the sum
!   of two quartic polynomial interpolation.
!
!   evaluate_quartic_linear evaluates the sum of a quartic and a linear 
!   polynomials on a given input point.
!
!   find_quartic_linear_extremum finds the extremum of the sum of
!   a quartic and a linear polynomials closest to a given input point.
!
!   evaluate_quartic_quadratic evaluates the sum of a quartic and a 
!   quadratic polynomial on a given input point.
!
!   find_quartic_quadratic_extremum finds the extremum of the sum of
!   a quartic and a quadratic polynomials closest to a given input point.
!   
!   print_chisq_quartic_quadratic writes on output the chi square of the sum
!   a quartic and quadratic polynomials interpolation.
!
!   evaluate_quartic_cubic evaluates the sum of a quartic and a cubic 
!   polynomials in a given input point.
!
!   find_quartic_cubic_extremum finds the extremum of the sum of
!   a quartic and a cubic polynomials closest to a given input point.
!   
!   print_chisq_quartic_cubic writes on output the chi square of the sum
!   of a quartic and a cubic polynomials interpolation
!
  USE kinds, ONLY : DP
  USE polynomial, ONLY : poly4
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: fit_multi_quartic, evaluate_fit_quartic,         &
            evaluate_quartic_grad, evaluate_quartic_hessian, &
            find_quartic_extremum,                           &
            print_quartic_polynomial, introduce_quartic_fit, &
            print_chisq_quartic,                             &
            evaluate_two_quartic,                            &
            find_two_quartic_extremum,                       &
            print_chisq_two_quartic,                         &
            evaluate_quartic_linear,                         &
            find_quartic_linear_extremum,                    &
            evaluate_quartic_quadratic,                      &
            find_quartic_quadratic_extremum,                 &
            print_chisq_quartic_quadratic,                   &
            evaluate_quartic_cubic,                          &
            find_quartic_cubic_extremum,                     &
            print_chisq_quartic_cubic                       

CONTAINS

SUBROUTINE fit_multi_quartic(ndata,nvar,lsolve,x,f,p4)
!
!  This routine receives as input a set of vectors x(nvar,ndata) and
!  function values f(ndata) and gives as output the coefficients of
!  a quartic interpolating polynomial p4.
!  In input ndata is the number of data points and nvar is the 
!  number of variables. 
!
!  The coefficients are organized as follows:
!  a0 + \sum_i v1(i) * x(i,idata) 
!     + \sum_i,j (j>=i) phi2(n(i,j)) * x(i,idata) * x(j,idata)  
!     + \sum_i,j,k (j>=i) (k>=j) phi3(m(i,j,k)) x(i,idata)*x(j,idata)*x(k,idata)
!     + \sum_i,j,k,l (j>=i) (k>=j) (l>=k) phi4(o(i,j,k,l)) *
!                             x(i,idata)*x(j,idata)*x(k,idata)*x(l,idata)
!
!  n(i,j), m(i,j,k), and o(i,j,k,l) are linear indeces that are computed as:
!
!  n=0
!  m=0
!  o=0
!  DO i=1,nvar
!     DO j=i,nvar
!        n=n+1
!        DO k=j,nvar
!           m=m+1
!           DO l=k,nvar
!              o=o+1
!           ENDDO
!        ENDDO
!     ENDDO
!  ENDDO
! 
!  The following number of coefficients are necessary depending on the
!  number of variables of the quartic polynomial.
!
!                                degree
!   number of variables   0   1    2    3    4 total number of coefficients
!
!          1              1   1    1    1    1         5            
!          2              1   2    3    4    5        15  
!          3              1   3    6   10   15        35                     
!          4              1   4   10   20   35        70  
!          5              1   5   15   35   70       126   
!          6              1   6   21   56  126       210   
!          ...                 ...                   ...
!
!  To interpolate the function it is better to give to fit_multi_quartic
!  a number of points equal or larger than the number of coefficients. 
!  The routine makes a least square fit of the data.
!
!  lsolve can be 1, 2 or 3. See the routine min_sqr_solve for an
!  explanation of its meaning. 
!
USE linear_solvers,     ONLY : min_sqr_solve
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ndata
INTEGER, INTENT(INOUT) :: lsolve
REAL(DP), INTENT(IN) :: x(nvar,ndata), f(ndata)
TYPE(poly4), INTENT(INOUT) :: p4
REAL(DP), ALLOCATABLE :: amat(:,:), coeff(:)

INTEGER :: i, j, k, l, m, n, o, idata, ncoeff

ncoeff=1+nvar+p4%ncoeff2+p4%ncoeff3+p4%ncoeff4

ALLOCATE(amat(ndata,ncoeff))
ALLOCATE(coeff(ncoeff))

IF (nvar<1) CALL errore('fit_multi_quartic','nvar must be larger than 0',1)
IF (ndata < ncoeff) &
   CALL errore('fit_multi_quartic','Too few sampling data',1)
!
!  prepare the auxiliary matrix
!
DO idata=1,ndata
   amat(idata,1) = 1.0_DP
   n=0
   m=0
   o=0
   DO i=1, nvar
      amat(idata,i+1)=x(i,idata)
      DO j=i,nvar
         n=n+1
         amat(idata,1+nvar+n)=x(i,idata)*x(j,idata)
         DO k=j,nvar
            m=m+1
            amat(idata,1+nvar+p4%ncoeff2+m)=x(i,idata)*x(j,idata)*x(k,idata)
            DO l=k, nvar
               o=o+1
               amat(idata,1+nvar+p4%ncoeff2+p4%ncoeff3+o)=&
                               x(i,idata)*x(j,idata)*x(k,idata)*x(l,idata)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

CALL min_sqr_solve(ndata, ncoeff, amat, f, coeff, lsolve)
!
!  assign the coefficients to the polynomial 
!
p4%a0=coeff(1)
n=0
m=0
o=0
DO i=1,nvar
   p4%phi1(i)=coeff(1+i)
   DO j=i,nvar
      n=n+1
      p4%phi2(n)=coeff(1+nvar+n)
      DO k=j,nvar
         m=m+1
         p4%phi3(m)=coeff(1+nvar+p4%ncoeff2+m)
         DO l=k,nvar
            o=o+1
            p4%phi4(o)=coeff(1+nvar+p4%ncoeff2+p4%ncoeff3+o)
         ENDDO 
      ENDDO
   ENDDO
ENDDO

DEALLOCATE(amat)
DEALLOCATE(coeff)

RETURN
END SUBROUTINE fit_multi_quartic

SUBROUTINE evaluate_fit_quartic(nvar,x,f,p4)
!
!  This routine evaluates the quartic polynomial at the point x.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
REAL(DP), INTENT(IN) :: x(nvar)
TYPE(poly4), INTENT(IN) :: p4
REAL(DP), INTENT(INOUT) :: f

INTEGER :: i, j, k, l, m, n, o

n=0
m=0
o=0
f=p4%a0
DO i=1,nvar
   f=f+p4%phi1(i)*x(i)
   DO j=i, nvar
      n=n+1
      f=f+p4%phi2(n)*x(i)*x(j)
      DO k=j,nvar
         m=m+1
         f=f+p4%phi3(m)*x(i)*x(j)*x(k)
         DO l=k,nvar
            o=o+1
            f=f+p4%phi4(o)*x(i)*x(j)*x(k)*x(l)
         ENDDO
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE evaluate_fit_quartic

SUBROUTINE evaluate_quartic_grad(nvar,x,f,p4)
!
!   This routine evaluates the gradient of the quartic polynomial at the 
!   point x.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar 
REAL(DP), INTENT(IN) :: x(nvar)
TYPE(poly4), INTENT(IN) :: p4
REAL(DP), INTENT(INOUT) :: f(nvar)

INTEGER :: i, j, k, l, n, m, o

n=0
m=0
o=0
f=p4%phi1
DO i=1, nvar
   DO j=i, nvar
      n=n+1
      f(i) = f(i) + p4%phi2(n) * x(j)  
      f(j) = f(j) + p4%phi2(n) * x(i)                     
      DO k=j, nvar
         m=m+1
         f(i) = f(i) + p4%phi3(m) * x(j) * x(k) 
         f(j) = f(j) + p4%phi3(m) * x(i) * x(k)    
         f(k) = f(k) + p4%phi3(m) * x(i) * x(j)    
         DO l=k, nvar
            o=o+1
            f(i)=f(i) + p4%phi4(o) * x(j) * x(k) * x(l)
            f(j)=f(j) + p4%phi4(o) * x(i) * x(k) * x(l)
            f(k)=f(k) + p4%phi4(o) * x(i) * x(j) * x(l)
            f(l)=f(l) + p4%phi4(o) * x(i) * x(j) * x(k)
         ENDDO
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE evaluate_quartic_grad

SUBROUTINE evaluate_quartic_hessian(nvar,x,f,p4)
!
!   This routine evaluates the hessian of the quartic polynomial at the 
!   point x.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
REAL(DP), INTENT(IN) :: x(nvar)
TYPE(poly4), INTENT(IN) :: p4
REAL(DP), INTENT(INOUT) :: f(nvar,nvar)

INTEGER :: i, j, k, l, m, n, o

n=0
m=0
o=0
f=0.0_DP
DO i=1, nvar
   DO j=i, nvar
      n=n+1
      f(i,j) = f(i,j) + p4%phi2(n) 
      f(j,i) = f(j,i) + p4%phi2(n) 
      DO k=j, nvar
         m=m+1
         f(i,j) = f(i,j) + p4%phi3(m) * x(k) 
         f(j,i) = f(j,i) + p4%phi3(m) * x(k) 
         f(i,k) = f(i,k) + p4%phi3(m) * x(j) 
         f(k,i) = f(k,i) + p4%phi3(m) * x(j) 
         f(j,k) = f(j,k) + p4%phi3(m) * x(i) 
         f(k,j) = f(k,j) + p4%phi3(m) * x(i)
         DO l=k, nvar
            o=o+1
            f(i,j) = f(i,j) + p4%phi4(o) * x(k) * x(l)
            f(j,i) = f(j,i) + p4%phi4(o) * x(k) * x(l)
            f(i,k) = f(i,k) + p4%phi4(o) * x(j) * x(l)
            f(k,i) = f(k,i) + p4%phi4(o) * x(j) * x(l)
            f(i,l) = f(i,l) + p4%phi4(o) * x(j) * x(k)
            f(l,i) = f(l,i) + p4%phi4(o) * x(j) * x(k)
            f(j,k) = f(j,k) + p4%phi4(o) * x(i) * x(l)
            f(k,j) = f(k,j) + p4%phi4(o) * x(i) * x(l)
            f(j,l) = f(j,l) + p4%phi4(o) * x(i) * x(k)
            f(l,j) = f(l,j) + p4%phi4(o) * x(i) * x(k)
            f(k,l) = f(k,l) + p4%phi4(o) * x(i) * x(j)
            f(l,k) = f(l,k) + p4%phi4(o) * x(i) * x(j)
         ENDDO
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE evaluate_quartic_hessian

SUBROUTINE find_quartic_extremum(nvar,x,f,p4)
!
!  This routine finds the extremum of the quartic polynomial p4 closest
!  to the input point x. The mimumum is written in x and the value of the
!  polynomial at the minimum in f.
!
USE linear_solvers, ONLY : linsolvx
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
REAL(DP), INTENT(INOUT) :: x(nvar), f
TYPE(poly4), INTENT(IN) :: p4

INTEGER, PARAMETER :: maxiter=300

INTEGER :: iter, ideg
REAL(DP), PARAMETER :: tol=2.D-11
REAL(DP) :: g(nvar), y(nvar), xold(nvar)
REAL(DP) :: j(nvar, nvar) 
REAL(DP) :: deltax, fmod

xold(:)=x(:)
DO iter=1,maxiter
   !
   CALL evaluate_quartic_grad(nvar,x,g,p4)
   !
   CALL evaluate_quartic_hessian(nvar,x,j,p4)
   !
   CALL linsolvx(j, nvar, g, y)
   !
   !  Use Newton's method to find the zero of the gradient
   !
   x(:)= x(:) - y(:)
   fmod=0.0_DP
   deltax=0.0_DP
   DO ideg=1,nvar
      fmod = fmod + g(ideg)**2
      deltax = deltax + (xold(ideg)-x(ideg))**2
   END DO
   !
!   WRITE(stdout,'(i5,2f20.12)') iter, SQRT(deltax), SQRT(fmod)
   IF (SQRT(fmod) < tol .OR. SQRT(deltax) < tol ) GOTO 100
   xold(:)=x(:)
   !
END DO
CALL errore('find_quartic_extremum','extremum not found',1)
100 CONTINUE
!
CALL evaluate_fit_quartic(nvar,x,f,p4)
!
RETURN
END SUBROUTINE find_quartic_extremum
!
SUBROUTINE print_quartic_polynomial(nvar, p4)
!
!  This subroutine writes on output the coefficients of a quartic
!  polynomial p4. 
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nvar
TYPE(poly4), INTENT(IN) ::  p4

INTEGER :: i, j, k, l, n, m, o
CHARACTER(LEN=6) :: int_to_char

WRITE(stdout,'(/,5x,"Quartic polynomial:")')
!
!  term of 0 degree
!
WRITE(stdout,'(5x,e15.7)') p4%a0
!
!  terms of degree 1
!
DO i=1,nvar
   WRITE(stdout,'(5x," +",e15.7,a)') p4%phi1(i), '  x'//TRIM(int_to_char(i))
ENDDO
!
!  terms of degree 2
!
n=0
DO i=1, nvar
   DO j=i,nvar
      n=n+1
      WRITE(stdout,'(5x," +",e15.7,a)') p4%phi2(n), &
               '  x'//TRIM(int_to_char(i))//' x'//TRIM(int_to_char(j))
   ENDDO
ENDDO
!
!  terms of degree 3
!
m=0
DO i=1, nvar
   DO j=i,nvar
      DO k=j,nvar
         m=m+1
         WRITE(stdout,'(5x," +",e15.7,a)') p4%phi3(m), &
                  '  x'//TRIM(int_to_char(i))//' x'//TRIM(int_to_char(j))// &
                   ' x'//TRIM(int_to_char(k))
      ENDDO
   ENDDO
ENDDO
!
!   terms of degree 4
!
o=0
DO i=1, nvar
   DO j=i,nvar
      DO k=j,nvar
         DO l=k,nvar
            o=o+1
            WRITE(stdout,'(5x," +",e15.7,a)') p4%phi4(o), &
                  '  x'//TRIM(int_to_char(i))//' x'//TRIM(int_to_char(j))// &
                   ' x'//TRIM(int_to_char(k))//' x'//TRIM(int_to_char(l))
         ENDDO
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE print_quartic_polynomial

SUBROUTINE introduce_quartic_fit(nvar, ncoeff, ndata)
!
!   This routine writes a small message with the information on the 
!   quartic polynomial fit.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ncoeff, ndata

WRITE(stdout,'(/,5x,"Fitting the data with a quartic polynomial:")')

WRITE(stdout,'(/,5x,"Number of variables:",10x,i5)') nvar
WRITE(stdout,'(5x,"Coefficients of the quartic polynomial:",2x,i5)')  ncoeff
WRITE(stdout,'(5x,"Number of fitting data:",7x,i5,/)') ndata

RETURN
END SUBROUTINE introduce_quartic_fit
!
SUBROUTINE print_chisq_quartic(ndata, nvar, x, f, p4)
!
!   This routine receives as input the values of a function f for ndata
!   values of the independent variables x, a set of ncoeff coefficients
!   of a quadratic interpolating polynomial and writes as output
!   the sum of the squares of the differences between the values of
!   the function and of the interpolating polynomial divided by the number
!   of data.
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: ndata, nvar
REAL(DP), INTENT(IN) :: x(nvar, ndata), f(ndata)
TYPE(poly4), INTENT(IN) :: p4

REAL(DP) :: chisq, perc, aux
INTEGER  :: idata

chisq=0.0_DP
perc=0.0_DP
DO idata=1,ndata
   CALL evaluate_fit_quartic(nvar,x(1,idata),aux,p4)
!     WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
   IF (ABS(f(idata))>1.D-12) perc= perc + ABS((f(idata)-aux) / f(idata))
ENDDO

WRITE(stdout,'(5x,"chi square quartic=",e18.5," relative error",e18.5,&
                                     &" %",/)') chisq/ndata, perc * 100 / ndata
RETURN
END SUBROUTINE print_chisq_quartic
!
SUBROUTINE evaluate_two_quartic(nvar,x,f,p4,p41)
!
!  This routine evaluates the sum of two quartic polynomials at the point x
!
USE polynomial, ONLY : init_poly, clean_poly
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
REAL(DP),INTENT(INOUT) :: x(nvar), f
TYPE(poly4), INTENT(IN) :: p4, p41

TYPE(poly4) :: p4s

CALL init_poly(nvar,p4s)

p4s%a0=p4%a0+p41%a0
p4s%phi1=p4%phi1+p41%phi1
p4s%phi2=p4%phi2+p41%phi2
p4s%phi3=p4%phi3+p41%phi3
p4s%phi4=p4%phi4+p41%phi4

CALL evaluate_fit_quartic(nvar,x,f,p4s)

CALL clean_poly(p4s)

RETURN
END SUBROUTINE evaluate_two_quartic
!
SUBROUTINE find_two_quartic_extremum(nvar,x,f,p4,p41)
!
!  This routine finds the extremum of the sum of two quartic polynomials 
!  p4 and p41 closest to the input point x. The extremum is written in x 
!  and the value of the sum of the two polynomials at the extremum in f.
!
USE polynomial, ONLY : init_poly, clean_poly
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
REAL(DP),INTENT(INOUT) :: x(nvar), f
TYPE(poly4),INTENT(IN) :: p4, p41

TYPE(poly4) :: p4s

CALL init_poly(nvar,p4s)

p4s%a0=p4%a0+p41%a0
p4s%phi1=p4%phi1+p41%phi1
p4s%phi2=p4%phi2+p41%phi2
p4s%phi3=p4%phi3+p41%phi3
p4s%phi4=p4%phi4+p41%phi4

CALL find_quartic_extremum(nvar,x,f,p4s)

CALL clean_poly(p4s)

RETURN
END SUBROUTINE find_two_quartic_extremum
!
SUBROUTINE print_chisq_two_quartic(ndata,nvar,x,f,p4,p41)
!
!  This routine writes on output the chi squared of the sum of two 
!  quartic polynomials that interpolate the function f in the ndata
!  points x.
!
IMPLICIT NONE
INTEGER  :: ndata, nvar
REAL(DP) :: x(nvar, ndata), f(ndata)
TYPE(poly4), INTENT(IN) :: p4, p41

REAL(DP) :: chisq, perc, aux
INTEGER  :: idata

chisq=0.0_DP
perc=0.0_DP
DO idata=1,ndata
   CALL evaluate_two_quartic(nvar,x(1,idata),aux,p4,p41)
!  WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
   IF (ABS(f(idata))>1.D-12) perc= perc + ABS((f(idata)-aux) / f(idata))
ENDDO
WRITE(stdout,'(5x,"chi square quartic=",e18.5," relative error",e18.5,&
                                     &" %",/)') chisq/ndata, perc * 100 / ndata
RETURN
END SUBROUTINE print_chisq_two_quartic
!
SUBROUTINE evaluate_quartic_linear(nvar,x,f,p4,p1) 
!
!  This routine evaluates the sum of a quartic and a linear polynomials 
!  at the point x.
!
USE polynomial, ONLY : poly1, init_poly, clean_poly
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
TYPE(poly4), INTENT(IN) :: p4
TYPE(poly1), INTENT(IN) :: p1
REAL(DP), INTENT(INOUT) :: x(nvar), f

TYPE(poly4) :: p4s

CALL init_poly(nvar,p4s)

p4s%a0=p4%a0+p1%a0
p4s%phi1=p4%phi1+p1%phi1
p4s%phi2=p4%phi2
p4s%phi3=p4%phi3
p4s%phi4=p4%phi4

CALL evaluate_fit_quartic(nvar,x,f,p4s)

CALL clean_poly(p4s)

RETURN
END SUBROUTINE evaluate_quartic_linear

SUBROUTINE find_quartic_linear_extremum(nvar,x,f,p4,p1) 
!
!  This routine finds the extremum of the sum of a quartic and a linear
!  polynomials p4 and p1 closest to the input point x. The extremum is 
!  written in x and the value of the sum of the two polynomials at 
!  the extremum in f.
!
USE polynomial, ONLY : poly1, init_poly, clean_poly
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
TYPE(poly4), INTENT(IN) :: p4
TYPE(poly1), INTENT(IN) :: p1
REAL(DP), INTENT(INOUT) :: x(nvar), f

TYPE(poly4) :: p4s

CALL init_poly(nvar,p4s)

p4s%a0=p4%a0+p1%a0
p4s%phi1=p4%phi1+p1%phi1
p4s%phi2=p4%phi2
p4s%phi3=p4%phi3
p4s%phi4=p4%phi4

CALL find_quartic_extremum(nvar,x,f,p4s)

CALL clean_poly(p4s)

RETURN
END SUBROUTINE find_quartic_linear_extremum
!
SUBROUTINE evaluate_quartic_quadratic(nvar,x,f,p4,p2)
!
!  This routine evaluates the sum of a quartic and a quadratic polynomials 
!  at the point x.
!
USE polynomial, ONLY : poly2, init_poly, clean_poly
IMPLICIT NONE

INTEGER, INTENT(IN) :: nvar
REAL(DP), INTENT(INOUT) :: x(nvar), f
TYPE(poly4), INTENT(IN) :: p4
TYPE(poly2), INTENT(IN) :: p2

TYPE(poly4) :: p4s

CALL init_poly(nvar,p4s)

p4s%a0=p4%a0+p2%a0
p4s%phi1=p4%phi1+p2%phi1
p4s%phi2=p4%phi2+p2%phi2
p4s%phi3=p4%phi3
p4s%phi4=p4%phi4

CALL evaluate_fit_quartic(nvar,x,f,p4s)

CALL clean_poly(p4s)

RETURN
END SUBROUTINE evaluate_quartic_quadratic
!
SUBROUTINE find_quartic_quadratic_extremum(nvar,x,f,p4,p2)
!
!  This routine finds the extremum of the sum of a quartic and a quadratic
!  polynomials p4 and p2 closest to the input point x. The extremum is 
!  written in x and the value of the sum of the two polynomials at 
!  the extremum in f.
!
USE polynomial, ONLY : poly2, init_poly, clean_poly
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
REAL(DP), INTENT(INOUT) :: x(nvar), f

TYPE(poly4), INTENT(IN) :: p4 
TYPE(poly2), INTENT(IN) :: p2

TYPE(poly4) :: p4s

CALL init_poly(nvar,p4s)

p4s%a0=p4%a0+p2%a0
p4s%phi1=p4%phi1+p2%phi1
p4s%phi2=p4%phi2+p2%phi2
p4s%phi3=p4%phi3
p4s%phi4=p4%phi4

CALL find_quartic_extremum(nvar,x,f,p4s)

CALL clean_poly(p4s)
RETURN
END SUBROUTINE find_quartic_quadratic_extremum
!
SUBROUTINE print_chisq_quartic_quadratic(ndata, nvar,x,f,p4,p2)
!
!  This routine writes on output the chi square of the sum of a
!  quadratic and a quartic polynomials that interpolate the function f in 
!  the ndata points x.
!
USE polynomial, ONLY : poly2
IMPLICIT NONE
INTEGER  :: ndata, nvar
REAL(DP) :: x(nvar, ndata), f(ndata)
TYPE(poly4), INTENT(IN) :: p4
TYPE(poly2), INTENT(IN) :: p2

REAL(DP) :: chisq, perc, aux
INTEGER  :: idata

chisq=0.0_DP
perc=0.0_DP
DO idata=1,ndata
   CALL evaluate_quartic_quadratic(nvar, x(1,idata), aux, p4, p2)
!  WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
   IF (ABS(f(idata))>1.D-12) perc= perc + ABS((f(idata)-aux) / f(idata))
ENDDO
WRITE(stdout,'(5x,"chi square quartic=",e18.5," relative error",e18.5,&
                                     &" %",/)') chisq/ndata, perc * 100 / ndata
RETURN
END SUBROUTINE print_chisq_quartic_quadratic
!
SUBROUTINE evaluate_quartic_cubic(nvar,x,f,p4,p3) 
!
!  This routine evaluates the sum of a quartic and a cubic polynomials 
!  at the point x.
!
USE polynomial, ONLY : poly3, init_poly, clean_poly
IMPLICIT NONE

INTEGER, INTENT(IN) :: nvar
TYPE(poly4), INTENT(IN) :: p4
TYPE(poly3), INTENT(IN) :: p3
REAL(DP), INTENT(INOUT) :: x(nvar), f

TYPE(poly4) :: p4s

CALL init_poly(nvar,p4s)

p4s%a0=p4%a0+p3%a0
p4s%phi1=p4%phi1+p3%phi1
p4s%phi2=p4%phi2+p3%phi2
p4s%phi3=p4%phi3+p3%phi3
p4s%phi4=p4%phi4

CALL evaluate_fit_quartic(nvar,x,f,p4s)

CALL clean_poly(p4s)

RETURN
END SUBROUTINE evaluate_quartic_cubic
!
SUBROUTINE find_quartic_cubic_extremum(nvar,x,f,p4,p3)
!
!  This routine finds the extremum of the sum of a quartic and a cubic
!  polynomials p4 and p3 closest to the input point x. The extremum is 
!  written in x and the value of the sum of the two polynomials at 
!  the extremum in f.
!
USE polynomial, ONLY : poly3, init_poly, clean_poly
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
TYPE(poly4), INTENT(IN) :: p4
TYPE(poly3), INTENT(IN) :: p3
REAL(DP), INTENT(INOUT) :: x(nvar), f

TYPE(poly4) :: p4s

CALL init_poly(nvar,p4s)

p4s%a0=p4%a0+p3%a0
p4s%phi1=p4%phi1+p3%phi1
p4s%phi2=p4%phi2+p3%phi2
p4s%phi3=p4%phi3+p3%phi3
p4s%phi4=p4%phi4

CALL find_quartic_extremum(nvar,x,f,p4s)

CALL clean_poly(p4s)

RETURN
END SUBROUTINE find_quartic_cubic_extremum
!
SUBROUTINE print_chisq_quartic_cubic(ndata,nvar,x,f,p4,p3)
!
!  This routine writes on output the chi square of the sum of a
!  cubic and a quartic polynomials that interpolate the function f in 
!  the ndata points x.
!
USE polynomial, ONLY : poly3
IMPLICIT NONE
INTEGER, INTENT(IN) :: ndata, nvar
TYPE(poly4), INTENT(IN) :: p4
TYPE(poly3), INTENT(IN) :: p3
REAL(DP), INTENT(INOUT) :: x(nvar,ndata), f(ndata)

REAL(DP) :: chisq, perc, aux
INTEGER  :: idata

chisq=0.0_DP
perc=0.0_DP
DO idata=1,ndata
   CALL evaluate_quartic_cubic(nvar,x(1,idata),aux,p4,p3)
!  WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
   IF (ABS(f(idata))>1.D-12) perc= perc + ABS((f(idata)-aux) / f(idata))
ENDDO
WRITE(stdout,'(5x,"chi square quartic=",e18.5," relative error",e18.5,&
                                     &" %",/)') chisq/ndata, perc * 100 / ndata
RETURN
END SUBROUTINE print_chisq_quartic_cubic

END MODULE quartic_surfaces
