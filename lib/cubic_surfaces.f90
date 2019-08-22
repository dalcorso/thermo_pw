!
! Copyright (C) 2019 A. Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE cubic_surfaces
!
!   this module contains the support routines for dealing with cubic
!   surfaces interpolation. 
!
!   It provides the following routines:
!
!   fit_multi_cubic receives as input the coordinates of some points and 
!   the values of a function in these points and finds the coefficients 
!   of the cubic polynomial that better fits the points.
!
!   evaluate_fit_cubic evaluates the cubic polynomial on a given input
!   point.
!
!   evaluate_cubic_grad evaluates the gradient of the cubic polynomial on 
!   a given input point.
!
!   evaluate_cubic_hessian evaluates the Hessian of the cubic polynomial
!   on a given input point.
!
!   find_cubic_extremum finds the extremum closest to the input point.
!
!   print_cubic_polynomial writes on output the coefficients of the
!   cubic polynomial.
!
!   introduce_cubic_fit writes a message with a few information on the
!   cubic polynomial and the number of data used to fit it.
!
!   print_chisq_cubic writes on output the chi square of a given cubic 
!   polynomial interpolation.
!
  USE kinds, ONLY : DP
  USE polynomial, ONLY : poly3
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: fit_multi_cubic, evaluate_fit_cubic, &
            evaluate_cubic_grad, evaluate_cubic_hessian, &
            find_cubic_extremum, &
            print_cubic_polynomial, introduce_cubic_fit, &
            print_chisq_cubic 

CONTAINS

SUBROUTINE fit_multi_cubic(ndata,nvar,lsolve,x,f,p3)
!
!  This routine receives as input a set of vectors x(nvar,ndata) and
!  function values f(ndata) and gives as output the coefficients of
!  a cubic interpolating polynomial p3. In input ndata is the number 
!  of data points and nvar is the number of independent variables. 
!
!  The coefficients are organized as follows:
!  a0 + \sum_i v1(i) * x(i,idata) 
!     + \sum_i,j (j>=i) phi2(n(i,j)) * x(i,idata) * x(j,idata)  
!     + \sum_i,j,k (j>=i) (k>=j) phi3(m(i,j,k)) x(i,idata)*x(j,idata)*x(k,idata)
!
!  n(i,j) and m(i,j,k) are linear indeces that are computed as:
!
!  n=0
!  m=0
!  DO i=1,nvar
!     DO j=i,nvar
!        n=n+1
!        DO k=j,nvar
!           m=m+1
!        ENDDO
!     ENDDO
!  ENDDO
!         
!  The total number of coefficients of a multivariate polynomial is 
!  (n+d)!/n!/d! where d is the degree of the polynomial (in the cubic
!  case d=3) and n is the number of variables. 
!
!  The following number of coefficients are necessary, depending on the
!  number of variables of the cubic polynomial.
!
!                          degree
!   number of variables   0   1   2   3  total number of coefficients
!          1              1   1   1   1             4	       
!          2              1   2   3   4             10
!          3              1   3   6  10             20                   
!          4              1   4  10  20             35
!          5              1   5  15  35             56
!          6              1   6  21  56             84
!          ...                 ...                 ...
!
!    To interpolate the function it is better to give to fit_multi_cubic 
!    a number of points equal or larger than the number of coefficients. 
!    The routine makes a least square fit of the data.
!
!    lsolve can be 1, 2 or 3. It chooses the method to compute the
!    polynomial coefficients. 
!    The problem can be written, in matrix form:
!    A X = B      where A is a matrix of dimension ndata x ncoeff
!                 Using 1 a matrix ncoeff x ncoeff is calculated as
!    A^T A X = A^T B     and the linear system provides X, the ncoeff 
!                 coefficients of the polynomial.
!                 Using 2 the overdetemined linear system AX=B is solved
!                 using QR or LQ factorization.
!                 Using 3 the overdetermined linear system AX=B is solved
!                 using SVD decomposition
!    If lsolve is not one of these values method 2 is used.
!
USE linear_solvers,     ONLY : linsolvx, linsolvx_sym, linsolvms, linsolvsvd
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ndata
INTEGER, INTENT(INOUT) :: lsolve
REAL(DP), INTENT(IN) :: x(nvar,ndata), f(ndata)
TYPE(poly3), INTENT(INOUT) :: p3

REAL(DP), ALLOCATABLE :: amat(:,:), aa(:,:), b(:), coeff(:) 

INTEGER :: ivar, jvar, idata, nv

INTEGER :: i, j, k, m, n, ncoeff

ncoeff=1+nvar+p3%ncoeff2+p3%ncoeff3

ALLOCATE(amat(ndata,ncoeff))
ALLOCATE(coeff(ncoeff))
ALLOCATE(aa(ncoeff,ncoeff))
ALLOCATE(b(ncoeff))

IF (nvar<1) CALL errore('fit_multi_cubic','nvar must be greater than 0',1)
IF (ndata < nvar) &
   WRITE(stdout,'(/,5x,"Be careful: there are too few sampling data")')
!
!  prepare the auxiliary matrix
!
amat=0.0_DP

DO idata=1,ndata
   amat(idata,1) = 1.0_DP
   n=0
   m=0
   DO i=1, nvar
      amat(idata,i+1)=x(i,idata)
      DO j=i,nvar
         n=n+1
         amat(idata,1+nvar+n)=x(i,idata)*x(j,idata)
         DO k=j,nvar
            m=m+1
            amat(idata,1+nvar+p3%ncoeff2+m)=x(i,idata)*x(j,idata)*x(k,idata)
         ENDDO
      ENDDO
   ENDDO
ENDDO

aa=0.0_DP
b =0.0_DP
DO ivar=1,ncoeff
   DO jvar=1,ncoeff
      DO idata=1,ndata
         aa(ivar,jvar)= aa(ivar,jvar) + amat(idata,ivar) * amat(idata,jvar)
      END DO
   END DO
   DO idata=1,ndata
      b(ivar) = b(ivar) + amat(idata,ivar) * f(idata)
   END DO
END DO
!
!   solve the linear system and find the coefficients
!
coeff=0.0_DP
IF (lsolve<1.OR.lsolve>3) lsolve=2
IF (lsolve==1) THEN
   WRITE(stdout,'(5x,"Finding the cubic polynomial using &
                                                   &nvar x nvar matrix")')  
   CALL linsolvx(aa,ncoeff,b,coeff)
!   CALL linsolvx_sym(aa,ncoeff,b,coeff)
ELSEIF(lsolve==2) THEN
   WRITE(stdout,'(5x,"Finding the cubic polynomial using &
                                                   &QR factorization")')  
   CALL linsolvms(amat,ndata,ncoeff,f,coeff)
ELSEIF(lsolve==3) THEN
   WRITE(stdout,'(5x,"Finding the cubic polynomial using SVD")')  
   CALL linsolvsvd(amat,ndata,ncoeff,f,coeff)
ENDIF

p3%a0=coeff(1)
n=0
m=0
DO i=1,nvar
   p3%phi1(i)=coeff(1+i)
   DO j=i,nvar
      n=n+1
      p3%phi2(n)=coeff(1+nvar+n)
      DO k=j,nvar
         m=m+1
         p3%phi3(m)=coeff(1+nvar+p3%ncoeff2+m)
      ENDDO 
   ENDDO
ENDDO

DEALLOCATE(amat)
DEALLOCATE(coeff)
DEALLOCATE(aa)
DEALLOCATE(b)

RETURN
END SUBROUTINE fit_multi_cubic

SUBROUTINE evaluate_fit_cubic(nvar,x,f,p3)
!
!  This routine evaluates the cubic polynomial at the point x
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
REAL(DP), INTENT(IN) :: x(nvar)
TYPE(poly3), INTENT(IN) :: p3
REAL(DP), INTENT(INOUT) :: f

INTEGER :: i, j, k, n, m
!
n=0
m=0
f = p3%a0
DO i=1,nvar
   f = f + p3%phi1(i)*x(i) 
   DO j=i, nvar
      n=n+1
      f=f+p3%phi2(n)*x(i)*x(j)
      DO k=j,nvar
         m=m+1
         f = f + p3%phi3(m)*x(i)*x(j)*x(k)
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE evaluate_fit_cubic

SUBROUTINE evaluate_cubic_grad(nvar,x,f,p3)
!
!  computes the gradient of the cubic polynomial at the point x.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
REAL(DP), INTENT(IN) :: x(nvar)
TYPE(poly3), INTENT(IN) :: p3
REAL(DP), INTENT(INOUT) :: f(nvar)

INTEGER :: i, j, k, m, n

n=0
m=0
f=p3%phi1
DO i=1, nvar
   DO j=i, nvar
      n=n+1
      f(i) = f(i) + p3%phi2(n) * x(j) 
      f(j) = f(j) + p3%phi2(n) * x(i) 
      DO k=j, nvar
         m=m+1
         f(i) = f(i) + p3%phi3(m) * x(j) * x(k)
         f(j) = f(j) + p3%phi3(m) * x(i) * x(k)
         f(k) = f(k) + p3%phi3(m) * x(i) * x(j)
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE evaluate_cubic_grad

SUBROUTINE evaluate_cubic_hessian(nvar,x,f,p3)
!
!  computes the Hessian of the cubic polynomial at the point x. 
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
REAL(DP), INTENT(IN) :: x(nvar)
TYPE(poly3), INTENT(IN) :: p3
REAL(DP), INTENT(INOUT) :: f(nvar,nvar)

INTEGER :: i, j, k, m, n

f=0.0_DP
n=0
m=0
DO i=1, nvar
   DO j=i, nvar
      n=n+1
      f(i,j) = f(i,j) + p3%phi2(n) 
      f(j,i) = f(i,j) + p3%phi2(n) 
      DO k=j+1, nvar
         m=m+1
         f(i,j) = f(i,j) + p3%phi3(m) * x(k)
         f(j,i) = f(j,i) + p3%phi3(m) * x(k)
         f(i,k) = f(i,k) + p3%phi3(m) * x(j)
         f(k,i) = f(k,i) + p3%phi3(m) * x(j)
         f(j,k) = f(j,k) + p3%phi3(m) * x(i)
         f(k,j) = f(k,j) + p3%phi3(m) * x(i)
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE evaluate_cubic_hessian

SUBROUTINE find_cubic_extremum(nvar,x,f,p3)
!
!  This routine starts from the point x and finds the extremum closest
!  to x. In output x are the coordinates of the extremum and f 
!  the value of the cubic function at the extremum
!
USE linear_solvers, ONLY : linsolvx
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
REAL(DP),INTENT(INOUT) :: x(nvar), f
TYPE(poly3), INTENT(IN) :: p3

INTEGER, PARAMETER :: maxiter=300

INTEGER :: iter, ideg
REAL(DP), PARAMETER :: tol=2.D-11
REAL(DP) :: g(nvar), y(nvar), xold(nvar)
REAL(DP) :: j(nvar, nvar) 
REAL(DP) :: deltax, fmod

xold(:)=x(:)
DO iter=1,maxiter
   !
   CALL evaluate_cubic_grad(nvar,x,g,p3)
   !
   CALL evaluate_cubic_hessian(nvar,x,j,p3)
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
CALL errore('find_cubic_extremum','extremum not found',1)
100 CONTINUE
CALL evaluate_fit_cubic(nvar,x,f,p3)

RETURN
END SUBROUTINE find_cubic_extremum

SUBROUTINE print_cubic_polynomial(nvar, p3)
!
!  This subroutine prints the coefficients of a cubic polynomial.
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nvar
TYPE(poly3), INTENT(IN) :: p3

INTEGER :: i, j, k, n, m
CHARACTER(LEN=6) :: int_to_char


WRITE(stdout,'(/,5x,"Cubic polynomial:")')

WRITE(stdout,'(5x,e15.7)') p3%a0

DO i=1,nvar
   WRITE(stdout,'(5x," +",e15.7,a)') p3%phi1(i), '  x'//TRIM(int_to_char(i))
ENDDO

n=0
DO i=1, nvar
   DO j=i,nvar
      n=n+1
      WRITE(stdout,'(5x," +",e15.7,a)') p3%phi2(n), &
                     '  x'//TRIM(int_to_char(i))//' x'//TRIM(int_to_char(j))
   ENDDO
ENDDO

m=0
DO i=1, nvar
   DO j=i,nvar
      DO k=j,nvar
         m=m+1
         WRITE(stdout,'(5x," +",e15.7,a)') p3%phi3(m), &
                  '  x'//TRIM(int_to_char(i))//' x'//TRIM(int_to_char(j))// &
                   ' x'//TRIM(int_to_char(k))
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE print_cubic_polynomial

SUBROUTINE introduce_cubic_fit(nvar, ncoeff, ndata)
!
!  This subroutine prints a few information on the cubic polynomial
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ncoeff, ndata

WRITE(stdout,'(/,5x,"Fitting the data with a cubic polynomial:")')

WRITE(stdout,'(/,5x,"Number of variables:",8x,i5)')  nvar
WRITE(stdout,'(5x,"Coefficients of the cubic polynomial:",2x,i5)')  ncoeff
WRITE(stdout,'(5x,"Number of fitting data:",5x,i5,/)')  ndata

RETURN
END SUBROUTINE introduce_cubic_fit

SUBROUTINE print_chisq_cubic(ndata, nvar, x, f, p3)
!
!   This routine receives as input the values of a function f for ndata
!   values of the independent variables x, a set of coefficients
!   of a cubic interpolating polynomial p3 and writes as output
!   the sum of the squares of the differences between the values of
!   the function and of the interpolating polynomial divided by the number
!   of data.
!
IMPLICIT NONE

INTEGER, INTENT(IN)  :: ndata, nvar
REAL(DP), INTENT(IN) :: x(nvar, ndata), f(ndata)
TYPE(poly3), INTENT(IN) :: p3

REAL(DP) :: chisq, perc, aux
INTEGER  :: idata

chisq=0.0_DP
perc=0.0_DP
DO idata=1,ndata
   CALL evaluate_fit_cubic(nvar,x(1,idata),aux,p3)
!   WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
   IF (ABS(f(idata))>1.D-12) perc= perc + ABS((f(idata)-aux) / f(idata))
ENDDO

WRITE(stdout,'(5x,"chi square cubic=",e18.5," relative error",e18.5,&
                                  &" %",/)') chisq / ndata, perc * 100 / ndata
RETURN
END SUBROUTINE print_chisq_cubic

END MODULE cubic_surfaces
