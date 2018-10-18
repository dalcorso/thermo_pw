!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE quadratic_surfaces
!
!   this module contains the support routines for dealing with quadratic
!   surfaces interpolation up to dimension 6. It is used to interpolate 
!   the total energy as a function of the celldm parameters. 
!   It provides basically two routines:
!
!   fit_multi_quadratic that receives as input the coordinates of some
!   points and the value of the quadratic function in these points and
!   finds the coefficients of the quadratic function that passes through
!   the points.
!
!   The following number of coefficients are necessary, depending on the
!   dimension of the quadratic function
!
!   dimension     number of coefficients
!      1                    3
!      2                    6
!      3                   10 
!      4                   15
!      5                   21
!      6                   28
!
!    To interpolate the function it is better to give to
!    multi_quadratic a number of points equal or larger than
!    to the number of coefficients. The routines makes a least square fit
!    of the data.
!   
!   find_fit_extremum receives as input the coefficients of the quadratic 
!   functions and finds the position of the extremum
!
!   find_two_fit_extremum receives as input the coefficients of two quadratic 
!   functions and finds the position of the extremum of their sum.
!
!   two additional routines are auxiliary.
!   evaluate_fit_quadratic, given the coordinates of a point, and the
!   coefficients of the quadratic function, evaluates the quadratic function 
!   on that point.
!
!   write_fit_hessian, given the coefficient of the quadratic function,
!   write them on output as the Hessian of the function. It also
!   diagonalizes the Hessian and find the eigenvalues and eigenvectors.
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: evaluate_fit_quadratic, write_fit_hessian, &
            fit_multi_quadratic, find_fit_extremum,  &
            find_two_fit_extremum, &
            evaluate_fit_grad_quadratic, polifit, write_poli, &
            print_quadratic_polynomial, summarize_fitting_data, &
            write_vector, introduce_quadratic_fit, print_chisq_quadratic

CONTAINS

SUBROUTINE fit_multi_quadratic(ndata,degree,nvar,x,f,coeff)
!
!  This routine receives as input a set of vectors x(degree,ndata) and
!  function values f(ndata) and gives as output the coefficients of
!  a quadratic interpolating polynomial coeff(nvar). In input
!  ndata is the number of data points. degree is the number of degrees of
!  freedom or equivalently the number of independent parameters (the
!  maximum is 6), and nvar is the number of coefficients of the
!  intepolating quadrating polynomial. nvar= (degree+1) * (degree+2) / 2. 
!  The coefficients are organized like this:
!  a_1 + a_2  x(1,i) + a_3  x(1,i)**2                                
!      + a_4  x(2,i) + a_5  x(2,i)**2 + a_6  x(1,i)*x(2,i) +         2
!      + a_7  x(3,i) + a_8  x(3,i)**2 + a_9  x(1,i)*x(3,i) 
!                                     + a_10 x(2,i)*x(3,i) +         3
!      + a_11 x(4,i) + a_12 x(4,i)**2 + a_13 x(1,i)*x(4,i) +
!                                     + a_14 x(2,i)*x(4,i) +
!                                     + a_15 x(3,i)*x(4,i) +         4
!      + a_16 x(5,i) + a_17 x(5,i)**2 + a_18 x(1,i)*x(5,i) +
!                                     + a_19 x(2,i)*x(5,i) +
!                                     + a_20 x(3,i)*x(5,i) +
!                                     + a_21 x(4,i)*x(5,i) +         5
!      + a_22 x(6,i) + a_23 x(6,i)**2 + a_24 x(1,i)*x(6,i) +
!                                     + a_25 x(2,i)*x(6,i) +
!                                     + a_26 x(3,i)*x(6,i) +
!                                     + a_27 x(4,i)*x(6,i) +
!                                     + a_28 x(5,i)*x(6,i)           6
!
!
USE kinds, ONLY : DP
USE linear_solvers, ONLY : linsolvx
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar, ndata
REAL(DP), INTENT(IN) :: x(degree,ndata), f(ndata)
REAL(DP), INTENT(OUT) :: coeff(nvar)

REAL(DP) :: amat(ndata,nvar), aa(nvar,nvar), b(nvar) 

INTEGER :: ivar, jvar, idata

IF (nvar /= (degree+1)*(degree + 2) / 2) &
   CALL errore('multi_quadratic','nvar and degree not compatible',1)
IF (degree>6.OR.degree<1) &
   CALL errore('multi_quadratic','degree must be from 1 to 6',1)
IF (ndata < 3) &
   CALL errore('multi_quadratic','Too few sampling data',1)
IF (ndata < nvar) &
   WRITE(stdout,'(/,5x,"Be careful: there are too few sampling data")')

!
!  prepare the auxiliary matrix
!
amat=0.0_DP

DO idata=1,ndata
   amat(idata,1) = 1.0_DP
   amat(idata,2) = x(1,idata)
   amat(idata,3) = x(1,idata)**2
   IF (degree>1) THEN
      amat(idata,4) = x(2,idata)
      amat(idata,5) = x(2,idata)**2
      amat(idata,6) = x(1,idata) * x(2,idata)
   ENDIF
   IF (degree>2) THEN
      amat(idata,7) = x(3,idata)
      amat(idata,8) = x(3,idata)**2
      amat(idata,9) = x(1,idata) * x(3,idata)
      amat(idata,10)= x(2,idata) * x(3,idata)
   ENDIF  
   IF (degree>3) THEN
      amat(idata,11) = x(4,idata)
      amat(idata,12) = x(4,idata)**2
      amat(idata,13) = x(1,idata) * x(4,idata)
      amat(idata,14) = x(2,idata) * x(4,idata)
      amat(idata,15) = x(3,idata) * x(4,idata)
   ENDIF
   IF (degree>4) THEN
      amat(idata,16) = x(5,idata)
      amat(idata,17) = x(5,idata)**2
      amat(idata,18) = x(1,idata) * x(5,idata)
      amat(idata,19) = x(2,idata) * x(5,idata)
      amat(idata,20) = x(3,idata) * x(5,idata)
      amat(idata,21) = x(4,idata) * x(5,idata)
   ENDIF
   IF (degree>5) THEN
      amat(idata,22) = x(6,idata)
      amat(idata,23) = x(6,idata)**2
      amat(idata,24) = x(1,idata) * x(6,idata)
      amat(idata,25) = x(2,idata) * x(6,idata)
      amat(idata,26) = x(3,idata) * x(6,idata)
      amat(idata,27) = x(4,idata) * x(6,idata)
      amat(idata,28) = x(5,idata) * x(6,idata)
   ENDIF
ENDDO
!
!   Now set the linear system
!
aa=0.0_DP
b =0.0_DP
DO ivar=1,nvar
   DO jvar=1,nvar 
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
CALL linsolvx(aa,nvar,b,coeff)

RETURN
END SUBROUTINE fit_multi_quadratic

SUBROUTINE find_fit_extremum(degree,nvar,x,f,coeff)
USE kinds, ONLY : DP
USE linear_solvers, ONLY : linsolvx
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: coeff(nvar)
REAL(DP), INTENT(OUT) :: x(degree), f

REAL(DP) :: amat(degree, degree), v(degree)

IF (degree > 6 .OR. degree < 1 ) &
   CALL errore('find_fit_extremum','incorrect degree',1)
IF (nvar /= (degree + 1)*(degree + 2) / 2) &
   CALL errore('find_fit_extremum','nvar and degree not compatible',1)

amat=0.0_DP
v=0.0_DP

v(1)= - coeff(2)
amat(1,1) = 2.0_DP * coeff(3)

IF (degree > 1) THEN
   v(2) = - coeff(4)
   amat(1,2) = coeff(6)
   amat(2,1) = coeff(6)
   amat(2,2) = 2.0_DP * coeff(5)
ENDIF

IF (degree > 2) THEN
   v(3) = - coeff(7)
   amat(1,3) = coeff(9)
   amat(3,1) = amat(1,3)
   amat(2,3) = coeff(10)
   amat(3,2) = amat(2,3)
   amat(3,3) = 2.0_DP * coeff(8)
ENDIF

IF (degree > 3) THEN
   v(4) = - coeff(11)
   amat(1,4) = coeff(13)
   amat(4,1) = amat(1,4)
   amat(2,4) = coeff(14)
   amat(4,2) = amat(2,4)
   amat(3,4) = coeff(15)
   amat(4,3) = amat(3,4)
   amat(4,4) = 2.0_DP * coeff(12)
ENDIF

IF (degree > 4) THEN
   v(5) = - coeff(16)
   amat(1,5) = coeff(18)
   amat(5,1) = amat(1,5)
   amat(2,5) = coeff(19)
   amat(5,2) = amat(2,5)
   amat(3,5) = coeff(20)
   amat(5,3) = amat(3,5)
   amat(4,5) = coeff(21)
   amat(5,4) = amat(4,5)
   amat(5,5) = 2.0_DP * coeff(17)
ENDIF

IF (degree > 5) THEN
   v(6) = - coeff(22)
   amat(1,6) = coeff(24)
   amat(6,1) = amat(1,6)
   amat(2,6) = coeff(25)
   amat(6,2) = amat(2,6)
   amat(3,6) = coeff(26)
   amat(6,3) = amat(3,6)
   amat(4,6) = coeff(27)
   amat(6,4) = amat(4,6)
   amat(5,6) = coeff(28)
   amat(6,5) = amat(5,6)
   amat(6,6) = 2.0_DP * coeff(23)
ENDIF

CALL linsolvx(amat,degree,v,x)

CALL evaluate_fit_quadratic(degree,nvar,x,f,coeff)

RETURN
END SUBROUTINE find_fit_extremum

SUBROUTINE find_two_fit_extremum(degree,nvar,x,f,coeff,coeff1)
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: coeff(nvar), coeff1(nvar)
REAL(DP), INTENT(OUT) :: x(degree), f

REAL(DP) :: coeff2(nvar)

coeff2(:) = coeff(:) + coeff1(:)

CALL find_fit_extremum(degree,nvar,x,f,coeff2)

RETURN
END SUBROUTINE find_two_fit_extremum

SUBROUTINE evaluate_fit_quadratic(degree,nvar,x,f,coeff)

USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: x(degree)
REAL(DP), INTENT(IN) :: coeff(nvar)
REAL(DP), INTENT(OUT) :: f

REAL(DP) :: aux

aux=coeff(1) + coeff(2) * x(1) + coeff(3) * x(1)**2
IF (degree>1) THEN
   aux = aux + coeff(4) * x(2) + coeff(5) * x(2)**2 + coeff(6) * x(1) * x(2)
ENDIF
IF (degree>2) THEN
   aux = aux + coeff(7)*x(3)+coeff(8)*x(3)**2+coeff(9)  * x(1) * x(3) + &
                                              coeff(10) * x(2) * x(3)
ENDIF
IF (degree>3) THEN
   aux = aux + coeff(11)*x(4)+coeff(12)*x(4)**2+coeff(13)*x(1)*x(4) + &
                                                coeff(14)*x(2)*x(4) + &
                                                coeff(15)*x(3)*x(4) 
ENDIF
IF (degree>4) THEN
   aux = aux + coeff(16)*x(5)+coeff(17)*x(5)**2+coeff(18)*x(1)*x(5) + &
                                                coeff(19)*x(2)*x(5) + &
                                                coeff(20)*x(3)*x(5) + &
                                                coeff(21)*x(4)*x(5) 
ENDIF
IF (degree>5) THEN
   aux = aux + coeff(22)*x(5)+coeff(23)*x(5)**2+coeff(24)*x(1)*x(6) + &
                                                coeff(25)*x(2)*x(6) + &
                                                coeff(26)*x(3)*x(6) + &
                                                coeff(27)*x(4)*x(6) + &
                                                coeff(28)*x(5)*x(6) 
ENDIF
f=aux

RETURN
END SUBROUTINE evaluate_fit_quadratic

SUBROUTINE evaluate_fit_grad_quadratic(degree,nvar,x,f,coeff)

USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: x(degree)
REAL(DP), INTENT(IN) :: coeff(nvar)
REAL(DP), INTENT(OUT) :: f(degree)

REAL(DP) :: aux(degree)

aux(1) = coeff(2) + 2.0_DP * coeff(3) * x(1) 

IF (degree>1) THEN
   aux(1) = aux(1) + coeff(6) * x(2) 
   aux(2) = coeff(4) + 2.0_DP * coeff(5) * x(2) + coeff(6) * x(1)
ENDIF

IF (degree>2) THEN
   aux(1) = aux(1) + coeff(9) * x(3)
   aux(2) = aux(2) + coeff(10) * x(3)
   aux(3) = coeff(7) + 2.0_DP * coeff(8) * x(3) + coeff(9)  * x(1)    & 
                                                + coeff(10) * x(2)
ENDIF

IF (degree>3) THEN
   aux(1) = aux(1) + coeff(13) * x(4) 
   aux(2) = aux(2) + coeff(14) * x(4) 
   aux(3) = aux(3) + coeff(15) * x(4) 
   aux(4) = coeff(11) + 2.0_DP * coeff(12) * x(4) + coeff(13) * x(1) &
                                                  + coeff(14) * x(2) &
                                                  + coeff(15) * x(3)
ENDIF
IF (degree>4) THEN
   aux(1) = aux(1) + coeff(18) * x(5)
   aux(2) = aux(2) + coeff(19) * x(5)
   aux(3) = aux(3) + coeff(20) * x(5)
   aux(4) = aux(4) + coeff(21) * x(5)
   aux(5) = coeff(16) + 2.0_DP * coeff(17) * x(5) + coeff(18) * x(1) &
                                                  + coeff(19) * x(2) &
                                                  + coeff(20) * x(3) &
                                                  + coeff(21) * x(4) 
ENDIF
IF (degree>5) THEN
   aux(1) = aux(1) + coeff(24) * x(6)
   aux(2) = aux(2) + coeff(25) * x(6)
   aux(3) = aux(3) + coeff(26) * x(6)
   aux(4) = aux(4) + coeff(27) * x(6)
   aux(5) = aux(5) + coeff(28) * x(6)
   aux(6) = coeff(22) + 2.0_DP * coeff(23) * x(6) + coeff(24) * x(1) &
                                                  + coeff(25) * x(2) &
                                                  + coeff(26) * x(3) &
                                                  + coeff(27) * x(4) &
                                                  + coeff(28) * x(5) 
ENDIF
f=aux

RETURN
END SUBROUTINE evaluate_fit_grad_quadratic

SUBROUTINE write_fit_hessian(degree,nvar,coeff, v, e)

USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: coeff(nvar)
REAL(DP), INTENT(INOUT) :: e(degree), v(degree,degree)
INTEGER ::  ideg

REAL(DP) :: amat(degree, degree)

amat(1,1)=coeff(3)
IF (degree > 1) THEN
   amat(1,2)=0.5_DP*coeff(6)
   amat(2,1)=amat(1,2)
   amat(2,2)=coeff(5)
ENDIF
IF (degree>2) THEN
   amat(1,3)=0.5_DP*coeff(9)
   amat(3,1)=amat(1,3)
   amat(2,3)=0.5_DP*coeff(10)
   amat(3,2)=amat(1,3)
   amat(3,3)=coeff(8)
ENDIF
IF (degree>3) THEN
   amat(1,4)=0.5_DP*coeff(13)
   amat(4,1)=amat(1,4)
   amat(2,4)=0.5_DP*coeff(14)
   amat(4,2)=amat(2,4)
   amat(3,4)=0.5_DP*coeff(15)
   amat(4,3)=amat(3,4)
   amat(4,4)=coeff(12)
ENDIF
IF (degree>4) THEN
   amat(1,5)=0.5_DP*coeff(18)
   amat(5,1)=amat(1,5)
   amat(2,5)=0.5_DP*coeff(19)
   amat(5,2)=amat(2,5)
   amat(3,5)=0.5_DP*coeff(20)
   amat(5,3)=amat(3,5)
   amat(4,5)=0.5_DP*coeff(21)
   amat(5,4)=amat(4,5)
   amat(5,5)=coeff(17)
ENDIF
IF (degree>5) THEN
   amat(1,6)=0.5_DP*coeff(24)
   amat(6,1)=amat(1,6)
   amat(2,6)=0.5_DP*coeff(25)
   amat(6,2)=amat(2,6)
   amat(3,6)=0.5_DP*coeff(26)
   amat(6,3)=amat(3,6)
   amat(4,6)=0.5_DP*coeff(27)
   amat(6,4)=amat(4,6)
   amat(5,6)=0.5_DP*coeff(28)
   amat(6,5)=amat(5,6)
   amat(6,6)=coeff(23)
ENDIF
amat(:,:) = 2.0_DP * amat(:,:) 

WRITE(stdout,'(/,5x,"Hessian:")')
DO ideg=1,degree
   WRITE(stdout,'(5x,6f18.8)') amat(ideg,:)
ENDDO

IF (degree > 1) THEN
   CALL diagonalize_r(degree,degree,amat,e,v)
   WRITE(stdout,'(/,5x,"Hessian eigenvalues:")')
   WRITE(stdout,'(5x,6f18.8)') e(1:degree)

   WRITE(stdout,'(/,5x,"Hessian eigenvectors (columns):")')
   DO ideg=1,degree
      WRITE(stdout,'(5x,6f18.8)') v(ideg,1:degree)
   ENDDO
ENDIF
WRITE(stdout,*)

RETURN
END SUBROUTINE write_fit_hessian

SUBROUTINE diagonalize_r (npw,nbnd,h,e,v)
!
! f90 interface to LAPACK routine ZHEEVX which calculates
! nbnd eigenvalues and eigenvectors of a complex hermitean matrix h. 
! 
use kinds, only : dp

IMPLICIT NONE

INTEGER, INTENT(IN) :: npw
INTEGER, INTENT(IN) :: nbnd

REAL(DP), INTENT(INOUT) :: h(npw,npw)   ! matrix to be diagonalized
REAL(DP), INTENT(OUT)   :: e(nbnd)      ! eigenvalues
REAL(DP), INTENT(OUT)   :: v(npw,nbnd)  ! eigenvectors (column-wise)

INTEGER  :: lwork,  &! auxiliary variable
            info,   &! flag saying if LAPACK execution was OK
            m        ! number of eigenvalues

CHARACTER(LEN=1) :: jobz, range, uplo  ! select the task in LAPACK

REAL(DP), ALLOCATABLE    :: work(:)      ! as above
INTEGER, ALLOCATABLE     :: iwork(:)     !    "
INTEGER, ALLOCATABLE     :: ifail(:)     !    "
REAL(DP)                 :: rdummy, zero ! dummy variable, zero
REAL(DP), ALLOCATABLE    :: ee(:)        ! axiliary space for eigenvalues
!
!   Initialize flags
!
jobz  = 'V' ! compute eigenvalues and eigenvectors
uplo  = 'U' ! LAPACK routines use the upper triangle of the input matrix
range = 'I' ! compute bands from 1 to nbnd

zero = 0.0_DP
v(:,:) = 0.0_DP
!
! allocate arrays of known size
!
ALLOCATE( ee(npw) )
ALLOCATE( iwork(5*npw) )
ALLOCATE( ifail(npw) )
ALLOCATE( work(16*npw) )
lwork=16*npw
!
! and diagonalize the matrix
!
CALL dsyevx(jobz, range, uplo, npw, h, npw, rdummy, rdummy, 1, nbnd, zero, &
            m, ee, v, npw, work, lwork, iwork, ifail, info)
!
IF (ABS(info) /= 0) THEN
   WRITE(stdout,'("Error in the diagonalization, info= ", i5)') info
   STOP 1
ENDIF
!
!
!  NB: h is overwritten by this routine. We save the eigenvalues on e
!      for plotting
!
e(1:nbnd)=ee(1:nbnd) 
       
DEALLOCATE(work)
DEALLOCATE(iwork)
DEALLOCATE(ifail)
DEALLOCATE(ee)

RETURN
END SUBROUTINE diagonalize_r
!
!
SUBROUTINE polifit(x,y,ndati,a,m1)
USE kinds, ONLY : DP
USE linear_solvers, ONLY : linsolvx
IMPLICIT NONE
INTEGER, INTENT(IN)   ::   ndati, &  ! number of data points
                           m1        ! polynomial coefficients 

REAL(DP), INTENT(IN)  :: x(ndati), y(ndati)
REAL(DP), INTENT(OUT) :: a(m1)

REAL(DP) :: amat(m1,m1), bvec(m1), eigv(m1,m1)
INTEGER  ::  i, j, k       ! counters

DO k=1,m1
   DO j=1,m1
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

CALL linsolvx(amat,m1,bvec,a)

RETURN
END SUBROUTINE polifit

SUBROUTINE write_poli(a,m1)
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: a(m1)

INTEGER :: i
CHARACTER(LEN=6) :: int_to_char

WRITE(stdout,'(/,5x,"Polynomial coefficients")')
DO i=1,m1
   WRITE(stdout,'(5x,a,e20.12)') "a"//TRIM(INT_TO_CHAR(i))//"=", a(i)
END DO

RETURN
END SUBROUTINE write_poli

SUBROUTINE print_quadratic_polynomial(degree, nvar, coeff)

USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE

INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: coeff(nvar)


  WRITE(stdout,'(/,5x,"Quadratic polynomial:")') 
  WRITE(stdout,'(5x,e15.7,"         +",e15.7," x1       +",e15.7," x1^2")') &
                                     coeff(1), coeff(2), coeff(3)
  IF (degree>1) THEN
     WRITE(stdout,'(4x,"+",e15.7," x2      +",e15.7," x2^2     +",e15.7,&
                         &" x1*x2")') coeff(4), coeff(5), coeff(6)
  ENDIF

  IF (degree>2) THEN
     WRITE(stdout,'(4x,"+",e15.7," x3 +",e15.7," x3^2 +",e15.7," x1*x3 +&
              &",e15.7," x2*x3")') coeff(7), coeff(8), coeff(9), coeff(10)
  ENDIF

  IF (degree>3) THEN
     WRITE(stdout,'(4x,"+",e15.7," x4 +",e15.7," x4^2 +",e15.7," x1*x4 +&
              &",e15.7," x2*x4",e15.7," x3*x4")') coeff(11), coeff(12), &
                                       coeff(13), coeff(14), coeff(15)
  ENDIF

  IF (degree>4) THEN
     WRITE(stdout,'(4x,"+",e15.7," x5 +",e15.7," x5^2 +",e15.7," x1*x5 +&
              &",e15.7," x2*x5",e15.7," x3*x5",e15.7, " x4*x5")') &
              coeff(16), coeff(17), coeff(18), coeff(19), coeff(20), coeff(21)
  ENDIF
  IF (degree>5) THEN
     WRITE(stdout,'(e15.7," x6 +",e15.7," x6^2 +",e15.7," x1*x6 +&
              &",e15.7," x2*x6",e15.7," x3*x6",e15.7, " x4*x6",&
              &e15.7," x5*x6")') coeff(22), coeff(23), coeff(24), &
                                 coeff(25), coeff(26), coeff(27), coeff(28)
  ENDIF

RETURN
END SUBROUTINE print_quadratic_polynomial

SUBROUTINE summarize_fitting_data(degree, ndata, x, f)

USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER :: degree, ndata
REAL(DP) :: x(degree,ndata), f(ndata)

INTEGER :: idata

  WRITE(stdout,'(5x," Fitting the following data")')
  IF (degree==1) THEN
     WRITE(stdout,'(10x,"x1",12x,"f")')
     DO idata=1,ndata
        WRITE(stdout,'(2f15.8)') x(1,idata), f(idata)
     ENDDO
  ELSEIF (degree==2) THEN
     WRITE(stdout,'(10x,"x1",12x,"x2",12x,"f")')
     DO idata=1,ndata
        WRITE(stdout,'(3f15.8)') x(1:2,idata), f(idata)
     ENDDO
  ELSEIF (degree==3) THEN
     WRITE(stdout,'(10x,"x1",12x,"x2",12x,"x3",12x,"f")')
     DO idata=1,ndata
        WRITE(stdout,'(4f15.8)') x(1:3,idata), f(idata)
     ENDDO
  ELSEIF (degree==4) THEN
     WRITE(stdout,'(10x,"x1",12x,"x2",12x,"x3",12x,12x,"x4",12x,"f")')
     DO idata=1,ndata
        WRITE(stdout,'(5f15.8)') x(1:4,idata), f(idata)
     ENDDO
  ELSEIF (degree==5) THEN
     WRITE(stdout,'(10x,"x1",12x,"x2",12x,"x3",12x,12x,"x4",12x,"x5",12x,"f")')
     DO idata=1,ndata
        WRITE(stdout,'(6f15.8)') x(1:5,idata), f(idata)
     ENDDO
  ELSEIF (degree==6) THEN
     WRITE(stdout,'(10x,"x1",12x,"x2",12x,"x3",12x,12x,"x4",12x,"x5",12x,&
                                                               &"x6",12x,"f")')
     DO idata=1,ndata
        WRITE(stdout,'(7f15.8)') x(1:6,idata), f(idata)
     ENDDO
  ENDIF

RETURN
END SUBROUTINE summarize_fitting_data

SUBROUTINE write_vector(degree, x)

IMPLICIT NONE
INTEGER, INTENT(IN) :: degree
REAL(DP) :: x(degree)

INTEGER :: i

  DO i=1, degree
     WRITE(stdout,'(23x,"x",i1,"=",f16.9)') i, x(i)
  END DO

RETURN
END SUBROUTINE write_vector

SUBROUTINE introduce_quadratic_fit(degree, nvar, ndata)

IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar, ndata

  WRITE(stdout,'(/,5x,"Fitting the data with a quadratic polynomial:")')

  WRITE(stdout,'(/,5x,"Number of variables:",10x,i5)')  degree
  WRITE(stdout,'(5x,"Coefficients of the quadratic:",i5)')  nvar
  WRITE(stdout,'(5x,"Number of fitting data:",7x,i5,/)')  ndata

RETURN
END SUBROUTINE introduce_quadratic_fit

SUBROUTINE print_chisq_quadratic(ndata, degree, nvar, x, f, coeff)

USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER :: ndata, degree, nvar
REAL(DP) :: x(degree, ndata), f(ndata), coeff(nvar)

REAL(DP) :: chisq, perc, aux
INTEGER :: idata

chisq=0.0_DP
perc=0.0_DP
DO idata=1,ndata
   CALL evaluate_fit_quadratic(degree,nvar,x(1,idata),aux,coeff)
!  WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
   IF (ABS(f(idata))>1.D-12) perc= perc + ABS((f(idata)-aux) / f(idata))
ENDDO
WRITE(stdout,'(5x,"chi square=",e18.5," relative error",e18.5,&
                                     &" %",/)') chisq, perc / ndata
RETURN
END SUBROUTINE print_chisq_quadratic

END MODULE quadratic_surfaces
