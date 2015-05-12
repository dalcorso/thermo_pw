!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE quartic_surfaces
!
!   this module contains the support routines for dealing with quartic
!   surfaces interpolation up to dimension 3. It is used to interpolate 
!   the total energy as a function of the celldm parameters. 
!   It provides the following routines:
!
!   fit_multi_quartic that receives as input the coordinates of some
!   points and the value of the quartic function in these points and
!   finds the coefficients of the quartic polynomial that passes through
!   the points.
!   evaluate_fit_quartic, given the coordinates of a point, and the
!   coefficients of the quartic polynomial, evaluates the quartic polynomial
!   at that point.
!   evaluate_fit_grad_quartic, given the coordinates of a point, and the
!   coefficients of the quartic polynomial, evaluates the gradient of the 
!   function at that point.
!   evaluate_fit_hess_quartic, given the coordinates of a point, and the
!   coefficients of the quartic polynomial, evaluates the hessian of the 
!   function at that point.
!
!
!   The following number of coefficients are necessary, depending on the
!   dimension of the quadratic function
!
!   dimension     number of coefficients
!      1                    5
!      2                    15
!      3                    35
!
!    To interpolate the function it is better to give to
!    multi_quartic a number of points equal or larger than
!    to the number of coefficients. The routines makes a least square fit
!    of the data.
!   
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: fit_multi_quartic,  evaluate_fit_quartic, &
            evaluate_fit_grad_quartic, evaluate_fit_hess_quartic, &
            compute_quartic_var

CONTAINS

SUBROUTINE fit_multi_quartic(ndata,degree,nvar,x,f,coeff)
!
!  This routine receives as input a set of vectors x(degree,ndata) and
!  function values f(ndata) and gives as output the coefficients of
!  a quadratic interpolating polynomial coeff(nvar). In input
!  ndata is the number of data points. degree is the number of degrees of
!  freedom or equivalently the number of independent parameters (the
!  maximum is 3), and nvar is the number of coefficients of the
!  intepolating quadrating polynomial. 
!        degree     nvar
!          1         5
!          2        15
!          3        35
!  The coefficients are organized like this:
!  a_1 + a_2  x(1,i) + a_3  x(1,i)**2 + a_4  x(1,i)**3 + a_5 x(1,i)**4        
!      + a_6  x(2,i) + a_7  x(2,i)**2 + a_8  x(2,i)**3 + a_9 x(2,i)**4        2
!      + a_10 x(1,i)*x(2,i) + a_11 x(1,i)*x(2,i)**2 + a_12  x(1,i)*x(2,i)**3
!                           + a_13 x(1,i)**2*x(2,i) + a_14  x(1,i)**2*x(2,i)**2 
!                           + a_15 x(1,i)**3*x(2,i) 
!      + a_16 x(3,i) + a_17 x(3,i)**2 + a_18 x(3,i)**3 + a_19 x(3,i)**4       3
!      + a_20 x(1,i)*x(3,i) + a_21 x(1,i)*x(3,i)**2 + a_22  x(1,i)*x(3,i)**3
!                           + a_23 x(1,i)**2*x(3,i) + a_24  x(1,i)**2*x(3,i)**2 
!                           + a_25 x(1,i)**3*x(3,i) 
!      + a_26 x(2,i)*x(3,i) + a_27 x(2,i)*x(3,i)**2 + a_28  x(2,i)*x(3,i)**3
!                           + a_29 x(2,i)**2*x(3,i) + a_30  x(2,i)**2*x(3,i)**2 
!                           + a_31 x(2,i)**3*x(3,i) 
!      + a_32 x(1,i) * x(2,i) * x(3,i) + a_33 x(1,i) * x(2,i)**2 * x(3,i)
!      + a_34 x(1,i) * x(2,i) * x(3,i)**2 + a_35 x(1,i)**2 * x(2,i) * x(3,i)
!
USE kinds, ONLY : DP
USE quadratic_surfaces, ONLY : linsolvx
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar, ndata
REAL(DP), INTENT(IN) :: x(degree,ndata), f(ndata)
REAL(DP), INTENT(OUT) :: coeff(nvar)

REAL(DP) :: amat(ndata,nvar), aa(nvar,nvar), b(nvar) 

INTEGER :: ivar, jvar, idata

IF (degree>3.OR.degree<1) &
   CALL errore('multi_quartic','degree must be from 1 to 3',1)
IF (ndata < 3) &
   CALL errore('multi_quartic','Too few sampling data',1)
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
   amat(idata,4) = x(1,idata)**3
   amat(idata,5) = x(1,idata)**4
   IF (degree>1) THEN
      amat(idata,6) = x(2,idata)
      amat(idata,7) = x(2,idata)**2
      amat(idata,8) = x(2,idata)**3
      amat(idata,9) = x(2,idata)**4
      amat(idata,10) = x(1,idata) * x(2,idata)
      amat(idata,11) = x(1,idata) * x(2,idata)**2
      amat(idata,12) = x(1,idata) * x(2,idata)**3
      amat(idata,13) = x(1,idata) **2 * x(2,idata)
      amat(idata,14) = x(1,idata) **2 * x(2,idata)**2
      amat(idata,15) = x(1,idata) **3 * x(2,idata)
   ENDIF
   IF (degree>2) THEN
      amat(idata,16) = x(3,idata)
      amat(idata,17) = x(3,idata)**2
      amat(idata,18) = x(3,idata)**3
      amat(idata,19) = x(3,idata)**4
      amat(idata,20) = x(1,idata) * x(3,idata)
      amat(idata,21) = x(1,idata) * x(3,idata)**2
      amat(idata,22) = x(1,idata) * x(3,idata)**3
      amat(idata,23) = x(1,idata) **2 * x(3,idata)
      amat(idata,24) = x(1,idata) **2 * x(3,idata)**2
      amat(idata,25) = x(1,idata) **3 * x(3,idata)
      amat(idata,26) = x(2,idata) * x(3,idata)
      amat(idata,27) = x(2,idata) * x(3,idata)**2
      amat(idata,28) = x(2,idata) * x(3,idata)**3
      amat(idata,29) = x(2,idata) **2 * x(3,idata)
      amat(idata,30) = x(2,idata) **2 * x(3,idata)**2
      amat(idata,31) = x(2,idata) **3 * x(3,idata)
      amat(idata,32) = x(1,idata) * x(2,idata) * x(3,idata)
      amat(idata,33) = x(1,idata) * x(2,idata)**2 * x(3,idata)
      amat(idata,34) = x(1,idata) * x(2,idata) * x(3,idata)**2
      amat(idata,35) = x(1,idata)**2 * x(2,idata) * x(3,idata)
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
END SUBROUTINE fit_multi_quartic

SUBROUTINE evaluate_fit_quartic(degree,nvar,x,f,coeff)

USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: x(degree)
REAL(DP), INTENT(IN) :: coeff(nvar)
REAL(DP), INTENT(INOUT) :: f

REAL(DP) :: aux

aux=coeff(1) + coeff(2) * x(1) + coeff(3) * x(1)**2 + coeff(4) * x(1)**3 + &
                                                      coeff(5) * x(1)**4
IF (degree>1) THEN
   aux = aux + coeff(6) * x(2) + coeff(7) * x(2)**2 + coeff(8) * x(2)**3 + &
                                                      coeff(9) * x(2)**4  &
             + coeff(10) * x(1) * x(2) + coeff(11) * x(1) * x(2)**2 &
             + coeff(12) * x(1) * x(2)**3 + coeff(13) * x(1)**2 * x(2) &
             + coeff(14) * x(1)**2 * x(2)**2 + coeff(15) * x(1)**3 * x(2)
ENDIF
IF (degree>2) THEN
   aux = aux + coeff(16) * x(3) + coeff(17) * x(3)**2 + coeff(18) * x(3)**3 + &
                                                        coeff(19) * x(3)**4  &
             + coeff(20) * x(1) * x(3) + coeff(21) * x(1) * x(3)**2 &
             + coeff(22) * x(1) * x(3)**3 + coeff(23) * x(1)**2 * x(3) &
             + coeff(24) * x(1)**2 * x(3)**2 + coeff(25) * x(1)**3 * x(3) &
             + coeff(26) * x(2) * x(3) + coeff(27) * x(2) * x(3)**2 &
             + coeff(28) * x(2) * x(3)**3 + coeff(29) * x(2)**2 * x(3) &
             + coeff(30) * x(2)**2 * x(3)**2 + coeff(31) * x(2)**3 * x(3) &
             + coeff(32) * x(1) * x(2) * x(3)     &
             + coeff(33) * x(1) * x(2)**2 * x(3)  &
             + coeff(34) * x(1) * x(2) * x(3)**2  &
             + coeff(35) * x(1)**2 * x(2) * x(3)  
ENDIF

f=aux

RETURN
END SUBROUTINE evaluate_fit_quartic

SUBROUTINE evaluate_fit_grad_quartic(degree,nvar,x,f,coeff)

USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: x(degree)
REAL(DP), INTENT(IN) :: coeff(nvar)
REAL(DP), INTENT(INOUT) :: f(degree)

REAL(DP) :: aux(degree)

aux(1) = coeff(2) + 2.0_DP * coeff(3) * x(1) + 3.0_DP * coeff(4)*x(1)**2 + &
                    4.0_DP * coeff(5) * x(1)**3

IF (degree>1) THEN
   aux(1) = aux(1) + coeff(10) * x(2) + coeff(11) * x(2)**2 + &
                     coeff(12) * x(2)**3 + 2.0_DP * coeff(13) * x(1) * x(2) + &
                     2.0_DP * coeff(14) * x(1) * x(2)**2 + & 
                     3.0_DP * coeff(15) * x(1)**2 * x(2)
   aux(2) = coeff(6) + 2.0_DP * coeff(7) * x(2) + 3.0_DP * coeff(8)*x(2)**2 + &
            4.0_DP*coeff(9)*x(2)**3 + coeff(10) * x(1) &
          + 2.0_DP*coeff(11)*x(1)*x(2) &
          + 3.0_DP*coeff(12)*x(1)*x(2)**2 + coeff(13)*x(1)**2 + &
           2.0_DP*coeff(14)*x(1)**2*x(2) + coeff(15)*x(1)**3 
ENDIF
IF (degree>2) THEN
   aux(1) = aux(1) + coeff(20) * x(3) + coeff(21) * x(3)**2 +               &
                     coeff(22) * x(3)**3 + 2.0_DP * coeff(23) * x(1)*x(3) + &
                     2.0_DP * coeff(24) * x(1)*x(3)**2 +                    & 
                     3.0_DP * coeff(25) * x(1)**2 * x(3) +                  &
                     coeff(32) * x(2) * x(3) + coeff(33) * x(2)**2 * x(3) + &
                     coeff(34) * x(2) * x(3)**2 +                           &
                     2.0_DP * coeff(35) * x(1) * x(2) * x(3) 
   aux(2) = aux(2) + coeff(26) * x(3) + coeff(27) * x(3)**2                 &
          + coeff(28) * x(3)**3 + 2.0_DP*coeff(29)*x(2)*x(3)                &
          + 2.0_DP*coeff(30)*x(2)*x(3)**2 + 3.0_DP*coeff(31)*x(2)**2*x(3)   &
          + coeff(32) * x(1) * x(3) + 2.0_DP * coeff(33) * x(1) * x(2) * x(3) &
          + coeff(34) * x(1) * x(3)**2 + coeff(35) * x(1)**2 * x(3) 
 
   aux(3) = coeff(16) + 2.0_DP*coeff(17)*x(3) + 3.0_DP*coeff(18)*x(3)**2  &
          + 4.0_DP*coeff(19)*x(3)**3 + coeff(20)*x(1)                     &
          + 2.0_DP*coeff(21)*x(1)*x(3)                                    &
          + 3.0_DP*coeff(22)*x(1)*x(3)**2 + coeff(23)*x(1)**2             &
          + 2.0_DP*coeff(24)*x(1)**2*x(3) + coeff(25)*x(1)**3             &
          +        coeff(26)*x(2)                                         &
          + 2.0_DP*coeff(27)*x(2)*x(3)                                    &
          + 3.0_DP*coeff(28)*x(2)*x(3)**2 + coeff(29)*x(2)**2             &
          + 2.0_DP*coeff(30)*x(2)**2*x(3) + coeff(31)*x(2)**3             &
          + coeff(32)*x(1)*x(2) + coeff(33)*x(1)*x(2)**2                  &
          + 2.0_DP*coeff(34)*x(1)*x(2)*x(3) + coeff(35)*x(1)**2*x(2) 
ENDIF

f=aux

RETURN
END SUBROUTINE evaluate_fit_grad_quartic

SUBROUTINE evaluate_fit_hess_quartic(degree,nvar,x,f,coeff)

USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: x(degree)
REAL(DP), INTENT(IN) :: coeff(nvar)
REAL(DP), INTENT(INOUT) :: f(degree,degree)

REAL(DP) :: aux(degree,degree)

aux(1,1) = 2.0_DP * coeff(3) + 6.0_DP * coeff(4) * x(1) + &
                              12.0_DP * coeff(5) * x(1)**2

IF (degree>1) THEN
   aux(1,1) = aux(1,1) + 2.0_DP* coeff(13) * x(2) + 2.0_DP*coeff(14)*x(2)**2+ &
                     6.0_DP * coeff(15) * x(1) * x(2)
   aux(1,2) = coeff(10) + 2.0_DP * coeff(11) * x(2) + &
                          3.0_DP * coeff(12) * x(2)**2 + &
                          2.0_DP * coeff(13) * x(1) +    &
                          4.0_DP * coeff(14) * x(1) * x(2) + &
                          3.0_DP * coeff(15) * x(1)**2
   aux(2,1) = aux(1,2)
   aux(2,2)=2.0_DP*coeff(7) + 6.0_DP*coeff(8)*x(2) + &
           12.0_DP*coeff(9)*x(2)**2 + 2.0_DP*coeff(11) * x(1) &
          + 6.0_DP*coeff(12)*x(1)*x(2) &
          + 2.0_DP*coeff(14)*x(1)**2 
ENDIF

IF (degree>2) THEN
   aux(1,1) = aux(1,1) + 2.0_DP* coeff(23) * x(3) + 2.0_DP*coeff(24)*x(3)**2  &
                   + 6.0_DP * coeff(25) * x(1) * x(3) &
                   + 2.0_DP * coeff(35) * x(2) * x(3)
   aux(1,2) = aux(1,2) + coeff(32) * x(3) +        &
                2.0_DP * coeff(33) * x(2) * x(3)   &
                       + coeff(34) * x(3)**2 +     &
                2.0_DP * coeff(35) * x(1) * x(3)    
   aux(2,1) = aux(1,2)
   aux(1,3) = coeff(20) + 2.0_DP*coeff(21)*x(3) + 3.0_DP*coeff(22)*x(3)**2 + &
              2.0_DP * coeff(23) * x(1) + 4.0_DP * coeff(24) * x(1) * x(3) + &
              3.0_DP * coeff(25) * x(1)**2 + coeff(32) * x(2) +             &
                       coeff(33) * x(2)**2 + 2.0_DP * coeff(34)*x(2)*x(3) + &
                       2.0_DP * coeff(35) * x(1)*x(2)
   aux(3,1) = aux(1,3)
   aux(2,2) = aux(2,2) + 2.0_DP*coeff(29) * x(3)        &
                       + 2.0_DP*coeff(30) * x(3)**2     &
                       + 6.0_DP*coeff(31) * x(2) * x(3) &
                       + 2.0_DP*coeff(33) * x(1) * x(3) 
   aux(2,3) = coeff(26) + 2.0_DP*coeff(27)*x(3) + 3.0_DP*coeff(28)*x(3)**2 + &
              2.0_DP * coeff(29) * x(2) + 4.0_DP * coeff(30) * x(2)*x(3) + &
              3.0_DP * coeff(31) * x(2)**2 + coeff(32) * x(1) +              &
                       2.0_DP * coeff(33) * x(1)*x(2) +                    &
                       2.0_DP * coeff(34) * x(1)*x(3) + coeff(35) * x(1)**2

   aux(3,2) = aux(2,3)
   aux(3,3) = 2.0_DP * coeff(17) + 6.0_DP * coeff(18)*x(3)      &
                                 + 12.0_DP * coeff(19)*x(3)**2  &
                     + 2.0_DP * coeff(21) * x(1)                &
                     + 6.0_DP * coeff(22) * x(1) * x(3)         &
                     + 2.0_DP * coeff(24) * x(1) ** 2           &
                     + 2.0_DP * coeff(27) * x(2)                &
                     + 6.0_DP * coeff(28) * x(2) * x(3)         &
                     + 2.0_DP * coeff(30) * x(2) ** 2           &    
                     + 2.0_DP * coeff(34) * x(1) * x(2)        
ENDIF

f(:,:)=aux(:,:)

RETURN
END SUBROUTINE evaluate_fit_hess_quartic

FUNCTION compute_quartic_var(degree)  

IMPLICIT NONE
INTEGER :: compute_quartic_var
INTEGER, INTENT(IN) :: degree

IF (degree==1) THEN
   compute_quartic_var=5
ELSEIF (degree==2) THEN
   compute_quartic_var=15
ELSEIF (degree==3) THEN
   compute_quartic_var=35
ELSE
   compute_quartic_var=0
ENDIF

RETURN
END FUNCTION compute_quartic_var

END MODULE quartic_surfaces
