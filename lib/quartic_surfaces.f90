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
!   find_quartic_extremum, find the extremum closest to the input point
!
!   The following number of coefficients are necessary, depending on the
!   dimension of the quadratic function
!
!   dimension     number of coefficients
!      1                    5
!      2                    15
!      3                    35
!      4                    70
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
            compute_quartic_var, find_quartic_extremum, &
            find_quartic_quadratic_extremum, evaluate_quartic_quadratic, &
            evaluate_two_quartic, find_two_quartic_extremum, &
            print_quartic_polynomial, introduce_quartic_fit

CONTAINS

SUBROUTINE fit_multi_quartic(ndata,degree,nvar,x,f,coeff)
!
!  This routine receives as input a set of vectors x(degree,ndata) and
!  function values f(ndata) and gives as output the coefficients of
!  a quadratic interpolating polynomial coeff(nvar). In input
!  ndata is the number of data points. degree is the number of degrees of
!  freedom or equivalently the number of independent parameters (the
!  maximum is 4), and nvar is the number of coefficients of the
!  intepolating quadrating polynomial. 
!        degree     nvar
!          1         5
!          2        15
!          3        35
!          4        70
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
!      + a_36 x(4,i) + a_37 x(4,i)**2 + a_38 x(4,i)**3 + a_39 x(4,i)**4       4
!      + a_40 x(1,i)*x(4,i) + a_41 x(1,i)*x(4,i)**2 + a_42  x(1,i)*x(4,i)**3
!                           + a_43 x(1,i)**2*x(4,i) + a_44  x(1,i)**2*x(4,i)**2 
!                           + a_45 x(1,i)**3*x(4,i) 
!      + a_46 x(2,i)*x(4,i) + a_47 x(2,i)*x(4,i)**2 + a_48  x(2,i)*x(4,i)**3
!                           + a_49 x(2,i)**2*x(4,i) + a_50  x(2,i)**2*x(4,i)**2 
!                           + a_51 x(2,i)**3*x(4,i) 
!      + a_52 x(3,i)*x(4,i) + a_53 x(3,i)*x(4,i)**2 + a_54  x(3,i)*x(4,i)**3
!                           + a_55 x(3,i)**2*x(4,i) + a_56  x(3,i)**2*x(4,i)**2 
!                           + a_57 x(3,i)**3*x(4,i) 
!      + a_58 x(1,i) * x(2,i) * x(4,i) + a_59 x(1,i) * x(2,i)**2 * x(4,i)
!      + a_60 x(1,i) * x(2,i) * x(4,i)**2 + a_61 x(1,i)**2 * x(2,i) * x(4,i)
!      + a_62 x(1,i) * x(3,i) * x(4,i) + a_63 x(1,i) * x(3,i)**2 * x(4,i)
!      + a_64 x(1,i) * x(3,i) * x(4,i)**2 + a_65 x(1,i)**2 * x(3,i) * x(4,i)
!      + a_66 x(2,i) * x(3,i) * x(4,i) + a_67 x(2,i) * x(3,i)**2 * x(4,i)
!      + a_68 x(2,i) * x(3,i) * x(4,i)**2 + a_69 x(2,i)**2 * x(3,i) * x(4,i)
!      + a_70 x(1,i) * x(2,i) * x(3,i) * x(4,i)
!
USE kinds, ONLY : DP
USE quadratic_surfaces, ONLY : linsolvx
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar, ndata
REAL(DP), INTENT(IN) :: x(degree,ndata), f(ndata)
REAL(DP), INTENT(OUT) :: coeff(nvar)

REAL(DP) :: amat(ndata,nvar), aa(nvar,nvar), b(nvar) 

INTEGER :: ivar, jvar, idata

IF (degree>4.OR.degree<1) &
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

   IF (degree>3) THEN
      amat(idata,36) = x(4,idata)
      amat(idata,37) = x(4,idata)**2
      amat(idata,38) = x(4,idata)**3
      amat(idata,39) = x(4,idata)**4
      amat(idata,40) = x(1,idata) * x(4,idata)
      amat(idata,41) = x(1,idata) * x(4,idata)**2
      amat(idata,42) = x(1,idata) * x(4,idata)**3
      amat(idata,43) = x(1,idata) **2 * x(4,idata)
      amat(idata,44) = x(1,idata) **2 * x(4,idata)**2
      amat(idata,45) = x(1,idata) **3 * x(4,idata)
      amat(idata,46) = x(2,idata) * x(4,idata)
      amat(idata,47) = x(2,idata) * x(4,idata)**2
      amat(idata,48) = x(2,idata) * x(4,idata)**3
      amat(idata,49) = x(2,idata) **2 * x(4,idata)
      amat(idata,50) = x(2,idata) **2 * x(4,idata)**2
      amat(idata,51) = x(2,idata) **3 * x(4,idata)
      amat(idata,52) = x(3,idata) * x(4,idata)
      amat(idata,53) = x(3,idata) * x(4,idata)**2
      amat(idata,54) = x(3,idata) * x(4,idata)**3
      amat(idata,55) = x(3,idata) **2 * x(4,idata)
      amat(idata,56) = x(3,idata) **2 * x(4,idata)**2
      amat(idata,57) = x(3,idata) **3 * x(4,idata)
      amat(idata,58) = x(1,idata) * x(2,idata) * x(4,idata)
      amat(idata,59) = x(1,idata) * x(2,idata)**2 * x(4,idata)
      amat(idata,60) = x(1,idata) * x(2,idata) * x(4,idata)**2
      amat(idata,61) = x(1,idata)**2 * x(2,idata) * x(4,idata)
      amat(idata,62) = x(1,idata) * x(3,idata) * x(4,idata)
      amat(idata,63) = x(1,idata) * x(3,idata)**2 * x(4,idata)
      amat(idata,64) = x(1,idata) * x(3,idata) * x(4,idata)**2
      amat(idata,65) = x(1,idata)**2 * x(3,idata) * x(4,idata)
      amat(idata,66) = x(2,idata) * x(3,idata) * x(4,idata)
      amat(idata,67) = x(2,idata) * x(3,idata)**2 * x(4,idata)
      amat(idata,68) = x(2,idata) * x(3,idata) * x(4,idata)**2
      amat(idata,69) = x(2,idata)**2 * x(3,idata) * x(4,idata)
      amat(idata,70) = x(1,idata) * x(2,idata) * x(3,idata) * x(4,idata)
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
IF (degree>3) THEN
   aux = aux + coeff(36) * x(4) + coeff(37) * x(4)**2 + coeff(38) * x(4)**3 + &
                                                        coeff(39) * x(4)**4   &
             + coeff(40) * x(1) * x(4) + coeff(41) * x(1) * x(4)**2           &
             + coeff(42) * x(1) * x(4)**3 + coeff(43) * x(1)**2 * x(4)        &
             + coeff(44) * x(1)**2 * x(4)**2 + coeff(45) * x(1)**3 * x(4)     &
             + coeff(46) * x(2) * x(4) + coeff(47) * x(2) * x(4)**2           &
             + coeff(48) * x(2) * x(4)**3 + coeff(49) * x(2)**2 * x(4)        &
             + coeff(50) * x(2)**2 * x(4)**2 + coeff(51) * x(2)**3 * x(4)     &
             + coeff(52) * x(3) * x(4) + coeff(53) * x(3) * x(4)**2           &
             + coeff(54) * x(3) * x(4)**3 + coeff(55) * x(3)**2 * x(4)        &
             + coeff(56) * x(3)**2 * x(4)**2 + coeff(57) * x(3)**3 * x(4)     &
             + coeff(58) * x(1) * x(2) * x(4)     &
             + coeff(59) * x(1) * x(2)**2 * x(4)  &
             + coeff(60) * x(1) * x(2) * x(4)**2  &
             + coeff(61) * x(1)**2 * x(2) * x(4)  &
             + coeff(62) * x(1) * x(3) * x(4)     &
             + coeff(63) * x(1) * x(3)**2 * x(4)  &
             + coeff(64) * x(1) * x(3) * x(4)**2  &
             + coeff(65) * x(1)**2 * x(3) * x(4)  &
             + coeff(66) * x(2) * x(3) * x(4)     &
             + coeff(67) * x(2) * x(3)**2 * x(4)  &
             + coeff(68) * x(2) * x(3) * x(4)**2  &
             + coeff(69) * x(2)**2 * x(3) * x(4)  & 
             + coeff(70) * x(1) * x(2) * x(3) * x(4)  
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

IF (degree>3) THEN
   aux(1) = aux(1) + coeff(40) * x(4) + coeff(41) * x(4)**2 +               &
                     coeff(42) * x(4)**3 + 2.0_DP * coeff(43) * x(1)*x(4) + &
                     2.0_DP * coeff(44) * x(1)*x(4)**2 +                    & 
                     3.0_DP * coeff(45) * x(1)**2 * x(4) +                  &
                     coeff(58) * x(2) * x(4) + coeff(59) * x(2)**2 * x(4) + &
                     coeff(60) * x(2) * x(4)**2 +                           &
                     2.0_DP * coeff(61) * x(1) * x(2) * x(4) +              &
                     coeff(62) * x(3) * x(4) + coeff(63) * x(3)**2 * x(4) + &
                     coeff(64) * x(3) * x(4)**2 +                           &
                     2.0_DP * coeff(65) * x(1) * x(3) * x(4) +              &
                     coeff(70) * x(2) * x(3) * x(4) 
   aux(2) = aux(2) + coeff(46) * x(4) + coeff(47) * x(4)**2 +               &
                     coeff(48) * x(4)**3 + 2.0_DP * coeff(49) * x(2)*x(4) + &
                     2.0_DP * coeff(50) * x(2)*x(4)**2 +                    & 
                     3.0_DP * coeff(51) * x(2)**2 * x(4) +                  &
                     coeff(58) * x(1) * x(4) +                              &
                     2.0_DP *coeff(59) * x(1)*x(2)*x(4) +                   &
                     coeff(60) * x(1) * x(4)**2 +                           &
                     coeff(61) * x(1) **2 * x(4) +                          &
                     coeff(66) * x(3) * x(4) + coeff(67) * x(3)**2 * x(4) + &
                     coeff(68) * x(3) * x(4)**2 +                           &
                     2.0_DP * coeff(69) * x(2) * x(3) * x(4) +              &
                     coeff(70) * x(1) * x(3) * x(4) 

   aux(3) = aux(3) + coeff(52) * x(4) + coeff(53) * x(4)**2 +               &
                    coeff(54) * x(4)**3 + 2.0_DP * coeff(55) * x(3)*x(4) + &
                    2.0_DP * coeff(56) * x(3)*x(4)**2 +                    & 
                    3.0_DP * coeff(57) * x(3)**2 * x(4) +                  &
                     coeff(62) * x(1) * x(4) +                              &
                     2.0_DP *coeff(63) * x(1)*x(3)*x(4) +                   &
                     coeff(64) * x(1) * x(4)**2 +                           &
                     coeff(65) * x(1) **2 * x(4) +                          &
                     coeff(66) * x(2) * x(4) +                              &
                     2.0_DP*coeff(67)*x(2)*x(3)*x(4) +                      &
                     coeff(68) * x(2) * x(4)**2 +                           &
                     coeff(69) * x(2) **2 * x(4) +                          &
                     coeff(70) * x(1) * x(2) * x(4) 
   aux(4) = coeff(36) + 2.0_DP * coeff(35) * x(4) +                         &
                        3.0_DP * coeff(38) * x(4)**2 +                      &
                        4.0_DP * coeff(39) * x(4)**3 +                      &
                                 coeff(40) * x(1) +                         &
                        2.0_DP * coeff(41) * x(1) * x(4) +                  &
                        3.0_DP * coeff(42) * x(1) * x(4) **2 +              &
                                 coeff(43) * x(1) ** 2 +                    &
                        2.0_DP * coeff(44) * x(1) **2 * x(4) +              &
                                 coeff(45) * x(1) ** 3 +                    &
                                 coeff(46) * x(2) +                         &
                        2.0_DP * coeff(47) * x(2) * x(4) +                  &
                        3.0_DP * coeff(48) * x(2) * x(4) **2 +              &
                                 coeff(49) * x(2) ** 2 +                    &
                        2.0_DP * coeff(50) * x(2) **2 * x(4) +              &
                                 coeff(51) * x(2) ** 3 +                    &
                                 coeff(52) * x(3) +                         &
                        2.0_DP * coeff(53) * x(3) * x(4) +                  &
                        3.0_DP * coeff(54) * x(3) * x(4) **2 +              &
                                 coeff(55) * x(3) ** 2 +                    &
                        2.0_DP * coeff(56) * x(3) **2 * x(4) +              &
                                 coeff(57) * x(3) ** 3 +                    &
                                 coeff(58) * x(1) * x(2) +                  &
                                 coeff(59) * x(1) * x(2) **2 +              &
                        2.0_DP * coeff(60) * x(1) * x(2) * x(4) +           &
                                 coeff(61) * x(1) **2 * x(2) +              &
                                 coeff(62) * x(1) * x(3) +                  &
                                 coeff(63) * x(1) * x(3) **2 +              &
                        2.0_DP * coeff(64) * x(1) * x(3) * x(4) +           &
                                 coeff(65) * x(1) **2 * x(3) +              &
                                 coeff(66) * x(2) * x(3) +                  &
                                 coeff(67) * x(2) * x(3) **2 +              &
                        2.0_DP * coeff(68) * x(2) * x(3) * x(4) +           &
                                 coeff(69) * x(2) **2 * x(3) +              &
                                 coeff(70) * x(1) * x(2) * x(3) 
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

IF (degree>3) THEN
   aux(1,1) = aux(1,1) + 2.0_DP* coeff(43) * x(4) + 2.0_DP*coeff(44)*x(4)**2  &
                   + 6.0_DP * coeff(45) * x(1) * x(4)                         &
                   + 2.0_DP * coeff(61) * x(2) * x(4)                         &
                   + 2.0_DP * coeff(65) * x(3) * x(4) 
   aux(1,2) = aux(1,2) + coeff(58) * x(4) +                     &
                2.0_DP * coeff(59) * x(2) * x(4)                &
                       + coeff(60) * x(4)**2 +                  &
                2.0_DP * coeff(61) * x(1) * x(4)                &
                       + coeff(70) * x(3) * x(4)    
   aux(2,1) = aux(1,2)
   aux(1,3) = aux(1,3) + coeff(62)*x(4)                                    &
              + 2.0_DP * coeff(63) * x(3) * x(4)                           &
              +          coeff(64) * x(4)**2                               &
              + 2.0_DP * coeff(65) * x(1) * x(4)                           &
              +          coeff(70) * x(2) * x(4)

   aux(3,1) = aux(1,3)
   aux(1,4) = coeff(40) + 2.0_DP*coeff(41)*x(4) + 3.0_DP*coeff(42)*x(4)**2 + &
              2.0_DP * coeff(43) * x(1) + 4.0_DP * coeff(44) * x(1) * x(4) + &
              3.0_DP * coeff(45) * x(1)**2 + coeff(58) * x(2) +              &
                       coeff(59) * x(2)**2 + 2.0_DP * coeff(60)*x(2)*x(4) +  &
                       2.0_DP * coeff(61) * x(1)*x(2) + coeff(62) * x(3) +   &
                       coeff(63) * x(3)**2 + 2.0_DP * coeff(64)*x(3)*x(4) +  &
                       2.0_DP * coeff(65) * x(1)*x(3) +                      &
                       2.0_DP * coeff(70) * x(2)*x(3)
   aux(4,1) = aux(1,4)
   aux(2,2) = aux(2,2) + 2.0_DP*coeff(49) * x(4)        &
                       + 2.0_DP*coeff(50) * x(4)**2     &
                       + 6.0_DP*coeff(51) * x(2) * x(4) &
                       + 2.0_DP*coeff(59) * x(1) * x(4) &
                       + 2.0_DP*coeff(69) * x(3) * x(4) 
   aux(2,3) = aux(2,3) + coeff(66)*x(4)                                    &
              + 2.0_DP * coeff(67) * x(3) * x(4)                           &
              +          coeff(68) * x(4)**2                               &
              + 2.0_DP * coeff(69) * x(2) * x(4)                           &
              +          coeff(70) * x(1) * x(4)
   aux(3,2) = aux(2,3)
   aux(2,4) = coeff(46) + 2.0_DP*coeff(47)*x(4) + 3.0_DP*coeff(48)*x(4)**2 +  &
              2.0_DP * coeff(49) * x(2) + 4.0_DP * coeff(50) * x(2) * x(4) +  &
              3.0_DP * coeff(51) * x(2)**2 + coeff(58) * x(1) +               &
              2.0_DP * coeff(59) * x(1) * x(2) + 2.0_DP*coeff(60)*x(1)*x(4) + &
                       coeff(61) * x(1)**2 + coeff(66) * x(3) +               &
                       coeff(67) * x(3)**2 + 2.0_DP * coeff(68)*x(3)*x(4) +   &
                       2.0_DP * coeff(69) * x(2) * x(3) +                     &
                                coeff(70) * x(1) * x(3)

   aux(4,2) = aux(2,4)
   aux(3,3) = aux(3,3) + 2.0_DP*coeff(55) * x(4)        &
                       + 2.0_DP*coeff(56) * x(4)**2     &
                       + 6.0_DP*coeff(57) * x(3) * x(4) &
                       + 2.0_DP*coeff(63) * x(1) * x(4) &
                       + 2.0_DP*coeff(67) * x(2) * x(4) 

   aux(3,4) = coeff(52) + 2.0_DP*coeff(53)*x(4) + 3.0_DP*coeff(54)*x(4)**2 +  &
              2.0_DP * coeff(55) * x(3) + 4.0_DP * coeff(56) * x(3) * x(4) +  &
              6.0_DP * coeff(57) * x(3)**2 + coeff(62) * x(1) +               &
              2.0_DP * coeff(63) * x(1) * x(3) + 2.0_DP*coeff(64)*x(1)*x(4) + &
                       coeff(65) * x(1)**2 + coeff(66) * x(2) +               &
              2.0_DP * coeff(67) * x(2)* x(3) + 2.0_DP*coeff(68)*x(2)*x(4) +  &
                       coeff(69) * x(2) **2  +                                &
                       coeff(70) * x(1) * x(2)
   aux(4,3) = aux(3,4)
   aux(4,4) = 2.0_DP * coeff(35) + 6.0_DP * coeff(38)*x(4)      &
                                 + 12.0_DP * coeff(39)*x(4)**2  &
                     + 2.0_DP * coeff(41) * x(1)                &
                     + 6.0_DP * coeff(42) * x(1) * x(4)         &
                     + 2.0_DP * coeff(44) * x(1) ** 2           &
                     + 2.0_DP * coeff(47) * x(2)                &
                     + 6.0_DP * coeff(48) * x(2) * x(4)         &
                     + 2.0_DP * coeff(50) * x(2) ** 2           &    
                     + 2.0_DP * coeff(53) * x(3)                &
                     + 6.0_DP * coeff(54) * x(3) * x(4)         &
                     + 2.0_DP * coeff(56) * x(3) ** 2           &    
                     + 2.0_DP * coeff(60) * x(1) * x(2)         &   
                     + 2.0_DP * coeff(64) * x(1) * x(3)         &
                     + 2.0_DP * coeff(68) * x(2) * x(3)        
ENDIF
f(:,:)=aux(:,:)

RETURN
END SUBROUTINE evaluate_fit_hess_quartic

SUBROUTINE find_quartic_extremum(degree,nvar,x,f,coeff)
!
!  This routine starts from the point x and finds the extremum closest
!  to x. In output x are the coordinates of the extremum and f 
!  the value of the quartic function at the minimum
!
USE quadratic_surfaces, ONLY : linsolvx
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar
REAL(DP),INTENT(INOUT) :: x(degree), f
REAL(DP),INTENT(IN) :: coeff(nvar)

INTEGER, PARAMETER :: maxiter=100

INTEGER :: iter, ideg
REAL(DP), PARAMETER :: tol=1.D-12
REAL(DP) :: g(degree), y(degree), xold(degree)
REAL(DP) :: j(degree, degree) 
REAL(DP) :: deltax, fmod

xold(:)=x(:)
DO iter=1,maxiter
   !
   CALL evaluate_fit_grad_quartic(degree,nvar,x,g,coeff)
   !
   CALL evaluate_fit_hess_quartic(degree,nvar,x,j,coeff)
   !
   CALL linsolvx(j, degree, g, y)
   !
   !  Use Newton's method to find the zero of the gradient
   !
   x(:)= x(:) - y(:)
   fmod=0.0_DP
   deltax=0.0_DP
   DO ideg=1,degree
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
CALL evaluate_fit_quartic(degree,nvar,x,f,coeff)

RETURN
END SUBROUTINE find_quartic_extremum

SUBROUTINE find_two_quartic_extremum(degree,nvar4,x,f,coeff,coeff1)

IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar4
REAL(DP),INTENT(INOUT) :: x(degree), f
REAL(DP),INTENT(IN) :: coeff(nvar4), coeff1(nvar4)

REAL(DP) :: coeffadd4(nvar4)

coeffadd4=coeff+coeff1

CALL find_quartic_extremum(degree,nvar4,x,f,coeffadd4)

RETURN
END SUBROUTINE find_two_quartic_extremum

SUBROUTINE evaluate_two_quartic(degree,nvar4,x,f,coeff,coeff1)

IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar4
REAL(DP),INTENT(INOUT) :: x(degree), f
REAL(DP),INTENT(IN) :: coeff(nvar4), coeff1(nvar4)

REAL(DP) :: coeffadd4(nvar4)

coeffadd4=coeff+coeff1

CALL evaluate_fit_quartic(degree,nvar4,x,f,coeffadd4)

RETURN
END SUBROUTINE evaluate_two_quartic

SUBROUTINE find_quartic_quadratic_extremum(degree,nvar4,nvar,x,f,coeff4,coeff)
!
!   This subroutines adds the coefficients of a quadratic polynomium
!   to those of a quartic polynomium and finds the minimum of the sum
!   of the two
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar, nvar4
REAL(DP), INTENT(IN) :: coeff(nvar), coeff4(nvar4) 
REAL(DP), INTENT(INOUT) :: x(degree), f

REAL(DP) :: coeffadd4(nvar4)

CALL set_quartic_coefficients(degree, nvar4, nvar, coeffadd4, coeff4, coeff)

CALL find_quartic_extremum(degree,nvar4,x,f,coeffadd4)

RETURN
END SUBROUTINE find_quartic_quadratic_extremum

SUBROUTINE evaluate_quartic_quadratic(degree,nvar4,nvar,x,f,coeff4,coeff)
IMPLICIT NONE

INTEGER, INTENT(IN) :: degree, nvar, nvar4
REAL(DP), INTENT(IN) :: coeff(nvar), coeff4(nvar4) 
REAL(DP), INTENT(INOUT) :: x(degree), f

REAL(DP) :: coeffadd4(nvar4)

CALL set_quartic_coefficients(degree, nvar4, nvar, coeffadd4, coeff4, coeff)

CALL evaluate_fit_quartic(degree,nvar4,x,f,coeffadd4)

RETURN
END SUBROUTINE evaluate_quartic_quadratic

SUBROUTINE set_quartic_coefficients(degree, nvar4, nvar, coeffadd4, coeff4, coeff)
IMPLICIT NONE

INTEGER, INTENT(IN) :: degree, nvar4, nvar
REAL(DP), INTENT(IN) :: coeff4(nvar4), coeff(nvar)
REAL(DP), INTENT(INOUT) :: coeffadd4(nvar4)

coeffadd4=coeff4
coeffadd4(1)=coeffadd4(1)+coeff(1)
coeffadd4(2)=coeffadd4(2)+coeff(2)
coeffadd4(3)=coeffadd4(3)+coeff(3)

IF (degree > 1) THEN
   coeffadd4(6)=coeffadd4(6)+coeff(4)
   coeffadd4(7)=coeffadd4(7)+coeff(5)
   coeffadd4(10)=coeffadd4(10)+coeff(6)
ENDIF

IF (degree > 2) THEN
   coeffadd4(16)=coeffadd4(16)+coeff(7)
   coeffadd4(17)=coeffadd4(17)+coeff(8)
   coeffadd4(20)=coeffadd4(20)+coeff(9)
   coeffadd4(26)=coeffadd4(26)+coeff(10)
END IF

IF (degree > 3) THEN
   coeffadd4(36)=coeffadd4(36)+coeff(11)
   coeffadd4(37)=coeffadd4(37)+coeff(12)
   coeffadd4(40)=coeffadd4(40)+coeff(13)
   coeffadd4(46)=coeffadd4(46)+coeff(14)
   coeffadd4(52)=coeffadd4(52)+coeff(15)
END IF

IF (degree > 4) THEN

END IF
IF (degree > 5) THEN

END IF

RETURN
END SUBROUTINE set_quartic_coefficients

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
ELSEIF (degree==4) THEN
   compute_quartic_var=70
ELSE
   compute_quartic_var=0
ENDIF

RETURN
END FUNCTION compute_quartic_var

SUBROUTINE print_quartic_polynomial(degree, nvar, coeff)

USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE

INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: coeff(nvar)

  WRITE(stdout,'(/,5x,"Quartic polynomial:")') 
  WRITE(stdout,'(5x,f13.7,"          +",f13.7," x1        +",f13.7," x1^2")') &
                                       coeff(1), coeff(2), coeff(3)
  WRITE(stdout,'(4x,"+",f13.7," x1^3     +",f13.7," x1^4")') coeff(4), &
                                                              coeff(5)
  IF (degree>1) THEN
     WRITE(stdout,'(4x,"+",f13.7," x2       +",f13.7," x2^2      +",&
                     &f13.7," x2^3")') coeff(6), coeff(7), coeff(8)
     WRITE(stdout,'(4x,"+",f13.7," x2^4     +",f13.7," x1 x2     +",f13.7,&
                              &" x1 x2^2")') coeff(9), coeff(10), coeff(11)
     WRITE(stdout,'(4x,"+",f13.7," x1 x2^3  +",f13.7," x1^2 x2   +",& 
                          &f13.7," x1^2 x2^2")') coeff(12), coeff(13), coeff(14)
     WRITE(stdout,'(4x,"+",f13.7," x1^3 x2")') coeff(15)
  ENDIF
  IF (degree>2) THEN
     WRITE(stdout,'(f15.8," x3 +",f15.8," x3^2 +",f15.8," x3^3")') coeff(16), &
                                                     coeff(17), coeff(18)
     WRITE(stdout,'(5x,"+",f15.8," x3^4 +",f15.8," x1 x3 +",f15.8,&
                              &" x1 x3^2")') coeff(19), coeff(20), coeff(21)
     WRITE(stdout,'(5x,"+",f15.8," x1 x3^3 +",f15.8," x1^2 x3 +",& 
                          &f15.8," x1^2 x3^2")') coeff(22), coeff(23), coeff(24)
     WRITE(stdout,'(5x,"+",f15.8," x1^3 x3 +",f15.8," x2 x3 +",& 
                          &f15.8," x2 x3^2")') coeff(25), coeff(26), coeff(27)
     WRITE(stdout,'(5x,"+",f15.8," x2 x3^3 +",f15.8," x2^2 x3 +",& 
                          &f15.8," x2^2 x3^2")') coeff(28), coeff(29), coeff(30)
     WRITE(stdout,'(5x,"+",f15.8," x2^2 x3 +",f15.8," x2^2 x3^2 +",& 
                          &f15.8," x2^3 x3")') coeff(29), coeff(30), coeff(31)
     WRITE(stdout,'(5x,"+",f15.8," x1 x2 x3 +",f15.8," x1 x2^2 x3 +",& 
                          &f15.8," x1 x2 x3^2")') coeff(32), coeff(33), &
                                                                     coeff(34)
     WRITE(stdout,'(5x,"+",f15.8," x1^2 x2 x3")') coeff(35)
  ENDIF

RETURN
END SUBROUTINE print_quartic_polynomial

SUBROUTINE introduce_quartic_fit(degree, nvar, ndata)
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar, ndata

WRITE(stdout,'(/,5x,"Fitting the data with a quartic polynomial:")')

WRITE(stdout,'(/,5x,"Number of variables:",10x,i5)')  degree
WRITE(stdout,'(5x,"Coefficients of the quartic:",2x,i5)')  nvar
WRITE(stdout,'(5x,"Number of fitting data:",7x,i5,/)')  ndata

RETURN
END SUBROUTINE introduce_quartic_fit

END MODULE quartic_surfaces
