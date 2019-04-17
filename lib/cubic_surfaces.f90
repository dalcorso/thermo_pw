!
! Copyright (C) 2018-2019 C. Malica
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE cubic_surfaces
!
!   this module contains the support routines for dealing with cubic
!   surfaces interpolation up to dimension 6. It is used to interpolate 
!   the total energy as a function of the celldm parameters. 
!   It provides the following routines:
!
!   fit_multi_cubic that receives as input the coordinates of some
!   points and the value of the cubic function in these points and
!   finds the coefficients of the cubic polynomial that passes through
!   the points.
!   evaluate_fit_cubic, given the coordinates of a point, and the
!   coefficients of the cubic polynomial, evaluates the cubic polynomial
!   at that point.
!   evaluate_fit_grad_cubic, given the coordinates of a point, and the
!   coefficients of the cubic polynomial, evaluates the gradient of the 
!   function at that point.
!   evaluate_fit_hess_cubic, given the coordinates of a point, and the
!   coefficients of the cubic polynomial, evaluates the hessian of the 
!   function at that point.
!
!   find_cubic_extremum, find the extremum closest to the input point
!
!   The number of coefficients of a multivariate polynomial is 
!   (n+d)!/n!/d! where d is the degree of the polynomial (in the cubic
!   case d=3) and n is the number of variables (from 1 till 6). 
!
!   The following number of coefficients are necessary, depending on the
!   dimension of the cubic function.
!
!   number of variables		number of coefficients
!          1				4	       
!          2	      		       10
!          3			       20                   
!          4                           35
!          5                           56
!          6                           84
!
!    To interpolate the function it is better to give to
!    multi_cubic a number of points equal or larger than
!    to the number of coefficients. The routines makes a 
!    least square fit of the data.
!
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: fit_multi_cubic,  evaluate_fit_cubic, &
            evaluate_fit_grad_cubic, evaluate_fit_hess_cubic, &
            cubic_var, find_cubic_extremum, &
            find_quartic_cubic_extremum, evaluate_quartic_cubic, &
            print_cubic_polynomial, introduce_cubic_fit, &
            print_chisq_cubic, & 
            print_chisq_quartic_cubic

CONTAINS

SUBROUTINE fit_multi_cubic(ndata,degree,nvar,lsolve,x,f,coeff)
!
!  This routine receives as input a set of vectors x(degree,ndata) and
!  function values f(ndata) and gives as output the coefficients of
!  a cubic interpolating polynomial coeff(nvar). In input
!  ndata is the number of data points. degree is the number of degrees of
!  freedom or equivalently the number of independent parameters (the
!  maximum is 6), and nvar is the number of coefficients of the
!  intepolating cubic polynomial. 
!        
!  lsolve can be 1, 2 or 3. It chooses the method to compute the
!  polynomial coefficients. Using 1 a matrix nvar x nvar is calculated
!                           Using 2 the overdetemined linear system is solved
!                           using QR or LQ factorization
!                           Using 3 the overdetermined linear system is solved
!                           using SVD decomposition
!                           If lsolve is not one of these values method 2 
!                           is used 
!
!  The coefficients are organized in the following manner:
!   
!  a_1 + a_2 x(1,i)  + a_3  x(1,i)**2 + a_4  x(1,i)**3					1
!        
!      + a_5 x(2,i)  + a_6  x(2,i)**2 + a_7  x(2,i)**3        				2
!      + a_8 x(1,i)*x(2,i) + a_9 x(1,i)*x(2,i)**2 + a_10 x(1,i)**2*x(2,i)
!
!      + a_11 x(3,i) + a_12 x(3,i)**2 + a_13 x(3,i)**3                   	     	3
!      + a_14 x(1,i)*x(3,i) + a_15 x(1,i)*x(3,i)**2 + a_16 x(1,i)**2*x(3,i) 
!      + a_17 x(2,i)*x(3,i) + a_18 x(2,i)*x(3,i)**2 + a_19 x(2,i)**2*x(3,i)  
!      + a_20 x(1,i) * x(2,i) * x(3,i)
!
!      + a_21 x(4,i) + a_22 x(4,i)**2 + a_23 x(4,i)**3			                4
!      + a_24 x(1,i)*x(4,i) + a_25 x(1,i)*x(4,i)**2 + a_26 x(1,i)**2*x(4,i) 
!      + a_27 x(2,i)*x(4,i) + a_28 x(2,i)*x(4,i)**2 + a_29 x(2,i)**2*x(4,i) 
!      + a_30 x(3,i)*x(4,i) + a_31 x(3,i)*x(4,i)**2 + a_32 x(3,i)**2*x(4,i) 
!      + a_33 x(1,i) * x(2,i) * x(4,i)
!      + a_34 x(1,i) * x(3,i) * x(4,i) 
!      + a_35 x(2,i) * x(3,i) * x(4,i) 
!        
!      + a_36 x(5,i) + a_37 x(5,i)**2 + a_38 x(5,i)**3					5
!      + a_39 x(1,i)*x(5,i) + a_40 x(1,i)*x(5,i)**2 + a_41 x(1,i)**2*x(5,i)         
!      + a_42 x(2,i)*x(5,i) + a_43 x(2,i)*x(5,i)**2 + a_44 x(2,i)**2*x(5,i) 
!      + a_45 x(3,i)*x(5,i) + a_46 x(3,i)*x(5,i)**2 + a_47 x(3,i)**2*x(5,i) 
!      + a_48 x(4,i)*x(5,i) + a_49 x(4,i)*x(5,i)**2 + a_50 x(4,i)**2*x(5,i)              
!      + a_51 x(1,i) * x(2,i) * x(5,i)
!      + a_52 x(1,i) * x(3,i) * x(5,i)
!      + a_53 x(1,i) * x(4,i) * x(5,i) 
!      + a_54 x(2,i) * x(3,i) * x(5,i) 
!      + a_55 x(2,i) * x(4,i) * x(5,i)
!      + a_56 x(3,i) * x(4,i) * x(5,i)
!
!      + a_57 x(6,i) + a_58 x(6,i)**2 + a_59 x(6,i)**3					6
!      + a_60 x(1,i)*x(6,i) + a_61 x(1,i)*x(6,i)**2 + a_62 x(1,i)**2*x(6,i)  
!      + a_63 x(2,i)*x(6,i) + a_64 x(2,i)*x(6,i)**2 + a_65 x(2,i)**2*x(6,i) 
!      + a_66 x(3,i)*x(6,i) + a_67 x(3,i)*x(6,i)**2 + a_68 x(3,i)**2*x(6,i)  
!      + a_69 x(4,i)*x(6,i) + a_70 x(4,i)*x(6,i)**2 + a_71 x(4,i)**2*x(6,i)  
!      + a_72 x(5,i)*x(6,i) + a_73 x(5,i)*x(6,i)**2 + a_74 x(5,i)**2*x(6,i)  
!      + a_75 x(1,i) * x(2,i) * x(6,i) 
!      + a_76 x(1,i) * x(3,i) * x(6,i) 
!      + a_77 x(1,i) * x(4,i) * x(6,i) 
!      + a_78 x(1,i) * x(5,i) * x(6,i) 
!      + a_79 x(2,i) * x(3,i) * x(6,i) 
!      + a_80 x(2,i) * x(4,i) * x(6,i) 
!      + a_81 x(2,i) * x(5,i) * x(6,i) 
!      + a_82 x(3,i) * x(4,i) * x(6,i) 
!      + a_83 x(3,i) * x(5,i) * x(6,i) 
!      + a_84 x(4,i) * x(5,i) * x(6,i) 
!
USE kinds, ONLY : DP
USE linear_solvers,     ONLY : linsolvx, linsolvms, linsolvsvd
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar, ndata
INTEGER, INTENT(INOUT) :: lsolve
REAL(DP), INTENT(IN) :: x(degree,ndata), f(ndata)
REAL(DP), INTENT(INOUT) :: coeff(nvar)

REAL(DP) :: amat(ndata,nvar), aa(nvar,nvar), b(nvar) 

INTEGER :: ivar, jvar, idata, nv

IF (degree>4.OR.degree<1) &
   CALL errore('multi_cubic','degree must be from 1 to 3',1)
IF (ndata < 3) &
   CALL errore('multi_cubic','Too few sampling data',1)
IF (ndata < nvar) &
   WRITE(stdout,'(/,5x,"Be careful: there are too few sampling data")')
!
!  prepare the auxiliary matrix
!
amat=0.0_DP

DO idata=1,ndata
   amat(idata,1) = 1.0_DP
   amat(idata,2) = x(1,idata)
   amat(idata,3) = x(1,idata)*x(1,idata)
   amat(idata,4) = x(1,idata)*x(1,idata)*x(1,idata)

   IF (degree>1) THEN
      amat(idata,5)  = x(2,idata)
      amat(idata,6)  = x(2,idata)*x(2,idata)
      amat(idata,7)  = x(2,idata)*x(2,idata)*x(2,idata)
      amat(idata,8)  = x(1,idata)*x(2,idata)
      amat(idata,9)  = x(1,idata)*x(2,idata)*x(2,idata)
      amat(idata,10) = x(1,idata)*x(1,idata)*x(2,idata)

   ENDIF

   IF (degree>2) THEN
      amat(idata,11) = x(3,idata)
      amat(idata,12) = x(3,idata)**2
      amat(idata,13) = x(3,idata)**3
      amat(idata,14) = x(1,idata)*x(3,idata)
      amat(idata,15) = x(1,idata)*x(3,idata)**2
      amat(idata,16) = x(1,idata)**2*x(3,idata)
      amat(idata,17) = x(2,idata)*x(3,idata)
      amat(idata,18) = x(2,idata)*x(3,idata)**2
      amat(idata,19) = x(2,idata)**2*x(3,idata)
      amat(idata,20) = x(1,idata)*x(2,idata)*x(3,idata)

   ENDIF

   IF (degree>3) THEN
      amat(idata,21) = x(4,idata)
      amat(idata,22) = x(4,idata)**2
      amat(idata,23) = x(4,idata)**3
      amat(idata,24) = x(1,idata)*x(4,idata)
      amat(idata,25) = x(1,idata)*x(4,idata)**2
      amat(idata,26) = x(1,idata)**2*x(4,idata)
      amat(idata,27) = x(2,idata)*x(4,idata)
      amat(idata,28) = x(2,idata)*x(4,idata)**2
      amat(idata,29) = x(2,idata)**2*x(4,idata)
      amat(idata,30) = x(3,idata)*x(4,idata)
      amat(idata,31) = x(3,idata)*x(4,idata)**2
      amat(idata,32) = x(3,idata)**2*x(4,idata)
      amat(idata,33) = x(1,idata)*x(2,idata)*x(4,idata)
      amat(idata,34) = x(1,idata)*x(3,idata)*x(4,idata)
      amat(idata,35) = x(2,idata)*x(3,idata)*x(4,idata)

   ENDIF

   IF (degree>4) THEN
      amat(idata,36) = x(5,idata)
      amat(idata,37) = x(5,idata)**2
      amat(idata,38) = x(5,idata)**3
      amat(idata,39) = x(1,idata)*x(5,idata)
      amat(idata,40) = x(1,idata)*x(5,idata)**2
      amat(idata,41) = x(1,idata)**2*x(5,idata)
      amat(idata,42) = x(2,idata)*x(5,idata)
      amat(idata,43) = x(2,idata)*x(5,idata)**2
      amat(idata,44) = x(2,idata)**2*x(5,idata)
      amat(idata,45) = x(3,idata)*x(5,idata)
      amat(idata,46) = x(3,idata)*x(5,idata)**2
      amat(idata,47) = x(3,idata)**2*x(5,idata)
      amat(idata,48) = x(4,idata)*x(5,idata)
      amat(idata,49) = x(4,idata)*x(5,idata)**2
      amat(idata,50) = x(4,idata)**2*x(5,idata)
      amat(idata,51) = x(1,idata)*x(2,idata)*x(5,idata)
      amat(idata,52) = x(1,idata)*x(3,idata)*x(5,idata)
      amat(idata,53) = x(1,idata)*x(4,idata)*x(5,idata)
      amat(idata,54) = x(2,idata)*x(3,idata)*x(5,idata)
      amat(idata,55) = x(2,idata)*x(4,idata)*x(5,idata)
      amat(idata,56) = x(3,idata)*x(4,idata)*x(5,idata)

   ENDIF

   IF (degree>5) THEN
      amat(idata,57) = x(6,idata)
      amat(idata,58) = x(6,idata)**2
      amat(idata,59) = x(6,idata)**3
      amat(idata,60) = x(1,idata)*x(6,idata)
      amat(idata,61) = x(1,idata)*x(6,idata)**2
      amat(idata,62) = x(1,idata)**2*x(6,idata)
      amat(idata,63) = x(2,idata)*x(6,idata)
      amat(idata,64) = x(2,idata)*x(6,idata)**2
      amat(idata,65) = x(2,idata)**2*x(6,idata)
      amat(idata,66) = x(3,idata)*x(6,idata)
      amat(idata,67) = x(3,idata)*x(6,idata)**2
      amat(idata,68) = x(3,idata)**2*x(6,idata)
      amat(idata,69) = x(4,idata)*x(6,idata)
      amat(idata,70) = x(4,idata)*x(6,idata)**2
      amat(idata,71) = x(4,idata)**2*x(6,idata)
      amat(idata,72) = x(5,idata)*x(6,idata)
      amat(idata,73) = x(5,idata)*x(6,idata)**2
      amat(idata,74) = x(5,idata)**2*x(6,idata)
      amat(idata,75) = x(1,idata)*x(2,idata)*x(6,idata)
      amat(idata,76) = x(1,idata)*x(3,idata)*x(6,idata)
      amat(idata,77) = x(1,idata)*x(4,idata)*x(6,idata)
      amat(idata,78) = x(1,idata)*x(5,idata)*x(6,idata)
      amat(idata,79) = x(2,idata)*x(3,idata)*x(6,idata)
      amat(idata,80) = x(2,idata)*x(4,idata)*x(6,idata)
      amat(idata,81) = x(2,idata)*x(5,idata)*x(6,idata)
      amat(idata,82) = x(3,idata)*x(4,idata)*x(6,idata)
      amat(idata,83) = x(3,idata)*x(5,idata)*x(6,idata)
      amat(idata,84) = x(4,idata)*x(5,idata)*x(6,idata)

   ENDIF

ENDDO

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
coeff=0.0_DP
IF (lsolve<1.OR.lsolve>3) lsolve=2
IF (lsolve==1) THEN
   WRITE(stdout,'(5x,"Finding the cubic polynomial using &
                                                   &nvar x nvar matrix")')  
   CALL linsolvx(aa,nvar,b,coeff)
ELSEIF(lsolve==2) THEN
   WRITE(stdout,'(5x,"Finding the cubic polynomial using &
                                                   &QR factorization")')  
   CALL linsolvms(amat,ndata,nvar,f,coeff)
ELSEIF(lsolve==3) THEN
   WRITE(stdout,'(5x,"Finding the cubic polynomial using SVD")')  
   CALL linsolvsvd(amat,ndata,nvar,f,coeff)
ENDIF

DO nv=1,nvar
   WRITE(stdout,'(5x, "coeff :",e15.7)') coeff(nv)
ENDDO

RETURN
END SUBROUTINE fit_multi_cubic

SUBROUTINE evaluate_fit_cubic(degree,nvar,x,f,coeff)

USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: x(degree)
REAL(DP), INTENT(IN) :: coeff(nvar)
REAL(DP), INTENT(INOUT) :: f

REAL(DP) :: aux
!
!  one variable
!
aux = coeff(1) + x(1)*(coeff(2)+x(1)*(coeff(3)+x(1)*coeff(4)))
!
!  two variables
!
IF (degree>1) THEN
   aux = aux + x(2)*(coeff(5)+x(2)*(coeff(6)+x(2)*coeff(7)+x(1)*coeff(9))  &
                     + x(1)*(coeff(8)+x(1)*coeff(10)))
ENDIF
!
!  three variabl
!
IF (degree>2) THEN
   aux = aux + coeff(11)*x(3) + coeff(12)*x(3)**2 + coeff(13)*x(3)**3                 &
             + coeff(14)*x(1)*x(3) + coeff(15)*x(1)*x(3)**2 + coeff(16)*x(1)**2*x(3)  &
             + coeff(17)*x(2)*x(3) + coeff(18)*x(2)*x(3)**2 + coeff(19)*x(2)**2*x(3)  &
             + coeff(20)*x(1)*x(2)*x(3)
ENDIF
!
!  four variables
!
IF (degree>3) THEN
   aux = aux + coeff(21)*x(4) + coeff(22)*x(4)**2 + coeff(23)*x(4)**3   &
             + coeff(24)*x(1)*x(4) + coeff(25)*x(1)*x(4)**2 + coeff(26)*x(1)**2*x(4)  &
             + coeff(27)*x(2)*x(4) + coeff(28)*x(2)*x(4)**2 + coeff(29)*x(2)**2*x(4)  &
             + coeff(30)*x(3)*x(4) + coeff(31)*x(3)*x(4)**2 + coeff(32)*x(3)**2*x(4)  &
             + coeff(33)*x(1)*x(2)*x(4)  &
             + coeff(34)*x(1)*x(3)*x(4)  &
             + coeff(35)*x(2)*x(3)*x(4)
ENDIF
!
!  five variables
!
IF (degree>4) THEN
   aux = aux + coeff(36)*x(5) + coeff(37)*x(5)**2 + coeff(38)*x(5)**3  &
             + coeff(39)*x(1)*x(5) + coeff(40)*x(1)*x(5)**2 + coeff(41)*x(1)**2*x(5)  &
             + coeff(42)*x(2)*x(5) + coeff(43)*x(2)*x(5)**2 + coeff(44)*x(2)**2*x(5)  &
             + coeff(45)*x(3)*x(5) + coeff(46)*x(3)*x(5)**2 + coeff(47)*x(3)**2*x(5)  &
             + coeff(48)*x(4)*x(5) + coeff(49)*x(4)*x(5)**2 + coeff(50)*x(4)**2*x(5)  &
             + coeff(51)*x(1)*x(2)*x(5)  &
             + coeff(52)*x(1)*x(3)*x(5)  &
             + coeff(53)*x(1)*x(4)*x(5)  &
             + coeff(54)*x(2)*x(3)*x(5)  &
             + coeff(55)*x(2)*x(4)*x(5)  &
             + coeff(56)*x(3)*x(4)*x(5)
ENDIF
!
!  six variables
!
IF (degree>5) THEN
   aux = aux + coeff(57)*x(6) + coeff(58)*x(6)**2 + coeff(59)*x(6)**3  &
             + coeff(60)*x(1)*x(6) + coeff(61)*x(1)*x(6)**2 + coeff(62)*x(1)**2*x(6)  &
             + coeff(63)*x(2)*x(6) + coeff(64)*x(2)*x(6)**2 + coeff(65)*x(2)**2*x(6)  &
             + coeff(66)*x(3)*x(6) + coeff(67)*x(3)*x(6)**2 + coeff(68)*x(3)**2*x(6)  &
             + coeff(69)*x(4)*x(6) + coeff(70)*x(4)*x(6)**2 + coeff(71)*x(4)**2*x(6)  &
             + coeff(72)*x(5)*x(6) + coeff(73)*x(5)*x(6)**2 + coeff(74)*x(5)**2*x(6)  &
             + coeff(75)*x(1)*x(2)*x(6)  &
             + coeff(76)*x(1)*x(3)*x(6)  &
             + coeff(77)*x(1)*x(4)*x(6)  &
             + coeff(78)*x(1)*x(5)*x(6)  &
             + coeff(79)*x(2)*x(3)*x(6)  &
             + coeff(80)*x(2)*x(4)*x(6)  &
             + coeff(81)*x(2)*x(5)*x(6)  &
             + coeff(82)*x(3)*x(4)*x(6)  &
             + coeff(83)*x(3)*x(5)*x(6)  &
             + coeff(84)*x(4)*x(5)*x(6)
          
ENDIF

f=aux

RETURN
END SUBROUTINE evaluate_fit_cubic

SUBROUTINE evaluate_fit_grad_cubic(degree,nvar,x,f,coeff)
!
!  computes the gradient of the cubic polynomial, up to 6 variables.
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: x(degree)
REAL(DP), INTENT(IN) :: coeff(nvar)
REAL(DP), INTENT(INOUT) :: f(degree)

REAL(DP) :: aux(degree)

IF (degree>6) CALL errore('evaluate_fit_grad_cubic','gradient not availble',1)

aux(1) = coeff(2) + 2.0_DP*coeff(3)*x(1) + 3.0_DP*coeff(4)*x(1)**2

IF (degree>1) THEN
   aux(1) = aux(1) + coeff(8)*x(2) + coeff(9)*x(2)**2 + 2.0_DP*coeff(10)*x(1)*x(2)
                     
   aux(2) = coeff(5) + 2.0_DP*coeff(6)*x(2) + 3.0_DP*coeff(7)*x(2)**2  &
            + coeff(8)*x(1) + 2.0_DP*coeff(9)*x(1)*x(2) + coeff(10)*x(1)**2 
ENDIF

IF (degree>2) THEN
   aux(1) = aux(1) + coeff(14)*x(3) + coeff(15)*x(3)**2 + 2.0_DP*coeff(16)*x(1)*x(3) &
                   + coeff(20)*x(2)*x(3) 
 
   aux(2) = aux(2) + coeff(17)*x(3) + coeff(18)*x(3)**2 + 2.0_DP*coeff(19)*x(2)*x(3) &
                   + coeff(20)*x(1)*x(3) 
 
   aux(3) = coeff(11) + 2.0_DP*coeff(12)*x(3) + 3.0_DP*coeff(13)*x(3)**2  &
            + coeff(14)*x(1) + 2.0_DP*coeff(15)*x(1)*x(3) + coeff(16)*x(1)**2  &
            + coeff(17)*x(2) + 2.0_DP*coeff(18)*x(2)*x(3) + coeff(19)*x(2)**2  &
            + coeff(20)*x(1)*x(2)
ENDIF


IF (degree>3) THEN
   aux(1) = aux(1) + coeff(24)*x(4) + coeff(25)*x(4)**2 + 2.0_DP*coeff(26)*x(1)*x(4)  &
                   + coeff(33)*x(2)*x(4)  &
                   + coeff(34)*x(3)*x(4)

   aux(2) = aux(2) + coeff(27)*x(4) + coeff(28)*x(4)**2 + 2.0_DP*coeff(29)*x(2)*x(4)  &
                   + coeff(33)*x(1)*x(4)  &
                   + coeff(35)*x(3)*x(4)

   aux(3) = aux(3) + coeff(30)*x(4) + coeff(31)*x(4)**2 + 2.0_DP*coeff(32)*x(3)*x(4)  &
                   + coeff(34)*x(1)*x(4)  &
                   + coeff(35)*x(2)*x(4) 

   aux(4) = coeff(21) + 2.0_DP*coeff(22)*x(4) + 3.0_DP*coeff(23)*x(4)**2  &
            + coeff(24)*x(1) + 2.0_DP*coeff(25)*x(1)*x(4) + coeff(26)*x(1)**2  &
            + coeff(27)*x(2) + 2.0_DP*coeff(28)*x(2)*x(4) + coeff(29)*x(2)**2  &
            + coeff(30)*x(3) + 2.0_DP*coeff(31)*x(3)*x(4) + coeff(32)*x(3)**2  &
            + coeff(33)*x(1)*x(2)  &
            + coeff(34)*x(1)*x(3)  &
            + coeff(35)*x(2)*x(3)

ENDIF

IF (degree>4) THEN
   aux(1) = aux(1) + coeff(39)*x(5) + coeff(40)*x(5)**2 + 2.0_DP*coeff(41)*x(1)*x(5)  &
                   + coeff(51)*x(2)*x(5)  &
                   + coeff(52)*x(3)*x(5)  &
                   + coeff(53)*x(4)*x(5)

   aux(2) = aux(2) + coeff(42)*x(5) + coeff(43)*x(5)**2 + 2.0_DP*coeff(44)*x(2)*x(5)  &
                   + coeff(51)*x(1)*x(5)  &
                   + coeff(54)*x(3)*x(5)  &
                   + coeff(55)*x(4)*x(5) 
                    
   aux(3) = aux(3) + coeff(45)*x(5) + coeff(46)*x(5)**2 + 2.0_DP*coeff(47)*x(3)*x(5)  &
                   + coeff(52)*x(1)*x(5)  &
                   + coeff(54)*x(2)*x(5)  &
                   + coeff(56)*x(4)*x(5)

   aux(4) = aux(4) + coeff(48)*x(5) + coeff(49)*x(5)**2 + 2.0_DP*coeff(50)*x(4)*x(5)  &
                   + coeff(53)*x(1)*x(5)  &
                   + coeff(55)*x(2)*x(5)  &
                   + coeff(56)*x(3)*x(5) 


   aux(5) = coeff(36) + 2.0_DP*coeff(37)*x(5) + 3.0_DP*coeff(38)*x(5)**2  &
            + coeff(39)*x(1) + 2.0_DP*coeff(40)*x(1)*x(5) + coeff(41)*x(1)**2  &
            + coeff(42)*x(2) + 2.0_DP*coeff(43)*x(2)*x(5) + coeff(44)*x(2)**2  &
            + coeff(45)*x(3) + 2.0_DP*coeff(46)*x(3)*x(5) + coeff(47)*x(3)**2  &
            + coeff(48)*x(4) + 2.0_DP*coeff(49)*x(4)*x(5) + coeff(50)*x(4)**2  &
            + coeff(51)*x(1)*x(2)  &
            + coeff(52)*x(1)*x(3)  &
            + coeff(53)*x(1)*x(4)  &
            + coeff(54)*x(2)*x(3)  &
            + coeff(55)*x(2)*x(4)  &
            + coeff(56)*x(3)*x(4)                
ENDIF

IF (degree>5) THEN
   aux(1) = aux(1) + coeff(60)*x(6) + coeff(61)*x(6)**2 + 2.0_DP*coeff(62)*x(1)*x(6)  &
                   + coeff(75)*x(2)*x(6)  & 
                   + coeff(76)*x(3)*x(6)  &
                   + coeff(77)*x(4)*x(6)  &
                   + coeff(78)*x(5)*x(6)

   aux(2) = aux(2) + coeff(63)*x(6) + coeff(64)*x(6)**2 + 2.0_DP*coeff(65)*x(2)*x(6)  &
                   + coeff(75)*x(1)*x(6)  &
                   + coeff(79)*x(3)*x(6)  &
                   + coeff(80)*x(4)*x(6)  &
                   + coeff(81)*x(5)*x(6)
                   
   aux(3) = aux(3) + coeff(66)*x(6) + coeff(67)*x(6)**2 + 2.0_DP*coeff(68)*x(3)*x(6)  &
                   + coeff(76)*x(1)*x(6)  &
                   + coeff(79)*x(2)*x(6)  &
                   + coeff(82)*x(4)*x(6)  &
                   + coeff(83)*x(5)*x(6) 

   aux(4) = aux(4) + coeff(69)*x(6) + coeff(70)*x(6)**2 + 2.0_DP*coeff(71)*x(4)*x(6)  &
                   + coeff(77)*x(1)*x(6)  &
                   + coeff(80)*x(2)*x(6)  &
                   + coeff(82)*x(3)*x(6)  &
                   + coeff(84)*x(5)*x(6)

   aux(5) = aux(5) + coeff(72)*x(6) + coeff(73)*x(6)**2 + 2.0_DP*coeff(74)*x(5)*x(6)  &
                   + coeff(78)*x(1)*x(6)  &
                   + coeff(81)*x(2)*x(6)  &
                   + coeff(83)*x(3)*x(6)  &
                   + coeff(84)*x(4)*x(6)


   aux(6) = coeff(57) + 2.0_DP*coeff(58)*x(6) + 3.0_DP*coeff(59)*x(6)**2  &
            + coeff(60)*x(1) + 2.0_DP*coeff(61)*x(1)*x(6) + coeff(62)*x(1)**2  &
            + coeff(63)*x(2) + 2.0_DP*coeff(64)*x(2)*x(6) + coeff(65)*x(2)**2  &
            + coeff(66)*x(3) + 2.0_DP*coeff(67)*x(3)*x(6) + coeff(68)*x(3)**2  &
            + coeff(69)*x(4) + 2.0_DP*coeff(70)*x(4)*x(6) + coeff(71)*x(4)**2  &
            + coeff(72)*x(5) + 2.0_DP*coeff(73)*x(5)*x(6) + coeff(74)*x(5)**2  &
            + coeff(75)*x(1)*x(2)  &
            + coeff(76)*x(1)*x(3)  &
            + coeff(77)*x(1)*x(4)  &
            + coeff(78)*x(1)*x(5)  &
            + coeff(79)*x(2)*x(3)  &
            + coeff(80)*x(2)*x(4)  &
            + coeff(81)*x(2)*x(5)  &
            + coeff(82)*x(3)*x(4)  &
            + coeff(83)*x(3)*x(5)  &
            + coeff(84)*x(4)*x(5)                     
ENDIF

f=aux

RETURN
END SUBROUTINE evaluate_fit_grad_cubic

SUBROUTINE evaluate_fit_hess_cubic(degree,nvar,x,f,coeff)

USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: x(degree)
REAL(DP), INTENT(IN) :: coeff(nvar)
REAL(DP), INTENT(INOUT) :: f(degree,degree)

REAL(DP) :: aux(degree,degree)

aux(1,1) = 2.0_DP*coeff(3) + 6.0_DP*coeff(4)*x(1)

IF (degree>1) THEN
   aux(1,1) = aux(1,1) + 2.0_DP*coeff(10)*x(2) 
   aux(1,2) = coeff(8) + 2.0_DP*coeff(9)*x(2) + 2.0_DP*coeff(10)*x(1)
   aux(2,1) = aux(1,2)
   aux(2,2) = 2.0_DP*coeff(6) + 6.0_DP*coeff(7)*x(2) + 2.0_DP*coeff(9)*x(1) 
ENDIF

IF (degree>2) THEN
   aux(1,1) = aux(1,1) + 2.0_DP*coeff(16)*x(3) 
   aux(1,2) = aux(1,2) + coeff(20)*x(3)     
   aux(2,1) = aux(1,2)
   aux(1,3) = coeff(14) + 2.0_DP*coeff(15)*x(3) + 2.0_DP*coeff(16)*x(1)  &
              + coeff(20)*x(2) 
   aux(3,1) = aux(1,3)
   aux(2,2) = aux(2,2) + 2.0_DP*coeff(19)*x(3) 
   aux(2,3) = coeff(17) + 2.0_DP*coeff(18)*x(3) + 2.0_DP*coeff(19)*x(2)  &
              + coeff(20)*x(1)
   aux(3,2) = aux(2,3)
   aux(3,3) = 2.0_DP*coeff(12) + 6.0_DP*coeff(13)*x(3)  &
              + 2.0_DP*coeff(15)*x(1)  &
              + 2.0_DP*coeff(18)*x(2)
ENDIF

IF (degree>3) THEN
   aux(1,1) = aux(1,1) + 2.0_DP*coeff(26)*x(4)
   aux(1,2) = aux(1,2) + coeff(33)*x(4)
   aux(2,1) = aux(1,2)
   aux(1,3) = aux(1,3) + coeff(34)*x(4)
   aux(3,1) = aux(1,3)
   aux(1,4) = coeff(24) + 2.0_DP*coeff(25)*x(4) + 2.0_DP*coeff(26)*x(1)  &
              + coeff(33)*x(2)  &
              + coeff(34)*x(3)
   aux(4,1) = aux(1,4)
   aux(2,2) = aux(2,2) + 2.0_DP*coeff(29)*x(4)
   aux(2,3) = aux(2,3) + coeff(35)*x(4)  
   aux(3,2) = aux(2,3)
   aux(2,4) = coeff(27) + 2.0_DP*coeff(28)*x(4) + 2.0_DP*coeff(29)*x(2)  &
              + coeff(33)*x(1)  &
              + coeff(35)*x(3) 
   aux(4,2) = aux(2,4)
   aux(3,3) = aux(3,3) + 2.0_DP*coeff(32)*x(4) 
   aux(3,4) = coeff(30) + 2.0_DP*coeff(31)*x(4) + 2.0_DP*coeff(32)*x(3)  &
              + coeff(34)*x(1)  &
              + coeff(35)*x(2) 
   aux(4,3) = aux(3,4)
   aux(4,4) = 2.0_DP*coeff(22) + 6.0_DP*coeff(23)*x(4)  &
              + 2.0_DP*coeff(25)*x(1)  &
              + 2.0_DP*coeff(28)*x(2)  &
              + 2.0_DP*coeff(31)*x(3)
ENDIF

IF (degree>4) THEN
   aux(1,1) = aux(1,1) + 2.0_DP*coeff(41)*x(5) 
   aux(1,2) = aux(1,2) + coeff(51)*x(5) 
   aux(2,1) = aux(1,2)
   aux(1,3) = aux(1,3) + coeff(52)*x(5)
   aux(3,1) = aux(1,3)
   aux(1,4) = aux(1,4) + coeff(53)*x(5)
   aux(4,1) = aux(1,4)
   aux(1,5) = coeff(39) + 2.0_DP*coeff(40)*x(5) + 2.0_DP*coeff(41)*x(1)  &
              + coeff(51)*x(2)  &
              + coeff(52)*x(3)  &
              + coeff(53)*x(4)
   aux(5,1) = aux(1,5)
   aux(2,2) = aux(2,2) + 2.0_DP*coeff(44)*x(5)  
   aux(2,3) = aux(2,3) + coeff(54)*x(5) 
   aux(3,2) = aux(2,3)
   aux(2,4) = aux(2,4) + coeff(55)*x(5) 
   aux(4,2) = aux(2,4)
   aux(2,5) = coeff(42) + 2.0_DP*coeff(43)*x(5) +  2.0_DP*coeff(44)*x(2)  &
              + coeff(51)*x(1)  &
              + coeff(54)*x(3)  &
              + coeff(55)*x(4)
   aux(5,2) = aux(2,5)
   aux(3,3) = aux(3,3) + 2.0_DP*coeff(47)*x(5)
   aux(3,4) = aux(3,4) + coeff(56)*x(5)
   aux(4,3) = aux(3,4)
   aux(3,5) = coeff(45) + 2.0_DP*coeff(46)*x(5) + 2.0_DP*coeff(47)*x(3)  &
              + coeff(52)*x(1)  &
              + coeff(54)*x(2)  &
              + coeff(56)*x(4) 
   aux(5,3) = aux(3,5)
   aux(4,4) = aux(4,4) + 2.0_DP*coeff(50)*x(5) 
   aux(4,5) = coeff(48) + 2.0_DP*coeff(49)*x(5) + 2.0_DP*coeff(50)*x(4)  &
              + coeff(53)*x(1)  &
              + coeff(55)*x(2)  &
              + coeff(56)*x(3)
   aux(5,4) = aux(4,5)
   aux(5,5) = 2.0_DP*coeff(37) + 6.0_DP*coeff(38)*x(5)  &
              + 2.0_DP*coeff(40)*x(1)  &
              + 2.0_DP*coeff(43)*x(2)  &
              + 2.0_DP*coeff(46)*x(3)  &
              + 2.0_DP*coeff(49)*x(4) 
ENDIF

IF (degree>5) THEN
   aux(1,1) = aux(1,1) + 2.0_DP*coeff(62)*x(6) 
   aux(1,2) = aux(1,2) + coeff(75)*x(6)
   aux(2,1) = aux(1,2)
   aux(1,3) = aux(1,3) + coeff(76)*x(6) 
   aux(3,1) = aux(1,3)
   aux(1,4) = aux(1,4) + coeff(77)*x(6)
   aux(4,1) = aux(1,4)
   aux(1,5) = aux(1,5) + coeff(78)*x(6) 
   aux(5,1) = aux(1,5)
   aux(1,6) = coeff(60) + 2.0_DP*coeff(61)*x(6) + 2.0_DP*coeff(62)*x(1)  &
              + coeff(75)*x(2)  &
              + coeff(76)*x(3)  &
              + coeff(77)*x(4)  &
              + coeff(78)*x(5) 
   aux(6,1) = aux(1,6)
   aux(2,2) = aux(2,2) + 2.0_DP*coeff(65)*x(6)
   aux(2,3) = aux(2,3) + coeff(79)*x(6)
   aux(3,2) = aux(2,3)
   aux(2,4) = aux(2,4) + coeff(80)*x(6)
   aux(4,2) = aux(2,4)
   aux(2,5) = aux(2,5) + coeff(81)*x(6) 
   aux(5,2) = aux(2,5)
   aux(2,6) = coeff(63) + 2.0_DP*coeff(64)*x(6) + 2.0_DP*coeff(65)*x(2)  &
              + coeff(75)*x(1)  &
              + coeff(79)*x(3)  &
              + coeff(80)*x(4)  &
              + coeff(81)*x(5)
   aux(6,2) = aux(2,6)
   aux(3,3) = aux(3,3) + 2.0_DP*coeff(68)*x(6) 
   aux(3,4) = aux(3,4) + coeff(82)*x(6)
   aux(4,3) = aux(3,4)
   aux(3,5) = aux(3,5) + coeff(83)*x(6)  
   aux(5,3) = aux(3,5)
   aux(3,6) = coeff(66) + 2.0_DP*coeff(67)*x(6) + 2.0_DP*coeff(68)*x(3)  &
              + coeff(76)*x(1)  &
              + coeff(79)*x(2)  &
              + coeff(82)*x(4)  &
              + coeff(83)*x(5)
   aux(6,3) = aux(3,6)
   aux(4,4) = aux(4,4) + 2.0_DP*coeff(71)*x(6)
   aux(4,5) = aux(4,5) + coeff(84)*x(6)
   aux(5,4) = aux(4,5)
   aux(4,6) = coeff(69) + 2.0_DP*coeff(70)*x(6) + 2.0_DP*coeff(71)*x(4)  &
              + coeff(77)*x(1)  &
              + coeff(80)*x(2)  &
              + coeff(82)*x(3)  &
              + coeff(84)*x(5) 
   aux(6,4) = aux(4,6)
   aux(5,5) = aux(5,5) + 2.0_DP*coeff(74)*x(6) 
   aux(5,6) = coeff(72) + 2.0_DP*coeff(73)*x(6) + 2.0_DP*coeff(74)*x(5)  &
              + coeff(78)*x(1)  &
              + coeff(81)*x(2)  &
              + coeff(83)*x(3)  &
              + coeff(84)*x(4) 
   aux(6,5) = aux(5,6)
   aux(6,6) = 2.0_DP*coeff(58) + 6.0_DP*coeff(59)*x(6)  &
              + 2.0_DP*coeff(61)*x(1)  &
              + 2.0_DP*coeff(64)*x(2)  &
              + 2.0_DP*coeff(67)*x(3)  &
              + 2.0_DP*coeff(70)*x(4)  &    
              + 2.0_DP*coeff(73)*x(3)
ENDIF

f(:,:)=aux(:,:)

RETURN
END SUBROUTINE evaluate_fit_hess_cubic

SUBROUTINE find_cubic_extremum(degree,nvar,x,f,coeff)
!
!  This routine starts from the point x and finds the extremum closest
!  to x. In output x are the coordinates of the extremum and f 
!  the value of the cubic function at the minimum
!
USE linear_solvers, ONLY : linsolvx
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar
REAL(DP),INTENT(INOUT) :: x(degree), f
REAL(DP),INTENT(IN) :: coeff(nvar)

INTEGER, PARAMETER :: maxiter=100

INTEGER :: iter, ideg
REAL(DP), PARAMETER :: tol=2.D-11
REAL(DP) :: g(degree), y(degree), xold(degree)
REAL(DP) :: j(degree, degree) 
REAL(DP) :: deltax, fmod

xold(:)=x(:)
DO iter=1,maxiter
   !
   CALL evaluate_fit_grad_cubic(degree,nvar,x,g,coeff)
   !
   CALL evaluate_fit_hess_cubic(degree,nvar,x,j,coeff)
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
CALL errore('find_cubic_extremum','extremum not found',1)
100 CONTINUE
CALL evaluate_fit_cubic(degree,nvar,x,f,coeff)

RETURN
END SUBROUTINE find_cubic_extremum

SUBROUTINE find_quartic_cubic_extremum(degree,nvar4,nvar,x,f,coeff4,coeff)
!
!   This subroutines adds the coefficients of a quartic polynomium
!   to those of a cubic polynomium and finds the minimum of the sum
!   of the two
USE quartic_surfaces, ONLY : find_quartic_extremum
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar, nvar4
REAL(DP), INTENT(IN) :: coeff(nvar), coeff4(nvar4) 
REAL(DP), INTENT(INOUT) :: x(degree), f

REAL(DP) :: coeffadd4(nvar4)

CALL set_quartic_cubic_coefficients(degree, nvar4, nvar, coeffadd4, &
                                                   coeff4, coeff)

CALL find_quartic_extremum(degree,nvar4,x,f,coeffadd4)

RETURN
END SUBROUTINE find_quartic_cubic_extremum

SUBROUTINE evaluate_quartic_cubic(degree,nvar4,nvar,x,f,coeff4,coeff)

USE quartic_surfaces, ONLY : evaluate_fit_quartic
IMPLICIT NONE

INTEGER, INTENT(IN) :: degree, nvar, nvar4
REAL(DP), INTENT(IN) :: coeff(nvar), coeff4(nvar4) 
REAL(DP), INTENT(INOUT) :: x(degree), f

REAL(DP) :: coeffadd4(nvar4)

CALL set_quartic_cubic_coefficients(degree, nvar4, nvar, coeffadd4, coeff4, &
                                                                    coeff)

CALL evaluate_fit_quartic(degree,nvar4,x,f,coeffadd4)

RETURN
END SUBROUTINE evaluate_quartic_cubic

SUBROUTINE set_quartic_cubic_coefficients(degree, nvar4, nvar, coeffadd4, &
                                                               coeff4, coeff)
IMPLICIT NONE

INTEGER, INTENT(IN) :: degree, nvar4, nvar
REAL(DP), INTENT(IN) :: coeff4(nvar4), coeff(nvar)
REAL(DP), INTENT(INOUT) :: coeffadd4(nvar4)

coeffadd4 = coeff4
coeffadd4(1) = coeffadd4(1) + coeff(1)
coeffadd4(2) = coeffadd4(2) + coeff(2)
coeffadd4(3) = coeffadd4(3) + coeff(3)
coeffadd4(4) = coeffadd4(4) + coeff(4)

IF (degree > 1) THEN
   coeffadd4(6)  =  coeffadd4(6) +  coeff(5)
   coeffadd4(7)  =  coeffadd4(7) +  coeff(6)
   coeffadd4(8)  =  coeffadd4(8) +  coeff(7)
   coeffadd4(10) = coeffadd4(10) +  coeff(8)
   coeffadd4(11) = coeffadd4(11) +  coeff(9)
   coeffadd4(13) = coeffadd4(13) + coeff(10)
ENDIF

IF (degree > 2) THEN
   coeffadd4(16) = coeffadd4(16) + coeff(11)
   coeffadd4(17) = coeffadd4(17) + coeff(12)
   coeffadd4(18) = coeffadd4(18) + coeff(13)
   coeffadd4(20) = coeffadd4(20) + coeff(14)
   coeffadd4(21) = coeffadd4(21) + coeff(15)
   coeffadd4(23) = coeffadd4(23) + coeff(16)
   coeffadd4(26) = coeffadd4(26) + coeff(17)
   coeffadd4(20) = coeffadd4(20) + coeff(14)
   coeffadd4(27) = coeffadd4(27) + coeff(18)
   coeffadd4(29) = coeffadd4(29) + coeff(19)
   coeffadd4(32) = coeffadd4(32) + coeff(20)
END IF

IF (degree > 3) THEN
   coeffadd4(36) = coeffadd4(36) + coeff(21)
   coeffadd4(37) = coeffadd4(37) + coeff(22)
   coeffadd4(38) = coeffadd4(38) + coeff(23)
   coeffadd4(40) = coeffadd4(40) + coeff(24)
   coeffadd4(41) = coeffadd4(41) + coeff(25)
   coeffadd4(43) = coeffadd4(43) + coeff(26)
   coeffadd4(46) = coeffadd4(46) + coeff(27)
   coeffadd4(47) = coeffadd4(47) + coeff(28)
   coeffadd4(49) = coeffadd4(49) + coeff(29)
   coeffadd4(52) = coeffadd4(52) + coeff(30)
   coeffadd4(53) = coeffadd4(53) + coeff(31)
   coeffadd4(55) = coeffadd4(55) + coeff(32)
   coeffadd4(58) = coeffadd4(58) + coeff(33)
   coeffadd4(62) = coeffadd4(62) + coeff(34)
   coeffadd4(66) = coeffadd4(66) + coeff(35)
END IF

IF (degree > 4) THEN
   coeffadd4(71) = coeffadd4(71) + coeff(36)
   coeffadd4(72) = coeffadd4(72) + coeff(37)
   coeffadd4(73) = coeffadd4(73) + coeff(38)
   coeffadd4(75) = coeffadd4(75) + coeff(39)
   coeffadd4(76) = coeffadd4(76) + coeff(40)
   coeffadd4(78) = coeffadd4(78) + coeff(41)
   coeffadd4(81) = coeffadd4(81) + coeff(42)
   coeffadd4(82) = coeffadd4(82) + coeff(43)
   coeffadd4(84) = coeffadd4(84) + coeff(44)
   coeffadd4(87) = coeffadd4(87) + coeff(45)
   coeffadd4(88) = coeffadd4(88) + coeff(46)
   coeffadd4(90) = coeffadd4(90) + coeff(47)
   coeffadd4(93) = coeffadd4(93) + coeff(48)
   coeffadd4(94) = coeffadd4(94) + coeff(49)
   coeffadd4(96) = coeffadd4(96) + coeff(50)
   coeffadd4(99) = coeffadd4(99) + coeff(51)
   coeffadd4(103) = coeffadd4(103) + coeff(52)
   coeffadd4(107) = coeffadd4(107) + coeff(53)
   coeffadd4(111) = coeffadd4(111) + coeff(54)
   coeffadd4(115) = coeffadd4(115) + coeff(55)
   coeffadd4(119) = coeffadd4(119) + coeff(56)
END IF

IF (degree > 5) THEN
   coeffadd4(127) = coeffadd4(127) + coeff(57)
   coeffadd4(128) = coeffadd4(128) + coeff(58)
   coeffadd4(129) = coeffadd4(129) + coeff(59)
   coeffadd4(131) = coeffadd4(131) + coeff(60)
   coeffadd4(132) = coeffadd4(132) + coeff(61)
   coeffadd4(134) = coeffadd4(134) + coeff(62)
   coeffadd4(137) = coeffadd4(137) + coeff(63)
   coeffadd4(138) = coeffadd4(138) + coeff(64)
   coeffadd4(140) = coeffadd4(140) + coeff(65)
   coeffadd4(143) = coeffadd4(143) + coeff(66)
   coeffadd4(144) = coeffadd4(144) + coeff(67)
   coeffadd4(146) = coeffadd4(146) + coeff(68)
   coeffadd4(149) = coeffadd4(149) + coeff(69)
   coeffadd4(150) = coeffadd4(150) + coeff(70)
   coeffadd4(152) = coeffadd4(152) + coeff(71)
   coeffadd4(155) = coeffadd4(155) + coeff(72)
   coeffadd4(156) = coeffadd4(156) + coeff(73)
   coeffadd4(158) = coeffadd4(158) + coeff(74)
   coeffadd4(161) = coeffadd4(161) + coeff(75)
   coeffadd4(165) = coeffadd4(165) + coeff(76)
   coeffadd4(169) = coeffadd4(169) + coeff(77)
   coeffadd4(173) = coeffadd4(173) + coeff(78)
   coeffadd4(177) = coeffadd4(177) + coeff(79)
   coeffadd4(181) = coeffadd4(181) + coeff(80)
   coeffadd4(185) = coeffadd4(185) + coeff(81)
   coeffadd4(189) = coeffadd4(189) + coeff(82)
   coeffadd4(193) = coeffadd4(193) + coeff(83)
   coeffadd4(197) = coeffadd4(197) + coeff(84)
END IF

RETURN
END SUBROUTINE set_quartic_cubic_coefficients

FUNCTION cubic_var(degree)  

IMPLICIT NONE
INTEGER :: cubic_var
INTEGER, INTENT(IN) :: degree

IF (degree==1) THEN
   cubic_var=4
ELSEIF (degree==2) THEN
   cubic_var=10
ELSEIF (degree==3) THEN
   cubic_var=20
ELSEIF (degree==4) THEN
   cubic_var=35
ELSEIF (degree==5) THEN
   cubic_var=56
ELSEIF (degree==6) THEN
   cubic_var=84
ELSE
   cubic_var=0
ENDIF

RETURN
END FUNCTION cubic_var

SUBROUTINE print_cubic_polynomial(degree, nvar, coeff)

USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE

INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: coeff(nvar)

  WRITE(stdout,'(/,5x,"Cubic polynomial:",/)') 
  WRITE(stdout,'(5x,    e20.7,11x,"+",e20.7," x1 ")') coeff(1), coeff(2) 
  WRITE(stdout,'(4x,"+",e20.7," x1^2",6x,"+",e20.7," x1^3")') coeff(3), coeff(4)
  
  IF (degree>1) THEN
     WRITE(stdout,'(4x,"+",e20.7," x2",8x,"+",e20.7," x2^2")')  coeff(5), &
                                                                coeff(6)
     WRITE(stdout,'(4x,"+",e20.7," x2^3",6x,"+",e20.7," x1 x2")') coeff(7), &
                                                                  coeff(8)
     WRITE(stdout,'(4x,"+",e20.7," x1 x2^2",6x,"+",e20.7," x1^2 x2")') coeff(9), &
                                                                       coeff(10)
  ENDIF

  IF (degree>2) THEN
     WRITE(stdout,'(4x,"+",e15.7," x3       +",e15.7," x3^2      +",&
                        &e15.7," x3^3")') coeff(11), coeff(12), coeff(13)
     WRITE(stdout,'(4x,"+",e15.7," x1 x3     +",e15.7," x1 x3^2   +",e13.7,&
                              &" x1^2 x3")') coeff(14), coeff(15), coeff(16)
     WRITE(stdout,'(4x,"+",e15.7," x2 x3    +",e15.7," x2 x3^2   +",& 
                          &e15.7," x2^2 x3  ")') coeff(17), coeff(18), coeff(19)
     WRITE(stdout,'(4x,"+",e15.7," x1 x2 x3  ")') coeff(20)
  ENDIF

  IF (degree>3) THEN
     WRITE(stdout,'(4x,"+",e15.7," x4       +",e15.7," x4^2      +",&
                        &e15.7," x4^3")') coeff(21), coeff(22), coeff(23)
     WRITE(stdout,'(4x,"+",e15.7," x1 x4    +",e15.7," x1 x4^2   +",e15.7,&
                        &" x1^2 x4")') coeff(24), coeff(25), coeff(26)
     WRITE(stdout,'(4x,"+",e15.7," x2 x4  +",e15.7," x2 x4^2   +",& 
                        &e15.7," x2^2 x4")') coeff(27), coeff(28), coeff(29)
     WRITE(stdout,'(4x,"+",e15.7," x3 x4    +",e15.7," x3 x4^2   +",& 
                        &e15.7," x3^2 x4")') coeff(30), coeff(31), coeff(32)
     WRITE(stdout,'(4x,"+",e15.7," x1 x2 x4 +",e15.7," x1 x3 x4  +",& 
                        &e15.7," x2 x3 x4")') coeff(33), coeff(34), coeff(35)
  ENDIF

  IF (degree>4) THEN
     WRITE(stdout,'(4x,"+",e15.7," x5       +",e15.7," x5^2      +",&
                          &e15.7," x5^3")') coeff(36), coeff(37), coeff(38)
     WRITE(stdout,'(4x,"+",e15.7," x1 x5    +",e15.7," x1 x5^2   +",e15.7,&
                          &" x1^2 x5")') coeff(39), coeff(40), coeff(41)
     WRITE(stdout,'(4x,"+",e15.7," x2 x5    +",e15.7," x2 x5^2   +",& 
                          &e15.7," x2^2 x5")') coeff(42), coeff(43), coeff(44)
     WRITE(stdout,'(4x,"+",e15.7," x3 x5    +",e15.7," x3 x5^2   +",&
               &e15.7," x3^2 x5")') coeff(45), coeff(46), coeff(47)
     WRITE(stdout,'(4x,"+",e15.7," x4 x5    +",e15.7," x4 x5^2   +",& 
                          &e15.7," x4^2 x5")') coeff(48), coeff(49), coeff(50)
     WRITE(stdout,'(4x,"+",e15.7," x1 x2 x5 +",e15.7," x1 x3 x5 +",& 
                          &e15.7," x1 x4 x5")') coeff(51), coeff(52), coeff(53)
     WRITE(stdout,'(4x,"+",e15.7," x2 x3 x5 +",e15.7," x2 x4 x5 +",& 
                          &e15.7," x3 x4 x5")') coeff(54), coeff(55), coeff(56)
  ENDIF

  IF (degree>5) THEN
     WRITE(stdout,'(4x,"+",e15.7," x6       +",e15.7," x6^2      +",      &
                        &e15.7," x6^3")') coeff(57), coeff(58), coeff(59)
     WRITE(stdout,'(4x,"+",e15.7," x1 x6    +",e15.7," x1 x6^2   +",e15.7,&
                      &" x1^2 x6")') coeff(60), coeff(61), coeff(62)
     WRITE(stdout,'(4x,"+",e15.7," x2 x6    +",e15.7," x2 x6^2   +",      & 
               &e15.7," x2^2 x6")') coeff(63), coeff(64), coeff(65)
     WRITE(stdout,'(4x,"+",e15.7," x3 x6    +",e15.7," x3 x6^2   +",      &
               &e15.7," x3^2 x6  ")') coeff(66), coeff(67), coeff(68)
     WRITE(stdout,'(4x,"+",e15.7,"  x4 x6  +",e15.7," x4 x6^2 +",      & 
               &e15.7," x4^2 x6")') coeff(69), coeff(70), coeff(71)
     WRITE(stdout,'(4x,"+",e15.7," x5 x6    +",e15.7," x5 x6^2   +",      &
               &e15.7," x5^2 x6  ")') coeff(72), coeff(73), coeff(74)
     WRITE(stdout,'(4x,"+",e15.7," x1 x2 x6  +",e15.7," x1 x3 x6 +",      & 
               &e15.7," x1 x4 x6")') coeff(75), coeff(76), coeff(77)
     WRITE(stdout,'(4x,"+",e15.7," x1 x5 x6   +",e15.7," x2 x3 x6  +",      &
               &e15.7," x2 x4 x6 ")') coeff(78), coeff(79), coeff(80)
     WRITE(stdout,'(4x,"+",e15.7," x2 x5 x6  +",e15.7," x3 x4 x6 +",      & 
               &e15.7," x3 x5 x6")') coeff(81), coeff(82), coeff(83)
     WRITE(stdout,'(4x,"+",e15.7," x4 x5 x6")') coeff(84)
  ENDIF

RETURN
END SUBROUTINE print_cubic_polynomial

SUBROUTINE introduce_cubic_fit(degree, nvar, ndata)
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar, ndata

WRITE(stdout,'(/,5x,"Fitting the data with a cubic polynomial:")')

WRITE(stdout,'(/,5x,"Number of variables:",8x,i5)')  degree
WRITE(stdout,'(5x,"Coefficients of the cubic:",2x,i5)')  nvar
WRITE(stdout,'(5x,"Number of fitting data:",5x,i5,/)')  ndata

RETURN
END SUBROUTINE introduce_cubic_fit


SUBROUTINE print_chisq_cubic(ndata, degree, nvar, x, f, coeff)
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER  :: ndata, degree, nvar
REAL(DP) :: x(degree, ndata), f(ndata), coeff(nvar)

REAL(DP) :: chisq, perc, aux
INTEGER  :: idata

chisq=0.0_DP
perc=0.0_DP
DO idata=1,ndata
   CALL evaluate_fit_cubic(degree,nvar,x(1,idata),aux,coeff)
   WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
   IF (ABS(f(idata))>1.D-12) perc= perc + ABS((f(idata)-aux) / f(idata))
ENDDO

WRITE(stdout,'(5x,"chi square cubic=",e18.5," relative error",e18.5,&
                                     &" %",/)') chisq, perc / ndata
RETURN
END SUBROUTINE print_chisq_cubic

SUBROUTINE print_chisq_quartic_cubic(ndata, degree, nvar4, nvar, x, &
                                                          f, coeff4, coeff)
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER  :: ndata, degree, nvar4, nvar
REAL(DP) :: x(degree, ndata), f(ndata), coeff4(nvar4), coeff(nvar)

REAL(DP) :: chisq, aux
INTEGER  :: idata

chisq=0.0_DP
DO idata=1,ndata
   CALL evaluate_quartic_cubic(degree, nvar4, nvar, x(1,idata), aux, &
                                                        coeff4, coeff)
!  WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
ENDDO
WRITE(stdout,'(5x,"chi square quartic cubic=",e18.5,/)') chisq

RETURN
END SUBROUTINE print_chisq_quartic_cubic

END MODULE cubic_surfaces
