!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE quartic_surfaces
!
!   This module contains the support routines for dealing with quartic
!   surfaces interpolation with 1 up to 6 variables. 
!
!   It provides the following routines:
!
!   fit_multi_quartic receives as input the coordinates of some
!   points and the value of a function in these points and
!   finds the coefficients of a quartic interpolating polynomial.
!
!   evaluate_fit_quartic, given the coordinates of a point, and the
!   coefficients of the quartic polynomial, evaluates the quartic polynomial
!   at that point.
!
!   evaluate_quartic_grad, given the coordinates of a point, and the
!   coefficients of the quartic polynomial, evaluates the gradient of the 
!   function at that point.
!
!   evaluate_quartic_hessian, given the coordinates of a point, and the
!   coefficients of the quartic polynomial, evaluates the hessian of the 
!   function at that point.
!
!   quartic_ncoeff, given the number of variables gives the number of
!                   coefficients of the quartic polynomial
!
!   find_quartic_extremum, find the extremum closest to the input point
!                   of the quartic polynomium
!
!   print_quartic_polynomium, write on output the coefficients of the
!                   polynomium
!
!   introduce_quartic_fit, gives a few information on output on the 
!                  number of data used for the fit and the number of
!                  coefficients to be found
!
!   print_chisq_quartic, writes on output the chi square of a given 
!                  quartic polynomial interpolation
!   
!   evaluate_two_quartic, adds the coefficients of two quartic before
!                  evaluating them at a input point x
!  
!   find_two_quartic_extremum, find the extremum of two quartic starting
!                  the search from the input x point 
!
!   print_chisq_two_quartic, writes on output the chi square of the sum
!                  of two quartic polynomial interpolation
!
!   set_quartic_linear_coefficients, adds the coefficients of a 
!                  quartic and a linear polynomial
!
!   find_quartic_linear_extremum, finds the extremum of the sum of
!                  a quartic and a linear polynomial starting from the
!                  point x given as input
!
!   set_quartic_quadratic_coefficients, adds the coefficients of a 
!                  quartic and a quadratic polynomial
!
!   evaluate_quartic_quadratic, evaluate the sum of a quartic and a 
!                  quadratic polynomial in a input x point
!
!   find_quartic_quadratic_extremum, finds the extremum of the sum of
!                  a quartic and a quadratic polynomial starting from the
!                  point x given as input
!   
!   print_chisq_quartic_quadratic, writes on output the chi square of the sum
!                  of a quartic and quadratic polynomial interpolation
!
!   set_quartic_cubic_coefficients, adds the coefficients of a 
!                  quartic and a cubic polynomial
!
!   evaluate_quartic_cubic, evaluate the sum of a quartic and a 
!                  cubic polynomial in an input x point
!
!   find_quartic_cubic_extremum, finds the extremum of the sum of
!                  a quartic and a cubic polynomial starting from the
!                  point x given as input
!   
!   print_chisq_quartic_cubic, writes on output the chi square of the sum
!                  of a quartic and cubic polynomial interpolation
!
!
!   The following number of coefficients are necessary, depending on the
!   dimension of the quadratic function
!
!   dimension     number of coefficients
!      1                    5
!      2                    15
!      3                    35
!      4                    70
!      5                   126
!      6                   210
!    To interpolate the function it is better to give to
!    multi_quartic a number of points equal or larger than
!    the number of coefficients. The routines makes a least square fit
!    of the data.
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: fit_multi_quartic, evaluate_fit_quartic,         &
            evaluate_quartic_grad, evaluate_quartic_hessian, &
            quartic_ncoeff, find_quartic_extremum,           &
            print_quartic_polynomial, introduce_quartic_fit, &
            print_chisq_quartic,                             &
            evaluate_two_quartic,                            &
            find_two_quartic_extremum,                       &
            print_chisq_two_quartic,                         &
            set_quartic_linear_coefficients,                 &
            find_quartic_linear_extremum,                    &
            set_quartic_quadratic_coefficients,              &
            evaluate_quartic_quadratic,                      &
            find_quartic_quadratic_extremum,                 &
            print_chisq_quartic_quadratic,                   &
            set_quartic_cubic_coefficients,                  &
            evaluate_quartic_cubic,                          &
            find_quartic_cubic_extremum,                     &
            print_chisq_quartic_cubic                       

CONTAINS

SUBROUTINE fit_multi_quartic(ndata,nvar,ncoeff,lsolve,x,f,coeff)
!
!  This routine receives as input a set of vectors x(nvar,ndata) and
!  function values f(ndata) and gives as output the coefficients of
!  a quartic interpolating polynomial coeff(ncoeff). In input
!  ndata is the number of data points. nvar is the number of variables
!  (the maximum is 6), and ncoeff is the number of coefficients of the
!  intepolating quartic polynomial. 
!  number of variables    number of coeffients
!          1                     5
!          2                    15
!          3                    35
!          4                    70
!          5                   126
!          6                   210
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
!  The coefficients are organized as follows:
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
!      + a_71 x(5,i) + a_72 x(5,i)**2 + a_73 x(5,i)**3 + a_74 x(5,i)**4       5
!      + a_75 x(1,i)*x(5,i) + a_76 x(1,i)*x(5,i)**2 + a_77 x(1,i)*x(5,i)**3
!                           + a_78 x(1,i)**2*x(5,i) + a_79 x(1,i)**2*x(5,i)**2 
!                           + a_80 x(1,i)**3*x(5,i) 
!      + a_81 x(2,i)*x(5,i) + a_82 x(2,i)*x(5,i)**2 + a_83 x(2,i)*x(5,i)**3
!                           + a_84 x(2,i)**2*x(5,i) + a_85 x(2,i)**2*x(5,i)**2 
!                           + a_86 x(2,i)**3*x(5,i) 
!      + a_87 x(3,i)*x(5,i) + a_88 x(3,i)*x(5,i)**2 + a_89 x(3,i)*x(5,i)**3
!                           + a_90 x(3,i)**2*x(5,i) + a_91 x(3,i)**2*x(5,i)**2 
!                           + a_92 x(3,i)**3*x(5,i) 
!      + a_93 x(4,i)*x(5,i) + a_94 x(4,i)*x(5,i)**2 + a_95 x(4,i)*x(5,i)**3
!                           + a_96 x(4,i)**2*x(5,i) + a_97 x(4,i)**2*x(5,i)**2 
!                           + a_98 x(4,i)**3*x(5,i) 
!      + a_99  x(1,i) * x(2,i) * x(5,i) + a_100 x(1,i) * x(2,i)**2 * x(5,i)
!      + a_101 x(1,i) * x(2,i) * x(5,i)**2 + a_102 x(1,i)**2 * x(2,i) * x(5,i)
!      + a_103 x(1,i) * x(3,i) * x(5,i) + a_104 x(1,i) * x(3,i)**2 * x(5,i)
!      + a_105 x(1,i) * x(3,i) * x(5,i)**2 + a_106 x(1,i)**2 * x(3,i) * x(5,i)
!      + a_107 x(1,i) * x(4,i) * x(5,i) + a_108 x(1,i) * x(4,i)**2 * x(5,i)
!      + a_109 x(1,i) * x(4,i) * x(5,i)**2 + a_110 x(1,i)**2 * x(4,i) * x(5,i)
!      + a_111 x(2,i) * x(3,i) * x(5,i) + a_112 x(2,i) * x(3,i)**2 * x(5,i)
!      + a_113 x(2,i) * x(3,i) * x(5,i)**2 + a_114 x(2,i)**2 * x(3,i) * x(5,i)
!      + a_115 x(2,i) * x(4,i) * x(5,i) + a_116 x(2,i) * x(4,i)**2 * x(5,i)
!      + a_117 x(2,i) * x(4,i) * x(5,i)**2 + a_118 x(2,i)**2 * x(4,i) * x(5,i)
!      + a_119 x(3,i) * x(4,i) * x(5,i) + a_120 x(3,i) * x(4,i)**2 * x(5,i)
!      + a_121 x(3,i) * x(4,i) * x(5,i)**2 + a_122 x(3,i)**2 * x(4,i) * x(5,i)
!      + a_123 x(1,i) * x(2,i) * x(3,i) * x(5,i)
!      + a_124 x(1,i) * x(2,i) * x(4,i) * x(5,i)
!      + a_125 x(1,i) * x(3,i) * x(4,i) * x(5,i)
!      + a_126 x(2,i) * x(3,i) * x(4,i) * x(5,i)
!      + a_127 x(6,i) + a_128 x(6,i)**2 + a_129 x(6,i)**3 + a_130 x(6,i)**4   6
!      + a_131 x(1,i)*x(6,i) + a_132 x(1,i)*x(6,i)**2 + a_133 x(1,i)*x(6,i)**3
!                         + a_134 x(1,i)**2*x(6,i) + a_135 x(1,i)**2*x(6,i)**2 
!                         + a_136 x(1,i)**3*x(6,i) 
!      + a_137 x(2,i)*x(6,i) + a_138 x(2,i)*x(6,i)**2 + a_139 x(2,i)*x(6,i)**3
!                         + a_140 x(2,i)**2*x(6,i) + a_141 x(2,i)**2*x(6,i)**2 
!                         + a_142 x(2,i)**3*x(6,i) 
!      + a_143 x(3,i)*x(6,i) + a_144 x(3,i)*x(6,i)**2 + a_145 x(3,i)*x(6,i)**3
!                         + a_146 x(3,i)**2*x(6,i) + a_147 x(3,i)**2*x(6,i)**2 
!                         + a_148 x(3,i)**3*x(6,i) 
!      + a_149 x(4,i)*x(6,i) + a_150 x(4,i)*x(6,i)**2 + a_151 x(4,i)*x(6,i)**3
!                         + a_152 x(4,i)**2*x(6,i) + a_153 x(4,i)**2*x(6,i)**2 
!                         + a_154 x(4,i)**3*x(6,i) 
!      + a_155 x(5,i)*x(6,i) + a_156 x(5,i)*x(6,i)**2 + a_157 x(5,i)*x(6,i)**3
!                         + a_158 x(5,i)**2*x(6,i) + a_159 x(5,i)**2*x(6,i)**2 
!                         + a_160 x(5,i)**3*x(6,i) 
!      + a_161 x(1,i) * x(2,i) * x(6,i) + a_162 x(1,i) * x(2,i)**2 * x(6,i)
!      + a_163 x(1,i) * x(2,i) * x(6,i)**2 + a_164 x(1,i)**2 * x(2,i) * x(6,i)
!      + a_165 x(1,i) * x(3,i) * x(6,i) + a_166 x(1,i) * x(3,i)**2 * x(6,i)
!      + a_167 x(1,i) * x(3,i) * x(6,i)**2 + a_168 x(1,i)**2 * x(3,i) * x(6,i)
!      + a_169 x(1,i) * x(4,i) * x(6,i) + a_170 x(1,i) * x(4,i)**2 * x(6,i)
!      + a_171 x(1,i) * x(4,i) * x(6,i)**2 + a_172 x(1,i)**2 * x(4,i) * x(6,i)
!      + a_173 x(1,i) * x(5,i) * x(6,i) + a_174 x(1,i) * x(5,i)**2 * x(6,i)
!      + a_175 x(1,i) * x(5,i) * x(6,i)**2 + a_176 x(1,i)**2 * x(5,i) * x(6,i)
!      + a_177 x(2,i) * x(3,i) * x(6,i) + a_178 x(2,i) * x(3,i)**2 * x(6,i)
!      + a_179 x(2,i) * x(3,i) * x(6,i)**2 + a_180 x(2,i)**2 * x(3,i) * x(6,i)
!      + a_181 x(2,i) * x(4,i) * x(6,i) + a_182 x(2,i) * x(4,i)**2 * x(6,i)
!      + a_183 x(2,i) * x(4,i) * x(6,i)**2 + a_184 x(2,i)**2 * x(4,i) * x(6,i)
!      + a_185 x(2,i) * x(5,i) * x(6,i) + a_186 x(2,i) * x(5,i)**2 * x(6,i)
!      + a_187 x(2,i) * x(5,i) * x(6,i)**2 + a_188 x(2,i)**2 * x(5,i) * x(6,i)
!      + a_189 x(3,i) * x(4,i) * x(6,i) + a_190 x(3,i) * x(4,i)**2 * x(6,i)
!      + a_191 x(3,i) * x(4,i) * x(6,i)**2 + a_192 x(3,i)**2 * x(4,i) * x(6,i)
!      + a_193 x(3,i) * x(5,i) * x(6,i) + a_194 x(3,i) * x(5,i)**2 * x(6,i)
!      + a_195 x(3,i) * x(5,i) * x(6,i)**2 + a_196 x(3,i)**2 * x(5,i) * x(6,i)
!      + a_197 x(4,i) * x(5,i) * x(6,i) + a_198 x(4,i) * x(5,i)**2 * x(6,i)
!      + a_199 x(4,i) * x(5,i) * x(6,i)**2 + a_200 x(4,i)**2 * x(5,i) * x(6,i)
!      + a_201 x(1,i) * x(2,i) * x(3,i) * x(6,i)
!      + a_202 x(1,i) * x(2,i) * x(4,i) * x(6,i)
!      + a_203 x(1,i) * x(2,i) * x(5,i) * x(6,i)
!      + a_204 x(1,i) * x(3,i) * x(4,i) * x(6,i)
!      + a_205 x(1,i) * x(3,i) * x(5,i) * x(6,i)
!      + a_206 x(1,i) * x(4,i) * x(5,i) * x(6,i)
!      + a_207 x(2,i) * x(3,i) * x(4,i) * x(6,i)
!      + a_208 x(2,i) * x(3,i) * x(5,i) * x(6,i)
!      + a_209 x(2,i) * x(4,i) * x(5,i) * x(6,i)
!      + a_210 x(3,i) * x(4,i) * x(5,i) * x(6,i)
!
USE linear_solvers,     ONLY : linsolvx, linsolvms, linsolvsvd
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ncoeff, ndata
INTEGER, INTENT(INOUT) :: lsolve
REAL(DP), INTENT(IN) :: x(nvar,ndata), f(ndata)
REAL(DP), INTENT(INOUT) :: coeff(ncoeff)

REAL(DP) :: amat(ndata,ncoeff), aa(ncoeff,ncoeff), b(ncoeff) 

INTEGER :: ivar, jvar, idata

IF (nvar>6.OR.nvar<1) &
   CALL errore('multi_quartic','nvar must be from 1 to 6',1)
IF (ndata < 3) &
   CALL errore('multi_quartic','Too few sampling data',1)
IF (ndata < ncoeff) &
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
   amat(idata,5) = x(1,idata)*x(1,idata)*x(1,idata)*x(1,idata)

   IF (nvar>1) THEN
      amat(idata,6)  = x(2,idata)
      amat(idata,7)  = x(2,idata) * x(2,idata)
      amat(idata,8)  = x(2,idata) * x(2,idata) * x(2,idata)
      amat(idata,9)  = x(2,idata) * x(2,idata) * x(2,idata) * x(2,idata)
      amat(idata,10) = x(1,idata) * x(2,idata)
      amat(idata,11) = x(1,idata) * x(2,idata) * x(2,idata)
      amat(idata,12) = x(1,idata) * x(2,idata) * x(2,idata) * x(2,idata)
      amat(idata,13) = x(1,idata) * x(1,idata) * x(2,idata)
      amat(idata,14) = x(1,idata) * x(1,idata) * x(2,idata) * x(2,idata)
      amat(idata,15) = x(1,idata) * x(1,idata) * x(1,idata) * x(2,idata)
   ENDIF

   IF (nvar>2) THEN
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

   IF (nvar>3) THEN
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

   IF (nvar>4) THEN
      amat(idata,71) = x(5,idata)
      amat(idata,72) = x(5,idata)**2
      amat(idata,73) = x(5,idata)**3
      amat(idata,74) = x(5,idata)**4
      amat(idata,75) = x(1,idata) * x(5,idata)
      amat(idata,76) = x(1,idata) * x(5,idata)**2
      amat(idata,77) = x(1,idata) * x(5,idata)**3
      amat(idata,78) = x(1,idata) **2 * x(5,idata)
      amat(idata,79) = x(1,idata) **2 * x(5,idata)**2
      amat(idata,80) = x(1,idata) **3 * x(5,idata)
      amat(idata,81) = x(2,idata) * x(5,idata)
      amat(idata,82) = x(2,idata) * x(5,idata)**2
      amat(idata,83) = x(2,idata) * x(5,idata)**3
      amat(idata,84) = x(2,idata) **2 * x(5,idata)
      amat(idata,85) = x(2,idata) **2 * x(5,idata)**2
      amat(idata,86) = x(2,idata) **3 * x(5,idata)
      amat(idata,87) = x(3,idata) * x(5,idata)
      amat(idata,88) = x(3,idata) * x(5,idata)**2
      amat(idata,89) = x(3,idata) * x(5,idata)**3
      amat(idata,90) = x(3,idata) **2 * x(5,idata)
      amat(idata,91) = x(3,idata) **2 * x(5,idata)**2
      amat(idata,92) = x(3,idata) **3 * x(5,idata)
      amat(idata,93) = x(4,idata) * x(5,idata)
      amat(idata,94) = x(4,idata) * x(5,idata)**2
      amat(idata,95) = x(4,idata) * x(5,idata)**3
      amat(idata,96) = x(4,idata) **2 * x(5,idata)
      amat(idata,97) = x(4,idata) **2 * x(5,idata)**2
      amat(idata,98) = x(4,idata) **3 * x(5,idata)
      amat(idata,99) = x(1,idata) * x(2,idata) * x(5,idata)
      amat(idata,100) = x(1,idata) * x(2,idata)**2 * x(5,idata)
      amat(idata,101) = x(1,idata) * x(2,idata) * x(5,idata)**2
      amat(idata,102) = x(1,idata)**2 * x(2,idata) * x(5,idata)
      amat(idata,103) = x(1,idata) * x(3,idata) * x(5,idata)
      amat(idata,104) = x(1,idata) * x(3,idata)**2 * x(5,idata)
      amat(idata,105) = x(1,idata) * x(3,idata) * x(5,idata)**2
      amat(idata,106) = x(1,idata)**2 * x(3,idata) * x(5,idata)
      amat(idata,107) = x(1,idata) * x(4,idata) * x(5,idata)
      amat(idata,108) = x(1,idata) * x(4,idata)**2 * x(5,idata)
      amat(idata,109) = x(1,idata) * x(4,idata) * x(5,idata)**2
      amat(idata,110) = x(1,idata)**2 * x(4,idata) * x(5,idata)
      amat(idata,111) = x(2,idata) * x(3,idata) * x(5,idata)
      amat(idata,112) = x(2,idata) * x(3,idata)**2 * x(5,idata)
      amat(idata,113) = x(2,idata) * x(3,idata) * x(5,idata)**2
      amat(idata,114) = x(2,idata)**2 * x(3,idata) * x(5,idata)
      amat(idata,115) = x(2,idata) * x(4,idata) * x(5,idata)
      amat(idata,116) = x(2,idata) * x(4,idata)**2 * x(5,idata)
      amat(idata,117) = x(2,idata) * x(4,idata) * x(5,idata)**2
      amat(idata,118) = x(2,idata)**2 * x(4,idata) * x(5,idata)
      amat(idata,119) = x(3,idata) * x(4,idata) * x(5,idata)
      amat(idata,120) = x(3,idata) * x(4,idata)**2 * x(5,idata)
      amat(idata,121) = x(3,idata) * x(4,idata) * x(5,idata)**2
      amat(idata,122) = x(3,idata)**2 * x(4,idata) * x(5,idata)
      amat(idata,123) = x(1,idata) * x(2,idata) * x(3,idata) * x(5,idata)
      amat(idata,124) = x(1,idata) * x(2,idata) * x(4,idata) * x(5,idata)
      amat(idata,125) = x(1,idata) * x(3,idata) * x(4,idata) * x(5,idata)
      amat(idata,126) = x(2,idata) * x(3,idata) * x(4,idata) * x(5,idata)
   ENDIF

   IF (nvar>5) THEN
      amat(idata,127) = x(6,idata)
      amat(idata,128) = x(6,idata)**2
      amat(idata,129) = x(6,idata)**3
      amat(idata,130) = x(6,idata)**4
      amat(idata,131) = x(1,idata) * x(6,idata)
      amat(idata,132) = x(1,idata) * x(6,idata)**2
      amat(idata,133) = x(1,idata) * x(6,idata)**3
      amat(idata,134) = x(1,idata) **2 * x(6,idata)
      amat(idata,135) = x(1,idata) **2 * x(6,idata)**2
      amat(idata,136) = x(1,idata) **3 * x(6,idata)
      amat(idata,137) = x(2,idata) * x(6,idata)
      amat(idata,138) = x(2,idata) * x(6,idata)**2
      amat(idata,139) = x(2,idata) * x(6,idata)**3
      amat(idata,140) = x(2,idata) **2 * x(6,idata)
      amat(idata,141) = x(2,idata) **2 * x(6,idata)**2
      amat(idata,142) = x(2,idata) **3 * x(6,idata)
      amat(idata,143) = x(3,idata) * x(6,idata)
      amat(idata,144) = x(3,idata) * x(6,idata)**2
      amat(idata,145) = x(3,idata) * x(6,idata)**3
      amat(idata,146) = x(3,idata) **2 * x(6,idata)
      amat(idata,147) = x(3,idata) **2 * x(6,idata)**2
      amat(idata,148) = x(3,idata) **3 * x(6,idata)
      amat(idata,149) = x(4,idata) * x(6,idata)
      amat(idata,150) = x(4,idata) * x(6,idata)**2
      amat(idata,151) = x(4,idata) * x(6,idata)**3
      amat(idata,152) = x(4,idata) **2 * x(6,idata)
      amat(idata,153) = x(4,idata) **2 * x(6,idata)**2
      amat(idata,154) = x(4,idata) **3 * x(6,idata)
      amat(idata,155) = x(5,idata) * x(6,idata)
      amat(idata,156) = x(5,idata) * x(6,idata)**2
      amat(idata,157) = x(5,idata) * x(6,idata)**3
      amat(idata,158) = x(5,idata) **2 * x(6,idata)
      amat(idata,159) = x(5,idata) **2 * x(6,idata)**2
      amat(idata,160) = x(5,idata) **3 * x(6,idata)
      amat(idata,161) = x(1,idata) * x(2,idata) * x(6,idata)
      amat(idata,162) = x(1,idata) * x(2,idata)**2 * x(6,idata)
      amat(idata,163) = x(1,idata) * x(2,idata) * x(6,idata)**2
      amat(idata,164) = x(1,idata)**2 * x(2,idata) * x(6,idata)
      amat(idata,165) = x(1,idata) * x(3,idata) * x(6,idata)
      amat(idata,166) = x(1,idata) * x(3,idata)**2 * x(6,idata)
      amat(idata,167) = x(1,idata) * x(3,idata) * x(6,idata)**2
      amat(idata,168) = x(1,idata)**2 * x(3,idata) * x(6,idata)
      amat(idata,169) = x(1,idata) * x(4,idata) * x(6,idata)
      amat(idata,170) = x(1,idata) * x(4,idata)**2 * x(6,idata)
      amat(idata,171) = x(1,idata) * x(4,idata) * x(6,idata)**2
      amat(idata,172) = x(1,idata)**2 * x(4,idata) * x(6,idata)
      amat(idata,173) = x(1,idata) * x(5,idata) * x(6,idata)
      amat(idata,174) = x(1,idata) * x(5,idata)**2 * x(6,idata)
      amat(idata,175) = x(1,idata) * x(5,idata) * x(6,idata)**2
      amat(idata,176) = x(1,idata)**2 * x(5,idata) * x(6,idata)
      amat(idata,177) = x(2,idata) * x(3,idata) * x(6,idata)
      amat(idata,178) = x(2,idata) * x(3,idata)**2 * x(6,idata)
      amat(idata,179) = x(2,idata) * x(3,idata) * x(6,idata)**2
      amat(idata,180) = x(2,idata)**2 * x(3,idata) * x(6,idata)
      amat(idata,181) = x(2,idata) * x(4,idata) * x(6,idata)
      amat(idata,182) = x(2,idata) * x(4,idata)**2 * x(6,idata)
      amat(idata,183) = x(2,idata) * x(4,idata) * x(6,idata)**2
      amat(idata,184) = x(2,idata)**2 * x(4,idata) * x(6,idata)
      amat(idata,185) = x(2,idata) * x(5,idata) * x(6,idata)
      amat(idata,186) = x(2,idata) * x(5,idata)**2 * x(6,idata)
      amat(idata,187) = x(2,idata) * x(5,idata) * x(6,idata)**2
      amat(idata,188) = x(2,idata)**2 * x(5,idata) * x(6,idata)
      amat(idata,189) = x(3,idata) * x(4,idata) * x(6,idata)
      amat(idata,190) = x(3,idata) * x(4,idata)**2 * x(6,idata)
      amat(idata,191) = x(3,idata) * x(4,idata) * x(6,idata)**2
      amat(idata,192) = x(3,idata)**2 * x(4,idata) * x(6,idata)
      amat(idata,193) = x(3,idata) * x(5,idata) * x(6,idata)
      amat(idata,194) = x(3,idata) * x(5,idata)**2 * x(6,idata)
      amat(idata,195) = x(3,idata) * x(5,idata) * x(6,idata)**2
      amat(idata,196) = x(3,idata)**2 * x(5,idata) * x(6,idata)
      amat(idata,197) = x(4,idata) * x(5,idata) * x(6,idata)
      amat(idata,198) = x(4,idata) * x(5,idata)**2 * x(6,idata)
      amat(idata,199) = x(4,idata) * x(5,idata) * x(6,idata)**2
      amat(idata,200) = x(4,idata)**2 * x(5,idata) * x(6,idata)
      amat(idata,201) = x(1,idata) * x(2,idata) * x(3,idata) * x(6,idata)
      amat(idata,202) = x(1,idata) * x(2,idata) * x(4,idata) * x(6,idata)
      amat(idata,203) = x(1,idata) * x(2,idata) * x(5,idata) * x(6,idata)
      amat(idata,204) = x(1,idata) * x(3,idata) * x(4,idata) * x(6,idata)
      amat(idata,205) = x(1,idata) * x(3,idata) * x(5,idata) * x(6,idata)
      amat(idata,206) = x(1,idata) * x(4,idata) * x(5,idata) * x(6,idata)
      amat(idata,207) = x(2,idata) * x(3,idata) * x(4,idata) * x(6,idata)
      amat(idata,208) = x(2,idata) * x(3,idata) * x(5,idata) * x(6,idata)
      amat(idata,209) = x(2,idata) * x(4,idata) * x(5,idata) * x(6,idata)
      amat(idata,210) = x(3,idata) * x(4,idata) * x(5,idata) * x(6,idata)
   ENDIF

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
   WRITE(stdout,'(5x,"Finding the quartic polynomial using &
                                              &ncoeff x ncoeff matrix")')  
   CALL linsolvx(aa,ncoeff,b,coeff)
ELSEIF(lsolve==2) THEN
   WRITE(stdout,'(5x,"Finding the quartic polynomial using &
                                                   &QR factorization")')  
   CALL linsolvms(amat,ndata,ncoeff,f,coeff)
ELSEIF(lsolve==3) THEN
   WRITE(stdout,'(5x,"Finding the quartic polynomial using SVD")')  
   CALL linsolvsvd(amat,ndata,ncoeff,f,coeff)
ENDIF

RETURN
END SUBROUTINE fit_multi_quartic

SUBROUTINE evaluate_fit_quartic(nvar,ncoeff,x,f,coeff)
!
!  This routine evaluates the quartic polynomial at the point x
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ncoeff
REAL(DP), INTENT(IN) :: x(nvar)
REAL(DP), INTENT(IN) :: coeff(ncoeff)
REAL(DP), INTENT(INOUT) :: f

REAL(DP) :: aux
!
!  one variable
!
aux = coeff(1) + x(1)*(coeff(2)+x(1)*(coeff(3)+x(1)*(coeff(4)+coeff(5)*x(1))))
!
!  two variables
!
IF (nvar>1) THEN
   aux=aux+x(2)*(coeff(6) + x(2)*( coeff(7) + coeff(11) * x(1) +     &
                                             coeff(14) * x(1)**2 +  &
              x(2)*(coeff(8) + coeff(9) * x(2) + coeff(12) * x(1)))  &
             + x(1) *( coeff(10)    &
             + x(1) * (coeff(13) + coeff(15) * x(1) ) ) )
ENDIF
!
!  three variables
!
IF (nvar>2) THEN
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
!
!  four variables
!
IF (nvar>3) THEN
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
!
!  five variables
!
IF (nvar>4) THEN
   aux = aux + coeff(71) * x(5) + coeff(72) * x(5)**2 + coeff(73) * x(5)**3 + &
                                                        coeff(74) * x(5)**4   &
             + coeff(75) * x(1) * x(5) + coeff(76) * x(1) * x(5)**2           &
             + coeff(77) * x(1) * x(5)**3 + coeff(78) * x(1)**2 * x(5)        &
             + coeff(79) * x(1)**2 * x(5)**2 + coeff(80) * x(1)**3 * x(5)     &
             + coeff(81) * x(2) * x(5) + coeff(82) * x(2) * x(5)**2           &
             + coeff(83) * x(2) * x(5)**3 + coeff(84) * x(2)**2 * x(5)        &
             + coeff(85) * x(2)**2 * x(5)**2 + coeff(86) * x(2)**3 * x(5)     &
             + coeff(87) * x(3) * x(5) + coeff(88) * x(3) * x(5)**2           &
             + coeff(89) * x(3) * x(5)**3 + coeff(90) * x(3)**2 * x(5)        &
             + coeff(91) * x(3)**2 * x(5)**2 + coeff(92) * x(3)**3 * x(5)     &
             + coeff(93) * x(4) * x(5) + coeff(94) * x(4) * x(5)**2           &
             + coeff(95) * x(4) * x(5)**3 + coeff(96) * x(4)**2 * x(5)        &
             + coeff(97) * x(4)**2 * x(5)**2 + coeff(98) * x(4)**3 * x(5)     &
             + coeff(99) * x(1) * x(2) * x(5)      &
             + coeff(100) * x(1) * x(2)**2 * x(5)  &
             + coeff(101) * x(1) * x(2) * x(5)**2  &
             + coeff(102) * x(1)**2 * x(2) * x(5)  &
             + coeff(103) * x(1) * x(3) * x(5)     &
             + coeff(104) * x(1) * x(3)**2 * x(5)  &
             + coeff(105) * x(1) * x(3) * x(5)**2  &
             + coeff(106) * x(1)**2 * x(3) * x(5)  &
             + coeff(107) * x(1) * x(4) * x(5)     &
             + coeff(108) * x(1) * x(4)**2 * x(5)  &
             + coeff(109) * x(1) * x(4) * x(5)**2  &
             + coeff(110) * x(1)**2 * x(4) * x(5)  &
             + coeff(111) * x(2) * x(3) * x(5)     &
             + coeff(112) * x(2) * x(3)**2 * x(5)  &
             + coeff(113) * x(2) * x(3) * x(5)**2  &
             + coeff(114) * x(2)**2 * x(3) * x(5)  & 
             + coeff(115) * x(2) * x(4) * x(5)     &
             + coeff(116) * x(2) * x(4)**2 * x(5)  &
             + coeff(117) * x(2) * x(4) * x(5)**2  &
             + coeff(118) * x(2)**2 * x(4) * x(5)  & 
             + coeff(119) * x(3) * x(4) * x(5)     &
             + coeff(120) * x(3) * x(4)**2 * x(5)  &
             + coeff(121) * x(3) * x(4) * x(5)**2  &
             + coeff(122) * x(3)**2 * x(4) * x(5)  & 
             + coeff(123) * x(1) * x(2) * x(3) * x(5) & 
             + coeff(124) * x(1) * x(2) * x(4) * x(5) & 
             + coeff(125) * x(1) * x(3) * x(4) * x(5) & 
             + coeff(126) * x(2) * x(3) * x(4) * x(5) 
ENDIF
!
!  six variables
!
IF (nvar>5) THEN
   aux = aux + coeff(127) * x(6) + coeff(128) * x(6)**2 +                     &
                      coeff(129) * x(6)**3 + coeff(130) * x(6)**4             &
             + coeff(131) * x(1) * x(6) + coeff(132) * x(1) * x(6)**2         &
             + coeff(133) * x(1) * x(6)**3 + coeff(134) * x(1)**2 * x(6)      &
             + coeff(135) * x(1)**2 * x(6)**2 + coeff(136) * x(1)**3 * x(6)   &
             + coeff(137) * x(2) * x(6) + coeff(138) * x(2) * x(6)**2         &
             + coeff(139) * x(2) * x(6)**3 + coeff(140) * x(2)**2 * x(6)      &
             + coeff(141) * x(2)**2 * x(6)**2 + coeff(142) * x(2)**3 * x(6)   &
             + coeff(143) * x(3) * x(6) + coeff(144) * x(3) * x(6)**2         &
             + coeff(145) * x(3) * x(6)**3 + coeff(146) * x(3)**2 * x(6)      &
             + coeff(147) * x(3)**2 * x(6)**2 + coeff(148) * x(3)**3 * x(6)   &
             + coeff(149) * x(4) * x(6) + coeff(150) * x(4) * x(6)**2         &
             + coeff(151) * x(4) * x(6)**3 + coeff(152) * x(4)**2 * x(6)      &
             + coeff(153) * x(4)**2 * x(6)**2 + coeff(154) * x(4)**3 * x(6)   &
             + coeff(155) * x(5) * x(6) + coeff(156) * x(5) * x(6)**2         &
             + coeff(157) * x(5) * x(6)**3 + coeff(158) * x(5)**2 * x(6)      &
             + coeff(159) * x(5)**2 * x(6)**2 + coeff(160) * x(5)**3 * x(6)   &
             + coeff(161) * x(1) * x(2) * x(6)        &
             + coeff(162) * x(1) * x(2)**2 * x(6)     &
             + coeff(163) * x(1) * x(2) * x(6)**2     &
             + coeff(164) * x(1)**2 * x(2) * x(6)     &
             + coeff(165) * x(1) * x(3) * x(6)        &
             + coeff(166) * x(1) * x(3)**2 * x(6)     &
             + coeff(167) * x(1) * x(3) * x(6)**2     &
             + coeff(168) * x(1)**2 * x(3) * x(6)     &
             + coeff(169) * x(1) * x(4) * x(6)        &
             + coeff(170) * x(1) * x(4)**2 * x(6)     &
             + coeff(171) * x(1) * x(4) * x(6)**2     &
             + coeff(172) * x(1)**2 * x(4) * x(6)     &
             + coeff(173) * x(1) * x(5) * x(6)        &
             + coeff(174) * x(1) * x(5)**2 * x(6)     & 
             + coeff(175) * x(1) * x(5) * x(6)**2     &
             + coeff(176) * x(1)**2 * x(5) * x(6)     &
             + coeff(177) * x(2) * x(3) * x(6)        &
             + coeff(178) * x(2) * x(3)**2 * x(6)     &
             + coeff(179) * x(2) * x(3) * x(6)**2     &
             + coeff(180) * x(2)**2 * x(3) * x(6)     & 
             + coeff(181) * x(2) * x(4) * x(6)        &
             + coeff(182) * x(2) * x(4)**2 * x(6)     &
             + coeff(183) * x(2) * x(4) * x(6)**2     &
             + coeff(184) * x(2)**2 * x(4) * x(6)     & 
             + coeff(185) * x(2) * x(5) * x(6)        &
             + coeff(186) * x(2) * x(5)**2 * x(6)     &
             + coeff(187) * x(2) * x(5) * x(6)**2     &
             + coeff(188) * x(2)**2 * x(5) * x(6)     & 
             + coeff(189) * x(3) * x(4) * x(6)        &
             + coeff(190) * x(3) * x(4)**2 * x(6)     &
             + coeff(191) * x(3) * x(4) * x(6)**2     &
             + coeff(192) * x(3)**2 * x(4) * x(6)     & 
             + coeff(193) * x(3) * x(5) * x(6)        &
             + coeff(194) * x(3) * x(5)**2 * x(6)     &
             + coeff(195) * x(3) * x(5) * x(6)**2     &
             + coeff(196) * x(3)**2 * x(5) * x(6)     & 
             + coeff(197) * x(4) * x(5) * x(6)        &
             + coeff(198) * x(4) * x(5)**2 * x(6)     &
             + coeff(199) * x(4) * x(5) * x(6)**2     &
             + coeff(200) * x(4)**2 * x(5) * x(6)     & 
             + coeff(201) * x(1) * x(2) * x(3) * x(6) & 
             + coeff(202) * x(1) * x(2) * x(4) * x(6) & 
             + coeff(203) * x(1) * x(2) * x(5) * x(6) & 
             + coeff(204) * x(1) * x(3) * x(4) * x(6) & 
             + coeff(205) * x(1) * x(3) * x(5) * x(6) & 
             + coeff(206) * x(1) * x(4) * x(5) * x(6) & 
             + coeff(207) * x(2) * x(3) * x(4) * x(6) &
             + coeff(208) * x(2) * x(3) * x(5) * x(6) &
             + coeff(209) * x(2) * x(4) * x(5) * x(6) &
             + coeff(210) * x(3) * x(4) * x(5) * x(6) 
ENDIF

f=aux

RETURN
END SUBROUTINE evaluate_fit_quartic

SUBROUTINE evaluate_quartic_grad(nvar,ncoeff,x,f,coeff)
!
!   This routine evaluates the gradient of the quartic polynomial
!   at the point x.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ncoeff
REAL(DP), INTENT(IN) :: x(nvar)
REAL(DP), INTENT(IN) :: coeff(ncoeff)
REAL(DP), INTENT(INOUT) :: f(nvar)

REAL(DP) :: aux(nvar)

IF (nvar>6) CALL errore('evaluate_quartic_grad','gradient not availble',1)

aux(1) = coeff(2) + 2.0_DP * coeff(3) * x(1) + 3.0_DP * coeff(4)*x(1)**2 + &
                    4.0_DP * coeff(5) * x(1)**3

IF (nvar>1) THEN
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

IF (nvar>2) THEN
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

IF (nvar>3) THEN
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
   aux(4) = coeff(36) + 2.0_DP * coeff(37) * x(4) +                         &
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

IF (nvar>4) THEN
   aux(1) = aux(1) + coeff(75) * x(5) + coeff(76) * x(5)**2 +               &
                     coeff(77) * x(5)**3 + 2.0_DP * coeff(78) * x(1)*x(5) + &
                     2.0_DP * coeff(79) * x(1)*x(5)**2 +                    & 
                     3.0_DP * coeff(80) * x(1)**2 * x(5) +                  &
                     coeff(99) * x(2) * x(5) + coeff(100) * x(2)**2 * x(5) + &
                     coeff(101) * x(2) * x(5)**2 +                           &
                     2.0_DP * coeff(102) * x(1) * x(2) * x(5) +              &
                     coeff(103) * x(3) * x(5) + coeff(104) * x(3)**2 * x(5) + &
                     coeff(105) * x(3) * x(5)**2 +                           &
                     2.0_DP * coeff(106) * x(1) * x(3) * x(5) +              &
                     coeff(107) * x(4) * x(5) + coeff(108) * x(4)**2 * x(5) + &
                     coeff(109) * x(4) * x(5)**2 +                           &
                     2.0_DP * coeff(110) * x(1) * x(4) * x(5) +              &
                     coeff(123) * x(2) * x(3) * x(5)  +                      &
                     coeff(124) * x(2) * x(4) * x(5)  +                      &
                     coeff(125) * x(3) * x(4) * x(5) 
   aux(2) = aux(2) + coeff(81) * x(5) + coeff(82) * x(5)**2 +                &
                     coeff(83) * x(5)**3 + 2.0_DP * coeff(84) * x(2)*x(5) +  &
                     2.0_DP * coeff(85) * x(2)*x(5)**2 +                     & 
                     3.0_DP * coeff(86) * x(2)**2 * x(5) +                   &
                     coeff(99) * x(1) * x(5) +                               &
                     2.0_DP *coeff(100) * x(1)*x(2)*x(5) +                   &
                     coeff(101) * x(1) * x(5)**2 +                           &
                     coeff(102) * x(1) **2 * x(5) +                          &
                     coeff(111) * x(3) * x(5) + coeff(112) * x(3)**2 * x(5) + &
                     coeff(113) * x(3) * x(5)**2 +                           &
                     2.0_DP * coeff(114) * x(2) * x(3) * x(5) +              &
                     coeff(115) * x(4) * x(5) + coeff(116) * x(4)**2 * x(5) + &
                     coeff(117) * x(4) * x(5)**2 +                           &
                     2.0_DP * coeff(118) * x(2) * x(4) * x(5) +              &
                     coeff(123) * x(1) * x(3) * x(5)  +                      &
                     coeff(124) * x(1) * x(4) * x(5)  +                      &
                     coeff(126) * x(3) * x(4) * x(5) 

   aux(3) = aux(3) + coeff(87) * x(5) + coeff(88) * x(5)**2 +               &
                    coeff(89) * x(5)**3 + 2.0_DP * coeff(90) * x(3)*x(5) + &
                    2.0_DP * coeff(91) * x(3)*x(5)**2 +                    & 
                    3.0_DP * coeff(92) * x(3)**2 * x(5) +                  &
                     coeff(103) * x(1) * x(5) +                            &
                     2.0_DP *coeff(104) * x(1)*x(3)*x(5) +                 &
                     coeff(105) * x(1) * x(5)**2 +                         &
                     coeff(106) * x(1) **2 * x(5) +                        &
                     coeff(111) * x(2) * x(5) +                            &
                     2.0_DP*coeff(112)*x(2)*x(3)*x(5) +                    &
                     coeff(113) * x(2) * x(5)**2 +                         &
                     coeff(114) * x(2) **2 * x(5) +                        &
                     coeff(119) * x(4) * x(5) +                            &
                     coeff(120) * x(4) **2 * x(5) +                        &
                     coeff(121) * x(4) * x(5)**2 +                         &
                     2.0_DP*coeff(122)*x(3)*x(4)*x(5) +                    &
                     coeff(123) * x(1) * x(2) * x(5) +                     &
                     coeff(125) * x(1) * x(4) * x(5) +                     &
                     coeff(126) * x(2) * x(4) * x(5) 

   aux(4) = aux(4) + coeff(93) * x(5) + coeff(94) * x(5)**2 +               &
                    coeff(95) * x(5)**3 + 2.0_DP * coeff(96) * x(4)*x(5) + &
                    2.0_DP * coeff(97) * x(4)*x(5)**2 +                    & 
                    3.0_DP * coeff(98) * x(4)**2 * x(5) +                  &
                     coeff(107) * x(1) * x(5) +                            &
                     2.0_DP *coeff(108) * x(1)*x(4)*x(5) +                 &
                     coeff(109) * x(1) * x(5)**2 +                         &
                     coeff(110) * x(1) **2 * x(5) +                        &
                     coeff(115) * x(2) * x(5) +                            &
                     2.0_DP*coeff(116)*x(2)*x(4)*x(5) +                    &
                     coeff(117) * x(2) * x(5)**2 +                         &
                     coeff(118) * x(2) **2 * x(5) +                        &
                     coeff(119) * x(3) * x(5) +                            &
                     2.0_DP*coeff(120)*x(3)*x(4)*x(5) +                    &
                     coeff(121) * x(3) * x(5)**2 +                         &
                     coeff(122) * x(3) **2 * x(5) +                        &
                     coeff(124) * x(1) * x(2) * x(5) +                     &
                     coeff(125) * x(1) * x(3) * x(5) +                     &
                     coeff(126) * x(2) * x(3) * x(5) 

   aux(5) = coeff(71) + 2.0_DP * coeff(72) * x(5) +                        &
                        3.0_DP * coeff(73) * x(5)**2 +                     &
                        4.0_DP * coeff(74) * x(5)**3 +                     &
                                 coeff(75) * x(1) +                        &
                        2.0_DP * coeff(76) * x(1) * x(5) +                 &
                        3.0_DP * coeff(77) * x(1) * x(5) **2 +             &
                                 coeff(78) * x(1) ** 2 +                   &
                        2.0_DP * coeff(79) * x(1) **2 * x(5) +             &
                                 coeff(80) * x(1) ** 3 +                   &
                                 coeff(81) * x(2) +                        &
                        2.0_DP * coeff(82) * x(2) * x(5) +                 &
                        3.0_DP * coeff(83) * x(2) * x(5) **2 +             &
                                 coeff(84) * x(2) ** 2 +                   &
                        2.0_DP * coeff(85) * x(2) **2 * x(5) +             &
                                 coeff(86) * x(2) ** 3 +                   &
                                 coeff(87) * x(3) +                        &
                        2.0_DP * coeff(88) * x(3) * x(5) +                 &
                        3.0_DP * coeff(89) * x(3) * x(5) **2 +             &
                                 coeff(90) * x(3) ** 2 +                   &
                        2.0_DP * coeff(91) * x(3) **2 * x(5) +             &
                                 coeff(92) * x(3) ** 3 +                   &
                                 coeff(93) * x(4) +                        &
                        2.0_DP * coeff(94) * x(4) * x(5) +                 &
                        3.0_DP * coeff(95) * x(4) * x(5) **2 +             &
                                 coeff(96) * x(4) ** 2 +                   &
                        2.0_DP * coeff(97) * x(4) **2 * x(5) +             &
                                 coeff(98) * x(4) ** 3 +                   &
                                 coeff(99) * x(1) * x(2) +                 &
                                 coeff(100) * x(1) * x(2) **2 +            &
                        2.0_DP * coeff(101) * x(1) * x(2) * x(5) +         &
                                 coeff(102) * x(1) **2 * x(2) +            &
                                 coeff(103) * x(1) * x(3) +                &
                                 coeff(104) * x(1) * x(3) **2 +            &
                        2.0_DP * coeff(105) * x(1) * x(3) * x(5) +         &
                                 coeff(106) * x(1) **2 * x(3) +            &
                                 coeff(107) * x(1) * x(4) +                &
                                 coeff(108) * x(1) * x(4) **2 +            &
                        2.0_DP * coeff(109) * x(1) * x(4) * x(5) +         &
                                 coeff(110) * x(1) **2 * x(4) +            &
                                 coeff(111) * x(2) * x(3) +                &
                                 coeff(112) * x(2) * x(3) **2 +            &
                        2.0_DP * coeff(113) * x(2) * x(3) * x(5) +         &
                                 coeff(114) * x(2) **2 * x(3) +            &
                                 coeff(115) * x(2) * x(4) +                &
                                 coeff(116) * x(2) * x(4) **2 +            &
                        2.0_DP * coeff(117) * x(2) * x(4) * x(5) +         &
                                 coeff(118) * x(2) **2 * x(4) +            &
                                 coeff(119) * x(3) * x(4) +                &
                                 coeff(120) * x(3) * x(4) **2 +            &
                        2.0_DP * coeff(121) * x(3) * x(4) * x(5) +         &
                                 coeff(122) * x(3) **2 * x(4) +            &
                                 coeff(123) * x(1) * x(2) * x(3) +         &
                                 coeff(124) * x(1) * x(2) * x(4) +         &
                                 coeff(125) * x(1) * x(3) * x(4) +         &
                                 coeff(126) * x(2) * x(3) * x(4) 
ENDIF

IF (nvar>5) THEN
   aux(1) = aux(1) + coeff(131) * x(6) + coeff(132) * x(6)**2 +              &
                     coeff(133) * x(6)**3 + 2.0_DP * coeff(134) * x(1)*x(6) + &
                     2.0_DP * coeff(135) * x(1)*x(6)**2 +                    & 
                     3.0_DP * coeff(136) * x(1)**2 * x(6) +                  &
                     coeff(161) * x(2) * x(6) + coeff(162) * x(2)**2 * x(6) + &
                     coeff(163) * x(2) * x(6)**2 +                           &
                     2.0_DP * coeff(164) * x(1) * x(2) * x(6) +              &
                     coeff(165) * x(3) * x(6) + coeff(166) * x(3)**2 * x(6) + &
                     coeff(167) * x(3) * x(6)**2 +                           &
                     2.0_DP * coeff(168) * x(1) * x(3) * x(6) +              &
                     coeff(169) * x(4) * x(6) + coeff(170) * x(4)**2 * x(6) + &
                     coeff(171) * x(4) * x(6)**2 +                           &
                     2.0_DP * coeff(172) * x(1) * x(4) * x(6) +              &
                     coeff(173) * x(5) * x(6) + coeff(174) * x(5)**2 * x(6) + &
                     coeff(175) * x(5) * x(6)**2 +                           &
                     2.0_DP * coeff(176) * x(1) * x(5) * x(6) +              &
                     coeff(201) * x(2) * x(3) * x(6)  +                      &
                     coeff(202) * x(2) * x(4) * x(6)  +                      &
                     coeff(203) * x(2) * x(5) * x(6)  +                      &
                     coeff(204) * x(3) * x(4) * x(6)  +                      &
                     coeff(205) * x(3) * x(5) * x(6)  +                      &
                     coeff(206) * x(4) * x(5) * x(6)  

   aux(2) = aux(2) + coeff(137) * x(6) + coeff(138) * x(6)**2 +              &
                     coeff(139) * x(6)**3 + 2.0_DP * coeff(140)*x(2)*x(6) +  &
                     2.0_DP * coeff(141) * x(2)*x(6)**2 +                    & 
                     3.0_DP * coeff(142) * x(2)**2 * x(6) +                  &
                     coeff(161) * x(1) * x(6) +                               &
                     2.0_DP *coeff(162) * x(1)*x(2)*x(6) +                   &
                     coeff(163) * x(1) * x(6)**2 +                           &
                     coeff(164) * x(1) **2 * x(6) +                          &
                     coeff(177) * x(3) * x(6) + coeff(178) * x(3)**2 * x(6) + &
                     coeff(179) * x(3) * x(6)**2 +                           &
                     2.0_DP * coeff(180) * x(2) * x(3) * x(6) +              &
                     coeff(181) * x(4) * x(6) + coeff(182) * x(4)**2 * x(6) + &
                     coeff(183) * x(4) * x(6)**2 +                           &
                     2.0_DP * coeff(184) * x(2) * x(4) * x(6) +              &
                     coeff(185) * x(5) * x(6) + coeff(186) * x(5)**2 * x(6) + &
                     coeff(187) * x(5) * x(6)**2 +                           &
                     2.0_DP * coeff(188) * x(2) * x(5) * x(6) +              &
                     coeff(201) * x(1) * x(3) * x(6) +                       &
                     coeff(202) * x(1) * x(4) * x(6) +                       &
                     coeff(203) * x(1) * x(5) * x(6) +                       &
                     coeff(207) * x(3) * x(4) * x(6) +                       &
                     coeff(208) * x(3) * x(5) * x(6) +                       &
                     coeff(209) * x(4) * x(5) * x(6)  

   aux(3) = aux(3) + coeff(143) * x(6) + coeff(144) * x(6)**2 +            &
                     coeff(145) * x(6)**3 + 2.0_DP * coeff(146) * x(3)*x(6) + &
                     2.0_DP * coeff(147) * x(3)*x(6)**2 +                  & 
                     3.0_DP * coeff(148) * x(3)**2 * x(6) +                &
                     coeff(165) * x(1) * x(6) +                            &
                     2.0_DP *coeff(166) * x(1)*x(3)*x(6) +                 &
                     coeff(167) * x(1) * x(6)**2 +                         &
                     coeff(168) * x(1) **2 * x(6) +                        &
                     coeff(177) * x(2) * x(6) +                            &
                     2.0_DP*coeff(178)*x(2)*x(3)*x(6) +                    &
                     coeff(179) * x(2) * x(6)**2 +                         &
                     coeff(180) * x(2) **2 * x(6) +                        &
                     coeff(189) * x(4) * x(6) +                            &
                     coeff(190) * x(4) **2 * x(6) +                        &
                     coeff(191) * x(4) * x(6)**2 +                         &
                     2.0_DP*coeff(192)*x(3)*x(4)*x(6) +                    &
                     coeff(193) * x(5) * x(6) +                            &
                     coeff(194) * x(5) **2 * x(6) +                        &
                     coeff(195) * x(5) * x(6)**2 +                         &
                     2.0_DP*coeff(196)*x(3)*x(5)*x(6) +                    &
                     coeff(201) * x(1) * x(2) * x(6) +                     &
                     coeff(204) * x(1) * x(4) * x(6) +                     &
                     coeff(205) * x(1) * x(5) * x(6) +                     &
                     coeff(207) * x(2) * x(4) * x(6) +                     &
                     coeff(208) * x(2) * x(5) * x(6) +                     &
                     coeff(210) * x(4) * x(5) * x(6)                      

   aux(4) = aux(4) + coeff(149) * x(6) + coeff(150) * x(6)**2 +            &
                     coeff(151) * x(6)**3 + 2.0_DP * coeff(152) * x(4)*x(6) + &
                     2.0_DP * coeff(153) * x(4)*x(6)**2 +                  & 
                     3.0_DP * coeff(154) * x(4)**2*x(6) +                  &
                     coeff(169) * x(1) * x(6) +                            &
                     2.0_DP *coeff(170) * x(1)*x(4)*x(6) +                 &
                     coeff(171) * x(1) * x(6)**2 +                         &
                     coeff(172) * x(1) **2 * x(6) +                        &
                     coeff(181) * x(2) * x(6) +                            &
                     2.0_DP*coeff(182)*x(2)*x(4)*x(6) +                    &
                     coeff(183) * x(2) * x(6)**2 +                         &
                     coeff(184) * x(2) **2 * x(6) +                        &
                     coeff(189) * x(3) * x(6) +                            &
                     2.0_DP*coeff(190)*x(3)*x(4)*x(6) +                    &
                     coeff(191) * x(3) * x(6)**2 +                         &
                     coeff(192) * x(3) **2 * x(6) +                        &
                     coeff(197) * x(5) * x(6) +                            &
                     coeff(198) * x(5) **2 * x(6) +                        &
                     coeff(199) * x(5) * x(6)**2 +                         &
                     2.0_DP*coeff(200)*x(4)*x(5)*x(6) +                    &
                     coeff(202) * x(1) * x(2) * x(6) +                     &
                     coeff(204) * x(1) * x(3) * x(6) +                     &
                     coeff(206) * x(1) * x(5) * x(6) +                     &
                     coeff(207) * x(2) * x(3) * x(6) +                     &
                     coeff(209) * x(2) * x(5) * x(6) +                     &
                     coeff(210) * x(3) * x(5) * x(6) 

   aux(5) = aux(5) + coeff(155) * x(6) + coeff(156) * x(6)**2 +            &
                     coeff(157) * x(6)**3 + 2.0_DP * coeff(158) * x(5)*x(6) + &
                     2.0_DP * coeff(159) * x(5)*x(6)**2 +                  & 
                     3.0_DP * coeff(160) * x(5)**2*x(6) +                  &
                     coeff(173) * x(1) * x(6) +                            &
                     2.0_DP *coeff(174) * x(1)*x(5)*x(6) +                 &
                     coeff(175) * x(1) * x(6)**2 +                         &
                     coeff(176) * x(1) **2 * x(6) +                        &
                     coeff(185) * x(2) * x(6) +                            &
                     2.0_DP*coeff(186)*x(2)*x(5)*x(6) +                    &
                     coeff(187) * x(2) * x(6)**2 +                         &
                     coeff(188) * x(2) **2 * x(6) +                        &
                     coeff(193) * x(3) * x(6) +                            &
                     2.0_DP*coeff(194)*x(3)*x(5)*x(6) +                    &
                     coeff(195) * x(3) * x(6)**2 +                         &
                     coeff(196) * x(3) **2 * x(6) +                        &
                     coeff(197) * x(4) * x(6) +                            &
                     2.0_DP*coeff(198)*x(4)*x(5)*x(6) +                    &
                     coeff(199) * x(4) * x(6)**2 +                         &
                     coeff(200) * x(4) **2 * x(6) +                        &
                     coeff(203) * x(1) * x(2) * x(6) +                     &
                     coeff(205) * x(1) * x(3) * x(6) +                     &
                     coeff(206) * x(1) * x(4) * x(6) +                     &
                     coeff(208) * x(2) * x(3) * x(6) +                     &
                     coeff(209) * x(2) * x(4) * x(6) +                     &
                     coeff(210) * x(3) * x(4) * x(6) 

   aux(6) = coeff(127) + 2.0_DP * coeff(128) * x(6) +                      &
                        3.0_DP * coeff(129) * x(6)**2 +                    &
                        4.0_DP * coeff(130) * x(6)**3 +                    &
                                 coeff(131) * x(1) +                       &
                        2.0_DP * coeff(132) * x(1) * x(6) +                &
                        3.0_DP * coeff(133) * x(1) * x(6) **2 +            &
                                 coeff(134) * x(1) ** 2 +                  &
                        2.0_DP * coeff(135) * x(1) **2 * x(6) +            &
                                 coeff(136) * x(1) ** 3 +                  &
                                 coeff(137) * x(2) +                       &
                        2.0_DP * coeff(138) * x(2) * x(6) +                &
                        3.0_DP * coeff(139) * x(2) * x(6) **2 +            &
                                 coeff(140) * x(2) ** 2 +                  &
                        2.0_DP * coeff(141) * x(2) **2 * x(6) +            &
                                 coeff(142) * x(2) ** 3 +                  &
                                 coeff(143) * x(3) +                       &
                        2.0_DP * coeff(144) * x(3) * x(6) +                &
                        3.0_DP * coeff(145) * x(3) * x(6) **2 +            &
                                 coeff(146) * x(3) ** 2 +                  &
                        2.0_DP * coeff(147) * x(3) **2 * x(6) +            &
                                 coeff(148) * x(3) ** 3 +                  &
                                 coeff(149) * x(4) +                       &
                        2.0_DP * coeff(150) * x(4) * x(6) +                &
                        3.0_DP * coeff(151) * x(4) * x(6) **2 +            &
                                 coeff(152) * x(4) ** 2 +                  &
                        2.0_DP * coeff(153) * x(4) **2 * x(6) +            &
                                 coeff(154) * x(4) ** 3 +                  &
                                 coeff(155) * x(5) +                       &
                        2.0_DP * coeff(156) * x(5) * x(6) +                &
                        3.0_DP * coeff(157) * x(5) * x(6) **2 +            &
                                 coeff(158) * x(5) ** 2 +                  &
                        2.0_DP * coeff(159) * x(5) **2 * x(6) +            &
                                 coeff(160) * x(5) ** 3 +                  &
                                 coeff(161) * x(1) * x(2) +                &
                                 coeff(162) * x(1) * x(2) **2 +            &
                        2.0_DP * coeff(163) * x(1) * x(2) * x(6) +         &
                                 coeff(164) * x(1) **2 * x(2) +            &
                                 coeff(165) * x(1) * x(3) +                &
                                 coeff(166) * x(1) * x(3) **2 +            &
                        2.0_DP * coeff(167) * x(1) * x(3) * x(6) +         &
                                 coeff(168) * x(1) **2 * x(3) +            &
                                 coeff(169) * x(1) * x(4) +                &
                                 coeff(170) * x(1) * x(4) **2 +            &
                        2.0_DP * coeff(171) * x(1) * x(4) * x(6) +         &
                                 coeff(172) * x(1) **2 * x(4) +            &
                                 coeff(173) * x(1) * x(5) +                &
                                 coeff(174) * x(1) * x(5) **2 +            &
                        2.0_DP * coeff(175) * x(1) * x(5) * x(6) +         &
                                 coeff(176) * x(1) **2 * x(5) +            &
                                 coeff(177) * x(2) * x(3) +                &
                                 coeff(178) * x(2) * x(3) **2 +            &
                        2.0_DP * coeff(179) * x(2) * x(3) * x(6) +         &
                                 coeff(180) * x(2) **2 * x(3) +            &
                                 coeff(181) * x(2) * x(4) +                &
                                 coeff(182) * x(2) * x(4) **2 +            &
                        2.0_DP * coeff(183) * x(2) * x(4) * x(6) +         &
                                 coeff(184) * x(2) **2 * x(4) +            &
                                 coeff(185) * x(2) * x(5) +                &
                                 coeff(186) * x(2) * x(5) **2 +            &
                        2.0_DP * coeff(187) * x(2) * x(5) * x(6) +         &
                                 coeff(188) * x(2) **2 * x(5) +            &
                                 coeff(189) * x(3) * x(4) +                &
                                 coeff(190) * x(3) * x(4) **2 +            &
                        2.0_DP * coeff(191) * x(3) * x(4) * x(6) +         &
                                 coeff(192) * x(3) **2 * x(4) +            &
                                 coeff(193) * x(3) * x(5) +                &
                                 coeff(194) * x(3) * x(5) **2 +            &
                        2.0_DP * coeff(195) * x(3) * x(5) * x(6) +         &
                                 coeff(196) * x(3) **2 * x(5) +            &
                                 coeff(197) * x(4) * x(5) +                &
                                 coeff(198) * x(4) * x(5) **2 +            &
                        2.0_DP * coeff(199) * x(4) * x(5) * x(6) +         &
                                 coeff(200) * x(4) **2 * x(5) +            &
                                 coeff(201) * x(1) * x(2) * x(3) +         &
                                 coeff(202) * x(1) * x(2) * x(4) +         &
                                 coeff(203) * x(1) * x(2) * x(5) +         &
                                 coeff(204) * x(1) * x(3) * x(4) +         &
                                 coeff(205) * x(1) * x(3) * x(5) +         &
                                 coeff(206) * x(1) * x(4) * x(5) +         &
                                 coeff(207) * x(2) * x(3) * x(4) +         & 
                                 coeff(208) * x(2) * x(3) * x(5) +         &
                                 coeff(209) * x(2) * x(4) * x(5) +         &
                                 coeff(210) * x(3) * x(4) * x(5)
ENDIF

f=aux

RETURN
END SUBROUTINE evaluate_quartic_grad

SUBROUTINE evaluate_quartic_hessian(nvar,ncoeff,x,f,coeff)
!
!   This routine evaluates the hessian of the quartic polynomial
!   at the point x.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ncoeff
REAL(DP), INTENT(IN) :: x(nvar)
REAL(DP), INTENT(IN) :: coeff(ncoeff)
REAL(DP), INTENT(INOUT) :: f(nvar,nvar)

REAL(DP) :: aux(nvar,nvar)

aux(1,1) = 2.0_DP * coeff(3) + 6.0_DP * coeff(4) * x(1) + &
                              12.0_DP * coeff(5) * x(1)**2

IF (nvar>1) THEN
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

IF (nvar>2) THEN
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

IF (nvar>3) THEN
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
                       coeff(70) * x(2)*x(3)
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
              3.0_DP * coeff(57) * x(3)**2 + coeff(62) * x(1) +               &
              2.0_DP * coeff(63) * x(1) * x(3) + 2.0_DP*coeff(64)*x(1)*x(4) + &
                       coeff(65) * x(1)**2 + coeff(66) * x(2) +               &
              2.0_DP * coeff(67) * x(2)* x(3) + 2.0_DP*coeff(68)*x(2)*x(4) +  &
                       coeff(69) * x(2) **2  +                                &
                       coeff(70) * x(1) * x(2)
   aux(4,3) = aux(3,4)
   aux(4,4) = 2.0_DP * coeff(37) + 6.0_DP * coeff(38)*x(4)      &
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

IF (nvar>4) THEN
   aux(1,1) = aux(1,1) + 2.0_DP* coeff(78) * x(5) + 2.0_DP*coeff(79)*x(5)**2  &
                   + 6.0_DP * coeff(80) * x(1) * x(5)                         &
                   + 2.0_DP * coeff(102) * x(2) * x(5)                        &
                   + 2.0_DP * coeff(106) * x(3) * x(5)                        & 
                   + 2.0_DP * coeff(110) * x(4) * x(5)                        

   aux(1,2) = aux(1,2) + coeff(99) * x(5) +                                 &
                2.0_DP * coeff(100) * x(2) * x(5)                           &
                       + coeff(101) * x(5)**2 +                             &
                2.0_DP * coeff(102) * x(1) * x(5)                           &
                       + coeff(123) * x(3) * x(5)                           &
                       + coeff(124) * x(4) * x(5)                
   aux(2,1) = aux(1,2)

   aux(1,3) = aux(1,3) + coeff(103)*x(5)                                    &
              + 2.0_DP * coeff(104) * x(3) * x(5)                           &
              +          coeff(105) * x(5)**2                               &
              + 2.0_DP * coeff(106) * x(1) * x(5)                           &
              +          coeff(123) * x(2) * x(5)                           &
              +          coeff(125) * x(4) * x(5)
   aux(3,1) = aux(1,3)

   aux(1,4) = aux(1,4) + coeff(107)*x(5)                                    &
              + 2.0_DP * coeff(108) * x(4) * x(5)                           &
              +          coeff(109) * x(5)**2                               &
              + 2.0_DP * coeff(110) * x(1) * x(5)                           &
              +          coeff(124) * x(2) * x(5)                           &
              +          coeff(125) * x(3) * x(5)
   aux(4,1) = aux(1,4)

   aux(1,5) = coeff(75) + 2.0_DP*coeff(76)*x(5) + 3.0_DP*coeff(77)*x(5)**2 + &
              2.0_DP * coeff(78) * x(1) + 4.0_DP * coeff(79) * x(1) * x(5) + &
              3.0_DP * coeff(80) * x(1)**2 + coeff(99) * x(2) +              &
                       coeff(100) * x(2)**2 + 2.0_DP * coeff(101)*x(2)*x(5) +  &
                       2.0_DP * coeff(102) * x(1)*x(2)                       &
                   + coeff(103) * x(3) +  coeff(104) * x(3)**2 +             &
              2.0_DP * coeff(105)*x(3)*x(5) +                                &
                       2.0_DP * coeff(106) * x(1)*x(3)                       &
                    + coeff(107) * x(4) +  coeff(108) * x(4)**2              &
                    +  2.0_DP * coeff(109)*x(4)*x(5) +                       &
                       2.0_DP * coeff(110) * x(1)*x(4) +                     &
                        coeff(123) * x(2)*x(3) +                             &
                        coeff(124) * x(2)*x(4) +                             &
                        coeff(125) * x(3)*x(4)

   aux(5,1) = aux(1,5)
   aux(2,2) = aux(2,2) + 2.0_DP*coeff(84) * x(5)                           &
                       + 2.0_DP*coeff(85) * x(5)**2                        &
                       + 6.0_DP*coeff(86) * x(2) * x(5)                    &
                       + 2.0_DP*coeff(100) * x(1) * x(5)                   &
                       + 2.0_DP*coeff(114) * x(3) * x(5)                   &
                       + 2.0_DP*coeff(118) * x(4) * x(5) 

   aux(2,3) = aux(2,3) + coeff(111)*x(5)                                    &
              + 2.0_DP * coeff(112) * x(3) * x(5)                           &
              +          coeff(113) * x(5)**2                               &
              + 2.0_DP * coeff(114) * x(2) * x(5)                           &
              +          coeff(123) * x(1) * x(5)                           &
              +          coeff(126) * x(4) * x(5)
   aux(3,2) = aux(2,3)

   aux(2,4) = aux(2,4) + coeff(115)*x(5)                                    &
              + 2.0_DP * coeff(116) * x(4) * x(5)                           &
              +          coeff(117) * x(5)**2                               &
              + 2.0_DP * coeff(118) * x(2) * x(5)                           &
              +          coeff(124) * x(1) * x(5)                           &
              +          coeff(126) * x(3) * x(5)
   aux(4,2) = aux(2,4)


   aux(2,5) = coeff(81) + 2.0_DP*coeff(82)*x(5) + 3.0_DP*coeff(83)*x(5)**2 +  &
              2.0_DP * coeff(84) * x(2) + 4.0_DP * coeff(85) * x(2) * x(5) +  &
              3.0_DP * coeff(86) * x(2)**2 + coeff(99) * x(1) +               &
              2.0_DP * coeff(100) *x(1)*x(2) + 2.0_DP*coeff(101)*x(1)*x(5) +  &
                       coeff(102) * x(1)**2 + coeff(111)*x(3) +               &
                       coeff(112) * x(3)**2 + 2.0_DP * coeff(113)*x(3)*x(5) + &
                       2.0_DP * coeff(114) * x(2) * x(3) +                    &
                       coeff(115)*x(4) +                                      &
                       coeff(116) * x(4)**2 + 2.0_DP * coeff(117)*x(4)*x(5) + &
                       2.0_DP * coeff(118) * x(2) * x(4) +                    &
                                coeff(123) * x(1) * x(3) +                    &
                                coeff(124) * x(1) * x(4) +                    &
                                coeff(126) * x(3) * x(4)

   aux(5,2) = aux(2,5)
   aux(3,3) = aux(3,3) + 2.0_DP*coeff(90) * x(5)         &
                       + 2.0_DP*coeff(91) * x(5)**2      &
                       + 6.0_DP*coeff(92) * x(3) * x(5)  &
                       + 2.0_DP*coeff(104) * x(1) * x(5) &
                       + 2.0_DP*coeff(112) * x(2) * x(5) &
                       + 2.0_DP*coeff(122) * x(4) * x(5) 

   aux(3,4) = aux(3,4) + coeff(119)*x(5)                                    &
              + 2.0_DP * coeff(120) * x(4) * x(5)                           &
              +          coeff(121) * x(5)**2                               &
              + 2.0_DP * coeff(122) * x(3) * x(5)                           &
              +          coeff(125) * x(1) * x(5)                           &
              +          coeff(126) * x(2) * x(5)
   aux(4,3) = aux(3,4)

   aux(3,5) = coeff(87) + 2.0_DP*coeff(88)*x(5) + 3.0_DP*coeff(89)*x(5)**2 +  &
              2.0_DP * coeff(90) * x(3) + 4.0_DP * coeff(91) * x(3) * x(5) +  &
              3.0_DP * coeff(92) * x(3)**2 + coeff(103) * x(1) +              &
              2.0_DP * coeff(104) *x(1)*x(3) + 2.0_DP*coeff(105)*x(1)*x(5) +  &
                       coeff(106) *x(1)**2 + coeff(111) * x(2) +              &
              2.0_DP * coeff(112) *x(2)*x(3) + 2.0_DP*coeff(113)*x(2)*x(5) +  &
                       coeff(114) * x(2) **2  +                               &
                       coeff(119) * x(4) +                                    &
                       coeff(120) * x(4) **2  +                               &
              2.0_DP * coeff(121) * x(4)*x(5) +                               & 
              2.0_DP * coeff(122) * x(3)*x(4) +                               &
                       coeff(123) * x(1) * x(2) +                             &
                       coeff(125) * x(1) * x(4) +                             &
                       coeff(126) * x(2) * x(4)

   aux(5,3) = aux(3,5)

   aux(4,4) = aux(4,4) + 2.0_DP*coeff(96) * x(5)         &
                       + 2.0_DP*coeff(97) * x(5)**2      &
                       + 6.0_DP*coeff(98) * x(4) * x(5)  &
                       + 2.0_DP*coeff(108) * x(1) * x(5) &
                       + 2.0_DP*coeff(116) * x(2) * x(5) &
                       + 2.0_DP*coeff(120) * x(3) * x(5) 

   aux(4,5) = coeff(93) + 2.0_DP*coeff(94)*x(5) + 3.0_DP*coeff(95)*x(5)**2 +  &
              2.0_DP * coeff(96) * x(4) + 4.0_DP * coeff(97) * x(4) * x(5) +  &
              3.0_DP * coeff(98) * x(4)**2 + coeff(107) * x(1) +              &
              2.0_DP * coeff(108) *x(1)*x(4) + 2.0_DP*coeff(109)*x(1)*x(5) +  &
                       coeff(110) *x(1)**2 + coeff(115) * x(2) +              &
              2.0_DP * coeff(116) *x(2)*x(4) + 2.0_DP*coeff(117)*x(2)*x(5) +  &
                       coeff(118) * x(2) **2  +                               &
                       coeff(119) * x(3) +                                    &
              2.0_DP * coeff(120) * x(3)*x(4) +                               &
              2.0_DP * coeff(121) * x(3)*x(5) +                               & 
                       coeff(122) * x(3) **2  +                               &
                       coeff(124) * x(1) * x(2) +                             &
                       coeff(125) * x(1) * x(3) +                             &
                       coeff(126) * x(2) * x(3)

   aux(5,4) = aux(4,5)
!
   aux(5,5) = 2.0_DP * coeff(72) + 6.0_DP * coeff(73)*x(5)      &
                                 + 12.0_DP * coeff(74)*x(5)**2  &
                     + 2.0_DP * coeff(76) * x(1)                &
                     + 6.0_DP * coeff(77) * x(1) * x(5)         &
                     + 2.0_DP * coeff(79) * x(1) ** 2           &
                     + 2.0_DP * coeff(82) * x(2)                &
                     + 6.0_DP * coeff(83) * x(2) * x(5)         &
                     + 2.0_DP * coeff(85) * x(2) ** 2           &    
                     + 2.0_DP * coeff(88) * x(3)                &
                     + 6.0_DP * coeff(89) * x(3) * x(5)         &
                     + 2.0_DP * coeff(91) * x(3) ** 2           &    
                     + 2.0_DP * coeff(94) * x(4)                &
                     + 6.0_DP * coeff(95) * x(4) * x(5)         &
                     + 2.0_DP * coeff(97) * x(4) ** 2           &    
                     + 2.0_DP * coeff(101) * x(1) * x(2)        &   
                     + 2.0_DP * coeff(105) * x(1) * x(3)        &
                     + 2.0_DP * coeff(109) * x(1) * x(4)        &
                     + 2.0_DP * coeff(113) * x(2) * x(3)        &
                     + 2.0_DP * coeff(117) * x(2) * x(4)        &
                     + 2.0_DP * coeff(121) * x(3) * x(4)         
ENDIF

IF (nvar>5) THEN
   aux(1,1) = aux(1,1) + 2.0_DP* coeff(134) * x(6) + 2.0_DP*coeff(135)*x(6)**2 &
                   + 6.0_DP * coeff(136) * x(1) * x(6)                       &
                   + 2.0_DP * coeff(164) * x(2) * x(6)                       &
                   + 2.0_DP * coeff(168) * x(3) * x(6)                       & 
                   + 2.0_DP * coeff(172) * x(4) * x(6)                       & 
                   + 2.0_DP * coeff(176) * x(5) * x(6)                        

   aux(1,2) = aux(1,2) + coeff(161) * x(6) +                                &
                2.0_DP * coeff(162) * x(2) * x(6)                           &
                       + coeff(163) * x(6)**2 +                             &
                2.0_DP * coeff(164) * x(1) * x(6)                           &
                       + coeff(201) * x(3) * x(6)                           &
                       + coeff(202) * x(4) * x(6)                           &
                       + coeff(203) * x(5) * x(6)                
   aux(2,1) = aux(1,2)

   aux(1,3) = aux(1,3) + coeff(165)*x(6)                                    &
              + 2.0_DP * coeff(166) * x(3) * x(6)                           &
              +          coeff(167) * x(6)**2                               &
              + 2.0_DP * coeff(168) * x(1) * x(6)                           &
              +          coeff(201) * x(2) * x(6)                           &
              +          coeff(204) * x(4) * x(6)                           &
              +          coeff(205) * x(5) * x(6)
   aux(3,1) = aux(1,3)

   aux(1,4) = aux(1,4) + coeff(169)*x(6)                                    &
              + 2.0_DP * coeff(170) * x(4) * x(6)                           &
              +          coeff(171) * x(6)**2                               &
              + 2.0_DP * coeff(172) * x(1) * x(6)                           &
              +          coeff(202) * x(2) * x(6)                           &
              +          coeff(204) * x(3) * x(6)                           &
              +          coeff(206) * x(5) * x(6)
   aux(4,1) = aux(1,4)

   aux(1,5) = aux(1,5) + coeff(173)*x(6)                                    &
              + 2.0_DP * coeff(174) * x(5) * x(6)                           &
              +          coeff(175) * x(6)**2                               &
              + 2.0_DP * coeff(176) * x(1) * x(6)                           &
              +          coeff(203) * x(2) * x(6)                           &
              +          coeff(205) * x(3) * x(6)                           &
              +          coeff(206) * x(4) * x(6)
   aux(5,1) = aux(1,5)

   aux(1,6) = coeff(131)+2.0_DP*coeff(132)*x(6)+3.0_DP*coeff(133)*x(6)**2 +   &
              2.0_DP * coeff(134) * x(1) + 4.0_DP * coeff(135) * x(1) * x(6) +&
              3.0_DP * coeff(136) * x(1)**2 + coeff(161) * x(2) +             &
                       coeff(162) * x(2)**2 + 2.0_DP * coeff(163)*x(2)*x(6) + &
                       2.0_DP * coeff(164) * x(1)*x(2)                        &
                   + coeff(165) * x(3) +  coeff(166) * x(3)**2 +              &
              2.0_DP * coeff(167)*x(3)*x(6) +                                 &
                       2.0_DP * coeff(168) * x(1)*x(3)                        &
                    + coeff(169) * x(4) +  coeff(170) * x(4)**2               &
                    +  2.0_DP * coeff(171)*x(4)*x(6) +                        &
                       2.0_DP * coeff(172) * x(1)*x(4)                        &
                    + coeff(173) * x(5) +  coeff(174) * x(5)**2               &
                    +  2.0_DP * coeff(175)*x(5)*x(6) +                        &
                       2.0_DP * coeff(176) * x(1)*x(5) +                      &
                        coeff(201) * x(2)*x(3) +                              &
                        coeff(202) * x(2)*x(4) +                              &
                        coeff(203) * x(2)*x(5) +                              &
                        coeff(204) * x(3)*x(4) +                              &
                        coeff(205) * x(3)*x(5) +                              &
                        coeff(206) * x(4)*x(5)

   aux(6,1) = aux(1,6)
   aux(2,2) = aux(2,2) + 2.0_DP*coeff(140) * x(6)                          &
                       + 2.0_DP*coeff(141) * x(6)**2                       &
                       + 6.0_DP*coeff(142) * x(2) * x(6)                   &
                       + 2.0_DP*coeff(162) * x(1) * x(6)                   &
                       + 2.0_DP*coeff(180) * x(3) * x(6)                   &
                       + 2.0_DP*coeff(184) * x(4) * x(6)                   &
                       + 2.0_DP*coeff(188) * x(5) * x(6) 

   aux(2,3) = aux(2,3) + coeff(177)*x(6)                                    &
              + 2.0_DP * coeff(178) * x(3) * x(6)                           &
              +          coeff(179) * x(6)**2                               &
              + 2.0_DP * coeff(180) * x(2) * x(6)                           &
              +          coeff(201) * x(1) * x(6)                           &
              +          coeff(207) * x(4) * x(6)                           &
              +          coeff(208) * x(5) * x(6)

   aux(3,2) = aux(2,3)

   aux(2,4) = aux(2,4) + coeff(181) * x(6)                                  &
              + 2.0_DP * coeff(182) * x(4) * x(6)                           &
              +          coeff(183) * x(6)**2                               &
              + 2.0_DP * coeff(184) * x(2) * x(6)                           &
              +          coeff(202) * x(1) * x(6)                           &
              +          coeff(207) * x(3) * x(6)                           &
              +          coeff(209) * x(5) * x(6)
   aux(4,2) = aux(2,4)

   aux(2,5) = aux(2,5) + coeff(185) * x(6)                                  &
              + 2.0_DP * coeff(186) * x(5) * x(6)                           &
              +          coeff(187) * x(6)**2                               &
              + 2.0_DP * coeff(188) * x(2) * x(6)                           &
              +          coeff(203) * x(1) * x(6)                           &
              +          coeff(208) * x(3) * x(6)                           &
              +          coeff(209) * x(4) * x(6)
   aux(5,2) = aux(2,5)

   aux(2,6) = coeff(137) + 2.0_DP*coeff(138)*x(6) + 3.0_DP*coeff(139)*x(6)**2 +&
              2.0_DP * coeff(140) * x(2) + 4.0_DP * coeff(141)*x(2)*x(6) +    &
              3.0_DP * coeff(142) * x(2)**2 + coeff(161) * x(1) +             &
              2.0_DP * coeff(162) *x(1)*x(2) + 2.0_DP*coeff(163)*x(1)*x(6) +  &
                       coeff(164) * x(1)**2 + coeff(177)*x(3) +               &
                       coeff(178) * x(3)**2 + 2.0_DP * coeff(179)*x(3)*x(6) + &
                       2.0_DP * coeff(180) * x(2) * x(3) +                    &
                       coeff(181)*x(4) +                                      &
                       coeff(182) * x(4)**2 + 2.0_DP * coeff(183)*x(4)*x(6) + &
                       2.0_DP * coeff(184) * x(2) * x(4) +                    &
                       coeff(185)*x(5) +                                      &
                       coeff(186) * x(5)**2 + 2.0_DP * coeff(187)*x(5)*x(6) + &
                       2.0_DP * coeff(188) * x(2) * x(5) +                    &
                                coeff(201) * x(1) * x(3) +                    &
                                coeff(202) * x(1) * x(4) +                    &
                                coeff(203) * x(1) * x(5) +                    &
                                coeff(207) * x(3) * x(4) +                    &
                                coeff(208) * x(3) * x(5) +                    &
                                coeff(209) * x(4) * x(5) 

   aux(6,2) = aux(2,6)
   aux(3,3) = aux(3,3) + 2.0_DP*coeff(146) * x(6)        &
                       + 2.0_DP*coeff(147) * x(6)**2     &
                       + 6.0_DP*coeff(148) * x(3) * x(6) &
                       + 2.0_DP*coeff(166) * x(1) * x(6) &
                       + 2.0_DP*coeff(178) * x(2) * x(6) &
                       + 2.0_DP*coeff(192) * x(4) * x(6) &
                       + 2.0_DP*coeff(196) * x(5) * x(6) 

   aux(3,4) = aux(3,4) + coeff(189) * x(6)                                  &
              + 2.0_DP * coeff(190) * x(4) * x(6)                           &
              +          coeff(191) * x(6)**2                               &
              + 2.0_DP * coeff(192) * x(3) * x(6)                           &
              +          coeff(204) * x(1) * x(6)                           &
              +          coeff(207) * x(2) * x(6)                           &
              +          coeff(210) * x(5) * x(6)

   aux(4,3) = aux(3,4)

   aux(3,5) = aux(3,5) + coeff(193) * x(6)                                  &
              + 2.0_DP * coeff(194) * x(5) * x(6)                           &
              +          coeff(195) * x(6)**2                               &
              + 2.0_DP * coeff(196) * x(3) * x(6)                           &
              +          coeff(205) * x(1) * x(6)                           &
              +          coeff(208) * x(2) * x(6)                           &
              +          coeff(210) * x(4) * x(6)

   aux(5,3) = aux(3,5)

   aux(3,6) = coeff(143)+2.0_DP*coeff(144)*x(6)+3.0_DP*coeff(145)*x(6)**2 +  &
              2.0_DP * coeff(146) * x(3) + 4.0_DP*coeff(147)*x(3)*x(6) +  &
              3.0_DP * coeff(148) * x(3)**2 + coeff(165) * x(1) +              &
              2.0_DP * coeff(166) *x(1)*x(3) + 2.0_DP*coeff(167)*x(1)*x(6) +  &
                       coeff(168) *x(1)**2 + coeff(177) * x(2) +              &
              2.0_DP * coeff(178) *x(2)*x(3) + 2.0_DP*coeff(179)*x(2)*x(6) +  &
                       coeff(180) * x(2) **2  +                               &
                       coeff(189) * x(4) +                                    &
                       coeff(190) * x(4) **2  +                               &
              2.0_DP * coeff(191) * x(4)*x(6) +                               & 
              2.0_DP * coeff(192) * x(3)*x(4) +                               &
                       coeff(193) * x(5) +                                    &
                       coeff(194) * x(5) **2  +                               &
              2.0_DP * coeff(195) * x(5)*x(6) +                               & 
              2.0_DP * coeff(196) * x(3)*x(5) +                               &
                       coeff(201) * x(1) * x(2) +                             &
                       coeff(204) * x(1) * x(4) +                             &
                       coeff(205) * x(1) * x(5) +                             &
                       coeff(207) * x(2) * x(4) +                             &
                       coeff(208) * x(2) * x(5) +                             &
                       coeff(210) * x(4) * x(5) 

   aux(6,3) = aux(3,6)

   aux(4,4) = aux(4,4) + 2.0_DP*coeff(152) * x(6)         &
                       + 2.0_DP*coeff(153) * x(6)**2      &
                       + 6.0_DP*coeff(154) * x(4) * x(6)  &
                       + 2.0_DP*coeff(170) * x(1) * x(6)  &
                       + 2.0_DP*coeff(182) * x(2) * x(6)  &
                       + 2.0_DP*coeff(190) * x(3) * x(6)  &
                       + 2.0_DP*coeff(200) * x(5) * x(6) 

   aux(4,5) = aux(4,5) + coeff(197) * x(6)                                  &
              + 2.0_DP * coeff(198) * x(5) * x(6)                           &
              +          coeff(199) * x(6)**2                               &
              + 2.0_DP * coeff(200) * x(4) * x(6)                           &
              +          coeff(206) * x(1) * x(6)                           &
              +          coeff(209) * x(2) * x(6)                           &
              +          coeff(210) * x(3) * x(6)

   aux(5,4) = aux(4,5)

   aux(4,6) = coeff(149) + 2.0_DP*coeff(150)*x(6)+3.0_DP*coeff(151)*x(6)**2 + &
              2.0_DP * coeff(152) * x(4) + 4.0_DP*coeff(153)*x(4)*x(6) +      &
              3.0_DP * coeff(154) * x(4)**2 + coeff(169) * x(1) +             &
              2.0_DP * coeff(170) *x(1)*x(4) + 2.0_DP*coeff(171)*x(1)*x(6) +  &
                       coeff(172) *x(1)**2 + coeff(181) * x(2) +              &
              2.0_DP * coeff(182) *x(2)*x(4) + 2.0_DP*coeff(183)*x(2)*x(6) +  &
                       coeff(184) * x(2) **2  +                               &
                       coeff(189) * x(3) +                                    &
              2.0_DP * coeff(190) * x(3)*x(4) +                               &
              2.0_DP * coeff(191) * x(3)*x(6) +                               & 
                       coeff(192) * x(3) **2  +                               &
                       coeff(197) * x(5) +                                    &
                       coeff(198) * x(5) **2  +                               &
              2.0_DP * coeff(199) * x(5)*x(6) +                               & 
              2.0_DP * coeff(200) * x(4)*x(5) +                               &
                       coeff(202) * x(1) * x(2) +                             &
                       coeff(204) * x(1) * x(3) +                             &
                       coeff(206) * x(1) * x(5) +                             &
                       coeff(207) * x(2) * x(3) +                             &
                       coeff(209) * x(2) * x(5) +                             &
                       coeff(210) * x(3) * x(5)
   aux(6,4) = aux(4,6)

   aux(5,5) = aux(5,5) + 2.0_DP*coeff(158) * x(6)         &
                       + 2.0_DP*coeff(159) * x(6)**2      &
                       + 6.0_DP*coeff(160) * x(5) * x(6)  &
                       + 2.0_DP*coeff(174) * x(1) * x(6)  &
                       + 2.0_DP*coeff(186) * x(2) * x(6)  &
                       + 2.0_DP*coeff(194) * x(3) * x(6)  &
                       + 2.0_DP*coeff(198) * x(4) * x(6) 

   aux(5,6) = coeff(155) + 2.0_DP*coeff(156)*x(6)+3.0_DP*coeff(157)*x(6)**2 + &
              2.0_DP * coeff(158) * x(5) + 4.0_DP*coeff(159)*x(5)*x(6) +      &
              3.0_DP * coeff(160) * x(5)**2 + coeff(173) * x(1) +             &
              2.0_DP * coeff(174) *x(1)*x(5) + 2.0_DP*coeff(175)*x(1)*x(6) +  &
                       coeff(176) *x(1)**2 + coeff(185) * x(2) +              &
              2.0_DP * coeff(186) *x(2)*x(5) + 2.0_DP*coeff(187)*x(2)*x(6) +  &
                       coeff(188) * x(2) **2  +                               &
                       coeff(193) * x(3) +                                    &
              2.0_DP * coeff(194) * x(3)*x(5) +                               &
              2.0_DP * coeff(195) * x(3)*x(6) +                               & 
                       coeff(196) * x(3) **2  +                               &
                       coeff(197) * x(4) +                                    &
              2.0_DP * coeff(198) * x(4)*x(5) +                               &
              2.0_DP * coeff(199) * x(4)*x(6) +                               & 
                       coeff(200) * x(4) **2  +                               &
                       coeff(203) * x(1) * x(2) +                             &
                       coeff(205) * x(1) * x(3) +                             &
                       coeff(206) * x(1) * x(4) +                             &
                       coeff(208) * x(2) * x(3) +                             &
                       coeff(209) * x(2) * x(4) +                             &
                       coeff(210) * x(3) * x(4)
   aux(6,5) = aux(5,6)
!
   aux(6,6) = 2.0_DP * coeff(128) + 6.0_DP * coeff(129)*x(6)    &
                                 + 12.0_DP * coeff(130)*x(6)**2 &
                     + 2.0_DP * coeff(132) * x(1)                &
                     + 6.0_DP * coeff(133) * x(1) * x(6)         &
                     + 2.0_DP * coeff(135) * x(1) ** 2           &
                     + 2.0_DP * coeff(138) * x(2)                &
                     + 6.0_DP * coeff(139) * x(2) * x(6)         &
                     + 2.0_DP * coeff(141) * x(2) ** 2           &    
                     + 2.0_DP * coeff(144) * x(3)                &
                     + 6.0_DP * coeff(145) * x(3) * x(6)         &
                     + 2.0_DP * coeff(147) * x(3) ** 2           &    
                     + 2.0_DP * coeff(150) * x(4)                &
                     + 6.0_DP * coeff(151) * x(4) * x(6)         &
                     + 2.0_DP * coeff(153) * x(4) ** 2           &    
                     + 2.0_DP * coeff(156) * x(5)                &
                     + 6.0_DP * coeff(157) * x(5) * x(6)         &
                     + 2.0_DP * coeff(159) * x(5) ** 2           &    
                     + 2.0_DP * coeff(163) * x(1) * x(2)         &   
                     + 2.0_DP * coeff(167) * x(1) * x(3)         &
                     + 2.0_DP * coeff(171) * x(1) * x(4)         &
                     + 2.0_DP * coeff(175) * x(1) * x(5)         &
                     + 2.0_DP * coeff(179) * x(2) * x(3)         &
                     + 2.0_DP * coeff(183) * x(2) * x(4)         &
                     + 2.0_DP * coeff(187) * x(2) * x(5)         &
                     + 2.0_DP * coeff(191) * x(3) * x(4)         &  
                     + 2.0_DP * coeff(195) * x(3) * x(5)         &
                     + 2.0_DP * coeff(199) * x(4) * x(5)         
ENDIF

f(:,:)=aux(:,:)

RETURN
END SUBROUTINE evaluate_quartic_hessian
!
FUNCTION quartic_ncoeff(nvar)  
!
!  This function gives the number of coeffiecients of a quartic 
!  polynomial with nvar variables
!
IMPLICIT NONE
INTEGER :: quartic_ncoeff
INTEGER, INTENT(IN) :: nvar

IF (nvar==1) THEN
   quartic_ncoeff=5
ELSEIF (nvar==2) THEN
   quartic_ncoeff=15
ELSEIF (nvar==3) THEN
   quartic_ncoeff=35
ELSEIF (nvar==4) THEN
   quartic_ncoeff=70
ELSEIF (nvar==5) THEN
   quartic_ncoeff=126
ELSEIF (nvar==6) THEN
   quartic_ncoeff=210
ELSE
   quartic_ncoeff=0
ENDIF

RETURN
END FUNCTION quartic_ncoeff

SUBROUTINE find_quartic_extremum(nvar,ncoeff,x,f,coeff)
!
!  This routine starts from the point x and finds the extremum closest
!  to x of the quartic polynomial. In output x are the coordinates of 
!  the extremum and f the value of the quartic polynomial at the extremum.
!
USE linear_solvers, ONLY : linsolvx
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ncoeff
REAL(DP),INTENT(INOUT) :: x(nvar), f
REAL(DP),INTENT(IN) :: coeff(ncoeff)

INTEGER, PARAMETER :: maxiter=300

INTEGER :: iter, ideg
REAL(DP), PARAMETER :: tol=2.D-11
REAL(DP) :: g(nvar), y(nvar), xold(nvar)
REAL(DP) :: j(nvar, nvar) 
REAL(DP) :: deltax, fmod

xold(:)=x(:)
DO iter=1,maxiter
   !
   CALL evaluate_quartic_grad(nvar,ncoeff,x,g,coeff)
   !
   CALL evaluate_quartic_hessian(nvar,ncoeff,x,j,coeff)
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
CALL evaluate_fit_quartic(nvar,ncoeff,x,f,coeff)

RETURN
END SUBROUTINE find_quartic_extremum
!
SUBROUTINE print_quartic_polynomial(nvar, ncoeff, coeff)
!
!  This subroutine writes on output the coefficients of a quartic
!  polynomial with nvar variables
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nvar, ncoeff
REAL(DP), INTENT(IN) :: coeff(ncoeff)

  WRITE(stdout,'(/,5x,"Quartic polynomial:",/)') 
  WRITE(stdout,'(5x,    e20.7,11x,"+",e20.7," x1 ")') coeff(1), coeff(2) 
  WRITE(stdout,'(4x,"+",e20.7," x1^2",6x,"+",e20.7," x1^3")') coeff(3), &
                                                              coeff(4)
  WRITE(stdout,'(4x,"+",e20.7," x1^4")') coeff(5)
  IF (nvar>1) THEN
     WRITE(stdout,'(4x,"+",e20.7," x2",8x,"+",e20.7," x2^2")')  coeff(6), &
                                                                coeff(7)
     WRITE(stdout,'(4x,"+",e20.7," x2^3",6x,"+",e20.7," x2^4")') coeff(8), &
                                                                coeff(9)
     WRITE(stdout,'(4x,"+",e20.7," x1 x2",5x,"+",e20.7," x1 x2^2")') &
                                                          coeff(10), coeff(11)
     WRITE(stdout,'(4x,"+",e20.7," x1 x2^3",3x,"+",e20.7," x1^2 x2")') &
                                                          coeff(12), coeff(13)
     WRITE(stdout,'(4x,"+",e20.7," x1^2 x2^2 +", e20.7," x1^3 x2")') &
                                             coeff(14), coeff(15)
  ENDIF

  IF (nvar>2) THEN
     WRITE(stdout,'(4x,"+",e15.7," x3       +",e15.7," x3^2      +",&
                        &e15.7," x3^3")') coeff(16), coeff(17), coeff(18)
     WRITE(stdout,'(4x,"+",e15.7," x3^4     +",e15.7," x1 x3     +",e13.7,&
                              &" x1 x3^2")') coeff(19), coeff(20), coeff(21)
     WRITE(stdout,'(4x,"+",e15.7," x1 x3^3  +",e15.7," x1^2 x3   +",& 
                          &e15.7," x1^2 x3^2")') coeff(22), coeff(23), &
                                                            coeff(24)
     WRITE(stdout,'(4x,"+",e15.7," x1^3 x3   ")') coeff(25)
     WRITE(stdout,'(4x,"+",e15.7," x2 x3    +",e15.7," x2 x3^2   +",&
               &e15.7," x2 x3^3  ")') coeff(26), coeff(27), coeff(28)
     WRITE(stdout,'(4x,"+",e15.7," x2^2 x3  +",e15.7," x2^2 x3^2 +",& 
                          &e15.7," x2^3 x3")') coeff(29), coeff(30), coeff(31)
     WRITE(stdout,'(4x,"+",e15.7," x1 x2 x3 +",e15.7," x1 x2^2 x3+",& 
                          &e15.7," x1 x2 x3^2")') coeff(32), coeff(33), &
                                                                     coeff(34)
     WRITE(stdout,'(4x,"+",e15.7," x1^2 x2 x3")') coeff(35)
  ENDIF

  IF (nvar>3) THEN
     WRITE(stdout,'(4x,"+",e15.7," x4       +",e15.7," x4^2      +",&
                        &e15.7," x4^3")') coeff(36), coeff(37), coeff(38)
     WRITE(stdout,'(4x,"+",e15.7," x4^4     +",e15.7," x1 x4     +",e15.7,&
                              &" x1 x4^2")') coeff(39), coeff(40), coeff(41)
     WRITE(stdout,'(4x,"+",e15.7," x1 x4^3  +",e15.7," x1^2 x4   +",& 
                          &e15.7," x1^2 x4^2")') coeff(42), coeff(43), &
                                                            coeff(44)
     WRITE(stdout,'(4x,"+",e15.7," x1^3 x4   ")') coeff(45)
     WRITE(stdout,'(4x,"+",e15.7," x2 x4    +",e15.7," x2 x4^2   +",&
               &e15.7," x2 x4^3  ")') coeff(46), coeff(47), coeff(48)
     WRITE(stdout,'(4x,"+",e15.7," x2^2 x4  +",e15.7," x2^2 x4^2 +",& 
                          &e15.7," x2^3 x4")') coeff(49), coeff(50), coeff(51)
     WRITE(stdout,'(4x,"+",e15.7," x3 x4    +",e15.7," x3 x4^2   +",&
               &e15.7," x3 x4^3  ")') coeff(52), coeff(53), coeff(54)
     WRITE(stdout,'(4x,"+",e15.7," x3^2 x4  +",e15.7," x3^2 x4^2 +",& 
                          &e15.7," x3^3 x4")') coeff(55), coeff(56), coeff(57)
     WRITE(stdout,'(4x,"+",e15.7," x1 x2 x4 +",e15.7," x1 x2^2 x4+",& 
                          &e15.7," x1 x2 x4^2")') coeff(58), coeff(59), &
                                                              coeff(60)
     WRITE(stdout,'(4x,"+",e15.7," x1^2 x2 x4")') coeff(61)
     WRITE(stdout,'(4x,"+",e15.7," x1 x3 x4 +",e15.7," x1 x3^2 x4+",& 
                          &e15.7," x1 x3 x4^2")') coeff(62), coeff(63), &
                                                             coeff(64)
     WRITE(stdout,'(4x,"+",e15.7," x1^2 x3 x4")') coeff(65)
     WRITE(stdout,'(4x,"+",e15.7," x2 x3 x4 +",e15.7," x2 x3^2 x4+",& 
                          &e15.7," x2 x3 x4^2")') coeff(66), coeff(67), &
                                                             coeff(68)
     WRITE(stdout,'(4x,"+",e15.7," x2^2 x3 x4")') coeff(69)
     WRITE(stdout,'(4x,"+",e15.7," x1 x2 x3 x4")') coeff(70)
  ENDIF

  IF (nvar>4) THEN
     WRITE(stdout,'(4x,"+",e15.7," x5       +",e15.7," x5^2      +",&
                        &e15.7," x5^3")') coeff(71), coeff(72), coeff(73)
     WRITE(stdout,'(4x,"+",e15.7," x5^4     +",e15.7," x1 x5     +",e15.7,&
                              &" x1 x5^2")') coeff(74), coeff(75), coeff(76)
     WRITE(stdout,'(4x,"+",e15.7," x1 x5^3  +",e15.7," x1^2 x5   +",& 
                          &e15.7," x1^2 x5^2")') coeff(77), coeff(78), &
                                                            coeff(79)
     WRITE(stdout,'(4x,"+",e15.7," x1^3 x5   ")') coeff(80)
     WRITE(stdout,'(4x,"+",e15.7," x2 x5    +",e15.7," x2 x5^2   +",&
               &e15.7," x2 x5^3  ")') coeff(81), coeff(82), coeff(83)
     WRITE(stdout,'(4x,"+",e15.7," x2^2 x5  +",e15.7," x2^2 x5^2 +",& 
                          &e15.7," x2^3 x5")') coeff(84), coeff(85), coeff(86)
     WRITE(stdout,'(4x,"+",e15.7," x3 x5    +",e15.7," x3 x5^2   +",&
               &e15.7," x3 x5^3  ")') coeff(87), coeff(88), coeff(89)
     WRITE(stdout,'(4x,"+",e15.7," x3^2 x5  +",e15.7," x3^2 x5^2 +",& 
                          &e15.7," x3^3 x5")') coeff(90), coeff(91), coeff(92)
     WRITE(stdout,'(4x,"+",e15.7," x4 x5    +",e15.7," x4 x5^2   +",&
               &e15.7," x4 x5^3  ")') coeff(93), coeff(94), coeff(95)
     WRITE(stdout,'(4x,"+",e15.7," x4^2 x5  +",e15.7," x4^2 x5^2 +",& 
                          &e15.7," x4^3 x5")') coeff(96), coeff(97), coeff(98)
     WRITE(stdout,'(4x,"+",e15.7," x1 x2 x5 +",e15.7," x1 x2^2 x5+",& 
                          &e15.7," x1 x2 x5^2")') coeff(99), coeff(100), &
                                                              coeff(101)
     WRITE(stdout,'(4x,"+",e15.7," x1^2 x2 x5")') coeff(102)
     WRITE(stdout,'(4x,"+",e15.7," x1 x3 x5 +",e15.7," x1 x3^2 x5+",& 
                          &e15.7," x1 x3 x5^2")') coeff(103), coeff(104), &
                                                              coeff(105)
     WRITE(stdout,'(4x,"+",e15.7," x1^2 x3 x5")') coeff(106)
     WRITE(stdout,'(4x,"+",e15.7," x1 x4 x5 +",e15.7," x1 x4^2 x5+",& 
                          &e15.7," x1 x4 x5^2")') coeff(107), coeff(108), &
                                                              coeff(109)
     WRITE(stdout,'(4x,"+",e15.7," x1^2 x4 x5")') coeff(110)
     WRITE(stdout,'(4x,"+",e15.7," x2 x3 x5 +",e15.7," x2 x3^2 x5+",& 
                          &e15.7," x2 x3 x5^2")') coeff(111), coeff(112), &
                                                              coeff(113)
     WRITE(stdout,'(4x,"+",e15.7," x2^2 x3 x5")') coeff(114)
     WRITE(stdout,'(4x,"+",e15.7," x2 x4 x5 +",e15.7," x2 x4^2 x5+",& 
                          &e15.7," x2 x4 x5^2")') coeff(115), coeff(116), &
                                                              coeff(117)
     WRITE(stdout,'(4x,"+",e15.7," x2^2 x4 x5")') coeff(118)
     WRITE(stdout,'(4x,"+",e15.7," x3 x4 x5 +",e15.7," x3 x4^2 x5+",& 
                          &e15.7," x3 x4 x5^2")') coeff(119), coeff(120), &
                                                              coeff(121)
     WRITE(stdout,'(4x,"+",e15.7," x3^2 x4 x5")') coeff(122)

     WRITE(stdout,'(4x,"+",e15.7," x1x2x3x5 +",e15.7," x1x2x4x5  +",& 
                          &e15.7," x1x3x4x5")') coeff(123), coeff(124), &
                                                              coeff(125)
     WRITE(stdout,'(4x,"+",e15.7," x2x3x4x5")') coeff(126)
  ENDIF

  IF (nvar>5) THEN
     WRITE(stdout,'(4x,"+",e15.7," x6       +",e15.7," x6^2      +",      &
                        &e15.7," x6^3")') coeff(127), coeff(128), coeff(129)
     WRITE(stdout,'(4x,"+",e15.7," x6^4     +",e15.7," x1 x6     +",e15.7,&
                      &" x1 x6^2")') coeff(130), coeff(131), coeff(132)
     WRITE(stdout,'(4x,"+",e15.7," x1 x6^3  +",e15.7," x1^2 x6   +",      & 
               &e15.7," x1^2 x6^2")') coeff(133), coeff(134), coeff(135)
     WRITE(stdout,'(4x,"+",e15.7," x1^3 x6   ")') coeff(136)
     WRITE(stdout,'(4x,"+",e15.7," x2 x6    +",e15.7," x2 x6^2   +",      &
               &e15.7," x2 x6^3  ")') coeff(137), coeff(138), coeff(139)
     WRITE(stdout,'(4x,"+",e15.7," x2^2 x6  +",e15.7," x2^2 x6^2 +",      & 
               &e15.7," x2^3 x6")') coeff(140), coeff(141), coeff(142)
     WRITE(stdout,'(4x,"+",e15.7," x3 x6    +",e15.7," x3 x6^2   +",      &
               &e15.7," x3 x6^3  ")') coeff(143), coeff(144), coeff(145)
     WRITE(stdout,'(4x,"+",e15.7," x3^2 x6  +",e15.7," x3^2 x6^2 +",      & 
               &e15.7," x3^3 x6")') coeff(146), coeff(147), coeff(148)
     WRITE(stdout,'(4x,"+",e15.7," x4 x6    +",e15.7," x4 x6^2   +",      &
               &e15.7," x4 x6^3  ")') coeff(149), coeff(150), coeff(151)
     WRITE(stdout,'(4x,"+",e15.7," x4^2 x6  +",e15.7," x4^2 x6^2 +",      & 
               &e15.7," x4^3 x6")') coeff(152), coeff(153), coeff(154)
     WRITE(stdout,'(4x,"+",e15.7," x5 x6    +",e15.7," x5 x6^2   +",      &
               &e15.7," x5 x6^3  ")') coeff(155), coeff(156), coeff(157)
     WRITE(stdout,'(4x,"+",e15.7," x5^2 x6  +",e15.7," x5^2 x6^2 +",      & 
               &e15.7," x5^3 x6")') coeff(158), coeff(159), coeff(160)
     WRITE(stdout,'(4x,"+",e15.7," x1 x2 x6 +",e15.7," x1 x2^2 x6+",      & 
                          &e15.7," x1 x2 x6^2")') coeff(161), coeff(162), &
                                                              coeff(163)
     WRITE(stdout,'(4x,"+",e15.7," x1^2 x2 x6")') coeff(164)
     WRITE(stdout,'(4x,"+",e15.7," x1 x3 x6 +",e15.7," x1 x3^2 x6+",      & 
                          &e15.7," x1 x3 x6^2")') coeff(165), coeff(166), &
                                                              coeff(167)
     WRITE(stdout,'(4x,"+",e15.7," x1^2 x3 x6")') coeff(168)
     WRITE(stdout,'(4x,"+",e15.7," x1 x4 x6 +",e15.7," x1 x4^2 x6+",      & 
                          &e15.7," x1 x4 x6^2")') coeff(169), coeff(170), &
                                                              coeff(171)
     WRITE(stdout,'(4x,"+",e15.7," x1^2 x4 x6")') coeff(172)
     WRITE(stdout,'(4x,"+",e15.7," x1 x5 x6 +",e15.7," x1 x5^2 x6+",      & 
                          &e15.7," x1 x5 x6^2")') coeff(173), coeff(174), &
                                                              coeff(175)
     WRITE(stdout,'(4x,"+",e15.7," x1^2 x5 x6")') coeff(176)
     WRITE(stdout,'(4x,"+",e15.7," x2 x3 x6 +",e15.7," x2 x3^2 x6+",      & 
                          &e15.7," x2 x3 x6^2")') coeff(177), coeff(178), &
                                                              coeff(179)
     WRITE(stdout,'(4x,"+",e15.7," x2^2 x3 x6")') coeff(180)
     WRITE(stdout,'(4x,"+",e15.7," x2 x4 x6 +",e15.7," x2 x4^2 x6+",      & 
                          &e15.7," x2 x4 x6^2")') coeff(181), coeff(182), &
                                                              coeff(183)
     WRITE(stdout,'(4x,"+",e15.7," x2^2 x4 x6")') coeff(184)
     WRITE(stdout,'(4x,"+",e15.7," x2 x5 x6 +",e15.7," x2 x5^2 x6+",      & 
                          &e15.7," x2 x5 x6^2")') coeff(185), coeff(186), &
                                                              coeff(187)
     WRITE(stdout,'(4x,"+",e15.7," x2^2 x5 x6")') coeff(188)

     WRITE(stdout,'(4x,"+",e15.7," x3 x4 x6 +",e15.7," x3 x4^2 x6+",      & 
                          &e15.7," x3 x4 x6^2")') coeff(189), coeff(190), &
                                                              coeff(191)
     WRITE(stdout,'(4x,"+",e15.7," x3^2 x4 x6")') coeff(192)
     WRITE(stdout,'(4x,"+",e15.7," x3 x5 x6 +",e15.7," x3 x5^2 x6+",      & 
                          &e15.7," x3 x5 x6^2")') coeff(193), coeff(194), &
                                                              coeff(195)
     WRITE(stdout,'(4x,"+",e15.7," x3^2 x5 x6")') coeff(196)
     WRITE(stdout,'(4x,"+",e15.7," x4 x5 x6 +",e15.7," x4 x5^2 x6+",      & 
                          &e15.7," x4 x5 x6^2")') coeff(197), coeff(198), &
                                                              coeff(199)
     WRITE(stdout,'(4x,"+",e15.7," x4^2 x5 x6")') coeff(200)
     WRITE(stdout,'(4x,"+",e15.7," x1x2x3x6 +",e15.7," x1x2x4x6  +",      & 
                          &e15.7," x1x2x5x6")') coeff(201), coeff(202),   &
                                                              coeff(203)
     WRITE(stdout,'(4x,"+",e15.7," x1x3x4x6 +",e15.7," x1x3x5x6  +",      & 
                          &e15.7," x1x4x5x6")') coeff(204), coeff(205),   &
                                                              coeff(206)
     WRITE(stdout,'(4x,"+",e15.7," x2x3x4x6 +",e15.7," x2x3x5x6  +",      & 
                          &e15.7," x2x4x5x6")') coeff(207), coeff(208),   &
                                                              coeff(209)
     WRITE(stdout,'(4x,"+",e15.7," x3x4x5x6")') coeff(210)
  ENDIF

RETURN
END SUBROUTINE print_quartic_polynomial

SUBROUTINE introduce_quartic_fit(nvar, ncoeff, ndata)
!
!   This routine writes a small message with the information on the 
!   quartic polynomial
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
SUBROUTINE print_chisq_quartic(ndata, nvar, ncoeff, x, f, coeff)
!
!   This routine receives as input the values of a function f for ndata
!   values of the independent variables x, a set of ncoeff coefficients
!   of a quadratic interpolating polynomial and writes as output
!   the sum of the squares of the differences between the values of
!   the function and of the interpolating polynomial 
!
IMPLICIT NONE

INTEGER  :: ndata, nvar, ncoeff
REAL(DP) :: x(nvar, ndata), f(ndata), coeff(ncoeff)

REAL(DP) :: chisq, perc, aux
INTEGER  :: idata

chisq=0.0_DP
perc=0.0_DP
DO idata=1,ndata
   CALL evaluate_fit_quartic(nvar,ncoeff,x(1,idata),aux,coeff)
!     WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
   IF (ABS(f(idata))>1.D-12) perc= perc + ABS((f(idata)-aux) / f(idata))
ENDDO

WRITE(stdout,'(5x,"chi square quartic=",e18.5," relative error",e18.5,&
                                     &" %",/)') chisq, perc / ndata
RETURN
END SUBROUTINE print_chisq_quartic
!
SUBROUTINE evaluate_two_quartic(nvar,ncoeff,x,f,coeff,coeff1)
!
!  This routine evaluates the sum of two quartic polynomials at the point x
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ncoeff
REAL(DP),INTENT(INOUT) :: x(nvar), f
REAL(DP),INTENT(IN) :: coeff(ncoeff), coeff1(ncoeff)

REAL(DP) :: coeffadd4(ncoeff)

coeffadd4=coeff+coeff1
CALL evaluate_fit_quartic(nvar,ncoeff,x,f,coeffadd4)

RETURN
END SUBROUTINE evaluate_two_quartic
!
SUBROUTINE find_two_quartic_extremum(nvar,ncoeff,x,f,coeff,coeff1)
!
!  This routine starts from the point x and finds the extremum closest
!  to x of the sum of two quartic polynomial. In output x are the 
!  coordinates of the extremum and f the value of the sum of the
!  two quartic polynomial at the extremum.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ncoeff
REAL(DP),INTENT(INOUT) :: x(nvar), f
REAL(DP),INTENT(IN) :: coeff(ncoeff), coeff1(ncoeff)

REAL(DP) :: coeffadd4(ncoeff)

coeffadd4=coeff+coeff1
CALL find_quartic_extremum(nvar,ncoeff,x,f,coeffadd4)

RETURN
END SUBROUTINE find_two_quartic_extremum
!
SUBROUTINE print_chisq_two_quartic(ndata, nvar, ncoeff, x, f, coeff, coeff1)
!
!  This routine writes on output the chi squared of the sum of two 
!  quartic polynomials that interpolate the function f in the ndata
!  points x.
!
IMPLICIT NONE
INTEGER  :: ndata, nvar, ncoeff
REAL(DP) :: x(nvar, ndata), f(ndata), coeff(ncoeff), coeff1(ncoeff)

REAL(DP) :: chisq, aux
INTEGER  :: idata

chisq=0.0_DP
DO idata=1,ndata
   CALL evaluate_two_quartic(nvar,ncoeff,x(1,idata),aux,coeff,coeff1)
!  WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
ENDDO
WRITE(stdout,'(5x,"chi square two quartic=",e18.5,/)') chisq

RETURN
END SUBROUTINE print_chisq_two_quartic
!
SUBROUTINE set_quartic_linear_coefficients(nvar, ncoeff4, ncoeff, &
                                        coeffadd4, coeff4, coeff)
!
!   This subroutine adds the coefficients of a linear polynomial
!   to those of a quartic polynomial and evaluates the sum at the point x.
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nvar, ncoeff4, ncoeff
REAL(DP), INTENT(IN) :: coeff4(ncoeff4), coeff(ncoeff)
REAL(DP), INTENT(INOUT) :: coeffadd4(ncoeff4)

coeffadd4 = coeff4
coeffadd4(1) = coeffadd4(1) + coeff(1)
coeffadd4(2) = coeffadd4(2) + coeff(2)

IF (nvar > 1) coeffadd4(6)  = coeffadd4(6) + coeff(3)
IF (nvar > 2) coeffadd4(16) = coeffadd4(16) + coeff(4)
IF (nvar > 3) coeffadd4(36) = coeffadd4(36) + coeff(5)
IF (nvar > 4) coeffadd4(71) = coeffadd4(71) + coeff(6)
IF (nvar > 5) coeffadd4(127) = coeffadd4(127) + coeff(7)

RETURN
END SUBROUTINE set_quartic_linear_coefficients
!
SUBROUTINE find_quartic_linear_extremum(nvar,ncoeff4,ncoeff,x,f, &
                                                          coeff4,coeff)
!
!   This subroutine adds the coefficients of a linear polynomial
!   to those of a quartic polynomial and finds the minimum of the sum
!   of the two.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ncoeff, ncoeff4
REAL(DP), INTENT(IN) :: coeff(ncoeff), coeff4(ncoeff4) 
REAL(DP), INTENT(INOUT) :: x(nvar), f

REAL(DP) :: coeffadd4(ncoeff4)

CALL set_quartic_linear_coefficients(nvar, ncoeff4, ncoeff, coeffadd4, &
                                                           coeff4, coeff)
CALL find_quartic_extremum(nvar,ncoeff4,x,f,coeffadd4)

RETURN
END SUBROUTINE find_quartic_linear_extremum
!
SUBROUTINE set_quartic_quadratic_coefficients(nvar, ncoeff4, ncoeff, &
                                        coeffadd4, coeff4, coeff)
!
!   This subroutine adds the coefficients of a quadratic polynomial
!   to those of a quartic polynomial and evaluates the sum at the point x.
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nvar, ncoeff4, ncoeff
REAL(DP), INTENT(IN) :: coeff4(ncoeff4), coeff(ncoeff)
REAL(DP), INTENT(INOUT) :: coeffadd4(ncoeff4)

coeffadd4 = coeff4
coeffadd4(1) = coeffadd4(1) + coeff(1)
coeffadd4(2) = coeffadd4(2) + coeff(2)
coeffadd4(3) = coeffadd4(3) + coeff(3)

IF (nvar > 1) THEN
   coeffadd4(6)  = coeffadd4(6) + coeff(4)
   coeffadd4(7)  = coeffadd4(7) + coeff(5)
   coeffadd4(10) = coeffadd4(10) + coeff(6)
ENDIF

IF (nvar > 2) THEN
   coeffadd4(16) = coeffadd4(16) + coeff(7)
   coeffadd4(17) = coeffadd4(17) + coeff(8)
   coeffadd4(20) = coeffadd4(20) + coeff(9)
   coeffadd4(26) = coeffadd4(26) + coeff(10)
END IF

IF (nvar > 3) THEN
   coeffadd4(36) = coeffadd4(36) + coeff(11)
   coeffadd4(37) = coeffadd4(37) + coeff(12)
   coeffadd4(40) = coeffadd4(40) + coeff(13)
   coeffadd4(46) = coeffadd4(46) + coeff(14)
   coeffadd4(52) = coeffadd4(52) + coeff(15)
END IF

IF (nvar > 4) THEN
   coeffadd4(71) = coeffadd4(71) + coeff(16)
   coeffadd4(72) = coeffadd4(72) + coeff(17)
   coeffadd4(75) = coeffadd4(75) + coeff(18)
   coeffadd4(81) = coeffadd4(81) + coeff(19)
   coeffadd4(87) = coeffadd4(87) + coeff(20)
   coeffadd4(93) = coeffadd4(93) + coeff(21)
END IF

IF (nvar > 5) THEN
   coeffadd4(127) = coeffadd4(127) + coeff(22)
   coeffadd4(128) = coeffadd4(128) + coeff(23)
   coeffadd4(131) = coeffadd4(131) + coeff(24)
   coeffadd4(137) = coeffadd4(137) + coeff(25)
   coeffadd4(143) = coeffadd4(143) + coeff(26)
   coeffadd4(149) = coeffadd4(149) + coeff(27)
   coeffadd4(155) = coeffadd4(155) + coeff(28)
END IF

RETURN
END SUBROUTINE set_quartic_quadratic_coefficients
!
SUBROUTINE evaluate_quartic_quadratic(nvar,ncoeff4,ncoeff,x,f,coeff4,coeff)
!
!   This subroutine adds the coefficients of a quadratic polynomial
!   to those of a quartic polynomial and evaluates the sum at the point x.
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nvar, ncoeff, ncoeff4
REAL(DP), INTENT(IN) :: coeff(nvar), coeff4(ncoeff4) 
REAL(DP), INTENT(INOUT) :: x(nvar), f

REAL(DP) :: coeffadd4(ncoeff4)

CALL set_quartic_quadratic_coefficients(nvar, ncoeff4, ncoeff, coeffadd4, &
                                                          coeff4, coeff)
CALL evaluate_fit_quartic(nvar,ncoeff4,x,f,coeffadd4)

RETURN
END SUBROUTINE evaluate_quartic_quadratic
!
SUBROUTINE find_quartic_quadratic_extremum(nvar,ncoeff4,ncoeff,x,f, &
                                                          coeff4,coeff)
!
!   This subroutine adds the coefficients of a quadratic polynomial
!   to those of a quartic polynomial and finds the minimum of the sum
!   of the two.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ncoeff, ncoeff4
REAL(DP), INTENT(IN) :: coeff(ncoeff), coeff4(ncoeff4) 
REAL(DP), INTENT(INOUT) :: x(nvar), f

REAL(DP) :: coeffadd4(ncoeff4)

CALL set_quartic_quadratic_coefficients(nvar, ncoeff4, ncoeff, coeffadd4, &
                                                               coeff4, coeff)
CALL find_quartic_extremum(nvar,ncoeff4,x,f,coeffadd4)

RETURN
END SUBROUTINE find_quartic_quadratic_extremum
!
SUBROUTINE print_chisq_quartic_quadratic(ndata, nvar, ncoeff4, ncoeff, x, &
                                                          f, coeff4, coeff)
!
!  This routine writes on output the chi square of the sum of a
!  quadratic and a quartic polynomials that interpolate the function f in 
!  the ndata points x.
!
IMPLICIT NONE
INTEGER  :: ndata, nvar, ncoeff4, ncoeff
REAL(DP) :: x(nvar, ndata), f(ndata), coeff4(ncoeff4), coeff(ncoeff)

REAL(DP) :: chisq, aux
INTEGER  :: idata

chisq=0.0_DP
DO idata=1,ndata
   CALL evaluate_quartic_quadratic(nvar, ncoeff4, ncoeff, x(1,idata), aux, &
                                                        coeff4, coeff)
!  WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
ENDDO
WRITE(stdout,'(5x,"chi square quartic quadratic=",e18.5,/)') chisq

RETURN
END SUBROUTINE print_chisq_quartic_quadratic
!
! Copyright (C) 2018-19 Cristiano Malica
!
SUBROUTINE set_quartic_cubic_coefficients(nvar, ncoeff4, ncoeff, coeffadd4, &
                                                               coeff4, coeff)
!
!   This subroutines adds the coefficients of a quartic polynomial
!   to those of a cubic polynomial,
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nvar, ncoeff4, ncoeff
REAL(DP), INTENT(IN) :: coeff4(ncoeff4), coeff(ncoeff)
REAL(DP), INTENT(INOUT) :: coeffadd4(ncoeff4)

coeffadd4 = coeff4
coeffadd4(1) = coeffadd4(1) + coeff(1)
coeffadd4(2) = coeffadd4(2) + coeff(2)
coeffadd4(3) = coeffadd4(3) + coeff(3)
coeffadd4(4) = coeffadd4(4) + coeff(4)

IF (nvar > 1) THEN
   coeffadd4(6)  =  coeffadd4(6) +  coeff(5)
   coeffadd4(7)  =  coeffadd4(7) +  coeff(6)
   coeffadd4(8)  =  coeffadd4(8) +  coeff(7)
   coeffadd4(10) = coeffadd4(10) +  coeff(8)
   coeffadd4(11) = coeffadd4(11) +  coeff(9)
   coeffadd4(13) = coeffadd4(13) + coeff(10)
ENDIF

IF (nvar > 2) THEN
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

IF (nvar > 3) THEN
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

IF (nvar > 4) THEN
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

IF (nvar > 5) THEN
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
!
SUBROUTINE evaluate_quartic_cubic(nvar,ncoeff4,ncoeff,x,f,coeff4,coeff)
!
!   This subroutines adds the coefficients of a quartic polynomial
!   to those of a cubic polynomial and evaluates the resulting polynomial
!   at the point x.
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nvar, ncoeff, ncoeff4
REAL(DP), INTENT(IN) :: coeff(ncoeff), coeff4(ncoeff4) 
REAL(DP), INTENT(INOUT) :: x(nvar), f

REAL(DP) :: coeffadd4(ncoeff4)

CALL set_quartic_cubic_coefficients(nvar, ncoeff4, ncoeff, coeffadd4, &
                                                             coeff4, coeff)

CALL evaluate_fit_quartic(nvar,ncoeff4,x,f,coeffadd4)

RETURN
END SUBROUTINE evaluate_quartic_cubic
!
SUBROUTINE find_quartic_cubic_extremum(nvar,ncoeff4,ncoeff,x,f,coeff4,coeff)
!
!   This subroutines adds the coefficients of a quartic polynomial
!   to those of a cubic polynomial and finds the extremum of the sum
!   of the two closest to the input value of x.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar, ncoeff, ncoeff4
REAL(DP), INTENT(IN) :: coeff(ncoeff), coeff4(ncoeff4) 
REAL(DP), INTENT(INOUT) :: x(nvar), f

REAL(DP) :: coeffadd4(ncoeff4)

CALL set_quartic_cubic_coefficients(nvar, ncoeff4, ncoeff, coeffadd4, &
                                                   coeff4, coeff)
CALL find_quartic_extremum(nvar,ncoeff4,x,f,coeffadd4)

RETURN
END SUBROUTINE find_quartic_cubic_extremum
!
SUBROUTINE print_chisq_quartic_cubic(ndata, nvar, ncoeff4, ncoeff, x, &
                                                         f, coeff4, coeff)
!
!  This routine writes on output the chi square of the sum of a
!  cubic and a quartic polynomials that interpolate the function f in 
!  the ndata points x.
!
IMPLICIT NONE
INTEGER  :: ndata, nvar, ncoeff4, ncoeff
REAL(DP) :: x(nvar, ndata), f(ndata), coeff4(ncoeff4), coeff(ncoeff)

REAL(DP) :: chisq, aux
INTEGER  :: idata

chisq=0.0_DP
DO idata=1,ndata
   CALL evaluate_quartic_cubic(nvar, ncoeff4, ncoeff, x(1,idata), aux, &
                                                        coeff4, coeff)
!  WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
ENDDO
WRITE(stdout,'(5x,"chi square quartic cubic=",e18.5,/)') chisq

RETURN
END SUBROUTINE print_chisq_quartic_cubic

END MODULE quartic_surfaces
