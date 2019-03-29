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
!   surfaces interpolation up to dimension 6. It is used to interpolate 
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
!      5                   126
!      6                   210
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
            quartic_var, find_quartic_extremum, &
            find_quartic_quadratic_extremum, evaluate_quartic_quadratic, &
            evaluate_two_quartic, find_two_quartic_extremum, &
            print_quartic_polynomial, introduce_quartic_fit, &
            print_chisq_quartic, print_chisq_two_quartic,    &
            print_chisq_quartic_quadratic

CONTAINS

SUBROUTINE fit_multi_quartic(ndata,degree,nvar,lsolve,x,f,coeff)
!
!  This routine receives as input a set of vectors x(degree,ndata) and
!  function values f(ndata) and gives as output the coefficients of
!  a quartic interpolating polynomial coeff(nvar). In input
!  ndata is the number of data points. degree is the number of degrees of
!  freedom or equivalently the number of independent parameters (the
!  maximum is 6), and nvar is the number of coefficients of the
!  intepolating quadrating polynomial. 
!        degree     nvar
!          1         5
!          2        15
!          3        35
!          4        70
!          5       126
!          6       210
!  lsolve can be 1, 2 or 3. It chooses the method to compute the
!  polynomial coefficients. Using 1 a matrix nvar x nvar is calculated
!                           Using 2 the overdetemined linear system is solved
!                           using QR or LQ factorization
!                           Using 3 the overdetermined linear system is solved
!                           using SVD decomposition
!                           If lsolve is not one of these values method 2 
!                           is used 
!
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
USE kinds, ONLY : DP
USE linear_solvers,     ONLY : linsolvx, linsolvms, linsolvsvd
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar, ndata
INTEGER, INTENT(INOUT) :: lsolve
REAL(DP), INTENT(IN) :: x(degree,ndata), f(ndata)
REAL(DP), INTENT(INOUT) :: coeff(nvar)

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
   amat(idata,3) = x(1,idata)*x(1,idata)
   amat(idata,4) = x(1,idata)*x(1,idata)*x(1,idata)
   amat(idata,5) = x(1,idata)*x(1,idata)*x(1,idata)*x(1,idata)

   IF (degree>1) THEN
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

   IF (degree>4) THEN
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

   IF (degree>5) THEN
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
   WRITE(stdout,'(5x,"Finding the quartic polynomial using &
                                                   &nvar x nvar matrix")')  
   CALL linsolvx(aa,nvar,b,coeff)
ELSEIF(lsolve==2) THEN
   WRITE(stdout,'(5x,"Finding the quartic polynomial using &
                                                   &QR factorization")')  
   CALL linsolvms(amat,ndata,nvar,f,coeff)
ELSEIF(lsolve==3) THEN
   WRITE(stdout,'(5x,"Finding the quartic polynomial using SVD")')  
   CALL linsolvsvd(amat,ndata,nvar,f,coeff)
ENDIF

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
!
!  one variable
!
aux = coeff(1) + x(1)*(coeff(2)+x(1)*(coeff(3)+x(1)*(coeff(4)+coeff(5)*x(1))))
!
!  two variables
!
IF (degree>1) THEN
   aux=aux+x(2)*(coeff(6) + x(2)*( coeff(7) + coeff(11) * x(1) +     &
                                             coeff(14) * x(1)**2 +  &
              x(2)*(coeff(8) + coeff(9) * x(2) + coeff(12) * x(1)))  &
             + x(1) *( coeff(10)    &
             + x(1) * (coeff(13) + coeff(15) * x(1) ) ) )
ENDIF
!
!  three variables
!
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
!
!  four variables
!
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
!
!  five variables
!
IF (degree>4) THEN
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
IF (degree>5) THEN
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

SUBROUTINE evaluate_fit_grad_quartic(degree,nvar,x,f,coeff)
!
!  computes the gradient of the quartic polynomial, up to 6 variables.
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: x(degree)
REAL(DP), INTENT(IN) :: coeff(nvar)
REAL(DP), INTENT(INOUT) :: f(degree)

REAL(DP) :: aux(degree)

IF (degree>6) CALL errore('evaluate_fit_grad_quartic','gradient not availble',1)

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

IF (degree>4) THEN
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

IF (degree>5) THEN
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

IF (degree>4) THEN
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

IF (degree>5) THEN
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
END SUBROUTINE evaluate_fit_hess_quartic

SUBROUTINE find_quartic_extremum(degree,nvar,x,f,coeff)
!
!  This routine starts from the point x and finds the extremum closest
!  to x. In output x are the coordinates of the extremum and f 
!  the value of the quartic function at the minimum
!
USE linear_solvers, ONLY : linsolvx
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar
REAL(DP),INTENT(INOUT) :: x(degree), f
REAL(DP),INTENT(IN) :: coeff(nvar)

INTEGER, PARAMETER :: maxiter=300

INTEGER :: iter, ideg
REAL(DP), PARAMETER :: tol=2.D-11
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

coeffadd4 = coeff4
coeffadd4(1) = coeffadd4(1) + coeff(1)
coeffadd4(2) = coeffadd4(2) + coeff(2)
coeffadd4(3) = coeffadd4(3) + coeff(3)

IF (degree > 1) THEN
   coeffadd4(6)  = coeffadd4(6) + coeff(4)
   coeffadd4(7)  = coeffadd4(7) + coeff(5)
   coeffadd4(10) = coeffadd4(10) + coeff(6)
ENDIF

IF (degree > 2) THEN
   coeffadd4(16) = coeffadd4(16) + coeff(7)
   coeffadd4(17) = coeffadd4(17) + coeff(8)
   coeffadd4(20) = coeffadd4(20) + coeff(9)
   coeffadd4(26) = coeffadd4(26) + coeff(10)
END IF

IF (degree > 3) THEN
   coeffadd4(36) = coeffadd4(36) + coeff(11)
   coeffadd4(37) = coeffadd4(37) + coeff(12)
   coeffadd4(40) = coeffadd4(40) + coeff(13)
   coeffadd4(46) = coeffadd4(46) + coeff(14)
   coeffadd4(52) = coeffadd4(52) + coeff(15)
END IF

IF (degree > 4) THEN
   coeffadd4(71) = coeffadd4(71) + coeff(16)
   coeffadd4(72) = coeffadd4(72) + coeff(17)
   coeffadd4(75) = coeffadd4(75) + coeff(18)
   coeffadd4(81) = coeffadd4(81) + coeff(19)
   coeffadd4(87) = coeffadd4(87) + coeff(20)
   coeffadd4(93) = coeffadd4(93) + coeff(21)
END IF

IF (degree > 5) THEN
   coeffadd4(127) = coeffadd4(127) + coeff(22)
   coeffadd4(128) = coeffadd4(128) + coeff(23)
   coeffadd4(131) = coeffadd4(131) + coeff(24)
   coeffadd4(137) = coeffadd4(137) + coeff(25)
   coeffadd4(143) = coeffadd4(143) + coeff(26)
   coeffadd4(149) = coeffadd4(149) + coeff(27)
   coeffadd4(155) = coeffadd4(155) + coeff(28)
END IF

RETURN
END SUBROUTINE set_quartic_coefficients

FUNCTION quartic_var(degree)  

IMPLICIT NONE
INTEGER :: quartic_var
INTEGER, INTENT(IN) :: degree

IF (degree==1) THEN
   quartic_var=5
ELSEIF (degree==2) THEN
   quartic_var=15
ELSEIF (degree==3) THEN
   quartic_var=35
ELSEIF (degree==4) THEN
   quartic_var=70
ELSEIF (degree==5) THEN
   quartic_var=126
ELSEIF (degree==6) THEN
   quartic_var=210
ELSE
   quartic_var=0
ENDIF

RETURN
END FUNCTION quartic_var

SUBROUTINE print_quartic_polynomial(degree, nvar, coeff)

USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE

INTEGER, INTENT(IN) :: degree, nvar
REAL(DP), INTENT(IN) :: coeff(nvar)

  WRITE(stdout,'(/,5x,"Quartic polynomial:",/)') 
  WRITE(stdout,'(5x,    e20.7,11x,"+",e20.7," x1 ")') coeff(1), coeff(2) 
  WRITE(stdout,'(4x,"+",e20.7," x1^2",6x,"+",e20.7," x1^3")') coeff(3), coeff(4)
  WRITE(stdout,'(4x,"+",e20.7," x1^4")') coeff(5)
  IF (degree>1) THEN
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

  IF (degree>2) THEN
     WRITE(stdout,'(4x,"+",e15.7," x3       +",e15.7," x3^2      +",&
                        &e15.7," x3^3")') coeff(16), coeff(17), coeff(18)
     WRITE(stdout,'(4x,"+",e15.7," x3^4     +",e15.7," x1 x3     +",e13.7,&
                              &" x1 x3^2")') coeff(19), coeff(20), coeff(21)
     WRITE(stdout,'(4x,"+",e15.7," x1 x3^3  +",e15.7," x1^2 x3   +",& 
                          &e15.7," x1^2 x3^2")') coeff(22), coeff(23), coeff(24)
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

  IF (degree>3) THEN
     WRITE(stdout,'(4x,"+",e15.7," x4       +",e15.7," x4^2      +",&
                        &e15.7," x4^3")') coeff(36), coeff(37), coeff(38)
     WRITE(stdout,'(4x,"+",e15.7," x4^4     +",e15.7," x1 x4     +",e15.7,&
                              &" x1 x4^2")') coeff(39), coeff(40), coeff(41)
     WRITE(stdout,'(4x,"+",e15.7," x1 x4^3  +",e15.7," x1^2 x4   +",& 
                          &e15.7," x1^2 x4^2")') coeff(42), coeff(43), coeff(44)
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

  IF (degree>4) THEN
     WRITE(stdout,'(4x,"+",e15.7," x5       +",e15.7," x5^2      +",&
                        &e15.7," x5^3")') coeff(71), coeff(72), coeff(73)
     WRITE(stdout,'(4x,"+",e15.7," x5^4     +",e15.7," x1 x5     +",e15.7,&
                              &" x1 x5^2")') coeff(74), coeff(75), coeff(76)
     WRITE(stdout,'(4x,"+",e15.7," x1 x5^3  +",e15.7," x1^2 x5   +",& 
                          &e15.7," x1^2 x5^2")') coeff(77), coeff(78), coeff(79)
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

  IF (degree>5) THEN
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

SUBROUTINE introduce_quartic_fit(degree, nvar, ndata)
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, nvar, ndata

WRITE(stdout,'(/,5x,"Fitting the data with a quartic polynomial:")')

WRITE(stdout,'(/,5x,"Number of variables:",10x,i5)')  degree
WRITE(stdout,'(5x,"Coefficients of the quartic:",2x,i5)')  nvar
WRITE(stdout,'(5x,"Number of fitting data:",7x,i5,/)')  ndata

RETURN
END SUBROUTINE introduce_quartic_fit


SUBROUTINE print_chisq_quartic(ndata, degree, nvar, x, f, coeff)
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER  :: ndata, degree, nvar
REAL(DP) :: x(degree, ndata), f(ndata), coeff(nvar)

REAL(DP) :: chisq, perc, aux
INTEGER  :: idata

chisq=0.0_DP
perc=0.0_DP
DO idata=1,ndata
   CALL evaluate_fit_quartic(degree,nvar,x(1,idata),aux,coeff)
!     WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
   IF (ABS(f(idata))>1.D-12) perc= perc + ABS((f(idata)-aux) / f(idata))
ENDDO

WRITE(stdout,'(5x,"chi square quartic=",e18.5," relative error",e18.5,&
                                     &" %",/)') chisq, perc / ndata
RETURN
END SUBROUTINE print_chisq_quartic

SUBROUTINE print_chisq_two_quartic(ndata, degree, nvar, x, f, coeff, coeff1)

USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER  :: ndata, degree, nvar
REAL(DP) :: x(degree, ndata), f(ndata), coeff(nvar), coeff1(nvar)

REAL(DP) :: chisq, aux
INTEGER  :: idata

chisq=0.0_DP
DO idata=1,ndata
   CALL evaluate_two_quartic(degree,nvar,x(1,idata),aux,coeff,coeff1)
!  WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
ENDDO
WRITE(stdout,'(5x,"chi square two quartic=",e18.5,/)') chisq

RETURN
END SUBROUTINE print_chisq_two_quartic

SUBROUTINE print_chisq_quartic_quadratic(ndata, degree, nvar4, nvar, x, &
                                                          f, coeff4, coeff)
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER  :: ndata, degree, nvar4, nvar
REAL(DP) :: x(degree, ndata), f(ndata), coeff4(nvar4), coeff(nvar)

REAL(DP) :: chisq, aux
INTEGER  :: idata

chisq=0.0_DP
DO idata=1,ndata
   CALL evaluate_quartic_quadratic(degree, nvar4, nvar, x(1,idata), aux, &
                                                        coeff4, coeff)
!  WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
   chisq = chisq + (aux - f(idata))**2
ENDDO
WRITE(stdout,'(5x,"chi square quartic quadratic=",e18.5,/)') chisq

RETURN
END SUBROUTINE print_chisq_quartic_quadratic

END MODULE quartic_surfaces
