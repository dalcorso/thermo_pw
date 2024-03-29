!> \brief \b DCABS1
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DCABS1(Z)
!
!       .. Scalar Arguments ..
!       COMPLEX*16 Z
!       ..
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DCABS1 computes |Re(.)| + |Im(.)| of a double complex number
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] Z
!> \verbatim
!>          Z is COMPLEX*16
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup double_blas_level1
!
!  =====================================================================
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION DCABS1_XG(Z)
#else
      FUNCTION DCABS1_XG(Z)
#endif
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16 Z
!     ..
!     ..
!  =====================================================================
      DOUBLE PRECISION  DCABS1_XG
!     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,DIMAG

#if defined(__CUDA)
      ATTRIBUTES(DEVICE) :: Z
#endif

      DCABS1_XG = ABS(DBLE(Z)) + ABS(DIMAG(Z))
      RETURN
!
!     End of DCABS1
!
      END

