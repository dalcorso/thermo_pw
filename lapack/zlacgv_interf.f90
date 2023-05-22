INTERFACE zlacgv_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZLACGV_XG( N, X, INCX )
#else
      SUBROUTINE ZLACGV_XG( N, X, INCX )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         X( * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  N, INCX 
      ATTRIBUTES(DEVICE) ::  X
#endif

      END SUBROUTINE ZLACGV_XG
END INTERFACE zlacgv_interf

