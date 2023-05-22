INTERFACE zladiv_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION ZLADIV_XG( X, Y )
#else
      FUNCTION ZLADIV_XG( X, Y )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      COMPLEX*16         ZLADIV_XG
      COMPLEX*16         X, Y
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) ::  X, Y
#endif


      END FUNCTION ZLADIV_XG

END INTERFACE zladiv_interf

