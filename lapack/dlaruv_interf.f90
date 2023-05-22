INTERFACE zhegvx_subs
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE DLARUV_XG( ISEED, N, X )
#else
      SUBROUTINE DLARUV_XG( ISEED, N, X )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            N
!     ..
!     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   X( N )
!     ..
!
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: N
      ATTRIBUTES(DEVICE) :: ISEED, X
#endif

!
!     End of DLARUV
!
      END SUBROUTINE DLARUV_XG
END INTERFACE zhegvx_subs
