INTERFACE zung2l_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZUNG2L_XG( M, N, K, A, LDA, TAU, WORK, INFO )
#else
      SUBROUTINE ZUNG2L_XG( M, N, K, A, LDA, TAU, WORK, INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  M, N, K, LDA
      ATTRIBUTES(DEVICE) ::  A, TAU, WORK, INFO
#endif


      END SUBROUTINE ZUNG2L_XG
END INTERFACE zung2l_interf
