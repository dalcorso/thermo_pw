INTERFACE zunm2r_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZUNM2R_XG( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, INFO )
#else
      SUBROUTINE ZUNM2R_XG( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  SIDE, TRANS, M, N, K, LDA, LDC
      ATTRIBUTES(DEVICE) ::  A, TAU, C, WORK, INFO
#endif

!
!     End of ZUNM2R
!
      END SUBROUTINE zunm2r_xg

END INTERFACE zunm2r_interf
