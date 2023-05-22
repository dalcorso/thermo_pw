INTERFACE zunmql_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZUNMQL_XG( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
#else
      SUBROUTINE ZUNMQL_XG( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!
!     Test the input arguments
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  SIDE, TRANS, M, N, K, LDA, LDC, LWORK 
      ATTRIBUTES(DEVICE) ::  A, C, TAU, WORK, INFO
#endif
!
!     End of ZUNMQL
!
      END SUBROUTINE ZUNMQL_XG
END INTERFACE zunmql_interf
