INTERFACE zungql_interf

#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZUNGQL_XG( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
#else
      SUBROUTINE ZUNGQL_XG( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: M, N, K, LDA, LWORK
      ATTRIBUTES(DEVICE) ::  A, TAU, WORK, INFO
#endif
!
      END  SUBROUTINE ZUNGQL_XG
END INTERFACE zungql_interf

