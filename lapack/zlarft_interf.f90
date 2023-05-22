INTERFACE zlarft_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZLARFT_XG( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
#else
      SUBROUTINE ZLARFT_XG( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         T( LDT, * ), TAU( * ), V( LDV, * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: DIRECT, STOREV, N, K, LDV, LDT
      ATTRIBUTES(DEVICE) ::  T, TAU, V
#endif
!
      END  SUBROUTINE ZLARFT_XG
END INTERFACE zlarft_interf
