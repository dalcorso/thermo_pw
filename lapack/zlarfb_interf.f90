INTERFACE zlarfb_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZLARFB_XG( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,T, LDT, C, LDC, WORK, LDWORK )
#else
      SUBROUTINE ZLARFB_XG( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,T, LDT, C, LDC, WORK, LDWORK )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 C( LDC, * ), T( LDT, * ), V( LDV, * ),   & 
                         WORK( LDWORK, * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::   SIDE, TRANS, DIRECT, STOREV, M, N, K, LDV, LDT, LDC, LDWORK
      ATTRIBUTES(DEVICE) ::  C, T, V, WORK
#endif

      END SUBROUTINE ZLARFB_XG
END INTERFACE zlarfb_interf

