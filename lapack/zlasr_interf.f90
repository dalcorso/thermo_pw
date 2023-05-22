INTERFACE zlasr_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZLASR_XG( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
#else
      SUBROUTINE ZLASR_XG( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( * ), S( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: SIDE, PIVOT, DIRECT, M, N, LDA
      ATTRIBUTES(DEVICE) ::  C, S, A
#endif
!
      END SUBROUTINE ZLASR_XG
END INTERFACE zlasr_interf
