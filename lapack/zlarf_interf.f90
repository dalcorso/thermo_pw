INTERFACE zlarf_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZLARF_XG( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
#else
      SUBROUTINE ZLARF_XG( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX*16         TAU
!     ..
!     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: SIDE, M, N, INCV, LDC
      ATTRIBUTES(DEVICE) ::  C, V, WORK,TAU
#endif
!
      END  SUBROUTINE ZLARF_XG
END INTERFACE zlarf_interf
