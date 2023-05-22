INTERFACE zsteqr_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZSTEQR_XG( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
#else
      SUBROUTINE ZSTEQR_XG( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * )
      COMPLEX*16         Z( LDZ, * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  COMPZ, N, LDZ
      ATTRIBUTES(DEVICE) ::  D, E, WORK, Z, INFO
#endif

      END SUBROUTINE ZSTEQR_XG
END INTERFACE zsteqr_interf
