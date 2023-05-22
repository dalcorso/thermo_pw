INTERFACE zstein_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZSTEIN_XG( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, &
                         IWORK, IFAIL, INFO )
#else
      SUBROUTINE ZSTEIN_XG( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, &
                         IWORK, IFAIL, INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDZ, M, N
!     ..
!     .. Array Arguments ..
      INTEGER IBLOCK( * ), IFAIL( * ), ISPLIT( * ),   & 
                         IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
      COMPLEX*16         Z( LDZ, * )
!     ..
!
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  N, M, LDZ
      ATTRIBUTES(DEVICE) :: IBLOCK, IFAIL, ISPLIT, IWORK, D, E, W, WORK, Z, INFO
#endif
!
!     End of ZSTEIN
!
      END SUBROUTINE ZSTEIN_XG
END INTERFACE zstein_interf

