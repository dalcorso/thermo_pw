INTERFACE zstebz_interf
!
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE DSTEBZ_XG( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E, &
                         M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK, &
                         INFO )
#else
      SUBROUTINE DSTEBZ_XG( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E, &
                         M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK, &
                         INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          ORDER, RANGE
      INTEGER            IL, INFO, IU, M, N, NSPLIT
      DOUBLE PRECISION   ABSTOL, VL, VU
!     ..
!     .. Array Arguments ..
      INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
!     ..
!

#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL 
      ATTRIBUTES(DEVICE) :: IBLOCK, ISPLIT, IWORK, D, E, W, WORK, M, INFO, NSPLIT
#endif
!     ..
!     .. Executable Statements ..
!
!
!     End of DSTEBZ
!
      END SUBROUTINE DSTEBZ_XG
END INTERFACE zstebz_interf
