INTERFACE dlaebz_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE DLAEBZ_XG( IJOB, NITMAX, N, MMAX, MINP, NBMIN, ABSTOL, &
                         RELTOL, PIVMIN, D, E, E2, NVAL, AB, C, MOUT, &
                         NAB, WORK, IWORK, INFO )
#else
      SUBROUTINE DLAEBZ_XG( IJOB, NITMAX, N, MMAX, MINP, NBMIN, ABSTOL, &
                         RELTOL, PIVMIN, D, E, E2, NVAL, AB, C, MOUT, &
                         NAB, WORK, IWORK, INFO )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IJOB, INFO, MINP, MMAX, MOUT, N, NBMIN, NITMAX
      DOUBLE PRECISION   ABSTOL, PIVMIN, RELTOL
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * ), NAB( MMAX, * ), NVAL( * )
      DOUBLE PRECISION AB( MMAX, * ), C( * ), D( * ), E( * ), E2( * ),   & 
                         WORK( * )
!     ..
!
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  IJOB, NITMAX, N, MMAX, MINP, NBMIN, ABSTOL, &
                         RELTOL, PIVMIN
      ATTRIBUTES(DEVICE) :: IWORK, NAB, NVAL, AB, C, D, E, E2, WORK, MOUT, INFO
#endif
      END SUBROUTINE DLAEBZ_XG
END INTERFACE dlaebz_interf
