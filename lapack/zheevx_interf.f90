INTERFACE zheevx_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZHEEVX_XG( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
                         ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, &
                         IWORK, IFAIL, INFO )
#else
      SUBROUTINE ZHEEVX_XG( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
                         ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, &
                         IWORK, IFAIL, INFO )
#endif


!
!  -- LAPACK driver routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
!     ..
!     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( LDA, * ), WORK( * ), Z( LDZ, * )
!     ..
!
!     .. Executable Statements ..
!
!     Test the input parameters.
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  JOBZ, RANGE, UPLO, N, LDA, VL, VU, IL, IU, &
                             ABSTOL, LDZ, LWORK
      ATTRIBUTES(DEVICE) :: IFAIL, IWORK, RWORK, W, A, WORK, Z, M, INFO
#endif

!
!     End of ZHEEVX
!
      END SUBROUTINE ZHEEVX_XG
END INTERFACE zheevx_interf
