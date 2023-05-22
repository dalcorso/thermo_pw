INTERFACE zheev_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZHEEV_XG( JOBZ, UPLO, N, A, LDA, W, &
                                           WORK, LWORK, RWORK, INFO )
#else
      SUBROUTINE ZHEEV_XG( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
#endif
!
!  -- LAPACK driver routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
!     ..

#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  JOBZ, UPLO, N, LDA, LWORK
      ATTRIBUTES(DEVICE) ::  RWORK, W, A, WORK, INFO
#endif


      END SUBROUTINE ZHEEV_XG

END INTERFACE zheev_interf
