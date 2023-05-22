INTERFACE zhegv_interf
!
!  =====================================================================
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZHEGV_XG( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, LWORK, RWORK, INFO )
#else
      SUBROUTINE ZHEGV_XG( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, LWORK, RWORK, INFO )
#endif

!
!  -- LAPACK driver routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
!     ..
!

#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: ITYPE, JOBZ, UPLO, N, LDA, LDB, LWORK
      ATTRIBUTES(DEVICE) ::  RWORK, W, A, B, WORK, INFO
#endif

      END SUBROUTINE ZHEGV_XG

END INTERFACE zhegv_interf
