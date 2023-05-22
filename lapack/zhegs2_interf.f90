INTERFACE zhegs2_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZHEGS2_XG( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
#else
      SUBROUTINE ZHEGS2_XG( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
#endif
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  ITYPE, UPLO, N, LDA, LDB
      ATTRIBUTES(DEVICE) ::  A, B, INFO
#endif

      END SUBROUTINE ZHEGS2_XG
END INTERFACE zhegs2_interf
