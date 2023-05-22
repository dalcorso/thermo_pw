INTERFACE zhegvx_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZHEGVX_XG( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, &
                         VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, &
                         LWORK, RWORK, IWORK, IFAIL, INFO )
#else
      SUBROUTINE ZHEGVX_XG( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, &
                         VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, &
                         LWORK, RWORK, IWORK, IFAIL, INFO )
#endif
!
!  -- LAPACK driver routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
!     ..
!     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16 A( LDA, * ), B( LDB, * ), WORK( * ),   & 
                         Z( LDZ, * )
!     ..
!
!     .. Executable Statements ..
!
!     Test the input parameters.
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  ITYPE, JOBZ, RANGE, UPLO, N, LDA, LDB, &
                         VL, VU, IL, IU, ABSTOL, LDZ, LWORK
      ATTRIBUTES(DEVICE) :: IFAIL, IWORK, RWORK, W, A, B, WORK, Z, M, INFO
#endif

!
!     End of ZHEGVX
!
      END SUBROUTINE ZHEGVX_XG
END INTERFACE zhegvx_interf
