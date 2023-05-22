INTERFACE zlascl_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZLASCL_XG( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
#else
      SUBROUTINE ZLASCL_XG( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: TYPE, KL, KU, CFROM, CTO, M, N, LDA
      ATTRIBUTES(DEVICE) ::  A, INFO
#endif

      END  SUBROUTINE ZLASCL_XG
END INTERFACE zlascl_interf
