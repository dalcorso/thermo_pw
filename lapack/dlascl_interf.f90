INTERFACE dlascl_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE DLASCL_XG( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
#else
      SUBROUTINE DLASCL_XG( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
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
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: TYPE, KL, KU, CFROM, CTO, M, N, LDA
      ATTRIBUTES(DEVICE) ::  A, INFO
#endif


      END SUBROUTINE DLASCL_XG
END INTERFACE dlascl_interf
