INTERFACE ilazlr_interf

#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION ILAZLR_XG( M, N, A, LDA )
#else
      FUNCTION ILAZLR_XG( M, N, A, LDA )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            ILAZLR_XG
      INTEGER            M, N, LDA
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
!     ..
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: M, N, LDA
      ATTRIBUTES(DEVICE) ::  A
#endif

      END FUNCTION ILAZLR_XG
END INTERFACE ilazlr_interf
