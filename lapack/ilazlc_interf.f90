INTERFACE ilazlc_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION ILAZLC_XG( M, N, A, LDA )
#else
      FUNCTION ILAZLC_XG( M, N, A, LDA )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            ILAZLC_XG
      INTEGER            M, N, LDA
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
!     ..
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: M, N, LDA
      ATTRIBUTES(DEVICE) ::  A
#endif


      END FUNCTION ILAZLC_XG
END INTERFACE ilazlc_interf
