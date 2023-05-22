INTERFACE lsame_interf
#if defined(__CUDA)
ATTRIBUTES(DEVICE) FUNCTION LSAME_XG(CA,CB)
#else
      FUNCTION LSAME_XG(CA,CB)
#endif
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      LOGICAL   LSAME_XG
      CHARACTER CA,CB
!     .
#if defined(__CUDA)   
      ATTRIBUTES(VALUE) ::  CA,CB 
#endif

      END FUNCTION LSAME_XG
END INTERFACE lsame_interf
