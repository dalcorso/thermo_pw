INTERFACE dlaev2_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE DLAEV2_XG( A, B, C, RT1, RT2, CS1, SN1 )
#else
      SUBROUTINE DLAEV2_XG( A, B, C, RT1, RT2, CS1, SN1 )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
!     ..
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  A, B, C
      ATTRIBUTES(DEVICE) ::  CS1, RT1, RT2, SN1
#endif
!
      END SUBROUTINE DLAEV2_XG
END INTERFACE dlaev2_interf

