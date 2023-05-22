#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION DLAMC3_XG( A, B )
#else
      FUNCTION DLAMC3_XG( A, B )
#endif
!
!  -- LAPACK auxiliary routine --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   DLAMC3_XG
      DOUBLE PRECISION   A, B
!     ..
! =====================================================================
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  A, B
#endif
!     .. Executable Statements ..
!
      DLAMC3_XG = A + B
!
      RETURN
!
!     End of DLAMC3_XG
!
      END
!

!  =====================================================================
