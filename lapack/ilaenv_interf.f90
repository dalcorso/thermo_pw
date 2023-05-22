INTERFACE ilaenv_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION ILAENV_XG( ISPEC, NAME_, OPTS, N1, N2, N3, N4 )
#else
      FUNCTION ILAENV_XG( ISPEC, NAME_, OPTS, N1, N2, N3, N4 )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            ILAENV_XG
      CHARACTER          NAME_, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  ISPEC, NAME_, OPTS, N1, N2, N3, N4
#endif
!
      END FUNCTION ILAENV_XG
END INTERFACE ilaenv_interf
