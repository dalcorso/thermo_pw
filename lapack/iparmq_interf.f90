INTERFACE iparmq_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) FUNCTION IPARMQ_XG( ISPEC, NAME_, OPTS, N, ILO, IHI, LWORK )
#else
      FUNCTION IPARMQ_XG( ISPEC, NAME_, OPTS, N, ILO, IHI, LWORK )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IPARMQ_XG
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER         NAME_, OPTS
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: ISPEC, NAME_, OPTS, N, ILO, IHI, LWORK
#endif
! 
      END FUNCTION IPARMQ_XG
END INTERFACE iparmq_interf

