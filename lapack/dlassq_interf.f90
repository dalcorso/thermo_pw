INTERFACE dlassq_interf
#if defined(__CUDA)
ATTRIBUTES(DEVICE) subroutine DLASSQ_XG( n, x, incx, scl, sumsq )
#else
subroutine DLASSQ_XG( n, x, incx, scl, sumsq )
#endif

 
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
   use LA_CONSTANTS_XG, only: wp=>dp
!  .. Scalar Arguments ..
   integer :: incx, n
   real(wp) :: scl, sumsq
!  ..
!  .. Array Arguments ..
   real(wp) :: x(*)

#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: n, incx
      ATTRIBUTES(DEVICE) :: x, scl, sumsq
#endif   
!  ..
end subroutine DLASSQ_XG
END INTERFACE dlassq_interf

