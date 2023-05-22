INTERFACE zlassq_interf
#if defined(__CUDA)
ATTRIBUTES(DEVICE) subroutine ZLASSQ_XG( n, x, incx, scl, sumsq )
#else
subroutine ZLASSQ_XG( n, x, incx, scl, sumsq )
#endif

use LA_CONSTANTS_XG, only: wp=>dp
 
!  .. Scalar Arguments ..
   integer :: incx, n
   real(wp) :: scl, sumsq
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(*)
!  ..
#if defined(__CUDA)
   ATTRIBUTES(VALUE) ::  n, incx
   ATTRIBUTES(DEVICE) ::  x, scl, sumsq
#endif

end subroutine ZLASSQ_XG

END INTERFACE zlassq_interf
