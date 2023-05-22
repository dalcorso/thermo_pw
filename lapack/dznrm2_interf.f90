INTERFACE dznrm2_interf
#if defined(__CUDA)
    ATTRIBUTES(DEVICE) function DZNRM2_XG( n, x, incx ) 
#else
    function DZNRM2_XG( n, x, incx ) 
#endif

   integer, parameter :: wp = kind(1.d0)
   real(wp) :: DZNRM2_XG
!  -- Reference BLAS level1 routine (version 3.9.1) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     March 2021
!
!  ..
!  .. Scalar Arguments ..
   integer :: incx, n
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(*)
!  ..
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: n, incx
      ATTRIBUTES(DEVICE) ::  x
#endif

end function DZNRM2_XG
END INTERFACE dznrm2_interf

