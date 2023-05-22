INTERFACE dnrm2_interf
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) function DNRM2_XG( n, x, incx ) 
#else
      function DNRM2_XG( n, x, incx ) 
#endif

   integer, parameter :: wp = kind(1.d0)
   real(wp) :: DNRM2_XG
!  .. Scalar Arguments ..
   integer :: incx, n
!  ..
!  .. Array Arguments ..
   real(wp) :: x(*)
   
!
!  -- Reference BLAS level1 routine (version 3.9.1) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     March 2021
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) ::  n, incx
      ATTRIBUTES(DEVICE) ::  x 
#endif
end function dnrm2_xg
END INTERFACE dnrm2_interf

