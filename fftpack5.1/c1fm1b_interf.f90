INTERFACE fftpack51_c1fm1b
#if defined(__CUDA)
ATTRIBUTES(DEVICE) subroutine c1fm1b ( n, inc, c, lenc, ch, wa, fnf, fac )
#else
   subroutine c1fm1b ( n, inc, c, lenc, ch, wa, fnf, fac )
#endif

!*****************************************************************************80
!
#if defined(__CUDA)
  USE cudafor
#endif
  implicit none

!  complex ( kind = 8 ) c(*)
  real ( kind = 8 ) c(2,lenc)
  real ( kind = 8 ) ch(*)
  real ( kind = 8 ) fac(*)
  real ( kind = 8 ) fnf
  integer ( kind = 4 ) inc,lenc
  integer ( kind = 4 ) n
  real ( kind = 8 ) wa(*)
#if defined(__CUDA)
   ATTRIBUTES(DEVICE) :: c, ch, fac, wa
   ATTRIBUTES(VALUE) :: lenc, n, inc, fnf
#endif


end subroutine c1fm1b
END INTERFACE fftpack51_c1fm1b
