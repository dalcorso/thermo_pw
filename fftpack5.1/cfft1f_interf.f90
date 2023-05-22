INTERFACE fftpack51_cfft1f 
#if defined(__CUDA)
ATTRIBUTES(DEVICE) subroutine cfft1f ( n, inc, c, lenc, wsave, lensav, work, &
                                                               lenwrk, ier )
#else
subroutine cfft1f ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )
#endif

!*****************************************************************************80
!
!
#if defined(__CUDA)
  USE cudafor
#endif
  implicit none

  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

!  complex ( kind = 8 ) cc(lenc)
  real ( kind = 8 ) c(2,lenc)
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
#if defined(__CUDA)
   ATTRIBUTES(DEVICE) :: c, work, wsave
   ATTRIBUTES(VALUE) :: lenc, lensav, lenwrk, n, inc
#endif

end subroutine cfft1f

END INTERFACE fftpack51_cfft1f 
