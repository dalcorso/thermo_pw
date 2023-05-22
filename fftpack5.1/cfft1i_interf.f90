INTERFACE fftpack51_cfft1i
#if defined(__CUDA)
ATTRIBUTES(DEVICE) subroutine cfft1i ( n, wsave, lensav, ier )
#else
subroutine cfft1i ( n, wsave, lensav, ier )
#endif

  implicit none

  integer ( kind = 4 ) lensav

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) n
  real ( kind = 8 ) wsave(lensav)

end subroutine cfft1i
END INTERFACE fftpack51_cfft1i

