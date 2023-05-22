INTERFACE fftpack51_r8_factor
#if defined(__CUDA)
ATTRIBUTES(DEVICE) subroutine r8_factor ( n, nf, fac )
#else
subroutine r8_factor ( n, nf, fac )
#endif
!
  implicit none

  real ( kind = 8 ) fac(*)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nf

end subroutine r8_factor
END INTERFACE fftpack51_r8_factor
