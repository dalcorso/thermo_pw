INTERFACE fftpack51_r8_mcfti1
#if defined(__CUDA)
ATTRIBUTES(DEVICE) subroutine r8_mcfti1 ( n, wa, fnf, fac )
#else
subroutine r8_mcfti1 ( n, wa, fnf, fac )
#endif

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) fac(*)
  real ( kind = 8 ) fnf
  integer ( kind = 4 ) n
  real ( kind = 8 ) wa(*)
end subroutine r8_mcfti1
END INTERFACE fftpack51_r8_mcfti1

