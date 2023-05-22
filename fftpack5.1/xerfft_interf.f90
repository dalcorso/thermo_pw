INTERFACE fftpack51_xerfft
#if defined(__CUDA)
ATTRIBUTES(DEVICE) subroutine xerfft ( info )
!ATTRIBUTES(DEVICE) subroutine xerfft ( srname, info )
#else
!subroutine xerfft ( srname, info )
subroutine xerfft ( info )
#endif

  implicit none

  integer ( kind = 4 ) info
!  character ( len = * ) srname

end subroutine xerfft
END INTERFACE fftpack51_xerfft
