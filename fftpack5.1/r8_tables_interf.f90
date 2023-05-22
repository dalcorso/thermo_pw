INTERFACE fftpack51_r8_tables
#if defined(__CUDA)
ATTRIBUTES(DEVICE) subroutine r8_tables ( ido, ip, wa )
#else
subroutine r8_tables ( ido, ip, wa )
#endif
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  real ( kind = 8 ) wa(ido,ip-1,2)

end subroutine r8_tables
END INTERFACE fftpack51_r8_tables
