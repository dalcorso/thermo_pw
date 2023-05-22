INTERFACE fftpack51_c1fgkf
#if defined(__CUDA) 
ATTRIBUTES(DEVICE) subroutine c1fgkf ( ido, ip, l1, lid, na, cc, cc1, in1, &
                                                          ch, ch1, in2, wa )
#else
   subroutine c1fgkf ( ido, ip, l1, lid, na, cc, cc1, in1, ch, ch1, in2, wa )
#endif

!*****************************************************************************80
!
#if defined(__CUDA)
  USE cudafor
#endif
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) lid

  real ( kind = 8 ) cc(in1,l1,ip,ido)
  real ( kind = 8 ) cc1(in1,lid,ip)
  real ( kind = 8 ) ch(in2,l1,ido,ip)
  real ( kind = 8 ) ch1(in2,lid,ip)
  integer ( kind = 4 ) na
  real ( kind = 8 ) wa(ido,ip-1,2)

#if defined(__CUDA)
   ATTRIBUTES(DEVICE) :: cc, cc1, ch, ch1, wa
   ATTRIBUTES(VALUE) :: ido, in1, in2, ip, l1, lid, na
#endif

end subroutine c1fgkf
END INTERFACE fftpack51_c1fgkf

