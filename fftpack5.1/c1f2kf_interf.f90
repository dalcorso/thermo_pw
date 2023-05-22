INTERFACE fftpack51_c1f2kf
#if defined(__CUDA) 
ATTRIBUTES(DEVICE) subroutine c1f2kf ( ido, l1, na, cc, in1, ch, in2, wa )
#else
subroutine c1f2kf ( ido, l1, na, cc, in1, ch, in2, wa )
#endif
!*****************************************************************************80
!
!
#if defined(__CUDA)
  USE cudafor
#endif
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,2)
  real ( kind = 8 ) ch(in2,l1,2,ido)
  integer ( kind = 4 ) na
  real ( kind = 8 ) wa(ido,1,2)
#if defined(__CUDA)
   ATTRIBUTES(DEVICE) :: cc, ch, wa
   ATTRIBUTES(VALUE) :: ido, in1, in2, l1, na
#endif

end subroutine c1f2kf
END INTERFACE fftpack51_c1f2kf

