INTERFACE dlartg_interf
#if defined(__CUDA)
ATTRIBUTES(DEVICE) subroutine DLARTG_XG( f, g, c, s, r )
#else
subroutine DLARTG_XG( f, g, c, s, r )
#endif

   use LA_CONSTANTS_XG, only: wp=>dp

!  .. Scalar Arguments ..
   real(wp) :: c, f, g, r, s
!  ..
!
#if defined(__CUDA)
    ATTRIBUTES(VALUE) ::  f, g
    ATTRIBUTES(DEVICE) ::  c,s, r
#endif

end subroutine DLARTG_XG
END INTERFACE dlartg_interf

