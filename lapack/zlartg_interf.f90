INTERFACE zlartg_interf
#if defined(__CUDA)
ATTRIBUTES(DEVICE) subroutine ZLARTG_XG( f, g, c, s, r )
#else
subroutine ZLARTG_XG( f, g, c, s, r )
#endif

   use LA_CONSTANTS_XG, &
   only: wp=>dp, zero=>dzero, one=>done, two=>dtwo, czero=>zzero, &
         safmin=>dsafmin, safmax=>dsafmax
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     February 2021
!
!  .. Scalar Arguments ..
   real(wp)           c
   complex(wp)        f, g, r, s
!  ..
#if defined(__CUDA)
   ATTRIBUTES(VALUE) :: f, g, s, r 
   ATTRIBUTES(DEVICE) :: c
#endif

end subroutine ZLARTG_XG
END INTERFACE zlartg_interf
