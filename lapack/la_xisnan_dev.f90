!module LA_XISNAN_XG
!   interface LA_ISNAN_XG
!
!   module procedure SISNAN_XG
!   module procedure DISNAN_XG
!
!   end interface
!
!contains

#if defined(__CUDA)
     ATTRIBUTES(DEVICE) logical function SISNAN_XG( x )
#else
     logical function SISNAN_XG( x )
#endif

   use LA_CONSTANTS_XG, only: wp=>sp
#ifdef USE_IEEE_INTRINSIC
   use, intrinsic :: ieee_arithmetic
#elif USE_ISNAN
   intrinsic :: isnan
#endif
   real(wp) :: x
#ifdef USE_IEEE_INTRINSIC
   sisnan_XG = ieee_is_nan(x)
#elif USE_ISNAN
   sisnan_XG = isnan(x)
#else
!   sisnan_XG = SLAISNAN_XG(x,x)
   sisnan_XG = (x.ne.x)

!   contains
!   logical function SLAISNAN_XG( x, y )
!   use LA_CONSTANTS_XG, only: wp=>sp
!   real(wp) :: x, y
!   SLAISNAN_XG = ( x.ne.y )
!   end function SLAISNAN_XG
#endif
   end function SISNAN_XG

#if defined(__CUDA)
     ATTRIBUTES(DEVICE) logical function DISNA1_XG ( x )
#else
     logical function DISNA1_XG( x )
#endif
! DISNAN -> DISNA1   
      use LA_CONSTANTS_XG, only: wp=>dp
#ifdef USE_IEEE_INTRINSIC
   use, intrinsic :: ieee_arithmetic
#elif USE_ISNAN
   intrinsic :: isnan
#endif
   real(wp) :: x
#ifdef USE_IEEE_INTRINSIC
   DISNA1_XG = ieee_is_nan(x)
#elif USE_ISNAN
   DISNA1_XG = isnan(x)
#else
!   DISNAN_XG = DLAISNAN_XG(x,x)
   DISNA1_XG = (x.ne.x)

!   contains
!   logical function DLAISNAN_XG( x, y )
!   use LA_CONSTANTS_XG, only: wp=>dp
!   real(wp) :: x, y
!   DLAISNAN_XG = ( x.ne.y )
!   end function DLAISNAN_XG
#endif
   end function DISNA1_XG

! end module LA_XISNAN_XG
