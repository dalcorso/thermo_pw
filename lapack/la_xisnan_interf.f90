INTERFACE laxinsan_xg_interf

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
   end function SISNAN_XG

#if defined(__CUDA)
     ATTRIBUTES(DEVICE) logical function DISNA1_XG( x )
#else
     logical function DISNA1_XG( x )
#endif
   
      use LA_CONSTANTS_XG, only: wp=>dp
#ifdef USE_IEEE_INTRINSIC
   use, intrinsic :: ieee_arithmetic
#elif USE_ISNAN
   intrinsic :: isnan
#endif
   real(wp) :: x
   end function DISNA1_XG

END INTERFACE laxinsan_xg_interf
