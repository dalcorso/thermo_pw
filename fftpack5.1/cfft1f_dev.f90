#if defined(__CUDA)
ATTRIBUTES(DEVICE) subroutine cfft1f ( n, inc, c, lenc, wsave, lensav, work, &
                                                               lenwrk, ier )
#else
subroutine cfft1f ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )
#endif

!*****************************************************************************80
!
!! CFFT1F: complex double precision forward fast Fourier transform, 1D.
!
!  Discussion:
!
!    CFFT1F computes the one-dimensional Fourier transform of a single 
!    periodic sequence within a complex array.  This transform is referred 
!    to as the forward transform or Fourier analysis, transforming the 
!    sequence from physical to spectral space.
!
!    This transform is normalized since a call to CFFT1F followed
!    by a call to CFFT1B (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in 
!    array C, of two consecutive elements within the sequence to be transformed.
!
!    Input/output, complex ( kind = 8 ) C(LENC) containing the sequence to 
!    be transformed.
!
!    Input, integer ( kind = 4 ) LENC, the dimension of the C array.  
!    LENC must be at least INC*(N-1) + 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to CFFT1I before the first call to routine CFFT1F 
!    or CFFT1B for a given transform length N.  WSAVE's contents may be re-used
!    for subsequent calls to CFFT1F and CFFT1B with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least 2*N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENC   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    20, input error returned by lower level routine.
!
#if defined(__CUDA)
  USE cudafor
#endif
  implicit none

#if defined(__CUDA)
#include<c1fm1f_interf.f90>
#endif

  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

!  complex ( kind = 8 ) cc(lenc)
  real ( kind = 8 ) :: c(2,lenc)
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) iw1
  integer ( kind = 4 ) n
  integer ( kind = 4 ) i, l
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
#if defined(__CUDA)
   ATTRIBUTES(DEVICE) :: c, work, wsave
   ATTRIBUTES(VALUE) :: lenc, lensav, lenwrk, n, inc
#endif
  ier = 0

!  if (lenc < inc*(n-1) + 1) then
!    ier = 1
!    call xerfft ('cfft1f ', 4)
!    call xerfft ( 4)
!  else if (lensav < 2*n + int(log( real ( n, kind = 8 ) ) &
!    /log( 2.0D+00 )) + 4) then
!    ier = 2
!!    call xerfft ('cfft1f ', 6)
!    call xerfft (6)
!  else if (lenwrk < 2*n) then
!    ier = 3
!    call xerfft (8)
!!    call xerfft ('cfft1f ', 8)
!  end if
!
!!  if (n == 1) then
!!    return
!!  end if

  iw1 = n + n + 1
  call c1fm1f ( n, inc, c, lenc, work, wsave, wsave(iw1), wsave(iw1+1) )

  return
end

