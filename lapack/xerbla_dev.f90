!> \brief \b XERBLA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download XERBLA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/xerbla.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/xerbla.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/xerbla.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE XERBLA( SRNAME, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER*(*)      SRNAME
!       INTEGER            INFO
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> XERBLA  is an error handler for the LAPACK routines.
!> It is called by an LAPACK routine if an input parameter has an
!> invalid value.  A message is printed and execution stops.
!>
!> Installers may consider modifying the STOP statement in order to
!> call system-specific exception-handling facilities.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SRNAME
!> \verbatim
!>          SRNAME is CHARACTER*(*)
!>          The name of the routine which called XERBLA.
!> \endverbatim
!>
!> \param[in] INFO
!> \verbatim
!>          INFO is INTEGER
!>          The position of the invalid parameter in the parameter list
!>          of the calling routine.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup OTHERauxiliary
!
!  =====================================================================
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE XERBLA_XG( SRNAME, INFO )
#else
      SUBROUTINE XERBLA_XG( SRNAME, INFO )
#endif
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!A ---- ZGEMV           
!B ---- ZGERC          
!C ---- ZHEEV          
!D ---- ZHEGS2          
!E ---- ZHEGST          
!F ---- ZHEGV          
!G ---- ZHEEM          
!H ---- ZHEMV          
!I ---- ZHER2          
!J ---- ZHER2K          
!K ---- ZHERK          
!L ---- ZHETD2          
!M ---- ZHETRD          
!N ---- ZLASCL          
!O ---- ZLASR          
!P ---- ZPOTRF2          
!Q ---- ZPOTRF          
!R ---- ZSTEQR          
!S ---- ZTRMM          
!T ---- ZTRMV          
!U ---- ZTRSM          
!V ---- ZTRSV          
!W ---- ZUNG2L          
!X ---- ZUNG2R          
!Y ---- ZUNGQL          
!Z ---- ZGEMM          
!1 ---- ZUNGQR          
!2 ---- ZUNGTR          
!3 ---- DLASCL          
!4 ---- DLASRT          
!5 ---- DSTERF     
!6 ---- ZUNM2R
!7 ---- ZUNM2L
!8 ---- ZUNMQR
!9 ---- ZUNMQL
!0 ---- ZSTEIN
!! ---- ZHEGVX
!@ ---- ZHEEVX
!# ---- DSTEBZ
!$ ---- DLAGTS
!% ---- DLAGTF
!^ ---- ZUNMTR
!& ---- DLAGTF
!* ---- DLAGTF
!( ---- DLAGTF
!
!     .. Scalar Arguments ..
!     CHARACTER*(*)      SRNAME
      CHARACTER          SRNAME
      INTEGER            INFO
!     ..
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          LEN_TRIM
!     ..
!     .. Executable Statements ..
!

#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: SRNAME
      ATTRIBUTES(DEVICE) :: INFO
#endif

!     WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
      PRINT *, SRNAME, INFO
!
      STOP
!
!      9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ',   & 
!            'an illegal value' )
!
!     End of XERBLA_XG
!
      END

!  =====================================================================
