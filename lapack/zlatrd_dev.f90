!> \brief \b ZLATRD reduces the first nb rows and columns of a symmetric/Hermitian matrix A to real tridiagonal form by an unitary similarity transformation.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLATRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlatrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlatrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlatrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDW, N, NB
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   E( * )
!       COMPLEX*16         A( LDA, * ), TAU( * ), W( LDW, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLATRD reduces NB rows and columns of a complex Hermitian matrix A to
!> Hermitian tridiagonal form by a unitary similarity
!> transformation Q**H * A * Q, and returns the matrices V and W which are
!> needed to apply the transformation to the unreduced part of A.
!>
!> If UPLO = 'U', ZLATRD reduces the last NB rows and columns of a
!> matrix, of which the upper triangle is supplied;
!> if UPLO = 'L', ZLATRD reduces the first NB rows and columns of a
!> matrix, of which the lower triangle is supplied.
!>
!> This is an auxiliary routine called by ZHETRD.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          Hermitian matrix A is stored:
!>          = 'U': Upper triangular
!>          = 'L': Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The number of rows and columns to be reduced.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!>          n-by-n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n-by-n lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>          On exit:
!>          if UPLO = 'U', the last NB columns have been reduced to
!>            tridiagonal form, with the diagonal elements overwriting
!>            the diagonal elements of A; the elements above the diagonal
!>            with the array TAU, represent the unitary matrix Q as a
!>            product of elementary reflectors;
!>          if UPLO = 'L', the first NB columns have been reduced to
!>            tridiagonal form, with the diagonal elements overwriting
!>            the diagonal elements of A; the elements below the diagonal
!>            with the array TAU, represent the  unitary matrix Q as a
!>            product of elementary reflectors.
!>          See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
!>          elements of the last NB columns of the reduced matrix;
!>          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
!>          the first NB columns of the reduced matrix.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (N-1)
!>          The scalar factors of the elementary reflectors, stored in
!>          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
!>          See Further Details.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX*16 array, dimension (LDW,NB)
!>          The n-by-nb matrix W required to update the unreduced part
!>          of A.
!> \endverbatim
!>
!> \param[in] LDW
!> \verbatim
!>          LDW is INTEGER
!>          The leading dimension of the array W. LDW >= max(1,N).
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
!> \ingroup complex16OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  If UPLO = 'U', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(n) H(n-1) . . . H(n-nb+1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
!>  and tau in TAU(i-1).
!>
!>  If UPLO = 'L', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(1) H(2) . . . H(nb).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
!>  and tau in TAU(i).
!>
!>  The elements of the vectors v together form the n-by-nb matrix V
!>  which is needed, with W, to apply the transformation to the unreduced
!>  part of the matrix, using a Hermitian rank-2k update of the form:
!>  A := A - V*W**H - W*V**H.
!>
!>  The contents of A on exit are illustrated by the following examples
!>  with n = 5 and nb = 2:
!>
!>  if UPLO = 'U':                       if UPLO = 'L':
!>
!>    (  a   a   a   v4  v5 )              (  d                  )
!>    (      a   a   v4  v5 )              (  1   d              )
!>    (          a   1   v5 )              (  v1  1   a          )
!>    (              d   1  )              (  v1  v2  a   a      )
!>    (                  d  )              (  v1  v2  a   a   a  )
!>
!>  where d denotes a diagonal element of the reduced matrix, a denotes
!>  an element of the original matrix that is unchanged, and vi denotes
!>  an element of the vector defining H(i).
!> \endverbatim
!>
!  =====================================================================
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZLATRD_XG( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
#else
      SUBROUTINE ZLATRD_XG( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
#endif

#include<lsame_interf.f90>
#include<zaxpy_interf.f90>
#include<zdotc_interf.f90>
#include<zgemv_interf.f90>
#include<zhemv_interf.f90>
#include<zlacgv_interf.f90>
#include<zlarfg_interf.f90>
#include<zscal_interf.f90>
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   E( * )
      COMPLEX*16         A( LDA, * ), TAU( * ), W( LDW, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ZERO, ONE, HALF
      PARAMETER ( ZERO = ( 0.0D+0, 0.0D+0 ),   & 
        ONE = ( 1.0D+0, 0.0D+0 ),   & 
                         HALF = ( 0.5D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IW
      COMPLEX*16         ALPHA
!     ..
!     .. External Subroutines ..
!      EXTERNAL           ZAXPY_XG, ZGEMV_XG, ZHEMV_XG, ZLACGV_XG, ZLARFG_XG, ZSCAL_XG
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME_XG
!      COMPLEX*16         ZDOTC_XG
!      EXTERNAL           LSAME_XG, ZDOTC_XG
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!

#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: UPLO, N, NB, LDA, LDW
      ATTRIBUTES(DEVICE) ::  E, A, TAU, W
#endif


      IF( N.LE.0 )   & 
         RETURN
!
      IF( LSAME_XG( UPLO, 'U' ) ) THEN
!
!        Reduce last NB columns of upper triangle
!
         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF( I.LT.N ) THEN
!
!              Update A(1:i,i)
!
               A( I, I ) = DBLE( A( I, I ) )
               CALL ZLACGV_XG( N-I, W( I, IW+1 ), LDW )
!      CALL ZGEMV_XG( 'No transpose', I, N-I, -ONE, A( 1, I+1 ),   & 
      CALL ZGEMV_XG( 'N', I, N-I, -ONE, A( 1, I+1 ),   & 
                           LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL ZLACGV_XG( N-I, W( I, IW+1 ), LDW )
               CALL ZLACGV_XG( N-I, A( I, I+1 ), LDA )
!      CALL ZGEMV_XG( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ),   & 
      CALL ZGEMV_XG( 'N', I, N-I, -ONE, W( 1, IW+1 ),   & 
                           LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
               CALL ZLACGV_XG( N-I, A( I, I+1 ), LDA )
               A( I, I ) = DBLE( A( I, I ) )
            END IF
            IF( I.GT.1 ) THEN
!
!              Generate elementary reflector H(i) to annihilate
!              A(1:i-2,i)
!
               ALPHA = A( I-1, I )
               CALL ZLARFG_XG( I-1, ALPHA, A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = DBLE( ALPHA )
               A( I-1, I ) = ONE
!
!              Compute W(1:i-1,i)
!
!      CALL ZHEMV_XG( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1,   & 
      CALL ZHEMV_XG( 'U', I-1, ONE, A, LDA, A( 1, I ), 1,   & 
                           ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
!      CALL ZGEMV_XG( 'Conjugate transpose', I-1, N-I, ONE,   & 
      CALL ZGEMV_XG( 'C', I-1, N-I, ONE,   & 
        W( 1, IW+1 ), LDW, A( 1, I ), 1, ZERO,   & 
                              W( I+1, IW ), 1 )
!      CALL ZGEMV_XG( 'No transpose', I-1, N-I, -ONE,   & 
      CALL ZGEMV_XG( 'N', I-1, N-I, -ONE,   & 
        A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE,   & 
                              W( 1, IW ), 1 )
!      CALL ZGEMV_XG( 'Conjugate transpose', I-1, N-I, ONE,   & 
      CALL ZGEMV_XG( 'C', I-1, N-I, ONE,   & 
        A( 1, I+1 ), LDA, A( 1, I ), 1, ZERO,   & 
                              W( I+1, IW ), 1 )
!      CALL ZGEMV_XG( 'No transpose', I-1, N-I, -ONE,   & 
      CALL ZGEMV_XG( 'N', I-1, N-I, -ONE,   & 
        W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE,   & 
                              W( 1, IW ), 1 )
               END IF
               CALL ZSCAL_XG( I-1, TAU( I-1 ), W( 1, IW ), 1 )
      ALPHA = -HALF*TAU( I-1 )*ZDOTC_XG( I-1, W( 1, IW ), 1,   & 
                       A( 1, I ), 1 )
               CALL ZAXPY_XG( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF
!
   10    CONTINUE
      ELSE
!
!        Reduce first NB columns of lower triangle
!
         DO 20 I = 1, NB
!
!           Update A(i:n,i)
!
            A( I, I ) = DBLE( A( I, I ) )
            CALL ZLACGV_XG( I-1, W( I, 1 ), LDW )
!      CALL ZGEMV_XG( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ),   & 
      CALL ZGEMV_XG( 'N', N-I+1, I-1, -ONE, A( I, 1 ),   & 
                        LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL ZLACGV_XG( I-1, W( I, 1 ), LDW )
            CALL ZLACGV_XG( I-1, A( I, 1 ), LDA )
!      CALL ZGEMV_XG( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ),   & 
      CALL ZGEMV_XG( 'N', N-I+1, I-1, -ONE, W( I, 1 ),   & 
                        LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            CALL ZLACGV_XG( I-1, A( I, 1 ), LDA )
            A( I, I ) = DBLE( A( I, I ) )
            IF( I.LT.N ) THEN
!
!              Generate elementary reflector H(i) to annihilate
!              A(i+2:n,i)
!
               ALPHA = A( I+1, I )
      CALL ZLARFG_XG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1,   & 
                            TAU( I ) )
               E( I ) = DBLE( ALPHA )
               A( I+1, I ) = ONE
!
!              Compute W(i+1:n,i)
!
!      CALL ZHEMV_XG( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA,   & 
      CALL ZHEMV_XG( 'L', N-I, ONE, A( I+1, I+1 ), LDA,   & 
                           A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
!      CALL ZGEMV_XG( 'Conjugate transpose', N-I, I-1, ONE,   & 
      CALL ZGEMV_XG( 'C', N-I, I-1, ONE,   & 
        W( I+1, 1 ), LDW, A( I+1, I ), 1, ZERO,   & 
                           W( 1, I ), 1 )
!      CALL ZGEMV_XG( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ),   & 
      CALL ZGEMV_XG( 'N', N-I, I-1, -ONE, A( I+1, 1 ),   & 
                           LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
!      CALL ZGEMV_XG( 'Conjugate transpose', N-I, I-1, ONE,   & 
      CALL ZGEMV_XG( 'C', N-I, I-1, ONE,   & 
        A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO,   & 
                           W( 1, I ), 1 )
!      CALL ZGEMV_XG( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ),   & 
      CALL ZGEMV_XG( 'N', N-I, I-1, -ONE, W( I+1, 1 ),   & 
                           LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL ZSCAL_XG( N-I, TAU( I ), W( I+1, I ), 1 )
      ALPHA = -HALF*TAU( I )*ZDOTC_XG( N-I, W( I+1, I ), 1,   & 
                       A( I+1, I ), 1 )
               CALL ZAXPY_XG( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF
!
   20    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLATRD_XG
!
      END

