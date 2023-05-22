!> \brief \b ZPOTRF2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       RECURSIVE SUBROUTINE ZPOTRF2( UPLO, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZPOTRF2 computes the Cholesky factorization of a Hermitian
!> positive definite matrix A using the recursive algorithm.
!>
!> The factorization has the form
!>    A = U**H * U,  if UPLO = 'U', or
!>    A = L  * L**H,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is lower triangular.
!>
!> This is the recursive version of the algorithm. It divides
!> the matrix into four submatrices:
!>
!>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
!>    A = [ -----|----- ]  with n1 = n/2
!>        [  A21 | A22  ]       n2 = n-n1
!>
!> The subroutine calls itself to factor A11. Update and scale A21
!> or A12, update A22 then call itself to factor A22.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, if INFO = 0, the factor U or L from the Cholesky
!>          factorization A = U**H*U or A = L*L**H.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the leading minor of order i is not
!>                positive definite, and the factorization could not be
!>                completed.
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
!> \ingroup complex16POcomputational
!
!  =====================================================================
!
!   In the original lapack routine this subroutine is recursive.
!   Recursive routines are not allowed in device, so we substitute
!   it with a non recursive Cholesky algorithm (from wikipedia)
!
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) SUBROUTINE ZPOTRF2_X_XG( UPLO, N, A, LDA, INFO )
#else
      SUBROUTINE ZPOTRF2_X_XG( UPLO, N, A, LDA, INFO )
#endif

#include<lsame_interf.f90>
#include<xerbla_interf.f90>
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER, PARAMETER :: DP=KIND(1.D0)
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
!     ..
!
!  =====================================================================
!
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, J, K
      COMPLEX*16         ASUM
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME_XG
!      EXTERNAL           LSAME_XG
!     ..
!     .. External Subroutines ..
!      EXTERNAL           XERBLA_XG !     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, DBLE, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
#if defined(__CUDA)
      ATTRIBUTES(VALUE) :: UPLO, N, LDA
      ATTRIBUTES(DEVICE) ::  A, INFO
#endif
      
      INFO = 0
      UPPER = LSAME_XG( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME_XG( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         ! CALL XERBLA_XG( 'ZPOTRF2', -INFO )
         CALL XERBLA_XG( 'P', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 )   & 
         RETURN
!
      IF (UPPER) THEN
         DO i = 1, N 
            DO j = i+1, N
               A(j,i) = CONJG(A(i,j))
            END DO
         END DO
      END IF
!
!   The algorithm assumes that the matrix is in the lower part of the 
!   matrix and saves the output triangular matrix on the upper part.
!
      DO i=1, N
         DO j=1, i
            asum=(0.0_DP,0.0_DP)
            DO k=1, j-1
               asum=asum + CONJG(A(k,i)) * A(k,j)
            ENDDO
            IF (i==j) THEN
               A(i,i) = CMPLX(SQRT(DBLE(A(i,i)-asum)), 0.0_DP, KIND=DP)
            ELSE
               A(j,i) = CONJG(A(i,j) - asum) / A(j,j)
            ENDIF
         ENDDO
      ENDDO
!
!   If the lower part is requested in output, reconstruct it from the
!   upper part.
!
      IF (.NOT.UPPER) THEN
         DO i = 1, N 
           DO j = 1, i-1
               A(i,j) = CONJG(A(j,i))
            END DO
         END DO
      END IF   
!
!     End of ZPOTRF2_X_XG
!
      END
