*> \brief \b ZPOTRF2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       RECURSIVE SUBROUTINE ZPOTRF2( UPLO, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZPOTRF2 computes the Cholesky factorization of a real symmetric
*> positive definite matrix A using the recursive algorithm.
*>
*> The factorization has the form
*>    A = U**H * U,  if UPLO = 'U', or
*>    A = L  * L**H,  if UPLO = 'L',
*> where U is an upper triangular matrix and L is lower triangular.
*>
*> This is the recursive version of the algorithm. It divides
*> the matrix into four submatrices:
*>
*>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
*>    A = [ -----|----- ]  with n1 = n/2
*>        [  A21 | A22  ]       n2 = n-n1
*>
*> The subroutine calls itself to factor A11. Update and scale A21
*> or A12, update A22 then call itself to factor A22.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangle of A is stored;
*>          = 'L':  Lower triangle of A is stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*>          N-by-N upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading N-by-N lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the strictly upper
*>          triangular part of A is not referenced.
*>
*>          On exit, if INFO = 0, the factor U or L from the Cholesky
*>          factorization A = U**H*U or A = L*L**H.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, the leading minor of order i is not
*>                positive definite, and the factorization could not be
*>                completed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16POcomputational
*
*  =====================================================================
      RECURSIVE SUBROUTINE zpotrf2( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( lda, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter                ( one = 1.0d+0, zero = 0.0d+0 )
      COMPLEX*16         CONE
      parameter                ( cone = (1.0d+0, 0.0d+0) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            N1, N2, IINFO
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           lsame, disnan
*     ..
*     .. External Subroutines ..
      EXTERNAL           zherk, ztrsm, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, dble, sqrt
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      info = 0
      upper = lsame( uplo, 'U' )
      IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZPOTRF2', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
*     N=1 case
*
      IF( n.EQ.1 ) THEN
*
*        Test for non-positive-definiteness
*
         ajj = dble( a( 1, 1 ) )
         IF( ajj.LE.zero.OR.disnan( ajj ) ) THEN
            info = 1
            RETURN
         END IF
*
*        Factor
*
         a( 1, 1 ) = sqrt( ajj )
*
*     Use recursive code
*
      ELSE
         n1 = n/2
         n2 = n-n1
*
*        Factor A11
*
         CALL zpotrf2( uplo, n1, a( 1, 1 ), lda, iinfo )
         IF ( iinfo.NE.0 ) THEN
            info = iinfo
            RETURN
         END IF
*
*        Compute the Cholesky factorization A = U**H*U
*
         IF( upper ) THEN
*
*           Update and scale A12
*
            CALL ztrsm( 'L', 'U', 'C', 'N', n1, n2, cone,
     $                  a( 1, 1 ), lda, a( 1, n1+1 ), lda )
*
*           Update and factor A22
*
            CALL zherk( uplo, 'C', n2, n1, -one, a( 1, n1+1 ), lda,
     $                  one, a( n1+1, n1+1 ), lda )
            CALL zpotrf2( uplo, n2, a( n1+1, n1+1 ), lda, iinfo )
            IF ( iinfo.NE.0 ) THEN
               info = iinfo + n1
               RETURN
            END IF
*
*        Compute the Cholesky factorization A = L*L**H
*
         ELSE
*
*           Update and scale A21
*
            CALL ztrsm( 'R', 'L', 'C', 'N', n2, n1, cone,
     $                  a( 1, 1 ), lda, a( n1+1, 1 ), lda )
*
*           Update and factor A22
*
            CALL zherk( uplo, 'N', n2, n1, -one, a( n1+1, 1 ), lda,
     $                  one, a( n1+1, n1+1 ), lda )
            CALL zpotrf2( uplo, n2, a( n1+1, n1+1 ), lda, iinfo )
            IF ( iinfo.NE.0 ) THEN
               info = iinfo + n1
               RETURN
            END IF
         END IF
      END IF
      RETURN
*
*     End of ZPOTRF2
*
      END