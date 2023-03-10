//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_tbmv.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?tbmv performs the matrix-vector operation:
 *    x=a*x or x=a'*x.
 *
 * \synopsis
 *    void ?tbmv(uplo, trans, diag, n, k, a, lda, x, incx)
 *    int           incx, k, lda, n
 *    char          diag, trans, uplo
 *    real-floating a(lda, *), x(*)
 *
 * \purpose
 *    ?tbmv performs the matrix-vector operation: x=a*x or x=a'*x where `x` is an n element vector and
 *    `a` is an n by n unit or non-unit, upper or lower triangular band matrix, with (k + 1) diagonals.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the matrix is an upper or lower triangular matrix as follows:
 *    uplo='U' or 'u', `a` is an upper triangular matrix.
 *    uplo='L' or 'l', `a` is a lower triangular matrix.
 *
 *    [in] trans - char. Specifies the operation to be performed as follows:
 *    trans='N' or 'n' x=a*x.
 *    trans='T' or 't' x=a'*x.
 *    trans='C' or 'c' x=a'*x.
 *
 *    [in] diag  - char. Specifies whether or not `a` is unit triangular as follows:
 *    diag='U' or 'u' `a` is assumed to be unit triangular.
 *    diag='N' or 'n' `a` is not assumed to be unit triangular.
 *
 *    [in] n     - int. Specifies the order of the matrix `a`, n must be at least zero.
 *    [in] k     - int. With uplo='U' or 'u', k specifies the number of super-diagonals of the matrix `a`.
 *    i.e the half-bandwidth of matrix `a`. With uplo='L' or 'l', k specifies the number of sub-diagonals
 *    of the matrix `a`.k must satisfy 0 < k.
 *
 *    [in] a     - real-floating array of dimension (lda, n).
 *    With uplo='U' or 'u', the leading (k + 1) by n part of the array `a` must contain the upper triangular band part of
 *    the matrix `a`. This matrix must be supplied column-by-column with the main diagonal of the matrix in row (k) of the
 *    array, the first super-diagonal starting at position 1 in row (k - 1), and so on. Elements in the array `a` which do
 *    not correspond to elements in the triangular band matrix (such as the top left k by k triangle) are not referenced.
 *
 *    @c-layout: the leading (k + 1) by n part of array `a` must contain the upper triangular band part of the matrix of
 *    coefficients. The matrix must be supplied row-by-row, with the leading diagonal of the matrix in column 0 of the array,
 *    the first super-diagonal starting at position 0 in column 1, and so on. The bottom right k by k triangle of array `a`
 *    is not referenced.
 *
 *    With uplo='L' or 'l', the leading (k + 1) by n part of the array `a` must contain the upper triangular band part of
 *    the matrix `a`. This matrix must be supplied column by column with the main diagonal of the matrix in row 0 of the
 *    array, the first sub-diagonal starting at position 0 in row 1, and so on. Elements in the array `a` that do not
 *    correspond to elements in the triangular band matrix (such as the bottom right k by k triangle) are not referenced.
 *
 *    @c-layout: the leading (k + 1) by n part of array `a` must contain the lower triangular band part of the matrix of
 *    coefficients, supplied row by row, with the leading diagonal of the matrix in column k of the array, the first sub-
 *    diagonal starting at position 1 in column (k - 1), and so on. The top left k by k triangle of array `a` is not
 *    referenced.
 *
 *    @see: \examples section about Triangular Band Matrix storage.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. lda must be at least (k + 1).
 *
 *    [out] x    - real-floating array of size at least (1+(n-1)*abs(incx)). The incremented array `x` must contain the
 *    n element vector `x`. `x` is overwritten with the transformed vector `x`.
 *
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 * \examples
 *     - A triangular band matrix `a` of n rows and n columns with k sub/super-diagonals and leading dimension lda. @note: using
 *       @c-layout, be aware that the true-memory leading or split dimension becomes the number of columns (The leading row syntax
 *       is only kept for the understanding of the transposed applied operations in reference to Fortran). For `a` a m by n matrix,
 *       the rational subscript operation should be as the following:
 *
 *       int row, col;
 *       for (row = 0; row < m; row++) {
 *          for (col = 0; col < n; col++) {
 *             a[(n * row) + col] = ...;
 *          }
 *       }
 *
 *     - Given uplo='U' or 'u': the following program segment transfers a band matrix from
 *       conventional full matrix storage b(ldb, *) to band storage a(lda, *):
 *
 *       real-floating b[ldb * n] = { ... };
 *       real-floating a[lda * n] = { 0 };
 *
 *       int i, j, m, k;
 *
 *       @c-fortan-layout: the following code is internal memory-storage independent.
 *       for (j = 1; j <= n; ++j) {
 *          m = k + 1 - j
 *          for (i = mc_maxmag(1, j - k); i <= j; ++i) {
 *             mc_blas_matrix_at(a, lda, n, m + i, j) = mc_blas_matrix_at(b, ldb, n, i, j);
 *          }
 *       }
 *       @fortan-layout: column-major + indexation starting at 0.
 *       for (j = 0; j < n; j++) {
 *          m = k - j;
 *          for (i = mc_maxmag(0, j - k); i <= j; i++) {
 *             a[(m + i) + j * lda] = b[i + j * ldb];
 *          }
 *       }
 *       @c-layout: row-major + indexation starting at 0.
 *       for (i = 0; i < n; i++) {
 *          m = -i;
 *          for (j = i; j < mc_minmag(n, i + k + 1); j++) {
 *             a[(m + j) + i * lda] = b[j + i * ldb];
 *          }
 *       }
 *
 *     - Given uplo='L' or 'l': the following program segment transfers a band matrix from
 *       conventional full matrix storage b(ldb, *) to band storage a(lda, *):
 *
 *       real-floating b[ldb * n] = { ... };
 *       real-floating a[lda * n] = { 0 };
 *
 *       int i, j, m, k;
 *
 *       @c-fortan-layout: the following code is internal memory-storage independent.
 *       for (j = 1; j <= n; ++j) {
 *          m = 1 - j
 *          for (i = j; i <= mc_minmag(n, j + k); ++i) {
 *             mc_blas_matrix_at(a, lda, n, m + i, j) = mc_blas_matrix_at(b, ldb, n, i, j);
 *          }
 *       }
 *       @fortan-layout: column-major + indexation starting at 0.
 *       for (j = 0; j < n; j++) {
 *          m = -j;
 *          for (i = j; i < mc_minmag(n, j + k + 1); i++) {
 *             a[(m + i) + j * lda] = b[i + j * ldb];
 *          }
 *       }
 *       @c-layout: row-major + indexation starting at 0.
 *       for (i = 0; i < n; i++) {
 *          m = k - i;
 *          for (j = mc_maxmag(0, i - k); j <= i; j++) {
 *             a[(m + j) + i * lda] = b[j + i * ldb];
 *          }
 *       }
 *
 * \level 2 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Argonne National Lab.
 *     \author Jeremy Du Croz, Nag Central Office.
 *     \author Sven Hammarling, Nag Central Office.
 *     \author Richard Hanson, Sandia National Labs.
 */

#include <macadam/lapack/blas/mc_blas_access.h>
#include <macadam/lapack/blas/mc_blas_lsame.h>

#ifndef MC_BLAS_NATIVE_TBMV_H
#define MC_BLAS_NATIVE_TBMV_H

#pragma mark - mc_blas_native_stbmv -

MC_TARGET_FUNC void mc_blas_native_stbmv(const char uplo, const char trans, const char diag, const int n, const int k, const float * a, const int lda, float * x, const int incx)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul      = mc_blas_lsame(uplo, 'U')  ? CblasUpper   : CblasLower;
	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);
	const enum CBLAS_DIAG di      = mc_blas_lsame(diag, 'U')  ? CblasUnit    : CblasNonUnit;

#	if MC_TARGET_CPP98
	::cblas_stbmv(ord, ul, ta, di, n, k, a, lda, x, incx);
#	else
	cblas_stbmv(ord, ul, ta, di, n, k, a, lda, x, incx);
#	endif
}

#pragma mark - mc_blas_native_dtbmv -

MC_TARGET_FUNC void mc_blas_native_dtbmv(const char uplo, const char trans, const char diag, const int n, const int k, const double * a, const int lda, double * x, const int incx)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul      = mc_blas_lsame(uplo, 'U')  ? CblasUpper   : CblasLower;
	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);
	const enum CBLAS_DIAG di      = mc_blas_lsame(diag, 'U')  ? CblasUnit    : CblasNonUnit;

#	if MC_TARGET_CPP98
	::cblas_dtbmv(ord, ul, ta, di, n, k, a, lda, x, incx);
#	else
	cblas_dtbmv(ord, ul, ta, di, n, k, a, lda, x, incx);
#	endif
}

/* \name
 *    ?tbmv performs the matrix-vector operation:
 *    x=a*x or x=a'*x or x=a_*x.
 *
 * \synopsis
 *    void ?tbmv(uplo, trans, diag, n, k, a, lda, x, incx)
 *    int     incx, k, lda, n
 *    char    diag, trans, uplo
 *    complex a(lda, *), x(*)
 *
 * \purpose
 *    ?tbmv performs the matrix-vector operation: x=a*x or x=a'*x or x=a_*x where `x` is an n element vector
 *    and `a` is an n by n unit or non-unit, upper or lower triangular band matrix, with (k + 1) diagonals.
 *
 * \parameters
 *    [in] uplo  - char. Specifies whether the matrix is an upper or lower triangular matrix as follows:
 *    uplo='U' or 'u', `a` is an upper triangular matrix.
 *    uplo='L' or 'l', `a` is a lower triangular matrix.
 *
 *    [in] trans - char. Specifies the operation to be performed as follows:
 *    trans='N' or 'n' x=a*x.
 *    trans='T' or 't' x=a'*x.
 *    trans='C' or 'c' x=a_*x.
 *
 *    [in] diag  - char. Specifies whether or not `a` is unit triangular as follows:
 *    diag='U' or 'u' `a` is assumed to be unit triangular.
 *    diag='N' or 'n' `a` is not assumed to be unit triangular.
 *
 *    [in] n     - int. Specifies the order of the matrix `a`, n must be at least zero.
 *    [in] k     - int. With uplo='U' or 'u', k specifies the number of super-diagonals of the matrix `a`.
 *    i.e the half-bandwidth of matrix `a`. With uplo='L' or 'l', k specifies the number of sub-diagonals
 *    of the matrix `a`.k must satisfy 0 < k.
 *
 *    [in] a     - complex array of dimension (lda, n).
 *    With uplo='U' or 'u', the leading (k + 1) by n part of the array `a` must contain the upper triangular band part of
 *    the matrix `a`. This matrix must be supplied column-by-column with the main diagonal of the matrix in row (k) of the
 *    array, the first super-diagonal starting at position 1 in row (k - 1), and so on. Elements in the array `a` which do
 *    not correspond to elements in the triangular band matrix (such as the top left k by k triangle) are not referenced.
 *
 *    @c-layout: the leading (k + 1) by n part of array `a` must contain the upper triangular band part of the matrix of
 *    coefficients. The matrix must be supplied row-by-row, with the leading diagonal of the matrix in column 0 of the array,
 *    the first super-diagonal starting at position 0 in column 1, and so on. The bottom right k by k triangle of array `a`
 *    is not referenced.
 *
 *    With uplo='L' or 'l', the leading (k + 1) by n part of the array `a` must contain the upper triangular band part of
 *    the matrix `a`. This matrix must be supplied column by column with the main diagonal of the matrix in row 0 of the
 *    array, the first sub-diagonal starting at position 0 in row 1, and so on. Elements in the array `a` that do not
 *    correspond to elements in the triangular band matrix (such as the bottom right k by k triangle) are not referenced.
 *
 *    @c-layout: the leading (k + 1) by n part of array `a` must contain the lower triangular band part of the matrix of
 *    coefficients, supplied row by row, with the leading diagonal of the matrix in column k of the array, the first sub-
 *    diagonal starting at position 1 in column (k - 1), and so on. The top left k by k triangle of array `a` is not
 *    referenced.
 *
 *    @see: \examples section about Triangular Band Matrix storage.
 *
 *    [in] lda   - int. Specifies the first dimension of `a`. lda must be at least (k + 1).
 *
 *    [out] x    - complex array of size at least (1+(n-1)*abs(incx)). The incremented array `x` must contain the
 *    n element vector `x`. `x` is overwritten with the transformed vector `x`.
 *
 *    [in]  incx - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 * \examples
 *     - A triangular band matrix `a` of n rows and n columns with k sub/super-diagonals and leading dimension lda. @note: using
 *       @c-layout, be aware that the true-memory leading or split dimension becomes the number of columns (The leading row syntax
 *       is only kept for the understanding of the transposed applied operations in reference to Fortran). For `a` a m by n matrix,
 *       the rational subscript operation should be as the following:
 *
 *       int row, col;
 *       for (row = 0; row < m; row++) {
 *          for (col = 0; col < n; col++) {
 *             a[(n * row) + col] = ...;
 *          }
 *       }
 *
 *     - Given uplo='U' or 'u': the following program segment transfers a band matrix from
 *       conventional full matrix storage b(ldb, *) to band storage a(lda, *):
 *
 *       complex b[ldb * n] = { ... };
 *       complex a[lda * n] = { 0 };
 *
 *       int i, j, m, k = lda - 1; // must statisfy ldb >= lda >= n > k.
 *
 *       @c-fortan-layout: the following code is internal memory-storage independent.
 *       for (j = 1; j <= n; ++j) {
 *          m = k + 1 - j
 *          for (i = mc_maxmag(1, j - k); i <= j; ++i) {
 *             mc_blas_matrix_at(a, lda, n, m + i, j) = mc_blas_matrix_at(b, ldb, n, i, j);
 *          }
 *       }
 *       @fortan-layout: column-major + indexation starting at 0.
 *       for (j = 0; j < n; j++) {
 *          m = k - j;
 *          for (i = mc_maxmag(0, j - k); i <= j; i++) {
 *             a[(m + i) + j * lda] = b[i + j * ldb];
 *          }
 *       }
 *       @c-layout: row-major + indexation starting at 0.
 *       for (i = 0; i < n; i++) {
 *          m = -i;
 *          for (j = i; j < mc_minmag(n, i + k + 1); j++) {
 *             a[(m + j) + i * lda] = b[j + i * ldb];
 *          }
 *       }
 *
 *     - Given uplo='L' or 'l': the following program segment transfers a band matrix from
 *       conventional full matrix storage b(ldb, *) to band storage a(lda, *):
 *
 *       complex b[ldb * n] = { ... };
 *       complex a[lda * n] = { 0 };
 *
 *       int i, j, m, k = lda - 1; // must statisfy ldb >= lda >= n > k.
 *
 *       @c-fortan-layout: the following code is internal memory-storage independent.
 *       for (j = 1; j <= n; ++j) {
 *          m = 1 - j
 *          for (i = j; i <= mc_minmag(n, j + k); ++i) {
 *             mc_blas_matrix_at(a, lda, n, m + i, j) = mc_blas_matrix_at(b, ldb, n, i, j);
 *          }
 *       }
 *       @fortan-layout: column-major + indexation starting at 0.
 *       for (j = 0; j < n; j++) {
 *          m = -j;
 *          for (i = j; i < mc_minmag(n, j + k + 1); i++) {
 *             a[(m + i) + j * lda] = b[i + j * ldb];
 *          }
 *       }
 *       @c-layout: row-major + indexation starting at 0.
 *       for (i = 0; i < n; i++) {
 *          m = k - i;
 *          for (j = mc_maxmag(0, i - k); j <= i; j++) {
 *             a[(m + j) + i * lda] = b[j + i * ldb];
 *          }
 *       }
 *
 * \level 2 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Jack Dongarra, Argonne National Lab.
 *     \author Jeremy Du Croz, Nag Central Office.
 *     \author Sven Hammarling, Nag Central Office.
 *     \author Richard Hanson, Sandia National Labs.
 */

#pragma mark - mc_blas_native_ctbmv -

MC_TARGET_FUNC void mc_blas_native_ctbmv(const char uplo, const char trans, const char diag, const int n, const int k, const mc_complex_float_t * a, const int lda, mc_complex_float_t * x, const int incx)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul      = mc_blas_lsame(uplo, 'U')  ? CblasUpper   : CblasLower;
	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);
	const enum CBLAS_DIAG di      = mc_blas_lsame(diag, 'U')  ? CblasUnit    : CblasNonUnit;

#	if MC_TARGET_CPP98
	::cblas_ctbmv(ord, ul, ta, di, n, k, a, lda, x, incx);
#	else
	cblas_ctbmv(ord, ul, ta, di, n, k, a, lda, x, incx);
#	endif
}

#pragma mark - mc_blas_native_ztbmv -

MC_TARGET_FUNC void mc_blas_native_ztbmv(const char uplo, const char trans, const char diag, const int n, const int k, const mc_complex_double_t * a, const int lda, mc_complex_double_t * x, const int incx)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_UPLO ul      = mc_blas_lsame(uplo, 'U')  ? CblasUpper   : CblasLower;
	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);
	const enum CBLAS_DIAG di      = mc_blas_lsame(diag, 'U')  ? CblasUnit    : CblasNonUnit;

#	if MC_TARGET_CPP98
	::cblas_ztbmv(ord, ul, ta, di, n, k, a, lda, x, incx);
#	else
	cblas_ztbmv(ord, ul, ta, di, n, k, a, lda, x, incx);
#	endif
}

#endif /* !MC_BLAS_NATIVE_TBMV_H */

/* EOF */