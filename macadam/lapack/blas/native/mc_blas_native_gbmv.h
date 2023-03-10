//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_native_gbmv.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?gbmv - performs one of the matrix-vector operations:
 *    y=alpha*a*x + beta*y or y=alpha*a'*x + beta*y.
 *
 * \synopsis
 *    void ?gbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
 *    real-floating alpha, beta
 *    int           incx, incy, kl, ku, lda, m, n
 *    char          trans
 *    real-floating a(lda,*), x(*), y(*)
 *
 * \purpose
 *   ?gbmv performs one of the matrix-vector operations: y=alpha*a*x + beta*y or y=alpha*a'*x + beta*y where alpha and beta
 *   are scalars, `x` and `y` are vectors and `a`is an m by n band matrix, with kl sub-diagonals and ku super-diagonals.
 *
 * \parameters
 *    [in] trans - char. Specifies the operation to be performed as follows:
 *    trans='N' or 'n' y=alpha*a*x  + beta*y.
 *    trans='T' or 't' y=alpha*a'*x + beta*y.
 *    trans='C' or 'c' y=alpha*a'*x + beta*y.
 *
 *    [in] m     - int. Specifies the number of rows of the matrix `a`, m must be at least zero.
 *    [in] n     - int. Specifies the number of columns of the matrix `a`, n must be at least zero.
 *    [in] kl    - int. Specifies the number of sub-diagonals of the matrix `a`. kl must satisfy 0 < kl.
 *    [in] ku    - int. Specifies the number of super-diagonals of the matrix `a`. ku must satisfy 0 < ku.
 *
 *    [in] alpha - real-floating. Specifies the scalar alpha.
 *
 *    [in] a     - real-floating array of dimension (lda, n). The leading (kl + ku + 1) by n part of
 *    the array `a` must contain the matrix of coefficients, supplied column by column, with the leading
 *    diagonal of the matrix in row (ku + 1) of the array, the first super-diagonal starting at position
 *    1 in row ku, the first sub-diagonal starting at position 0 in row (ku + 1), and so on. Elements in
 *    the array `a` that do not correspond to elements in the band matrix (such as the top left ku by ku
 *    triangle) are not referenced.
 *
 *    @c-layout: the leading (kl + ku + 1) by m part of the array `a` must contain the matrix of coefficients.
 *    This matrix must be supplied row-by-row, with the leading diagonal of the matrix in column (kl) of the
 *    array, the first super-diagonal starting at position 0 in column (kl + 1), the first sub-diagonal starting
 *    at position 1 in row (kl - 1), and so on. Elements in the array `a` that do not correspond to elements in
 *    the band matrix (such as the top left kl by kl triangle) are not referenced.
 *
 *    @see: \examples section about Band Matrix storage.
 *
 *    [in] lda   - int. Specifies the first dimension of a band matrix, lda must be at least (kl + ku + 1).
 *
 *    [in] x     - real-floating array of dimension (at least) (1+(n-1)*abs(incx)) when trans = 'N' or 'n'
 *    and at least (1+(m-1)*abs(incx)) otherwise. The incremented array `x` must contain the vector `x`.
 *
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in] beta  - real-floating. Specifies the scalar beta. when beta is supplied as zero then `y` need not
 *    be set on input.
 *
 *    [out] y    - real-floating array of dimension (at least) (1+(m-1)*abs(incy)) when trans = 'N' or 'n'
 *    and at least (1+(n-1)*abs(incy)) otherwise. The incremented array `y` must contain the vector `y`, y
 *    is overwritten by the updated vector `y`.
 *
 *    [in] incy  - int. Specifies the increment for the elements of `y`. incy must not be zero.
 *
 * \examples
 *     - A general band matrix `a` of m rows and n columns with kl sub-diagonals, ku super-diagonals, and leading dimension
 *       lda. @note: using @c-layout, be aware that the true-memory leading or split dimension becomes the number of columns
 *       (The leading row syntax is only kept for the understanding of the transposed applied operations in reference to
 *       Fortran). For `a` a m by n matrix, the rational subscript operation should be as the following:
 *
 *       int row, col;
 *       for (row = 0; row < m; row++) {
 *          for (col = 0; col < n; col++) {
 *             a[(n * row) + col] = ...;
 *          }
 *       }
 *
 *     - The following program segment transfers a band matrix from conventional full matrix storage b(ldb, *) to band
 *       storage a(lda, *):
 *
 *       real-floating b[ldb * n] = { ... };
 *       real-floating a[lda * n] = { 0 };
 *
 *       int i, j, ku, kl, m;
 *
 *       @c-fortan-layout: the following code is internal memory-storage independent.
 *       for (j = 1; j <= n; ++j) {
 *          k = ku + 1 - j;
 *          for (i = mc_maxmag(1, j - ku); i <= mc_minmag(m, j + kl); ++i) {
 *             mc_blas_matrix_at(a, lda, n, k + i, j) = mc_blas_matrix_at(b, ldb, n, i, j);
 *          }
 *       }
 *       @fortan-layout: column-major + indexation starting at 0.
 *       for (j = 0; j < n; j++) {
 *          k = ku - j;
 *          for (i = mc_maxmag(0, j - ku); i < mc_minmagin(m, j + kl + 1); i++) {
 *             a[(k + i) + j * lda] = b[i + j * ldb];
 *          }
 *       }
 *       @c-layout: row-major + indexation starting at 0.
 *       for (i = 0; i < m; i++) {
 *          k = kl - i;
 *          for (j = mc_maxmag(0, i - kl); j < mc_minmagin(n, i + ku + 1); j++) {
 *             a[(k + j) + i * lda] = b[j + i * ldb];
 *          }
 *       }
 *
 *              | 1 1 1 0 |
 *              | 2 2 2 2 |
 *     a[5x4] = | 3 3 3 3 |
 *              | 4 4 4 4 |
 *              | 0 5 5 5 |
 *
 *     const real-floating a_band[] = {
 *          0, 0, 1, 2
 *        , 0, 1, 2, 3
 *        , 1, 2, 3, 4
 *        , 2, 3, 4, 5
 *        , 3, 4, 5, 0
 *        , 4, 5, 0, 0
 *        , 0, 0, 0, 0
 *        , 0, 0, 0, 0
 *     };
 *     const real-floating x[] = { 1, 2, 3, 4 };
 *           real-floating y[] = { 1, 0 , 2, 0 , 3, 0 , 4, 0 , 5, 0 };
 *     mc_blas_?gbmv('N', 5, 4, 3, 2, 2, a, 8, x, 1, 10, y, 2);
 *     on output -> y = { 22, 0, 60, 0, 90, 0, 120, 0, 140, 0 }
 *
 *              | 1 1 1 1 1 |
 *     a[4x5] = | 2 2 2 2 2 |
 *              | 3 3 3 3 3 |
 *              | 4 4 4 4 4 |
 *
 *     const real-floating a_band[] = {
 *          0, 0, 0, 0, 0
 *        , 0, 0, 0, 0, 1
 *        , 0, 0, 0, 1, 2
 *        , 0, 0, 1, 2, 3
 *        , 0, 1, 2, 3, 4
 *        , 1, 2, 3, 4, 0
 *        , 2, 3, 4, 0, 0
 *        , 3, 4, 0, 0, 0
 *        , 4, 0, 0, 0, 0
 *        , 0, 0, 0, 0, 0
 *        , 0, 0, 0, 0, 0
 *        , 0, 0, 0, 0, 0
 *     };
 *     const real-floating x[] = { 1, 2, 3, 4, 5 };
 *           real-floating y[] = { 1, 0, 2, 0, 3, 0, 4, 0 };
 *     mc_blas_?gbmv('N', 4, 5, 6, 5, 2, a_band, 12, x, 1, 10, y, 2);
 *     on output -> y = { 40, 0 , 80, 0 , 120, 0 , 160, 0 }
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

#include <macadam/lapack/blas/mc_blas_lsame.h>

#ifndef MC_BLAS_NATIVE_GBMV_H
#define MC_BLAS_NATIVE_GBMV_H

#pragma mark - mc_blas_native_sgbmv -

MC_TARGET_FUNC void mc_blas_native_sgbmv(const char trans, const int m, const int n, const int kl, const int ku, const float alpha, const float * a, const int lda, const float * x, const int incx, const float beta, float * y, const int incy)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_sgbmv(ord, ta, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
#	else
	cblas_sgbmv(ord, ta, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
#	endif
}

#pragma mark - mc_blas_native_dgbmv -

MC_TARGET_FUNC void mc_blas_native_dgbmv(const char trans, const int m, const int n, const int kl, const int ku, const double alpha, const double * a, const int lda, const double * x, const int incx, const double beta, double * y, const int incy)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_dgbmv(ord, ta, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
#	else
	cblas_dgbmv(ord, ta, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
#	endif
}

/* \name
 *    ?gbmv - performs one of the matrix-vector operations:
 *    y=alpha*a*x + beta*y or y=alpha*a'*x + beta*y or y=alpha*a_*x + beta*y.
 *
 * \synopsis
 *    void ?gbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
 *    complex alpha, beta
 *    int     incx, incy, kl, ku, lda, m, n
 *    char    trans
 *    complex a(lda,*), x(*), y(*)
 *
 * \purpose
 *   ?gbmv performs one of the matrix-vector operations: y=alpha*a*x + beta*y or y=alpha*a'*x + beta*y or y=alpha*a_*x + beta*y where
 *   alpha and beta are scalars, `x` and `y` are vectors and `a`is an m by n band matrix, with kl sub-diagonals and ku super-diagonals.
 *
 * \parameters
 *    [in] trans - char. Specifies the operation to be performed as follows:
 *    trans='N' or 'n' y=alpha*a*x  + beta*y.
 *    trans='T' or 't' y=alpha*a'*x + beta*y.
 *    trans='C' or 'c' y=alpha*a_*x + beta*y.
 *
 *    [in] m     - int. Specifies the number of rows of the matrix `a`, m must be at least zero.
 *    [in] n     - int. Specifies the number of columns of the matrix `a`, n must be at least zero.
 *    [in] kl    - int. Specifies the number of sub-diagonals of the matrix `a`. kl must satisfy 0 < kl.
 *    [in] ku    - int. Specifies the number of super-diagonals of the matrix `a`. ku must satisfy 0 < ku.
 *
 *    [in] alpha - complex. Specifies the scalar alpha.
 *
 *    [in] a     - complex array of dimension (lda, n). The leading (kl + ku + 1) by n part of
 *    the array `a` must contain the matrix of coefficients, supplied column by column, with the leading
 *    diagonal of the matrix in row (ku + 1) of the array, the first super-diagonal starting at position
 *    1 in row ku, the first sub-diagonal starting at position 0 in row (ku + 1), and so on. Elements in
 *    the array `a` that do not correspond to elements in the band matrix (such as the top left ku by ku
 *    triangle) are not referenced.
 *
 *    @c-layout: the leading (kl + ku + 1) by m part of the array `a` must contain the matrix of coefficients.
 *    This matrix must be supplied row-by-row, with the leading diagonal of the matrix in column (kl) of the
 *    array, the first super-diagonal starting at position 0 in column (kl + 1), the first sub-diagonal starting
 *    at position 1 in row (kl - 1), and so on. Elements in the array `a` that do not correspond to elements in
 *    the band matrix (such as the top left kl by kl triangle) are not referenced.
 *
 *    @see: \examples section about Band Matrix storage.
 *
 *    [in] lda   - int. Specifies the first dimension of a band matrix, lda must be at least (kl + ku + 1).
 *
 *    [in] x     - complex array of dimension (at least) (1+(n-1)*abs(incx)) when trans = 'N' or 'n'
 *    and at least (1+(m-1)*abs(incx)) otherwise. The incremented array `x` must contain the vector `x`.
 *
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 *    [in] beta  - complex. Specifies the scalar beta. when beta is supplied as zero then `y` need not
 *    be set on input.
 *
 *    [out] y    - complex array of dimension (at least) (1+(m-1)*abs(incy)) when trans = 'N' or 'n'
 *    and at least (1+(n-1)*abs(incy)) otherwise. The incremented array `y` must contain the vector `y`, y
 *    is overwritten by the updated vector `y`.
 *
 *    [in] incy  - int. Specifies the increment for the elements of `y`. incy must not be zero.
 *
 * \examples
 *     - A general band matrix `a` of m rows and n columns with kl sub-diagonals, ku super-diagonals, and leading dimension
 *       lda. @note: using @c-layout, be aware that the true-memory leading or split dimension becomes the number of columns
 *       (The leading row syntax is only kept for the understanding of the transposed applied operations in reference to
 *       Fortran). For `a` a m by n matrix, the rational subscript operation should be as the following:
 *
 *       int row, col;
 *       for (row = 0; row < m; row++) {
 *          for (col = 0; col < n; col++) {
 *             a[(n * row) + col] = ...;
 *          }
 *       }
 *
 *     - The following program segment transfers a band matrix from conventional full matrix storage b(ldb, *) to band
 *       storage a(lda, *):
 *
 *       real-floating b[ldb * n] = { ... };
 *       real-floating a[lda * n] = { 0 };
 *
 *       int i, j, ku, kl, m;
 *
 *       @c-fortan-layout: the following code is internal memory-storage independent.
 *       for (j = 1; j <= n; ++j) {
 *          k = ku + 1 - j;
 *          for (i = mc_maxmag(1, j - ku); i <= mc_minmag(m, j + kl); ++i) {
 *             mc_blas_matrix_at(a, lda, n, k + i, j) = mc_blas_matrix_at(b, ldb, n, i, j);
 *          }
 *       }
 *       @fortan-layout: column-major + indexation starting at 0.
 *       for (j = 0; j < n; j++) {
 *          k = ku - j;
 *          for (i = mc_maxmag(0, j - ku); i < mc_minmagin(m, j + kl + 1); i++) {
 *             a[(k + i) + j * lda] = b[i + j * ldb];
 *          }
 *       }
 *       @c-layout: row-major + indexation starting at 0.
 *       for (i = 0; i < m; i++) {
 *          k = kl - i;
 *          for (j = mc_maxmag(0, i - kl); j < mc_minmagin(n, i + ku + 1); j++) {
 *             a[(k + j) + i * lda] = b[j + i * ldb];
 *          }
 *       }
 *
 *              | 1 1 1 0 |
 *              | 2 2 2 2 |
 *     a[5x4] = | 3 3 3 3 |
 *              | 4 4 4 4 |
 *              | 0 5 5 5 |
 *
 *     const real-floating a_band[] = {
 *          0, 0, 1, 2
 *        , 0, 1, 2, 3
 *        , 1, 2, 3, 4
 *        , 2, 3, 4, 5
 *        , 3, 4, 5, 0
 *        , 4, 5, 0, 0
 *        , 0, 0, 0, 0
 *        , 0, 0, 0, 0
 *     };
 *     const real-floating x[] = { 1, 2, 3, 4 };
 *           real-floating y[] = { 1, 0 , 2, 0 , 3, 0 , 4, 0 , 5, 0 };
 *     mc_blas_?gbmv('N', 5, 4, 3, 2, 2, a, 8, x, 1, 10, y, 2);
 *     on output -> y = { 22, 0, 60, 0, 90, 0, 120, 0, 140, 0 }
 *
 *              | 1 1 1 1 1 |
 *     a[4x5] = | 2 2 2 2 2 |
 *              | 3 3 3 3 3 |
 *              | 4 4 4 4 4 |
 *
 *     const real-floating a_band[] = {
 *          0, 0, 0, 0, 0
 *        , 0, 0, 0, 0, 1
 *        , 0, 0, 0, 1, 2
 *        , 0, 0, 1, 2, 3
 *        , 0, 1, 2, 3, 4
 *        , 1, 2, 3, 4, 0
 *        , 2, 3, 4, 0, 0
 *        , 3, 4, 0, 0, 0
 *        , 4, 0, 0, 0, 0
 *        , 0, 0, 0, 0, 0
 *        , 0, 0, 0, 0, 0
 *        , 0, 0, 0, 0, 0
 *     };
 *     const real-floating x[] = { 1, 2, 3, 4, 5 };
 *           real-floating y[] = { 1, 0, 2, 0, 3, 0, 4, 0 };
 *     mc_blas_?gbmv('N', 4, 5, 6, 5, 2, a_band, 12, x, 1, 10, y, 2);
 *     on output -> y = { 40, 0 , 80, 0 , 120, 0 , 160, 0 }
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

#pragma mark - mc_blas_native_cgbmv -

MC_TARGET_FUNC void mc_blas_native_cgbmv(const char trans, const int m, const int n, const int kl, const int ku, const mc_complex_float_t alpha, const mc_complex_float_t * a, const int lda, const mc_complex_float_t * x, const int incx, const mc_complex_float_t beta, mc_complex_float_t * y, const int incy)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_cgbmv(ord, ta, m, n, kl, ku, &alpha, a, lda, x, incx, &beta, y, incy);
#	else
	cblas_cgbmv(ord, ta, m, n, kl, ku, &alpha, a, lda, x, incx, &beta, y, incy);
#	endif
}

#pragma mark - mc_blas_native_cgbmv -

MC_TARGET_FUNC void mc_blas_native_zgbmv(const char trans, const int m, const int n, const int kl, const int ku, const mc_complex_double_t alpha, const mc_complex_double_t * a, const int lda, const mc_complex_double_t * x, const int incx, const mc_complex_double_t beta, mc_complex_double_t * y, const int incy)
{
#	if MC_TARGET_BLAS_USE_CLAYOUT
	const enum CBLAS_ORDER ord = CblasRowMajor;
#	else
	const enum CBLAS_ORDER ord = CblasColMajor;
#	endif

	const enum CBLAS_TRANSPOSE ta = mc_blas_lsame(trans, 'N') ? CblasNoTrans : (mc_blas_lsame(trans, 'T') ? CblasTrans : CblasConjTrans);

#	if MC_TARGET_CPP98
	::cblas_zgbmv(ord, ta, m, n, kl, ku, &alpha, a, lda, x, incx, &beta, y, incy);
#	else
	cblas_zgbmv(ord, ta, m, n, kl, ku, &alpha, a, lda, x, incx, &beta, y, incy);
#	endif
}

#endif /* !MC_BLAS_NATIVE_GBMV_H */

/* EOF */